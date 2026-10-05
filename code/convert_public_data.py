"""Format conversion only; all model fitting and statistical tests run in Rust.

Run with pinned dependencies in data/tools/format-env. Preserve missing methylation,
cell IDs and independent LARRY labels; never infer clones from CpGs in this converter.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import rdata
import scipy.io
import scipy.sparse
from rdata.parser import RObjectType


def resolve(obj: Any) -> Any:
    while obj.info.type == RObjectType.REF:
        obj = obj.referenced_object
    return obj


def slots(obj: Any) -> dict[str, Any]:
    result: dict[str, Any] = {}
    node = resolve(obj).attributes
    while node is not None:
        node = resolve(node)
        if node.info.type == RObjectType.NILVALUE:
            break
        key = resolve(resolve(node.tag).value).value.decode()
        result[key] = node.value[0]
        node = node.value[1]
    return result


def convert_matrix(obj: Any, converter: Any) -> tuple[Any, list[str], list[str]]:
    obj = resolve(obj)
    attributes = slots(obj)
    if "Dim" in attributes:
        dimensions = resolve(attributes["Dim"]).value
        matrix = scipy.sparse.csc_matrix(
            (
                resolve(attributes["x"]).value,
                resolve(attributes["i"]).value,
                resolve(attributes["p"]).value,
            ),
            shape=tuple(dimensions),
        )
        names = converter.convert(attributes["Dimnames"])
    else:
        dimensions = resolve(attributes["dim"]).value
        matrix = np.asarray(obj.value).reshape(dimensions, order="F")
        names = converter.convert(attributes["dimnames"])
    return matrix, list(names[0]), list(names[1])


def main() -> None:
    target = Path("data/derived")
    target.mkdir(parents=True, exist_ok=True)
    with (target / "maintenance_counts.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["chrom", "position", "m0", "u0", "m1", "u1", "m4", "u4", "m16", "u16"]
        )
        for chrom in [1, 22]:
            with tarfile.open(f"data/raw/maintenance_chr{chrom}.tar.gz") as archive:
                stream = archive.extractfile(f"AllDat_chr{chrom}.mat")
                assert stream is not None
                data = scipy.io.loadmat(stream)
            counts = data["AllDat"]
            positions = data["sites"].ravel()
            assert counts.shape == (len(positions), 4, 2)
            assert np.all(counts >= 0) and np.all(counts == np.floor(counts))
            # Fixed input-index sampling, before any coverage or outcome selection.
            for i in range(0, len(positions), 32):
                writer.writerow(
                    [
                        chrom,
                        int(positions[i]),
                        *counts[i].astype(np.int64).ravel().tolist(),
                    ]
                )

    parsed = rdata.parser.parse_file("data/raw/EPIClone_M123.rds")
    attributes = slots(parsed.object)
    converter = rdata.conversion.SimpleConverter()
    metadata = converter.convert(attributes["meta.data"])
    metadata.index.name = "cell"
    assays = dict(
        zip(
            converter.convert(slots(attributes["assays"])["names"]),
            resolve(attributes["assays"]).value,
        )
    )
    for name, filename in [
        ("DNAm", "clone_methylation.tsv"),
        ("AB", "clone_proteins.tsv"),
    ]:
        matrix, features, cells = convert_matrix(slots(assays[name])["data"], converter)
        assert len(cells) == matrix.shape[1] and len(features) == matrix.shape[0]
        assert cells == list(
            metadata.index
        ), "Seurat metadata/matrix cell order differs"
        matrix = matrix.toarray() if scipy.sparse.issparse(matrix) else matrix
        if name == "DNAm":
            observed = matrix[np.isfinite(matrix)]
            assert np.all(
                (observed == 0) | (observed == 1)
            ), "DNAm data are not binary calls"
        with (target / filename).open("w") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["cell", *features])
            for j, cell in enumerate(cells):
                writer.writerow(
                    [
                        cell,
                        *[
                            format(float(x), ".9g") if np.isfinite(x) else "NA"
                            for x in matrix[:, j]
                        ],
                    ]
                )
    metadata.to_csv(target / "clone_metadata.tsv", sep="\t", na_rep="NA")
    print(
        f"Converted observed kinetics and {len(metadata)} clone-study cells; missing calls preserved."
    )


def canonical_sequence(data: bytes) -> bytes:
    """Discard API retrieval timestamps, preserving reference sequence identity."""
    obj = json.loads(data)
    stable = {
        key: obj[key]
        for key in ["genome", "chrom", "start", "end", "dna"]
        if key in obj
    }
    return (json.dumps(stable, sort_keys=True) + "\n").encode()


def broad_inputs(download: bool) -> None:
    """Acquire frozen public inputs by recorded byte range and verify content."""
    for manifest in ["broad_data_manifest.csv", "broad_reference_manifest.csv"]:
        with Path(f"code/results/{manifest}").open() as handle:
            for row in csv.DictReader(handle):
                path = Path(row["path"])
                if not path.exists() and download:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    partial = path.with_suffix(path.suffix + ".partial")
                    span = row["archive_byte_range"]
                    chunks: list[tuple[int, int] | None] = []
                    if span:
                        start, end = map(int, span.split("-"))
                        chunks = [
                            (a, min(a + 8 * 1024 * 1024 - 1, end))
                            for a in range(start, end + 1, 8 * 1024 * 1024)
                        ]
                    with partial.open("wb") as output:
                        for chunk in chunks or [None]:
                            headers = (
                                {"Range": f"bytes={chunk[0]}-{chunk[1]}"}
                                if chunk
                                else {}
                            )
                            request = urllib.request.Request(
                                row["source_url"], headers=headers
                            )
                            for attempt in range(4):
                                try:
                                    try:
                                        with urllib.request.urlopen(
                                            request, timeout=90
                                        ) as response:
                                            data = response.read()
                                            if (
                                                chunk
                                                and response.headers.get(
                                                    "Content-Range", ""
                                                ).split("/")[0]
                                                != f"bytes {chunk[0]}-{chunk[1]}"
                                            ):
                                                raise ValueError(
                                                    "server did not honor requested archive range"
                                                )
                                    except urllib.error.HTTPError as error:
                                        # One hg38 numeric interval is past chr19's end;
                                        # preserve that negative assembly check as empty JSON.
                                        if (
                                            error.code != 400
                                            or path.name != "guide_region_25_hg38.json"
                                        ):
                                            raise
                                        data = error.read()
                                    if chunk and len(data) != chunk[1] - chunk[0] + 1:
                                        raise ValueError("incomplete archive range")
                                    if path.name.startswith("guide_region_"):
                                        data = canonical_sequence(data)
                                    output.write(data)
                                    break
                                except (OSError, ValueError):
                                    if attempt == 3:
                                        raise
                                    time.sleep(2)
                    if (
                        hashlib.sha256(partial.read_bytes()).hexdigest()
                        != row["sha256"]
                    ):
                        raise ValueError(f"download checksum mismatch: {path}")
                    partial.replace(path)
                if (
                    path.stat().st_size != int(row["bytes"])
                    or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]
                ):
                    raise ValueError(f"input checksum mismatch: {path}")
    print("Verified all frozen broad-validation inputs and reference sequences.")


def broad_formats() -> None:
    """Decode XLSX and translate intact reference intervals; no statistical fitting."""
    root = Path("data/raw/broad_validation")
    target = Path("data/derived")
    target.mkdir(parents=True, exist_ok=True)
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

    def sheets(path: Path) -> list[list[list[str]]]:
        result = []
        with zipfile.ZipFile(path) as archive:
            strings = [
                "".join(node.itertext())
                for node in ET.fromstring(archive.read("xl/sharedStrings.xml"))
            ]
            for name in sorted(archive.namelist()):
                if not name.startswith("xl/worksheets/sheet") or not name.endswith(
                    ".xml"
                ):
                    continue
                rows = []
                for node in ET.fromstring(archive.read(name)).findall(".//s:row", ns):
                    values: dict[int, str] = {}
                    for cell in node.findall("s:c", ns):
                        column = 0
                        for letter in cell.attrib["r"]:
                            if letter.isalpha():
                                column = column * 26 + ord(letter) - ord("A") + 1
                        value = cell.find("s:v", ns)
                        if value is not None:
                            values[column - 1] = (
                                strings[int(value.text or "0")]
                                if cell.attrib.get("t") == "s"
                                else value.text or ""
                            )
                    rows.append(
                        [values.get(i, "") for i in range(max(values, default=-1) + 1)]
                    )
                result.append(rows)
        return result

    with (target / "broad_screen.tsv").open("w") as handle:
        output = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for i, rows in enumerate(sheets(root / "screen_S9.xlsx")):
            output.writerows(rows if i == 0 else rows[1:])

    # UCSC chain target is hg19; query is hg38. Map only entire intervals
    # contained in one aligned block, retaining strand and unique mapping.
    blocks: dict[str, list[tuple[int, int, str, int, int, str]]] = {}
    with gzip.open(root / "hg19ToHg38.over.chain.gz", "rt") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "chain":
                chrom, pos = parts[2], int(parts[5])
                qchrom, qsize, strand, qpos = (
                    parts[7],
                    int(parts[8]),
                    parts[9],
                    int(parts[10]),
                )
                assert parts[4] == "+"
            else:
                size = int(parts[0])
                blocks.setdefault(chrom, []).append(
                    (pos, pos + size, qchrom, qpos, qsize, strand)
                )
                pos += size + (int(parts[1]) if len(parts) == 3 else 0)
                qpos += size + (int(parts[2]) if len(parts) == 3 else 0)

    def reverse(sequence: str) -> str:
        return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    with (target / "broad_regions.tsv").open("w") as handle:
        output = csv.writer(handle, delimiter="\t", lineterminator="\n")
        output.writerow(
            [
                "gene",
                "chrom",
                "start",
                "end",
                "assembly",
                "hg19_guides",
                "hg38_guides",
                "mapping",
                "cpg_density",
            ]
        )
        for i, row in enumerate(sheets(root / "screen_guides_S8.xlsx")[0][1:]):
            chrom, span = row[5].split(":")
            start, end = map(int, span.split("-"))
            start -= 1
            sequences = [row[10].upper(), row[15].upper()]
            dna = {}
            hits = {}
            for assembly in ["hg19", "hg38"]:
                dna[assembly] = (
                    json.loads((root / f"guide_region_{i}_{assembly}.json").read_text())
                    .get("dna", "")
                    .upper()
                )
                hits[assembly] = sum(
                    bool(s) and (s in dna[assembly] or reverse(s) in dna[assembly])
                    for s in sequences
                )
            assembly = (
                "hg19"
                if hits["hg19"] == 2 and hits["hg38"] != 2
                else "hg38" if hits["hg38"] == 2 and hits["hg19"] != 2 else "ambiguous"
            )
            mappings = []
            if assembly == "hg19":
                for a, b, qc, qp, qs, strand in blocks.get(chrom, []):
                    if a <= start and end <= b:
                        x, y = qp + start - a, qp + end - a
                        mappings.append(
                            (qc, x, y) if strand == "+" else (qc, qs - y, qs - x)
                        )
            elif assembly == "hg38":
                mappings = [(chrom, start, end)]
            mapped = mappings[0] if len(mappings) == 1 else ("NA", 0, 0)
            source = dna.get(assembly, "")
            density = source.count("CG") / max(1, len(source))
            output.writerow(
                [
                    row[1],
                    *mapped,
                    assembly,
                    hits["hg19"],
                    hits["hg38"],
                    "unique" if len(mappings) == 1 else "unresolved",
                    density,
                ]
            )
    print("Decoded screen wells and sequence-checked, translated edited regions.")


if __name__ == "__main__":
    if sys.argv[1:] == ["broad"]:
        broad_formats()
    elif sys.argv[1:] == ["broad-download"]:
        broad_inputs(True)
    elif sys.argv[1:] == ["broad-verify"]:
        broad_inputs(False)
    else:
        main()
