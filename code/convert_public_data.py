"""Format conversion only; all model fitting and statistical tests run in Rust.

Run with pinned dependencies in data/tools/format-env. Preserve missing methylation,
cell IDs and independent LARRY labels; never infer clones from CpGs in this converter.
"""

from __future__ import annotations

import csv
import tarfile
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


if __name__ == "__main__":
    main()
