"""Acquire/verify frozen public dog data and decode the authors' R metadata.

All matching, effect estimation, resampling and conditional power are Rust code.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path

import rdata

MANIFEST = Path("code/results/dog_aging_sources.csv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("download", "verify", "decode"))
    args = parser.parse_args()
    with MANIFEST.open(newline="") as source:
        sources = list(csv.DictReader(source))
    for row in sources:
        path = Path(row["path"])
        if args.command == "download" and not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix + ".partial")
            subprocess.run(
                [
                    "curl",
                    "--fail",
                    "--location",
                    "--retry",
                    "2",
                    "--max-time",
                    "600",
                    row["url"],
                    "--output",
                    str(temporary),
                ],
                check=True,
            )
            if sha256(temporary) != row["sha256"]:
                raise ValueError(f"Changed upstream source: {path}")
            temporary.replace(path)
        if not path.exists() or sha256(path) != row["sha256"]:
            raise ValueError(f"Missing or changed source: {path}")
    print(f"Verified {len(sources)} frozen sources")
    if args.command == "decode":
        metadata = rdata.read_rds("data/raw/dog_aging/dap_rrbs-metaData.rds")
        expected = (
            "lid_pid",
            "dog_id",
            "Cohort",
            "first_rrbs",
            "sid",
            "prep_date",
            "Sex",
            "Age_at_sample",
            "fixed",
            "Breed_Size_Class_at_HLES",
            "predicted_height",
            "Heterozyg",
        )
        if tuple(metadata.columns) != expected:
            raise ValueError("Unexpected author metadata schema")
        if metadata.lid_pid.duplicated().any():
            raise ValueError("Duplicate author sample identifiers")
        destination = Path("data/derived/dog_aging/author_metadata.csv")
        destination.parent.mkdir(parents=True, exist_ok=True)
        metadata.to_csv(destination, index=False)
        print(f"Decoded {len(metadata)} metadata rows into {destination}")


if __name__ == "__main__":
    main()
