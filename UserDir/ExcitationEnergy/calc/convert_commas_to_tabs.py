#!/usr/bin/env python3
"""Convert comma-delimited files in a folder to tab-delimited.

Defaults to writing new files alongside inputs with a .tsv suffix.
Use --inplace to overwrite originals.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert comma-delimited files in a folder to tab-delimited files."
    )
    parser.add_argument(
        "folder",
        type=Path,
        help="Folder containing comma-delimited files",
    )
    parser.add_argument(
        "--glob",
        default="*",
        help="Glob pattern for files to include (default: '*')",
    )
    parser.add_argument(
        "--ext",
        default=None,
        help="Only process files with this extension (e.g., .dat).",
    )
    parser.add_argument(
        "--suffix",
        default=".tsv",
        help="Suffix to append to output filenames when not using --inplace (default: .tsv)",
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Overwrite input files in place",
    )
    return parser.parse_args()


def output_path_for(input_path: Path, suffix: str) -> Path:
    if suffix and not suffix.startswith("."):
        suffix = f".{suffix}"
    if input_path.suffix:
        return input_path.with_suffix(input_path.suffix + suffix)
    return input_path.with_name(input_path.name + suffix)


def convert_file(input_path: Path, output_path: Path) -> bool:
    # Return True if conversion occurred, False if skipped.
    try:
        with input_path.open("r", newline="", encoding="utf-8") as src:
            sample = src.read(4096)
            src.seek(0)
            if "," not in sample:
                return False

            reader = csv.reader(src, delimiter=",")
            with output_path.open("w", newline="", encoding="utf-8") as dst:
                writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
                for row in reader:
                    writer.writerow(row)
        return True
    except UnicodeDecodeError:
        # Skip non-text/binary files.
        return False


def main() -> int:
    args = parse_args()
    folder = args.folder

    if not folder.exists() or not folder.is_dir():
        print(f"Folder not found or not a directory: {folder}", file=sys.stderr)
        return 2

    if args.ext and not args.ext.startswith("."):
        args.ext = f".{args.ext}"

    converted = 0
    skipped = 0

    for path in sorted(folder.glob(args.glob)):
        if not path.is_file():
            continue
        if args.ext and path.suffix != args.ext:
            continue

        if args.inplace:
            out_path = path
        else:
            out_path = output_path_for(path, args.suffix)
            if out_path == path or out_path.exists():
                # Avoid clobbering an existing output.
                skipped += 1
                continue

        did_convert = convert_file(path, out_path)
        if did_convert:
            converted += 1
        else:
            skipped += 1

    print(f"Converted: {converted}")
    print(f"Skipped: {skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
