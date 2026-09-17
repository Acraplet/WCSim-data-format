#!/usr/bin/env python3
"""Caller script for the validation module.

Loads one or more flattened files into DataFrames (module.validation_module
does the actual file I/O) and then runs whatever checks/plots we want on
them. Files can be given either directly on the command line, or listed in
a text file with --files-list (format TBD - written by a separate script).

TODO: fill in the actual validation logic once files are loading correctly.
"""

import argparse
from pathlib import Path

from module import load_runs, plot_histogram_by_particle_and_energy


def parse_files_list(path) -> list:
    """Read a text file of files to load into (name, file_path) pairs.

    TODO: decide on the exact file format (one path per line? "name path"
    per line? comments with '#'?) once the file-writing script exists -
    see runs_file.py for a similar convention already used elsewhere.
    """
    path = Path(path)
    runs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # placeholder: treat the whole line as a file path, name = stem
            file_path = line
            name = Path(file_path).stem
            runs.append((name, file_path))
    return runs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "files", nargs="*",
        help="Flattened file(s) to load, given directly on the command line.",
    )
    parser.add_argument(
        "--files-list", type=str, default=None,
        help="Text file listing files to load (one per line), as an alternative to passing them directly.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.files_list:
        runs = parse_files_list(args.files_list)
    else:
        runs = [(Path(f).stem, f) for f in args.files]

    if not runs:
        raise SystemExit("No files given - pass files as arguments or use --files-list.")

    df = load_runs(runs)

    for name, group in df.groupby("run_name"):
        print(f"{name}: loaded {len(group)} events")
    print(f"total: {len(df)} events from {df['run_name'].nunique()} run(s)")

    plot_histogram_by_particle_and_energy(df, "true_max_scatter_angle_deg")

    # TODO: actual validation logic (comparisons, plots, cuts, ...) goes here.


if __name__ == "__main__":
    main()
