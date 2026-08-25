#!/usr/bin/env python3
"""Strip cell outputs/execution counts from a Jupyter notebook.

Used as a git 'clean' filter (see .gitattributes) so notebooks are stripped
of outputs when staged/committed, while the working-tree copy on disk keeps
its outputs untouched for local viewing.

Usage:
    python3 strip_notebook_outputs.py < in.ipynb > out.ipynb
    python3 strip_notebook_outputs.py notebook.ipynb   # strips in place
"""
import json
import sys


def strip(nb):
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        cell["outputs"] = []
        cell["execution_count"] = None
        cell.get("metadata", {}).pop("execution", None)
    return nb


def main():
    if len(sys.argv) > 1:
        path = sys.argv[1]
        with open(path) as f:
            nb = json.load(f)
        nb = strip(nb)
        with open(path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")
    else:
        nb = json.load(sys.stdin)
        nb = strip(nb)
        json.dump(nb, sys.stdout, indent=1)
        sys.stdout.write("\n")


if __name__ == "__main__":
    main()
