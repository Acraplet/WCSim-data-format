from pathlib import Path
from typing import List

from .flattening import RunSpec


def parse_runs_file(path) -> List[RunSpec]:
    """Parse a runs list file into RunSpecs.

    One run per line, whitespace-separated:

        <name>  <raw_file>  [is_mdt]

    - raw_file may be absolute, or relative to the current working directory
      (resolved later by Flattener, not here).
    - is_mdt is optional: -1 (default, auto-detect), 0 (force plain WCSim),
      or 1 (force MDT) - see flatten_wcsim.C.
    - Blank lines and lines starting with '#' are ignored.

    Example:
        # name              raw_file                                          is_mdt
        pi-_293MeV          wcsim_wCDS_pi-_Beam_293MeV_5cm_0000.root
        pi-_293MeV_MDT      wcsim_wCDS_pi-_Beam_293MeV_5cm_0000_MDT.root       1
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"runs file not found: {path}")

    runs = []
    names_seen = set()
    with open(path) as f:
        for lineno, raw_line in enumerate(f, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue

            fields = line.split()
            if len(fields) not in (2, 3):
                raise ValueError(
                    f"{path}:{lineno}: expected '<name> <raw_file> [is_mdt]', got: {raw_line!r}"
                )

            name, raw_file = fields[0], fields[1]
            is_mdt = int(fields[2]) if len(fields) == 3 else -1

            if name in names_seen:
                raise ValueError(f"{path}:{lineno}: duplicate run name {name!r}")
            names_seen.add(name)

            runs.append(RunSpec(name=name, raw_file=raw_file, is_mdt=is_mdt))

    return runs
