import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from .paths import AnalysisPaths


@dataclass
class RunSpec:
    """One WCSim run to validate.

    is_mdt: -1 (default) auto-detect, 0 force plain-WCSim truth-matching,
    1 force MDT truth-matching - passed straight through to flatten_wcsim.C.
    flat_file: override the default '<raw stem>_flat.root' output path.
    """

    name: str
    raw_file: str
    is_mdt: int = -1
    flat_file: Optional[str] = None


class Flattener:
    """Runs flatten_single_file.sh (-> flatten_wcsim.C) for RunSpecs, skipping
    any run whose flattened output already exists.
    """

    def __init__(self, paths: AnalysisPaths):
        self.paths = paths

    def resolve_raw(self, run: RunSpec) -> Path:
        return self.paths.raw_path(run.raw_file)

    def resolve_flat(self, run: RunSpec) -> Path:
        if run.flat_file is not None:
            p = Path(run.flat_file)
            if p.is_absolute():
                return p
            base = Path(self.paths.flat_data_dir) if self.paths.flat_data_dir else self.resolve_raw(run).parent
            return base / p
        return self.paths.flat_path(self.resolve_raw(run))

    def ensure_flattened(self, run: RunSpec, force: bool = False) -> Path:
        """Flatten `run` if its flattened output doesn't exist yet (or force=True).
        Returns the path to the flattened file either way.
        """
        raw = self.resolve_raw(run)
        flat = self.resolve_flat(run)

        if not raw.exists():
            raise FileNotFoundError(f"[{run.name}] raw WCSim file not found: {raw}")

        if flat.exists() and not force:
            print(f"[{run.name}] flattened file already exists, skipping flatten: {flat}")
            return flat

        script = self.paths.flatten_script
        if not script.exists():
            raise FileNotFoundError(f"flatten script not found: {script}")

        flat.parent.mkdir(parents=True, exist_ok=True)

        env = os.environ.copy()
        env["WCSIM_BUILD_DIR"] = self.paths.wcsim_build_dir

        cmd = [str(script), str(raw), str(flat), str(run.is_mdt)]
        print(f"[{run.name}] flattening: {' '.join(cmd)}")
        result = subprocess.run(cmd, env=env)
        if result.returncode != 0:
            raise RuntimeError(f"[{run.name}] flattening failed (exit code {result.returncode})")

        return flat

    def ensure_all(self, runs: List[RunSpec], force: bool = False) -> Dict[str, Path]:
        """Flatten every run that needs it, returning {run name: flat file path}."""
        return {run.name: self.ensure_flattened(run, force=force) for run in runs}
