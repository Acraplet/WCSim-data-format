import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class AnalysisPaths:
    """Central place to configure every filesystem path the validation suite needs.

    wcsim_build_dir must contain lib/libWCSimRoot.so matching the ROOT version
    on PATH - flatten_wcsim.C loads it with R__LOAD_LIBRARY.
    """

    wcsim_build_dir: str
    analysis_tools_dir: str = "/eos/user/a/acraplet/WCSim/Analysis/codes/analysis_tools"
    # None -> flattened files are written next to their raw input file
    flat_data_dir: Optional[str] = None
    plots_dir: str = "/eos/user/a/acraplet/WCSim/Analysis/codes/validations/plots"

    @classmethod
    def from_env(cls) -> "AnalysisPaths":
        """Build paths from the environment variables set by setup_env.sh.

        Only WCSIM_BUILD_DIR is required; the rest fall back to this class's
        defaults if unset (or unset/empty for flat_data_dir).
        """
        wcsim_build_dir = os.environ.get("WCSIM_BUILD_DIR")
        if not wcsim_build_dir:
            raise EnvironmentError(
                "WCSIM_BUILD_DIR is not set - source setup_env.sh before running this script."
            )

        kwargs = {"wcsim_build_dir": wcsim_build_dir}
        env_to_field = {
            "VALIDATION_ANALYSIS_TOOLS_DIR": "analysis_tools_dir",
            "VALIDATION_FLAT_DATA_DIR": "flat_data_dir",
            "VALIDATION_PLOTS_DIR": "plots_dir",
        }
        for env_var, field in env_to_field.items():
            value = os.environ.get(env_var)
            if value:
                kwargs[field] = value

        return cls(**kwargs)

    @property
    def flatten_script(self) -> Path:
        return Path(self.analysis_tools_dir) / "flatten_single_file.sh"

    def flat_path(self, raw_file) -> Path:
        """Default flattened-file path for a given raw file: '<stem>_flat.root'."""
        raw = Path(raw_file)
        out_dir = Path(self.flat_data_dir) if self.flat_data_dir else raw.parent
        return out_dir / f"{raw.stem}_flat.root"
