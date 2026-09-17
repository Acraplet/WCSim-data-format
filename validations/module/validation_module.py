"""Core module for the validation scripts.

Holds the logic to open a single flattened file, pull out the branches we
care about, and append each event's data into a pandas DataFrame. The caller
script(s) in validations/ are responsible for argument parsing and looping
over the list of files; this module only knows how to deal with one file at
a time (plus a helper to combine several already-loaded DataFrames).

TODO: fill in the branch list and the per-event extraction logic together.
"""

from pathlib import Path
from typing import Iterable, List, Optional

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot


# TODO: branches to read from the flattened ROOT file (see run_loader.py's
# sim_branches for the full list of what's available).
DEFAULT_BRANCHES: List[str] = [
     # per-hit
    "hit_pmt_charges", "hit_pmt_calibrated_times", "hit_mpmt_slot_ids", "hit_pmt_position_ids",
    "hit_track_id", "hit_x", "hit_y", "hit_z",
    # per-event truth
    "n_digihits", "event_number", "true_pdg", "true_E", "true_ke",
    "true_vtx_x", "true_vtx_y", "true_vtx_z", "true_start_x", "true_start_y", "true_start_z",
    "true_stop_x", "true_stop_y", "true_stop_z",
    "stop_process", "had_inelastic", "had_elastic", "n_prim_daughters",
    "true_exit_ke", "true_stopvol", "n_elastic", "n_inelastic",
    # per-track (jagged, one entry per saved track in the event)
    "track_id", "track_parent_id", "track_pdg", "track_process", "track_ke",
    "track_start_x", "track_start_y", "track_start_z",
    "track_end_x", "track_end_y", "track_end_z", "track_nhits",
    #largest single-step deflection along the primary [deg] 
    "true_max_scatter_angle_deg", # largest single-step deflection along the primary [deg] (-1 = none)
    "true_max_scatter_process",   #process that defined that step
    "true_max_scatter_x",# // position of that step [cm]
    "true_max_scatter_y",
    "true_max_scatter_z", 
    "true_max_scatter_predir_x", # track direction just before that step (unit vector)
    "true_max_scatter_predir_y", 
    "true_max_scatter_predir_z", 
    "true_max_scatter_postdir_x", # track direction just after that step (unit vector)
    "true_max_scatter_postdir_y", 
    "true_max_scatter_postdir_z", 
    "true_max_scatter_ke_pre", # KE just before that step [MeV]
    "true_max_scatter_ke_post",# KE just after that step [MeV]
    "true_n_scatters_above_7deg",
]


def load_file(file_path, branches: Optional[List[str]] = None, tree_name: str = "hits") -> pd.DataFrame:
    """Open one flattened file and build a DataFrame with one row per event.

    Per-event scalar branches (e.g. true_pdg, true_E) become plain columns.
    Per-hit / per-track branches (jagged - multiple entries per event) are
    kept as list-valued columns, one Python list per event.

    TODO: some of the jagged branches may need summarizing (e.g. n_hits,
    max/mean) rather than being kept as raw lists - revisit once we know
    what the plots need.
    """
    file_path = Path(file_path)
    branches = branches or DEFAULT_BRANCHES

    with uproot.open(file_path) as f:
        data = f[tree_name].arrays(branches, library="ak")

    rows = []
    for event in data:
        row = {}
        for branch in branches:
            value = event[branch]
            row[branch] = ak.to_list(value) if isinstance(value, ak.Array) else value
        rows.append(row)

    return pd.DataFrame(rows)


def load_run(name: str, file_path, branches: Optional[List[str]] = None, tree_name: str = "hits") -> pd.DataFrame:
    """load_file, tagged with which run/file each row came from.

    Adds `run_name` and `file_path` columns so identity survives once this
    gets concatenated with other runs in load_runs().
    """
    df = load_file(file_path, branches=branches, tree_name=tree_name)
    df.insert(0, "run_name", name)
    df.insert(1, "file_path", str(file_path))
    return df


def load_runs(runs: Iterable, branches: Optional[List[str]] = None, tree_name: str = "hits") -> pd.DataFrame:
    """Load several (name, file_path) pairs and concatenate into one DataFrame.

    `runs` is an iterable of (name, file_path) tuples - see runs_file.py for
    where these come from when read from a text file. Rows carry `run_name`/
    `file_path` columns (see load_run) so you can still filter/group back
    down to a single file after everything's been combined.
    """
    dfs = [load_run(name, file_path, branches=branches, tree_name=tree_name) for name, file_path in runs]
    return pd.concat(dfs, ignore_index=True)


# Fixed-order categorical palette (CVD-safe adjacent ordering) - cycle only if
# there are more (pdg, energy) combinations than slots.
_CATEGORICAL_COLORS = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#4a3aa7",  # violet
    "#e34948",  # red
]


def plot_histogram_by_particle_and_energy(
    df: pd.DataFrame,
    variable: str,
    pdg_col: str = "true_pdg",
    energy_col: str = "true_E",
    bins=50,
    out_dir=None,
    out_name: Optional[str] = None,
) -> Path:
    """Histogram `variable`, overlaying one series per (pdg_col, energy_col)
    combination in df on a single plot.

    Groups the combined DataFrame by true PID and starting energy (each
    (pdg, energy) pair is assumed to be one simulated beam configuration).
    Saves one PDF containing all series and returns its path.
    """
    out_dir = Path(out_dir) if out_dir else Path(__file__).parent.parent / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    groups = list(df.groupby([pdg_col, energy_col]))

    # shared bin edges so every series is comparable on the same axis
    bin_edges = np.histogram_bin_edges(df[variable].dropna(), bins=bins)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, ((pdg, energy), group) in enumerate(groups):
        color = _CATEGORICAL_COLORS[i % len(_CATEGORICAL_COLORS)]
        values = group[variable].dropna()
        ax.hist(
            values, bins=bin_edges, histtype="step", linewidth=2, color=color,
            label=f"true_pdg={pdg}, {energy_col}={energy:.1f}",
        )

    ax.set_xlabel(variable)
    ax.set_ylabel("events")
    ax.set_title(variable)
    ax.legend(frameon=False)
    fig.tight_layout()

    out_path = out_dir / (out_name or f"{variable}_by_pdg_and_{energy_col}.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved plot to {out_path}")

    return out_path
