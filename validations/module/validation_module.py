"""Core module for the validation scripts.

Holds the logic to open a single flattened file, pull out the branches we
care about, and append each event's data into a pandas DataFrame. The caller
script(s) in validations/ are responsible for argument parsing and looping
over the list of files; this module only knows how to deal with one file at
a time (plus a helper to combine several already-loaded DataFrames).

TODO: fill in the branch list and the per-event extraction logic together.
"""

import copy
import os
import sys
from pathlib import Path
from typing import Iterable, List, Optional

import awkward as ak
import matplotlib.colors as mcolors
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

    Adds an `event_index` column: each row's position within THIS file
    (0-based). Unlike `event_number` (a truth branch, not guaranteed unique
    across different runs/files once several get concatenated), event_index
    is assigned here. Note this is only unique within a single file - once
    several files are combined via load_runs(), it gets overwritten there
    with the row's position in the combined DataFrame instead (see
    load_runs docstring).

    TODO: some of the jagged branches may need summarizing (e.g. n_hits,
    max/mean) rather than being kept as raw lists - revisit once we know
    what the plots need.
    """
    file_path = Path(file_path)
    branches = branches or DEFAULT_BRANCHES

    with uproot.open(file_path) as f:
        data = f[tree_name].arrays(branches, library="ak")

    rows = []
    for event_index, event in enumerate(data):
        row = {"event_index": event_index}
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

    `event_index` is reset here to the row's position in the COMBINED
    DataFrame (0-based), overwriting the per-file value load_file assigned.
    Otherwise every file would restart at 0 and event_index would collide
    across files once concatenated, defeating its purpose as a unique row id.
    """
    dfs = [load_run(name, file_path, branches=branches, tree_name=tree_name) for name, file_path in runs]
    df = pd.concat(dfs, ignore_index=True)
    df["event_index"] = df.index
    return df


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


# --- 2D event display ---------------------------------------------------------
# Draws hits projected onto the mPMT plane via the (separately maintained)
# EventDisplay class - see examples/compare_sim_data.ipynb section 7, which this
# is adapted from. Loading the mPMT position CSV is slow, so the EventDisplay
# instance is built once per event_display_dir and cached.

_EVENT_DISPLAY_CACHE = {}


def _get_event_display(event_display_dir):
    """Lazily construct (and cache) the EventDisplay for event_display_dir."""
    event_display_dir = str(event_display_dir)
    if event_display_dir in _EVENT_DISPLAY_CACHE:
        return _EVENT_DISPLAY_CACHE[event_display_dir]

    sys.path.insert(0, event_display_dir)
    from EventDisplay import EventDisplay  # noqa: E402

    cwd = os.getcwd()
    os.chdir(event_display_dir)  # load_mPMT_positions resolves its CSV relative to the module dir
    try:
        ed = EventDisplay()
        ed.load_mPMT_positions("mPMT_2D_projection_angles.csv")
    finally:
        os.chdir(cwd)

    _EVENT_DISPLAY_CACHE[event_display_dir] = ed
    return ed


def _channel_values(slots, positions, values, n_channels, how="sum"):
    """Aggregate per-hit values onto per-(slot, position) channel values.

    how="sum": for continuous values (e.g. charge) - unhit channels come back
    as NaN (not 0), so they're left blank rather than drawn as the bottom of
    the colour scale.
    how="first": for categorical values (e.g. a track id) where summing or
    averaging would be meaningless - picks the first hit seen per channel.
    """
    out = np.full(n_channels, np.nan)
    hit = np.zeros(n_channels, dtype=bool)
    for slot, pos, v in zip(slots, positions, values):
        ch = 19 * int(slot) + int(pos)
        if how == "sum":
            out[ch] = v if not hit[ch] else out[ch] + v
        elif how == "first":
            if not hit[ch]:
                out[ch] = v
        else:
            raise ValueError(f"how must be 'sum' or 'first', got {how!r}")
        hit[ch] = True
    return out


def _track_color_labels(row, track_ids_in_event):
    """One legend string per track id, in the same order as its colour index.

    Looks up each track id in row's per-track branches (track_id, track_pdg,
    track_process) to print what that colour actually corresponds to, instead
    of a bare, meaningless integer on the colorbar.
    """
    track_id_list, pdg_list, process_list = row["track_id"], row["track_pdg"], row["track_process"]

    labels = []
    for tid in track_ids_in_event:
        matches = [j for j, t in enumerate(track_id_list) if t == tid]
        if matches:
            j = matches[0]
            labels.append(f"track_id={tid}  pdg={pdg_list[j]}  {process_list[j]}")
        elif tid == -1:
            labels.append("track_id=-1  (dark noise / no truth match)")
        else:
            labels.append(f"track_id={tid}  (no truth match)")
    return labels


def _draw_event(ax, ed, channel_values, cmap, norm, title, color_label=None,
                 colorbar_ticks=None, colorbar_ticklabels=None):
    """Draw one event's mPMT-plane hit pattern into `ax`.

    colorbar_ticks/colorbar_ticklabels let a discrete/categorical colour scale
    (e.g. track id) spell out what each colour means right next to the bar,
    instead of a bare integer scale.
    """
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Circle

    coordinates = ed.coordinates_eachChannel()
    cmap = copy.copy(cmap)
    cmap.set_bad(color="white")

    pmt_circles = [Circle((x, y), radius=0.48) for x, y in ed.mPMT_2D_projection[:, 1:3]]
    ax.add_collection(PatchCollection(pmt_circles, facecolor="none", linewidths=1, edgecolors="0.85"))
    pts = ax.scatter(coordinates[:, 0], coordinates[:, 1], c=channel_values, s=25, cmap=cmap, norm=norm)

    ax.set_title(title)
    ax.set_aspect("equal")
    ax.axis("off")

    cbar = plt.colorbar(pts, ax=ax, pad=0.01, label=color_label)
    if colorbar_ticks is not None:
        cbar.set_ticks(colorbar_ticks)
        cbar.set_ticklabels(colorbar_ticklabels)
    return pts


def _event_summary_title(row: pd.Series) -> str:
    """Two-line MC-truth summary for the primary particle in one event.

    Line 1: what the primary is and how it stopped. Line 2: its largest
    recorded single-step scatter (angle, energy transferred, process, and
    where it happened) - or a note that no scatter was recorded, when
    true_max_scatter_angle_deg == -1 (see DEFAULT_BRANCHES).
    """
    line1 = (
        f"primary pdg={row['true_pdg']}  true_E={row['true_E']:.1f} MeV  "
        f"stop_process={row['stop_process'] or '(ranged out / ionization only)'}  "
        f"stopvol={row['true_stopvol']}"
    )

    angle = row["true_max_scatter_angle_deg"]
    if angle == -1:
        line2 = "max single-step scatter: none recorded"
    else:
        delta_ke = row["true_max_scatter_ke_pre"] - row["true_max_scatter_ke_post"]
        line2 = (
            f"max single-step scatter: {angle:.1f} deg, {delta_ke:.2f} MeV transferred, "
            f"process={row['true_max_scatter_process']}, "
            f"at (x,y,z)=({row['true_max_scatter_x']:.1f}, {row['true_max_scatter_y']:.1f}, {row['true_max_scatter_z']:.1f}) cm"
        )

    return f"{line1}\n{line2}"


def plot_event_display(
    row: pd.Series,
    event_display_dir,
    out_dir=None,
    out_name: Optional[str] = None,
) -> Path:
    """Side-by-side 2D event display for one event: charge (left) vs.
    hit-to-track assignment (right).

    `row` is one row of a DataFrame from load_file/load_run/load_runs - needs
    the per-hit branches (hit_mpmt_slot_ids, hit_pmt_position_ids,
    hit_pmt_charges, hit_track_id) and per-track truth (track_id, track_pdg,
    track_process) to label the right-hand colorbar.

    event_display_dir is the directory containing EventDisplay.py and
    mPMT_2D_projection_angles.csv (see run_loader.py's EVENT_DISPLAY_DIR).
    """
    ed = _get_event_display(event_display_dir)

    slots, positions = row["hit_mpmt_slot_ids"], row["hit_pmt_position_ids"]
    channel_charge = _channel_values(slots, positions, row["hit_pmt_charges"], ed.nChannels, how="sum")

    hit_track_ids = row["hit_track_id"]
    track_ids_in_event = sorted(set(hit_track_ids))
    track_id_to_index = {tid: i for i, tid in enumerate(track_ids_in_event)}
    hit_color_index = [track_id_to_index[tid] for tid in hit_track_ids]
    channel_track = _channel_values(slots, positions, hit_color_index, ed.nChannels, how="first")

    n_tracks = len(track_ids_in_event)
    track_labels = _track_color_labels(row, track_ids_in_event)

    event_number = int(row["event_number"])
    event_index = int(row["event_index"])
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    _draw_event(
        axes[0], ed, channel_charge, cmap=plt.cm.plasma,
        norm=mcolors.Normalize(vmin=0, vmax=np.nanmax(channel_charge) if np.any(~np.isnan(channel_charge)) else 1),
        title=f"event {event_number}: charge",
        color_label="charge [p.e.]",
    )
    # discrete colour scale (crisp per-track blocks, as in run_loader.py's tab20 display),
    # but sampled from turbo (not a fixed-size qualitative map) so it scales to any
    # number of tracks without two of them ever landing on the same colour
    track_cmap = plt.cm.get_cmap("turbo", max(n_tracks, 1))
    track_norm = mcolors.BoundaryNorm(np.arange(-0.5, n_tracks + 0.5, 1), track_cmap.N)
    _draw_event(
        axes[1], ed, channel_track, cmap=track_cmap, norm=track_norm,
        title=f"event {event_number}: Track information",
        colorbar_ticks=list(range(n_tracks)),
        colorbar_ticklabels=track_labels,
    )

    fig.suptitle(_event_summary_title(row), fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.88])

    out_dir = Path(out_dir) if out_dir else Path(__file__).parent.parent / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / (out_name or f"event_display_evt{event_index}.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved plot to {out_path}")

    return out_path
