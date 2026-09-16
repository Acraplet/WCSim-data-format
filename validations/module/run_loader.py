from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import awkward as ak
import uproot

# Same branch list used in examples/compare_sim_data.ipynb section 1, so a
# RunData.events plays nicely with code lifted from that notebook.
DEFAULT_SIM_BRANCHES = [
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
]


@dataclass
class RunData:
    name: str
    flat_file: Path
    events: ak.Array

    def __len__(self) -> int:
        return len(self.events)


class RunLoader:
    """Reads flatten_wcsim.C output ("hits" TTree) into awkward arrays, one
    RunData per run - the same load step as compare_sim_data.ipynb section 1.
    """

    def __init__(self, branches: Optional[List[str]] = None, tree_name: str = "hits"):
        self.branches = branches or DEFAULT_SIM_BRANCHES
        self.tree_name = tree_name

    def load(self, name: str, flat_file) -> RunData:
        flat_file = Path(flat_file)
        with uproot.open(flat_file) as f:
            events = f[self.tree_name].arrays(self.branches, library="ak")
        print(f"[{name}] loaded {len(events)} events from {flat_file}")
        return RunData(name=name, flat_file=flat_file, events=events)

    def load_all(self, flat_files: Dict[str, Path]) -> Dict[str, RunData]:
        return {name: self.load(name, path) for name, path in flat_files.items()}


def classify_outcome(stop_process: str, had_inelastic: int, had_elastic: int, true_exit_ke: float) -> str:
    """Human-readable summary of what happened to the primary, from MC truth.

    Same classification as examples/compare_sim_data.ipynb section 5.
    """
    if "Decay" in stop_process:
        return "decayed in flight"
    if "CaptureAtRest" in stop_process:
        return "captured at rest"
    if had_inelastic:
        return "hadronic inelastic interaction"
    if had_elastic:
        return "hadronic elastic scatter"
    if true_exit_ke >= 0:
        return "exited the tank"
    return "ranged out (ionization only)"
