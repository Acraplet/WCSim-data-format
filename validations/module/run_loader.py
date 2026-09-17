from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np

import sys

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import awkward as ak
import uproot

SIM_FILE = "/eos/user/m/mprincov/WCSim/WCSim-data-format/analysis_tools/flat_wCDS_mu-_Beam_350MeV_0cm_0000.root"

EVENT_DISPLAY_DIR = "/eos/user/m/mprincov/WCTE_event_display"

PARTICLE = "muon"  
# --- standard per-species definitions, keyed by PARTICLE. Add more rows here
#     (proton, kaon, electron, ...) as needed. ---
# PARTICLE_PDG: |true_pdg| of the MC primary for this species.
PARTICLE_PDG = {"muon": 13, "pion": 211}

# DECAY_DAUGHTER_PDG: |pdg| of the expected decay-in-flight daughter, used in
# section 6's decay-in-flight kinematics (pi- -> mu- nu_mu_bar, mu- -> e- nu nubar).
DECAY_DAUGHTER_PDG = {"pion": 13, "muon": 11}

print(f"PARTICLE = {PARTICLE!r}  (MC selects |true_pdg| == {PARTICLE_PDG[PARTICLE]})")


sim_branches = [
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

with uproot.open(SIM_FILE) as f:
    sim = f["hits"].arrays(sim_branches, library="ak")

print(f"Loaded {len(sim)} simulated events from {SIM_FILE}")

# --- all interaction modes in the simulated sample ---------------------------------
# flatten_wcsim.C records the primary's fate in `stop_process` (the discrete process
# at its stop point: "" = ranged out / ionization only) and every Geant4 creator
# process of every saved track in `track_process`. Print both.
from collections import Counter

n_ev = len(sim)

# 1. primary-particle stop process (one per event)
stop_counts = Counter(ak.to_list(sim["stop_process"]))
print(f"Primary stop_process over {n_ev} events:")
for proc, n in stop_counts.most_common():
    label = proc if proc else "(ranged out / ionization only)"
    print(f"  {label:<28s} {n:6d}  ({100*n/n_ev:5.1f}%)")

# 2. hadronic-interaction flags (elastic / inelastic anywhere along the primary lineage)
print("\nPrimary hadronic-interaction summary:")
print(f"  had_elastic == 1   : {int(ak.sum(sim['had_elastic'] == 1)):6d}")
print(f"  had_inelastic == 1 : {int(ak.sum(sim['had_inelastic'] == 1)):6d}")
print(f"  n_elastic distribution   : {dict(Counter(ak.to_list(sim['n_elastic'])))}")
print(f"  n_inelastic distribution : {dict(Counter(ak.to_list(sim['n_inelastic'])))}")

# 3. every creator process seen among all saved tracks in the sample
proc_counts = Counter(
    proc for ev in ak.to_list(sim["track_process"]) for proc in ev
)
print(f"\nAll track creator processes (across every saved track, {sum(proc_counts.values())} tracks):")
for proc, n in proc_counts.most_common():
    print(f"  {proc:<28s} {n:8d}")


species_pdg = PARTICLE_PDG[PARTICLE]
is_species  = np.abs(ak.to_numpy(sim["true_pdg"])) == species_pdg
has_decay   = np.array(["Decay" in s for s in ak.to_list(sim["stop_process"])])

dif_mask   = is_species & has_decay      # decayed in flight
nodif_mask = is_species & ~has_decay     # did NOT decay in flight

print("Number of decay-in-flight events:", np.flatnonzero(dif_mask))

sim_species = sim[is_species]
sim_dif     = sim[dif_mask]
sim_nodif   = sim[nodif_mask]

print(f"MC: {is_species.sum()} {PARTICLE} events "
      f"({dif_mask.sum()} decay-in-flight, {nodif_mask.sum()} no decay-in-flight)")
print("\nstop_process breakdown of the NO-decay-in-flight sample:")
for proc, n in Counter(np.array(ak.to_list(sim["stop_process"]))[nodif_mask].tolist()).most_common():
    print(f"  {proc or '(ranged out / ionization only)':<28s} {n:6d}")


#-----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(ak.to_numpy(sim["true_max_scatter_angle_deg"]), bins=50)
ax.set_xlabel("true_max_scatter_angle_deg")
ax.set_ylabel("events")
ax.set_title(f"Largest single-step scatter angle of the primary ({PARTICLE})")
fig.tight_layout()

out_pdf = Path(__file__).parent / "plots/true_max_scatter_angle_deg.pdf"
fig.savefig(out_pdf)
print(f"Saved plot to {out_pdf}")

#-----------------------------------------------------------------
# stop_process breakdown, restricted to events with a recorded single-step scatter
# (true_max_scatter_angle_deg != -1, i.e. a scatter actually happened)
has_scatter = ak.to_numpy(sim["true_max_scatter_angle_deg"]) != -1
scattered_stop_process = [
    sp if sp else "(ranged out / ionization only)"
    for sp, keep in zip(ak.to_list(sim["stop_process"]), has_scatter) if keep
]
stop_process_counts = Counter(scattered_stop_process).most_common()
labels, counts = zip(*stop_process_counts)

fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(labels, counts)
ax.set_xlabel("stop_process")
ax.set_ylabel("events")
ax.set_title(f"stop_process for events with a recorded single-step scatter ({PARTICLE})")
ax.tick_params(axis="x", rotation=45)
fig.tight_layout()

out_pdf = Path(__file__).parent / "plots/stop_process_scattered.pdf"
fig.savefig(out_pdf)
print(f"Saved plot to {out_pdf}")

#-----------------------------------------------------------------


sys.path.insert(0, EVENT_DISPLAY_DIR)
from EventDisplay import EventDisplay  # noqa: E402
import os
os.chdir(EVENT_DISPLAY_DIR)  # load_mPMT_positions resolves its CSV relative to the module dir

ed = EventDisplay()
ed.load_mPMT_positions("mPMT_2D_projection_angles.csv")


def channel_first_value(slots, positions, values, n_channels):
    """One representative value per (slot, position) channel - the first hit seen.

    Unlike EventDisplay.process_data (sum, or min with a 0-vs-unset ambiguity), this
    is what we want for a *categorical* value like a track id: pick a single
    representative hit per channel rather than summing/averaging track ids.
    """
    out = np.full(n_channels, np.nan)
    for slot, pos, v in zip(slots, positions, values):
        ch = 19 * int(slot) + int(pos)
        if np.isnan(out[ch]):
            out[ch] = v
    return out

def classify_mpmt_slots(good_encoded, n_mpmts=106, n_pos=19):
    """Classify each mPMT slot by how many of its 19 positions are in the
    good_wcte_pmts whitelist (slot*100+position encoding).

    Returns (fully_off_slots, partially_off_slots): slots with 0/19 good
    channels, and slots with some-but-not-all good channels, respectively.
    Fully-good slots (19/19) are not returned - nothing to highlight there.
    """
    good_set = set(int(v) for v in good_encoded)
    fully_off, partially_off = [], []
    for slot in range(n_mpmts):
        n_good = sum((slot * 100 + pos) in good_set for pos in range(n_pos))
        if n_good == 0:
            fully_off.append(slot)
        elif n_good < n_pos:
            partially_off.append(slot)
    return np.array(fully_off), np.array(partially_off)

def draw_event(ax, ed, data, cmap, norm, title, color_label=None,
                fully_off_slots=(), partially_off_slots=()):
    """Minimal reimplementation of EventDisplay.plotEventDisplay that draws into a
    given ax instead of creating its own figure, so two events can share one row.

    fully_off_slots / partially_off_slots (mPMT slot IDs, see classify_mpmt_slots)
    are drawn as extra highlighted circles on top of the normal per-mPMT outline,
    so it's visually obvious which blank regions are due to the good_wcte_pmts
    mask rather than just this one event happening to have no hits there.
    """
    import copy
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection

    coordinates = ed.coordinates_eachChannel()
    plot_data = data.copy()
    #plot_data[plot_data == 0] = np.nan  # keep unhit channels blank, as plotEventDisplay does

    cmap = copy.copy(cmap)
    cmap.set_bad(color="white")

    pmt_circles = [Circle((x, y), radius=0.48) for x, y in ed.mPMT_2D_projection[:, 1:3]]
    ax.add_collection(PatchCollection(pmt_circles, facecolor="none", linewidths=1, edgecolors="0.85"))
    pts = ax.scatter(coordinates[:, 0], coordinates[:, 1], c=plot_data, s=25, cmap=cmap, norm=norm)

    # highlight mPMTs excluded (fully or partially) by the good_wcte_pmts mask
    if len(fully_off_slots):
        xy = ed.mPMT_2D_projection[np.asarray(fully_off_slots), 1:3]
        ax.add_collection(PatchCollection(
            [Circle((x, y), radius=0.75) for x, y in xy],
            facecolor="none", edgecolors="red", linewidths=2, zorder=5))
    if len(partially_off_slots):
        xy = ed.mPMT_2D_projection[np.asarray(partially_off_slots), 1:3]
        ax.add_collection(PatchCollection(
            [Circle((x, y), radius=0.75) for x, y in xy],
            facecolor="none", edgecolors="orange", linewidths=2, linestyles="--", zorder=5))

    ax.set_title(title)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.colorbar(pts, ax=ax, pad=0.01, label=color_label)
    return pts

def classify_outcome(stop_process, had_inelastic, had_elastic, true_exit_ke):
    """Human-readable summary of what happened to the primary, from MC truth."""
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

def find_primary_index(ev):
    """Index into ev.track_* of the primary track (matched by start position)."""
    matches = np.flatnonzero(
        (np.asarray(ev["track_start_x"]) == ev["true_start_x"]) &
        (np.asarray(ev["track_start_y"]) == ev["true_start_y"]) &
        (np.asarray(ev["track_start_z"]) == ev["true_start_z"])
    )
    return int(matches[0]) if len(matches) else None

# --- per-hit view of the simulation, unmasked (no good_wcte_pmts filtering applied) ---
sim_hits = ak.zip({
    "slot": sim["hit_mpmt_slot_ids"],
    "pos":  sim["hit_pmt_position_ids"],
    "q":    sim["hit_pmt_charges"],
    "t":    sim["hit_pmt_calibrated_times"],
    "hit_track_id": sim["hit_track_id"],
})

# --- display every simulated event with true_max_scatter_angle_deg > 4 deg
#     (no good_wcte_pmts mask applied) ---
MAX_SCATTER_ANGLE_CUT = 4.0
large_scatter_event_indices = np.flatnonzero(
    ak.to_numpy(sim["true_max_scatter_angle_deg"]) > MAX_SCATTER_ANGLE_CUT
)
print(f"{len(large_scatter_event_indices)} events with true_max_scatter_angle_deg > {MAX_SCATTER_ANGLE_CUT}")

for sim_event_idx in large_scatter_event_indices:
    sim_event_idx = int(sim_event_idx)
    sim_ev, sim_ev_hits = sim[sim_event_idx], sim_hits[sim_event_idx]
    scatter_angle = float(sim_ev["true_max_scatter_angle_deg"])

    fig, ax = plt.subplots(figsize=(7, 6))

    # colored by parent TRACK ID (not charge)
    track_ids_in_event = sorted(set(ak.to_list(sim_ev_hits["hit_track_id"])))
    track_id_to_color_index = {tid: i for i, tid in enumerate(track_ids_in_event)}
    hit_color_index = [track_id_to_color_index[tid] for tid in ak.to_list(sim_ev_hits["hit_track_id"])]

    sim_channel_track = channel_first_value(
        ak.to_list(sim_ev_hits["slot"]), ak.to_list(sim_ev_hits["pos"]),
        hit_color_index, ed.nChannels,
    )
    draw_event(
        ax, ed, sim_channel_track, cmap=plt.cm.tab20,
        norm=mcolors.Normalize(vmin=0, vmax=max(len(track_ids_in_event) - 1, 1)),
        title=f"MC {PARTICLE} (event {int(sim_ev['event_number'])}): {classify_outcome(str(sim_ev['stop_process']), sim_ev['had_inelastic'], sim_ev['had_elastic'], sim_ev['true_exit_ke'])}, "
              f"max scatter {scatter_angle:.1f} deg",
        color_label="track index (see legend below)",
    )
    fig.tight_layout()

    out_pdf = Path(__file__).parent / "plots" / f"event_display_evt{sim_event_idx}_scatter{scatter_angle:.1f}deg.pdf"
    fig.savefig(out_pdf)
    print(f"Saved plot to {out_pdf}")
    plt.close(fig)

    sim_hit_tids = ak.to_list(sim_ev_hits["hit_track_id"])
    hits_per_track = Counter(sim_hit_tids)

    # legend: which track index corresponds to which particle
    print(f"Track color legend for event {sim_event_idx}:")
    for tid, idx in track_id_to_color_index.items():
        matches = np.flatnonzero(np.asarray(sim_ev["track_id"]) == tid)
        n = hits_per_track[tid]
        if len(matches):
            j = matches[0]
            tag = " <- PRIMARY" if tid == sim_ev["track_id"][find_primary_index(sim_ev)] else ""
            print(f"  color {idx}: track_id={tid}  pdg={sim_ev['track_pdg'][j]}  "
                  f"process={sim_ev['track_process'][j]}  ke={sim_ev['track_ke'][j]:.2f} MeV  ({n} hits){tag}")
        else:
            print(f"  color {idx}: track_id={tid}  (dark noise / not in truth track list)" if tid == -1
                  else f"  color {idx}: track_id={tid}  (no truth match)")

