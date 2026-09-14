# WCTE sim/data comparison tools (work in progress)

Early-stage code for flattening WCSim output and comparing it against real WCTE
data. **Not finished** - the plan is to merge this into the main
[`analysis_tools`](https://github.com/) repo once it's stable, so the layout
here loosely mirrors that repo's structure (`analysis_tools/` for reusable
code, `examples/` for notebooks) to make that merge easier later.

## Contents

- `analysis_tools/flatten_wcsim.C` - (Copied and adapated from Bruno's repo:  https://github.com/bferrazzi/wcsim-wcte-easy/tree/main) 
  ROOT macro that flattens a raw WCSim
  output file into a simple, uproot-readable tree (`hits`). Mirrors the real
  WCTE data schema where possible (branch names, WCTE-frame coordinates,
  WCTE mPMT slot/position convention) so downstream code can treat simulation
  and data uniformly. Must be run inside the WCSim container
  (needs `libWCSimRoot.so`):
  ```
  root -l -b -q 'flatten_wcsim.C("wcsim.root","flat.root")'
  ```
- `analysis_tools/read_flatten_root_file.C` - basic template macro
  demonstrating how to read every branch of the flattened tree.
- `analysis_tools/read_data_root_file.C` - basic template macro demonstrating
  how to read every branch of the real merged WCTE data file
  (`WCTEReadoutWindows` plus the single-entry summary trees).
- `examples/compare_sim_data.ipynb` - template notebook comparing a flattened
  simulation file against real merged WCTE data: loads/selects muons from data
  the standard way (`analysis_tools.DataLoader` + `BeamSelection`), applies
  the same `good_wcte_pmts` channel mask to both, compares hit multiplicities,
  summarizes simulated track/interaction outcomes, and produces a side-by-side
  event display (data vs. simulation, simulation hits colored by parent track).

## Status / known gaps

- `flatten_wcsim.C` needs to be re-run inside the WCSim container to
  regenerate flat files with the current branch schema - some existing flat
  files predate recent schema changes.
- No automated tests yet.
- Not yet merged into `analysis_tools`; branch/file names may still change.
- **MDT-processed files** (`*_MDT.root`, produced by the separate
  [MDT](https://github.com/hyperk/MDT) tool - Merging/Digitizing/Triggering,
  downstream of raw WCSim output) currently fail to flatten cleanly: some
  digihits carry a `TubeId` far outside the file's own geometry (e.g. index
  ~78865 into a 1644/1844-PMT array), which used to crash `flatten_wcsim.C`
  outright via an out-of-bounds `TClonesArray::At`. Investigated so far:
  MDT's own source passes `TubeId` through **unmodified** from the input
  WCSim file (`WCRootData.cc:89,363-368` - no merged/global re-encoding), and
  its `wcsimGeoT` output is a byte-for-byte copy of the input geometry (not
  trimmed) - so the mismatch isn't explained by either of those. One untested
  lead: `parameter/MDTParamenter_WCTE.txt` hardcodes the nominal full-WCTE
  `MaxTubeID` (2014) regardless of the actual run's PMT count, which could be
  involved in dark-noise hit generation. `flatten_wcsim.C` now skips any hit
  with `TubeId` outside `[1, geo->GetWCNumPMT()]` instead of crashing, and
  prints a count + sample of the bad values at the end - rerun inside the
  container to get real numbers for further diagnosis. Plain (non-MDT)
  WCSim output is unaffected.
- **`hit_track_id` on MDT-processed files was reading the wrong field
  entirely**, not just occasionally out-of-bounds. `flatten_wcsim.C` (like the
  official WCSim examples) assumed `dh->GetPhotonIds()[0]` is an index into
  `trig->GetCherenkovHitTimes()`, which is true for plain WCSim output but not
  for MDT: MDT's own digitizer
  ([`HitDigitizer.cc`](https://github.com/hyperk/MDT/blob/angular_response_gain/cpp/src/HitDigitizer.cc),
  `parent_composition.push_back(PEs[iPE]->GetParentId())`) writes the true
  parent **track ID directly** into that same field
  ([`WCRootData.cc`](https://github.com/hyperk/MDT/blob/angular_response_gain/app/utilities/WCRootData/src/WCRootData.cc),
  `true_pe_comp = aPH->GetParentCompositionDigi(i)` passed straight to
  `AddCherenkovDigiHit`). Treating a track ID as an array index either throws
  `Error in <TClonesArray::At>: index ... out of bounds` (large track IDs,
  e.g. ~79403 into a 1662-entry array) or "succeeds" by accident on small
  track IDs (e.g. 1, 2, 3 - common ones, since low IDs tend to be the primary
  and early secondaries) and returns an unrelated hit's parent track - which
  is why `hit_track_id` on MDT files looked collapsed onto ~1 dominant track
  per event instead of the real per-hit diversity seen on plain WCSim output.
  `flatten_wcsim.C` now takes an `isMDT` argument (default `-1` = auto-detect,
  via `DetectIsMDT()`: plain WCSim writes `wcsimrootevent2`/`wcsimrootevent_OD`
  branches in `wcsimT` alongside `wcsimrootevent`, MDT (run with its default
  branch list) writes only `wcsimrootevent` - checked on both plain and
  MDT-processed pi-/mu- samples here, but pass `isMDT=0`/`1` explicitly if a
  file doesn't match this pattern) and uses `photonIds[0]` directly as the
  track id in MDT mode, with no `CherenkovHitTimes` lookup at all. The
  `nBadPhotonId` bounds-check diagnostic below only applies to (and only fires
  in) plain-WCSim mode.

## Notebook outputs

Notebooks in `examples/` are committed **without cell outputs** to keep diffs
readable - only source changes when the code itself changes, not every time
someone re-runs a cell. This is enforced via a git clean filter
(`tools/strip_notebook_outputs.py`) rather than a separate package
(`nbstripout`) so no extra install is required.

If you clone this repo fresh, enable the filter once with:
```
git config filter.strip-notebook-output.clean "python3 tools/strip_notebook_outputs.py"
git config filter.strip-notebook-output.smudge cat
```
Your local working copy keeps whatever outputs you produce by running the
notebook - only what gets staged/committed is stripped.


## WCSim Truth information

For each event we store the hit information, some true information about the event and the 
track information where each track correspond to one particle produced. 

The tracks are stored as an array (similar to all the hits) and for each hit the ID of the 
track that produced it is stored in `hit_track_id`. The ID of the parent particle track is 
stored in `track_parent_id` which allows us to trace the lineage of a given track. 

### `true_deflection_angle`

The angle (degrees) between the primary's initial direction (`true_dir_*`) and the
straight-line chord from `true_start_*` to `true_stop_*`. This is a cheap, whole-track
proxy for "how much did this track bend end-to-end" - **not** per-scatter truth: it
carries no information about *where* along the track the bending happened, or whether
it came from one big kink or many small ones. `0` means dead straight; `-1` means
undefined (zero-length chord).

In particular this does **not** let you distinguish continuous electromagnetic multiple
scattering (`msc` - many tiny random-angle kicks smeared along the whole path, which in
Geant4 deflects the *same* track in place rather than creating a new one) from a single
large-angle interaction. Geant4 doesn't record per-step direction/process information
for charged particles anywhere in the current WCSim output, so genuine per-scatter EM
truth would require instrumenting WCSim's simulation source itself (extending
`WCSimTrajectory`'s `AppendStep`, which already records per-step *positions* in memory
for visualization but never writes them to ROOT) and rebuilding - out of scope here.

For hadronic beam particles, a large, localized deflection is more reliably diagnosed
via `had_elastic`/`n_elastic` (an actual `hadElastic` interaction, which *does* create a
new track - see below) than via this angle alone.

### Hadronic-interaction flags: `had_inelastic`, `had_elastic`, `n_elastic`, `n_inelastic`

These four fields all describe whether the *primary* particle underwent a hard hadronic
interaction, but they answer two different questions and aren't interchangeable.

- **`had_inelastic` / `had_elastic`** - simple 0/1 flags looking only at the primary's own
  **direct** daughters (one generation deep): is there a daughter created by a Geant4
  process whose name contains `"Inelastic"` (e.g. `protonInelastic`, `pi-Inelastic`), or
  by a process that is exactly `"hadElastic"`? They don't care what species that daughter
  is, and they say nothing about anything beyond the primary's first interaction vertex.

- **`n_elastic` / `n_inelastic`** - follow the primary's own **lineage** as far as it goes,
  vertex by vertex. Starting from the primary, at each step look among its direct
  daughters for: (a) any daughter created via a process containing `"Inelastic"`, or
  (b) a daughter with the **same PDG code as the primary** created via `hadElastic` (if
  several qualify, take the highest-KE one - the "surviving" scattered beam particle). If
  (a) is found the chain stops there and `n_inelastic` is set to 1. Otherwise, if (b) is
  found, `n_elastic` is incremented and the walk continues from that scattered particle.
  The walk stops (`n_inelastic` staying 0) as soon as neither is found, i.e. the particle
  ranged out, decayed, or was captured at rest with no further hard interaction.

  Because this follows the whole chain rather than stopping after one generation,
  `n_elastic` can be greater than 1 (several sequential elastic scatters), and
  `n_inelastic` can be 1 even when `had_inelastic` is 0 - if the eventual breakup happens
  after one or more elastic scatters rather than at the very first vertex.

  Example combinations, for a hadronic beam particle (proton, pion, kaon, ...):

  | `had_inelastic` | `had_elastic` | `n_elastic` | `n_inelastic` | Meaning |
  |---|---|---|---|---|
  | 0 | 0 | 0 | 0 | No hadronic interaction at all - ranged out, decayed, or captured at rest |
  | 1 | 0 | 0 | 1 | Broke up in a single inelastic interaction straight away |
  | 0 | 1 | 1 | 0 | Scattered elastically once, then ranged out/decayed/captured - no breakup |
  | 0 | 1 | 2 | 1 | Scattered elastically twice in a row, then broke up on the third interaction |

  This lineage-following logic is written around hadronic process names (`hadElastic` /
  `*Inelastic`), so for lepton beams (muons, electrons) it will normally read all zero -
  their dominant processes (`muIoni`, `muMinusCaptureAtRest`, `Decay`, `muBrems`, ...)
  don't match either pattern. `stop_process` (the creator-process name of whichever direct
  daughter sits closest to the primary's stop point) is usually the more informative field
  for those particles.