# WCTE sim/data comparison tools (work in progress)

Early-stage code for flattening WCSim output and comparing it against real WCTE
data. **Not finished** - the plan is to merge this into the main
[`analysis_tools`](https://github.com/) repo once it's stable, so the layout
here loosely mirrors that repo's structure (`analysis_tools/` for reusable
code, `examples/` for notebooks) to make that merge easier later.

## Contents

- `analysis_tools/flatten_wcsim.C` - ROOT macro that flattens a raw WCSim
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
