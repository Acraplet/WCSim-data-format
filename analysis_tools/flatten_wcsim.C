// flatten_wcsim.C — flatten a WCSim output file into a simple, uproot-readable tree.
//
// Reads the custom WCSim classes (needs libWCSimRoot) and writes a flat TTree
// "hits" where every digitized hit carries its charge, time and PMT (x,y,z),
// plus per-event truth (primary particle PDG, energy, vertex, direction).
//
// Run inside the container:
//   root -l -b -q 'flatten_wcsim.C("wcsim.root","flat.root")'
//
// isMDT: -1 (default) = auto-detect whether `fname` is raw WCSim output or has
// been through the MDT (Merging/Digitizing/Triggering, https://github.com/hyperk/MDT)
// tool; 0 = force raw-WCSim truth-matching; 1 = force MDT truth-matching. See
// DetectIsMDT() and the digihit truth-match block below for why this matters:
// MDT's own digitizer (HitDigitizer.cc) repurposes WCSimRootCherenkovDigiHit's
// photon-id field to carry the true parent TRACK ID directly, instead of an
// index into trig->GetCherenkovHitTimes() like raw WCSim - using the wrong
// interpretation silently collapses hit_track_id onto ~1 track per event
// (small track-id values misread as small, "successful" array indices).
R__LOAD_LIBRARY($WCSIM_BUILD_DIR/lib/libWCSimRoot.so)
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>
#include "TKey.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"

// Raw WCSim writes separate "wcsimrootevent2" / "wcsimrootevent_OD" branches in
// wcsimT alongside "wcsimrootevent" (second/OD detector copies); MDT's own
// WCRootData::CreateTree, run with its default (single-detector) branch list,
// writes only "wcsimrootevent". Checked on both plain and MDT-processed pi-/mu-
// samples in this analysis; not guaranteed for every possible MDT/WCSim
// configuration, hence the isMDT override above.
bool DetectIsMDT(TTree* t)
{
    return t->GetBranch("wcsimrootevent_OD") == nullptr;
}

void flatten_wcsim(const char* fname, const char* foutname = "flat.root", Int_t isMDT = -1)
{
    TFile* f = TFile::Open(fname);
    if (!f || !f->IsOpen()) { printf("ERROR: cannot open %s\n", fname); return; }

    TTree* t = (TTree*)f->Get("wcsimT");
    WCSimRootEvent* superevt = new WCSimRootEvent();
    t->SetBranchAddress("wcsimrootevent", &superevt);

    bool useMDTTruth = (isMDT < 0) ? DetectIsMDT(t) : (isMDT != 0);
    printf("flatten_wcsim: %s input as %s (pass isMDT=0/1 to override)\n",
           (isMDT < 0) ? "auto-detected" : "forced", useMDTTruth ? "MDT-processed" : "plain WCSim");

    WCSimRootGeom* geo = 0;
    TTree* geotree = (TTree*)f->Get("wcsimGeoT");
    geotree->SetBranchAddress("wcsimrootgeom", &geo);
    geotree->GetEntry(0);

    // Does THIS FILE's own WCSimRootTrack layout include the "largest
    // scatter" fields (added to WCSim after this flattener was first
    // written)? Checked from the file's embedded TStreamerInfo rather than
    // just calling the getters, so older files (written before that WCSim
    // change) are read fine even though the getters themselves always exist
    // in whatever (current) WCSim library this macro is compiled against -
    // ROOT's own schema evolution would already default them safely, but
    // this makes the fallback explicit and printed instead of silent.
    bool hasMaxScatterInfo = false;
    {
        TList* streamerInfoList = f->GetStreamerInfoList();
        TStreamerInfo* si = streamerInfoList ? (TStreamerInfo*)streamerInfoList->FindObject("WCSimRootTrack") : nullptr;
        if (si && si->GetElements() && si->GetElements()->FindObject("fMaxScatterAngleDeg")) {
            hasMaxScatterInfo = true;
        }
        if (hasMaxScatterInfo) {
            printf("flatten_wcsim: input file has the largest-scatter track fields - filling true_max_scatter_*\n");
        } else {
            printf("flatten_wcsim: NOTE - this information is not available because %s was produced with an "
                   "older WCSim container/build that predates the largest-scatter tracking feature "
                   "(no 'fMaxScatterAngleDeg' field in its WCSimRootTrack). "
                   "true_max_scatter_* branches will be filled with their default/sentinel values "
                   "(-1 / \"\" / 0) for every event in this file.\n", fname);
        }
    }

    TFile* fout = new TFile(foutname, "RECREATE");
    TTree* out  = new TTree("hits", "flattened WCSim digihits with geometry + truth");

    // WCSim's native geometry is in its OWN frame, offset from the WCTE
    // convention along y only (see analysis_examples/WCSim_coordinates_and_
    // mapping_example.ipynb -> WCSimCoordinateTransform, vertical_offset =
    // 424.7625 mm = 42.47625 cm). All positions written below (hit_*, vtx_*,
    // start_*, stop_*, track_start/end_*, cher_end_*) have this added to
    // their y-component so this file uses the WCTE convention directly,
    // like the rest of the analysis (x and z are unaffected).
    const float kWCSimToWCTE_Yoffset_cm = 42.47625f;

    // Diagnostic for hits whose TubeId falls outside this file's own geometry
    // (seen on MDT-processed inputs - see README "Status / known gaps"). Rather than
    // crash into an out-of-bounds TClonesArray lookup, skip the hit and keep
    // a small sample of the offending values to report at the end.
    int nBadTubeId = 0;
    std::set<int> badTubeIdSamples;
    const int kMaxBadTubeIdSamples = 20;

    // Diagnostic for the plain-WCSim truth-match path only (useMDTTruth==false):
    // dh->GetPhotonIds()[0] is supposed to index trig->GetCherenkovHitTimes(),
    // but can come back pointing outside it. Skip the truth match for that hit
    // (hit_track_id stays -999) instead of crashing into an out-of-bounds
    // TClonesArray lookup. Never triggers in MDT mode (see useMDTTruth below -
    // no array lookup happens there).
    int nBadPhotonId = 0;
    std::set<int> badPhotonIdSamples;
    const int kMaxBadPhotonIdSamples = 20;

    // hit_pmt_charges / hit_pmt_calibrated_times use the real WCTE data's
    // branch names + type (double) instead of WCSim's native float, so this
    // tree can be read with the same code as WCTE_merged_production_*.root.
    std::vector<double> hit_pmt_charges, hit_pmt_calibrated_times;
    std::vector<float> hit_x, hit_y, hit_z;
    // hit_mpmt_slot_ids / hit_pmt_position_ids follow the WCTE convention (see
    // analysis_examples/WCSim_coordinates_and_mapping_example.ipynb): the slot
    // number matches WCSim's native GetmPMTNo() directly, but the position
    // must be converted from WCSim's 1-indexed GetmPMT_PMTNo() (1-19) to the
    // WCTE 0-indexed convention (0-18) used by the real WCTE DAQ data
    // (hit_pmt_position_ids in WCTE_merged_production_*.root).
    std::vector<int>   hit_tube, hit_mpmt_slot_ids, hit_pmt_position_ids;
    // hit_track_id: track id of the parent of the (first) true Cherenkov photon
    // behind this digitized hit. A digihit can combine several true photon
    // hits within the DAQ time window (pileup/dark noise), so this stores the
    // parent track of the earliest-indexed contributor as a representative
    // value, not an exhaustive list. -1 = dark noise; -999 = no truth match found.
    std::vector<int>   hit_track_id;
    int   n_digihits = 0, true_pdg = 0, event_number = 0;
    float true_E = -1, true_p = -1, true_ke = -1, true_length = -1;
    // true_deflection_angle [degrees]: angle between the primary's INITIAL
    // direction (true_dir_*) and the straight-line chord from true_start_*
    // to true_stop_*. A cheap, whole-track proxy for "how much did this
    // track bend end-to-end" - NOT per-scatter truth (no timing/location of
    // individual scatters, since Geant4 multiple/Coulomb scattering doesn't
    // create separate tracks the way hadronic elastic/inelastic does; see
    // the README). 0 = dead straight; -1 = undefined (zero-length chord,
    // e.g. the primary stopped right where it started, or wasn't found).
    float true_deflection_angle = -1;
    float vtx_x = 0, vtx_y = 0, vtx_z = 0, dir_x = 0, dir_y = 0, dir_z = 0;
    float start_x = 0, start_y = 0, start_z = 0, stop_x = 0, stop_y = 0, stop_z = 0;

    // --- what happened to the primary (derived from its Geant4 daughters) ---
    // `stop_process` = the discrete fate at the stop point ("" = ranged out /
    // ionization stop): a hadronic process (*Inelastic / hadElastic), an in-flight
    // "Decay", or an at-rest capture (*CaptureAtRest, e.g. hBertiniCaptureAtRest
    // for stopped K-/pi-). `had_inelastic`/`had_elastic` flag hadronic interactions.
    std::string stop_process = "";
    int had_inelastic = 0, had_elastic = 0, n_prim_daughters = 0;

    // --- every saved track in the event (the primary + all its saved
    // descendants at any depth - i.e. every track WCSim kept because it, or
    // one of ITS descendants, produced a PMT-hit photon). Parallel arrays,
    // one entry per track. ---
    std::vector<int>         track_id, track_parent_id, track_pdg, track_nhits;
    std::vector<std::string> track_process;
    std::vector<float>       track_ke;
    std::vector<float>       track_start_x, track_start_y, track_start_z;
    std::vector<float>       track_end_x,   track_end_y,   track_end_z;
    std::vector<float>       track_dir_x,   track_dir_y,   track_dir_z;
    // --- where the primary itself stopped making Cherenkov light ---
    // Emission point of the farthest-along true Cherenkov photon whose parent is
    // the primary; arclength is measured along the primary direction (cm). This
    // is a lower bound: only photons that produced a PMT hit are stored.
    float cher_end_len = -1, cher_end_x = 0, cher_end_y = 0, cher_end_z = 0;
    int   cher_n_prim = 0;
    // --- did the primary leave the inner detector (vs range out in the water)? ---
    // `true_exit_ke` = kinetic energy [MeV] at the primary's blacksheet crossing
    // (the ID edge); -1 means it never crossed -> it ranged out / interacted inside.
    // `true_stopvol` is WCSim's stopping-volume code (for reference).
    int   true_stopvol = -999;
    float true_exit_ke = -1;
    // --- beam-proton lineage: # elastic scatters, and whether it ends inelastic ---
    // Follows the highest-KE hadElastic proton continuation from the primary;
    // `n_elastic` counts elastic scatters, `n_inelastic` = 1 if the chain ends in
    // a protonInelastic breakup (so n_elastic>0 & n_inelastic==1 = elastic-then-inelastic).
    int   n_elastic = 0, n_inelastic = 0;

    // --- largest single-step deflection ("scatter") along the primary track,
    // from WCSimRootTrack::GetMaxScatter* (see WCSimTrajectory::AppendStep in
    // the sim). -1 / "" = no primary found or the primary took no steps. This
    // is a genuine per-step (not whole-track chord) deflection: a discrete
    // process (hadElastic, *Inelastic, ...) shows up as one step with a large
    // angle, while continuous msc is spread over many small-angle steps.
    float true_max_scatter_angle_deg = -1;
    std::string true_max_scatter_process = "";
    float true_max_scatter_x = 0, true_max_scatter_y = 0, true_max_scatter_z = 0;
    float true_max_scatter_predir_x = 0, true_max_scatter_predir_y = 0, true_max_scatter_predir_z = 0;
    float true_max_scatter_postdir_x = 0, true_max_scatter_postdir_y = 0, true_max_scatter_postdir_z = 0;
    float true_max_scatter_ke_pre = -1, true_max_scatter_ke_post = -1;
    int   true_n_scatters_above_7deg = 0;

    out->Branch("hit_pmt_charges", &hit_pmt_charges);
    out->Branch("hit_pmt_calibrated_times", &hit_pmt_calibrated_times);
    out->Branch("hit_x", &hit_x);
    out->Branch("hit_y", &hit_y);
    out->Branch("hit_z", &hit_z);
    out->Branch("hit_tube", &hit_tube);             // WCSim tube id (1-based)
    out->Branch("hit_mpmt_slot_ids", &hit_mpmt_slot_ids);         // mPMT module slot number (WCTE convention)
    out->Branch("hit_pmt_position_ids", &hit_pmt_position_ids);   // PMT position within the module, 0-indexed (WCTE convention)
    out->Branch("hit_track_id", &hit_track_id);      // track id of the parent of this hit's (first) true Cherenkov photon
    out->Branch("n_digihits", &n_digihits);
    out->Branch("event_number", &event_number);      // placeholder = event loop index (real data's event_number has no MC equivalent)
    out->Branch("true_pdg", &true_pdg);
    out->Branch("true_E", &true_E);              // total energy [MeV]
    out->Branch("true_p", &true_p);              // momentum [MeV/c]
    out->Branch("true_ke", &true_ke);            // kinetic energy [MeV]
    out->Branch("true_length", &true_length);    // track length start->stop [cm]
    out->Branch("true_deflection_angle", &true_deflection_angle);  // angle [deg] between true_dir_* and the start->stop chord (-1 = undefined)
    out->Branch("true_vtx_x", &vtx_x);
    out->Branch("true_vtx_y", &vtx_y);
    out->Branch("true_vtx_z", &vtx_z);
    out->Branch("true_dir_x", &dir_x);
    out->Branch("true_dir_y", &dir_y);
    out->Branch("true_dir_z", &dir_z);
    out->Branch("true_start_x", &start_x);
    out->Branch("true_start_y", &start_y);
    out->Branch("true_start_z", &start_z);
    out->Branch("true_stop_x", &stop_x);
    out->Branch("true_stop_y", &stop_y);
    out->Branch("true_stop_z", &stop_z);
    out->Branch("stop_process", &stop_process);      // hadronic process at stop ("" = ranged out)
    out->Branch("had_inelastic", &had_inelastic);    // primary had a *Inelastic daughter
    out->Branch("had_elastic", &had_elastic);        // primary had a hadElastic daughter
    out->Branch("n_prim_daughters", &n_prim_daughters);
    out->Branch("track_id", &track_id);                    // every saved track (primary + all its saved descendants):
    out->Branch("track_parent_id", &track_parent_id);      //   track id of its parent (primary's parent is WCSim's placeholder incident-track entry, not itself in this list)
    out->Branch("track_pdg", &track_pdg);                  //   PDG code
    out->Branch("track_process", &track_process);          //   creator process name
    out->Branch("track_ke", &track_ke);                    //   kinetic energy [MeV]
    out->Branch("track_start_x", &track_start_x);          //   start position [cm]
    out->Branch("track_start_y", &track_start_y);
    out->Branch("track_start_z", &track_start_z);
    out->Branch("track_end_x", &track_end_x);              //   stop position [cm]
    out->Branch("track_end_y", &track_end_y);
    out->Branch("track_end_z", &track_end_z);
    out->Branch("track_dir_x", &track_dir_x);              //   initial direction (unit vector)
    out->Branch("track_dir_y", &track_dir_y);
    out->Branch("track_dir_z", &track_dir_z);
    out->Branch("track_nhits", &track_nhits);               //   # true PMT hits it produced
    out->Branch("cher_end_len", &cher_end_len);      // arclength of last primary Cherenkov emission [cm]
    out->Branch("cher_end_x", &cher_end_x);          // its position [cm]
    out->Branch("cher_end_y", &cher_end_y);
    out->Branch("cher_end_z", &cher_end_z);
    out->Branch("cher_n_prim", &cher_n_prim);        // # true Cherenkov photons from the primary
    out->Branch("true_stopvol", &true_stopvol);      // WCSim stopping-volume code
    out->Branch("true_exit_ke", &true_exit_ke);      // KE at ID-edge crossing [MeV] (-1 = ranged out inside)
    out->Branch("n_elastic", &n_elastic);            // # hadElastic scatters along the beam-proton lineage
    out->Branch("n_inelastic", &n_inelastic);        // 1 if that lineage ends in a protonInelastic breakup
    out->Branch("true_max_scatter_angle_deg", &true_max_scatter_angle_deg);  // largest single-step deflection along the primary [deg] (-1 = none)
    out->Branch("true_max_scatter_process", &true_max_scatter_process);     // process that defined that step
    out->Branch("true_max_scatter_x", &true_max_scatter_x);  // position of that step [cm]
    out->Branch("true_max_scatter_y", &true_max_scatter_y);
    out->Branch("true_max_scatter_z", &true_max_scatter_z);
    out->Branch("true_max_scatter_predir_x", &true_max_scatter_predir_x);   // track direction just before that step (unit vector)
    out->Branch("true_max_scatter_predir_y", &true_max_scatter_predir_y);
    out->Branch("true_max_scatter_predir_z", &true_max_scatter_predir_z);
    out->Branch("true_max_scatter_postdir_x", &true_max_scatter_postdir_x); // track direction just after that step (unit vector)
    out->Branch("true_max_scatter_postdir_y", &true_max_scatter_postdir_y);
    out->Branch("true_max_scatter_postdir_z", &true_max_scatter_postdir_z);
    out->Branch("true_max_scatter_ke_pre", &true_max_scatter_ke_pre);   // KE just before that step [MeV]
    out->Branch("true_max_scatter_ke_post", &true_max_scatter_ke_post); // KE just after that step [MeV]
    out->Branch("true_n_scatters_above_7deg", &true_n_scatters_above_7deg); // # steps along the primary with deflection > 7 deg

    Long64_t nev = t->GetEntries();
    for (Long64_t i = 0; i < nev; i++) {
        delete superevt; superevt = 0;          // EXTREMELY IMPORTANT (per WCSim examples)
        t->GetEntry(i);
        WCSimRootTrigger* trig = superevt->GetTrigger(0);

        hit_pmt_charges.clear(); hit_pmt_calibrated_times.clear();
        hit_x.clear(); hit_y.clear(); hit_z.clear();
        hit_tube.clear(); hit_mpmt_slot_ids.clear(); hit_pmt_position_ids.clear();
        hit_track_id.clear();
        event_number = (int)i;

        // --- truth: vertex + first primary track ---
        vtx_x = trig->GetVtx(0); vtx_y = trig->GetVtx(1) + kWCSimToWCTE_Yoffset_cm; vtx_z = trig->GetVtx(2);
        true_E = true_p = true_ke = true_length = -1; true_pdg = 0;
        true_deflection_angle = -1;
        dir_x = dir_y = dir_z = 0;
        start_x = start_y = start_z = stop_x = stop_y = stop_z = 0;
        true_stopvol = -999; true_exit_ke = -1;
        true_max_scatter_angle_deg = -1; true_max_scatter_process = "";
        true_max_scatter_x = true_max_scatter_y = true_max_scatter_z = 0;
        true_max_scatter_predir_x = true_max_scatter_predir_y = true_max_scatter_predir_z = 0;
        true_max_scatter_postdir_x = true_max_scatter_postdir_y = true_max_scatter_postdir_z = 0;
        true_max_scatter_ke_pre = true_max_scatter_ke_post = -1;
        true_n_scatters_above_7deg = 0;
        int primId = -999;
        int ntrack = trig->GetNtrack();
        for (int it = 0; it < ntrack; it++) {
            WCSimRootTrack* tr = dynamic_cast<WCSimRootTrack*>(trig->GetTracks()->At(it));
            if (!tr) continue;
            // beam primary: Parenttype 0 and Flag 0 (skip WCSim's flag -1/-2
            // incident/target pseudo-tracks, which have m=0 and a far-upstream start)
            if (tr->GetParenttype() == 0 && tr->GetFlag() == 0) {
                primId   = tr->GetId();
                true_E   = tr->GetE();
                true_p   = tr->GetP();
                true_ke  = tr->GetE() - tr->GetM();
                true_pdg = tr->GetIpnu();
                dir_x = tr->GetDir(0); dir_y = tr->GetDir(1); dir_z = tr->GetDir(2);
                start_x = tr->GetStart(0); start_y = tr->GetStart(1) + kWCSimToWCTE_Yoffset_cm; start_z = tr->GetStart(2);
                stop_x  = tr->GetStop(0);  stop_y  = tr->GetStop(1) + kWCSimToWCTE_Yoffset_cm;  stop_z  = tr->GetStop(2);
                true_length = std::sqrt((stop_x-start_x)*(stop_x-start_x) +
                                        (stop_y-start_y)*(stop_y-start_y) +
                                        (stop_z-start_z)*(stop_z-start_z));
                if (true_length > 0) {
                    float cx = (stop_x-start_x)/true_length;
                    float cy = (stop_y-start_y)/true_length;
                    float cz = (stop_z-start_z)/true_length;
                    float cosang = dir_x*cx + dir_y*cy + dir_z*cz;
                    cosang = std::max(-1.f, std::min(1.f, cosang));   // guard float rounding past [-1,1]
                    true_deflection_angle = std::acos(cosang) * (180.0f / 3.14159265358979323846f);
                }
                true_stopvol = tr->GetStopvol();
                true_exit_ke = -1;                    // KE at the ID-edge (blacksheet) crossing
                {
                    std::vector<int>   bt = tr->GetBoundaryTypes();
                    std::vector<float> bk = tr->GetBoundaryKEs();
                    for (size_t k = 0; k < bt.size() && k < bk.size(); k++)
                        if (bt[k] == 1) true_exit_ke = bk[k];   // keep the last blacksheet crossing
                }
                if (hasMaxScatterInfo) {
                    true_max_scatter_angle_deg = tr->GetMaxScatterAngleDeg();
                    true_max_scatter_process   = tr->GetMaxScatterProcess();
                    true_max_scatter_x = tr->GetMaxScatterPos(0);
                    true_max_scatter_y = tr->GetMaxScatterPos(1) + kWCSimToWCTE_Yoffset_cm;
                    true_max_scatter_z = tr->GetMaxScatterPos(2);
                    true_max_scatter_predir_x = tr->GetMaxScatterPreDir(0);
                    true_max_scatter_predir_y = tr->GetMaxScatterPreDir(1);
                    true_max_scatter_predir_z = tr->GetMaxScatterPreDir(2);
                    true_max_scatter_postdir_x = tr->GetMaxScatterPostDir(0);
                    true_max_scatter_postdir_y = tr->GetMaxScatterPostDir(1);
                    true_max_scatter_postdir_z = tr->GetMaxScatterPostDir(2);
                    true_max_scatter_ke_pre  = tr->GetMaxScatterKEPre();
                    true_max_scatter_ke_post = tr->GetMaxScatterKEPost();
                    true_n_scatters_above_7deg = tr->GetNScattersAbove7Deg();
                }
                break;
            }
        }

        // --- primary's direct daughters: process at/along the track ---
        stop_process = ""; had_inelastic = 0; had_elastic = 0; n_prim_daughters = 0;
        n_elastic = 0; n_inelastic = 0;
        // true PMT hits per creating track: parentSavedTrackID -> count
        std::map<int,int> hitsByParent;
        for (int ih = 0; ih < trig->GetNcherenkovhittimes(); ih++) {
            WCSimRootCherenkovHitTime* ht =
                dynamic_cast<WCSimRootCherenkovHitTime*>(trig->GetCherenkovHitTimes()->At(ih));
            if (ht) hitsByParent[ht->GetParentSavedTrackID()]++;
        }
        float best_d2 = 1e18f;   // pick the hadronic daughter nearest the primary stop
        if (primId != -999) {
            for (int it = 0; it < ntrack; it++) {
                WCSimRootTrack* tr = dynamic_cast<WCSimRootTrack*>(trig->GetTracks()->At(it));
                if (!tr || tr->GetParentId() != primId || tr->GetFlag() != 0) continue;
                std::string proc = tr->GetCreatorProcessName();
                float dsx = tr->GetStart(0), dsy = tr->GetStart(1) + kWCSimToWCTE_Yoffset_cm, dsz = tr->GetStart(2);
                n_prim_daughters++;
                bool inel = proc.find("Inelastic") != std::string::npos;
                bool elas = (proc == "hadElastic");
                bool decy = (proc == "Decay");
                bool capt = proc.find("CaptureAtRest") != std::string::npos;
                if (inel) had_inelastic = 1;
                if (elas) had_elastic = 1;
                if (inel || elas || decy || capt) {   // discrete fate: keep the one at the stop
                    float d2 = (dsx-stop_x)*(dsx-stop_x) + (dsy-stop_y)*(dsy-stop_y) +
                               (dsz-stop_z)*(dsz-stop_z);
                    if (d2 < best_d2) { best_d2 = d2; stop_process = proc; }
                }
            }
        }

        // --- every saved track in the event (primary + all its saved
        // descendants at any depth), independent of the direct-daughters
        // loop above which only looks at the primary's immediate fate ---
        track_id.clear(); track_parent_id.clear(); track_pdg.clear(); track_nhits.clear();
        track_process.clear(); track_ke.clear();
        track_start_x.clear(); track_start_y.clear(); track_start_z.clear();
        track_end_x.clear();   track_end_y.clear();   track_end_z.clear();
        track_dir_x.clear();   track_dir_y.clear();   track_dir_z.clear();
        for (int it = 0; it < ntrack; it++) {
            WCSimRootTrack* tr = dynamic_cast<WCSimRootTrack*>(trig->GetTracks()->At(it));
            if (!tr || tr->GetFlag() != 0) continue;   // skip WCSim's -1/-2 incident/target pseudo-tracks
            int tid = tr->GetId();
            track_id.push_back(tid);
            track_parent_id.push_back(tr->GetParentId());
            track_pdg.push_back(tr->GetIpnu());
            track_process.push_back(tr->GetCreatorProcessName());
            track_ke.push_back(tr->GetE() - tr->GetM());
            track_start_x.push_back(tr->GetStart(0)); track_start_y.push_back(tr->GetStart(1) + kWCSimToWCTE_Yoffset_cm); track_start_z.push_back(tr->GetStart(2));
            track_end_x.push_back(tr->GetStop(0));    track_end_y.push_back(tr->GetStop(1) + kWCSimToWCTE_Yoffset_cm);    track_end_z.push_back(tr->GetStop(2));
            track_dir_x.push_back(tr->GetDir(0)); track_dir_y.push_back(tr->GetDir(1)); track_dir_z.push_back(tr->GetDir(2));
            track_nhits.push_back(hitsByParent.count(tid) ? hitsByParent[tid] : 0);
        }

        // --- follow the beam-particle lineage: elastic scatters, then breakup? ---
        // Each track ends in ONE hard process; hadElastic spawns a scattered
        // same-species continuation (take the highest-KE one), *Inelastic ends the
        // chain. Follows the primary's own PDG (proton, kaon-, ...).
        if (primId != -999) {
            int cur = primId;
            std::set<int> visited;
            while (cur != -999 && !visited.count(cur)) {
                visited.insert(cur);
                int nextSame = -999; float bestKE = -1; bool curInelastic = false;
                for (int it = 0; it < ntrack; it++) {
                    WCSimRootTrack* tr = dynamic_cast<WCSimRootTrack*>(trig->GetTracks()->At(it));
                    if (!tr || tr->GetParentId() != cur || tr->GetFlag() != 0) continue;
                    std::string proc = tr->GetCreatorProcessName();
                    if (proc.find("Inelastic") != std::string::npos) {
                        curInelastic = true;
                    } else if (proc == "hadElastic" && tr->GetIpnu() == true_pdg) {
                        float ke = tr->GetE() - tr->GetM();     // scattered beam particle = highest KE
                        if (ke > bestKE) { bestKE = ke; nextSame = tr->GetId(); }
                    }
                }
                if (curInelastic) { n_inelastic++; break; }     // chain ends in breakup
                if (nextSame != -999) { n_elastic++; cur = nextSame; }  // elastic -> continue
                else break;                                     // ranged out / decayed / captured
            }
        }

        // --- Cherenkov emission endpoint of the primary (photon pos: mm -> cm) ---
        cher_end_len = -1; cher_end_x = cher_end_y = cher_end_z = 0; cher_n_prim = 0;
        if (primId != -999) {
            float best_s = -1e18f;
            int nht = trig->GetNcherenkovhittimes();
            for (int ih = 0; ih < nht; ih++) {
                WCSimRootCherenkovHitTime* ht =
                    dynamic_cast<WCSimRootCherenkovHitTime*>(trig->GetCherenkovHitTimes()->At(ih));
                if (!ht || ht->GetParentSavedTrackID() != primId) continue;
                float px = ht->GetPhotonStartPos(0) * 0.1f;   // mm -> cm
                float py = ht->GetPhotonStartPos(1) * 0.1f + kWCSimToWCTE_Yoffset_cm;
                float pz = ht->GetPhotonStartPos(2) * 0.1f;
                float s = (px-start_x)*dir_x + (py-start_y)*dir_y + (pz-start_z)*dir_z;
                cher_n_prim++;
                if (s > best_s) { best_s = s; cher_end_x = px; cher_end_y = py; cher_end_z = pz; }
            }
            if (cher_n_prim > 0) cher_end_len = best_s;
        }

        // --- digitized hits + PMT positions ---
        int ndigi = trig->GetNcherenkovdigihits();
        for (int idigi = 0; idigi < ndigi; idigi++) {
            WCSimRootCherenkovDigiHit* dh =
                dynamic_cast<WCSimRootCherenkovDigiHit*>(trig->GetCherenkovDigiHits()->At(idigi));
            if (!dh) continue;
            int tubeId = dh->GetTubeId();
            if (tubeId < 1 || tubeId > geo->GetWCNumPMT()) {
                nBadTubeId++;
                if ((int)badTubeIdSamples.size() < kMaxBadTubeIdSamples) badTubeIdSamples.insert(tubeId);
                continue;   // don't index geo->GetPMT() out of bounds
            }
            WCSimRootPMT pmt = geo->GetPMT(tubeId - 1);
            hit_pmt_charges.push_back(dh->GetQ());
            hit_pmt_calibrated_times.push_back(dh->GetT());
            hit_x.push_back(pmt.GetPosition(0));
            hit_y.push_back(pmt.GetPosition(1) + kWCSimToWCTE_Yoffset_cm);
            hit_z.push_back(pmt.GetPosition(2));
            hit_tube.push_back(tubeId);
            hit_mpmt_slot_ids.push_back(pmt.GetmPMTNo());
            hit_pmt_position_ids.push_back(pmt.GetmPMT_PMTNo() - 1); // 1-indexed -> WCTE 0-indexed

            // truth-match: which track's Cherenkov photon produced this digit?
            // a digihit can combine several true photon hits (pileup within the
            // digitization window / dark noise); take the parent track of the
            // first (earliest-indexed) contributing true hit as a representative value.
            //
            // The two producers disagree on what dh->GetPhotonIds() actually holds:
            //  - plain WCSim: indices into trig->GetCherenkovHitTimes() - look the
            //    entry up and read ITS GetParentSavedTrackID().
            //  - MDT (HitDigitizer.cc's parent_composition, built from
            //    TrueHit::GetParentId()): the true parent TRACK ID itself, already -
            //    no CherenkovHitTimes lookup needed (or valid: that array is laid
            //    out differently in MDT output, indexed per-tube-across-the-whole
            //    -event, not addressable by these ids at all).
            int trackId = -999;
            std::vector<int> photonIds = dh->GetPhotonIds();
            if (!photonIds.empty()) {
                if (useMDTTruth) {
                    trackId = photonIds[0];   // already a track id; -1 = dark noise (MDT's own convention)
                } else {
                    int photonId = photonIds[0];
                    if (photonId < 0 || photonId >= trig->GetNcherenkovhittimes()) {
                        nBadPhotonId++;
                        if ((int)badPhotonIdSamples.size() < kMaxBadPhotonIdSamples) badPhotonIdSamples.insert(photonId);
                    } else {
                        WCSimRootCherenkovHitTime* ht =
                            dynamic_cast<WCSimRootCherenkovHitTime*>(trig->GetCherenkovHitTimes()->At(photonId));
                        if (ht) trackId = ht->GetParentSavedTrackID();   // -1 = dark noise
                    }
                }
            }
            hit_track_id.push_back(trackId);
        }
        n_digihits = (int)hit_pmt_charges.size();
        out->Fill();
    }

    fout->cd();
    out->Write();

    // Copy every other TTree in the input file to the output unchanged
    // (byte-for-byte, no reprocessing) - e.g. "AllSecondaries",
    // "AllSecondaryPhotons", "Settings", "wcsimGeoT", "wcsimRootOptionsT" -
    // so the flat file also carries whatever else WCSim wrote, exactly as
    // it wrote it. Skips "wcsimT" itself, which is what "hits" above was
    // built from. Each tree name is only copied once even though ROOT
    // lists one key per write cycle (e.g. "AllSecondaries;1" and ";2");
    // f->Get(name) always returns the latest cycle.
    std::set<std::string> copiedTreeNames;
    TIter nextkey(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        if (std::string(key->GetClassName()) != "TTree") continue;
        std::string treeName = key->GetName();
        if (treeName == "wcsimT" || copiedTreeNames.count(treeName)) continue;
        copiedTreeNames.insert(treeName);
        TTree* srcTree = (TTree*)f->Get(treeName.c_str());
        if (!srcTree) continue;
        fout->cd();
        TTree* treeCopy = srcTree->CloneTree(-1, "fast");
        treeCopy->Write();
        printf("flatten_wcsim: copied tree '%s' unchanged (%lld entries)\n",
               treeName.c_str(), treeCopy->GetEntries());
    }

    fout->Close();
    f->Close();
    printf("Wrote %lld events to %s\n", nev, foutname);
    if (nBadTubeId > 0) {
        printf("WARNING: skipped %d digihit(s) with TubeId outside [1,%d] (this file's PMT count). Sample bad values: ",
               nBadTubeId, geo->GetWCNumPMT());
        for (int v : badTubeIdSamples) printf("%d ", v);
        printf("\n");
    }
    if (nBadPhotonId > 0) {
        printf("WARNING: skipped truth-match for %d digihit(s) with a photon id outside this trigger's "
               "CherenkovHitTimes array. Sample bad values: ", nBadPhotonId);
        for (int v : badPhotonIdSamples) printf("%d ", v);
        printf("\n");
    }
}
