// read_flatten_root_file.C
//
// Basic template macro to read the flattened WCSim digihit file:
//   /eos/user/a/acraplet/WCSim/wcsim-wcte-easy/notebooks/flat_100MeV_30cm_0000.root
//
// The file contains a single TTree called "hits" ("flattened WCSim digihits
// with geometry + truth"). Each entry in the tree corresponds to one event.
// Per-hit quantities (there can be many hits per event) are stored as
// std::vector<...> branches, while per-event quantities (true vertex,
// direction, energy, etc.) are stored as simple scalar branches.
//
// Run with:
//   root -l read_flatten_root_file.C
// or from inside a ROOT session:
//   .x read_flatten_root_file.C

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <vector>
#include <string>
#include <iostream>

void read_flatten_root_file()
{
    // ------------------------------------------------------------------
    // 1) Open the input file and grab the TTree
    // ------------------------------------------------------------------
    const char* file_path =
        "/eos/user/a/acraplet/WCSim/wcsim-wcte-easy/notebooks/flat_100MeV_30cm_0000.root";

    TFile* file = TFile::Open(file_path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: could not open file " << file_path << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("hits");
    if (!tree) {
        std::cerr << "ERROR: could not find TTree 'hits' in file" << std::endl;
        file->Close();
        return;
    }

    Long64_t n_entries = tree->GetEntries();
    std::cout << "Opened " << file_path << std::endl;
    std::cout << "Tree 'hits' has " << n_entries << " entries (events)" << std::endl;

    // ------------------------------------------------------------------
    // 2) Declare variables to hold branch data, and set branch addresses
    //
    //    Per-hit branches are vectors (one entry per digitized hit in the
    //    event). Per-event branches are simple scalars (one value per
    //    event). ROOT requires pointers for vector branches when using
    //    SetBranchAddress.
    // ------------------------------------------------------------------

    // --- Per-hit (vector) branches: hit charge, time, position, tube IDs ---
    // hit_pmt_charges / hit_pmt_calibrated_times use the real WCTE data's
    // branch names + type (double), so this file can be read with the same
    // code as WCTE_merged_production_*.root.
    std::vector<double>* hit_pmt_charges          = nullptr; // digitized hit charge (p.e.)
    std::vector<double>* hit_pmt_calibrated_times = nullptr; // digitized hit time (ns)
    std::vector<float>* hit_x          = nullptr; // hit PMT position, x (cm)
    std::vector<float>* hit_y          = nullptr; // hit PMT position, y (cm)
    std::vector<float>* hit_z          = nullptr; // hit PMT position, z (cm)
    std::vector<int>*   hit_tube              = nullptr; // WCSim tube (PMT) ID
    std::vector<int>*   hit_mpmt_slot_ids     = nullptr; // mPMT module slot number (WCTE convention)
    std::vector<int>*   hit_pmt_position_ids  = nullptr; // PMT position within its mPMT, 0-indexed (WCTE convention)
    std::vector<int>*   hit_track_id          = nullptr; // track id of the parent of this hit's (first) true Cherenkov photon

    // --- Per-event scalar branches: digihit count and MC truth ---
    Int_t   n_digihits    = -999; // number of digitized hits in this event
    Int_t   event_number  = -999; // event loop index (MC placeholder for the data's event_number)
    Int_t   true_pdg      = -999; // PDG code of the primary particle
    Float_t true_E        = -999; // true total energy (MeV)
    Float_t true_p        = -999; // true momentum magnitude (MeV/c)
    Float_t true_ke       = -999; // true kinetic energy (MeV)
    Float_t true_length   = -999; // true track length (cm)
    Float_t true_vtx_x    = -999; // true vertex position, x (cm)
    Float_t true_vtx_y    = -999; // true vertex position, y (cm)
    Float_t true_vtx_z    = -999; // true vertex position, z (cm)
    Float_t true_dir_x    = -999; // true initial direction, x
    Float_t true_dir_y    = -999; // true initial direction, y
    Float_t true_dir_z    = -999; // true initial direction, z
    Float_t true_start_x  = -999; // true track start position, x (cm)
    Float_t true_start_y  = -999; // true track start position, y (cm)
    Float_t true_start_z  = -999; // true track start position, z (cm)
    Float_t true_stop_x   = -999; // true track stop position, x (cm)
    Float_t true_stop_y   = -999; // true track stop position, y (cm)
    Float_t true_stop_z   = -999; // true track stop position, z (cm)
    std::string* stop_process = nullptr; // Geant4 process name that stopped the primary
    Int_t   had_inelastic = -999; // flag/count: hadronic inelastic interaction occurred
    Int_t   had_elastic   = -999; // flag/count: hadronic elastic interaction occurred
    Int_t   n_prim_daughters = -999; // number of secondary particles from the primary

    // --- Track (vector) branches: EVERY saved track in the event (the
    // primary + all its saved descendants at any depth), one entry per track ---
    std::vector<int>*    track_id        = nullptr; // track ID
    std::vector<int>*    track_parent_id = nullptr; // track ID of this track's parent
    std::vector<int>*    track_pdg       = nullptr; // PDG code of each track
    std::vector<std::string>* track_process = nullptr; // creation process of each track
    std::vector<float>*  track_ke      = nullptr; // kinetic energy of each track (MeV)
    std::vector<float>*  track_start_x = nullptr; // track start position, x (cm)
    std::vector<float>*  track_start_y = nullptr; // track start position, y (cm)
    std::vector<float>*  track_start_z = nullptr; // track start position, z (cm)
    std::vector<float>*  track_end_x   = nullptr; // track stop position, x (cm)
    std::vector<float>*  track_end_y   = nullptr; // track stop position, y (cm)
    std::vector<float>*  track_end_z   = nullptr; // track stop position, z (cm)
    std::vector<float>*  track_dir_x   = nullptr; // track initial direction, x
    std::vector<float>*  track_dir_y   = nullptr; // track initial direction, y
    std::vector<float>*  track_dir_z   = nullptr; // track initial direction, z
    std::vector<int>*    track_nhits   = nullptr; // number of true PMT hits produced by each track

    // --- Cherenkov / stopping-related scalar branches ---
    Float_t cher_end_len  = -999; // path length to end of Cherenkov emission (cm)
    Float_t cher_end_x    = -999; // Cherenkov emission end position, x (cm)
    Float_t cher_end_y    = -999; // Cherenkov emission end position, y (cm)
    Float_t cher_end_z    = -999; // Cherenkov emission end position, z (cm)
    Int_t   cher_n_prim   = -999; // number of primary Cherenkov photons/tracks
    Int_t   true_stopvol  = -999; // ID of the volume where the primary stopped
    Float_t true_exit_ke  = -999; // kinetic energy when exiting the fiducial volume (MeV)
    Int_t   n_elastic     = -999; // number of elastic scatters
    Int_t   n_inelastic   = -999; // number of inelastic scatters

    // Attach each variable above to its corresponding branch in the tree
    tree->SetBranchAddress("hit_pmt_charges",          &hit_pmt_charges);
    tree->SetBranchAddress("hit_pmt_calibrated_times", &hit_pmt_calibrated_times);
    tree->SetBranchAddress("hit_x",           &hit_x);
    tree->SetBranchAddress("hit_y",           &hit_y);
    tree->SetBranchAddress("hit_z",           &hit_z);
    tree->SetBranchAddress("hit_tube",        &hit_tube);
    tree->SetBranchAddress("hit_mpmt_slot_ids",    &hit_mpmt_slot_ids);
    tree->SetBranchAddress("hit_pmt_position_ids", &hit_pmt_position_ids);
    tree->SetBranchAddress("hit_track_id",         &hit_track_id);

    tree->SetBranchAddress("n_digihits",  &n_digihits);
    tree->SetBranchAddress("event_number", &event_number);
    tree->SetBranchAddress("true_pdg",    &true_pdg);
    tree->SetBranchAddress("true_E",      &true_E);
    tree->SetBranchAddress("true_p",      &true_p);
    tree->SetBranchAddress("true_ke",     &true_ke);
    tree->SetBranchAddress("true_length", &true_length);
    tree->SetBranchAddress("true_vtx_x",  &true_vtx_x);
    tree->SetBranchAddress("true_vtx_y",  &true_vtx_y);
    tree->SetBranchAddress("true_vtx_z",  &true_vtx_z);
    tree->SetBranchAddress("true_dir_x",  &true_dir_x);
    tree->SetBranchAddress("true_dir_y",  &true_dir_y);
    tree->SetBranchAddress("true_dir_z",  &true_dir_z);
    tree->SetBranchAddress("true_start_x", &true_start_x);
    tree->SetBranchAddress("true_start_y", &true_start_y);
    tree->SetBranchAddress("true_start_z", &true_start_z);
    tree->SetBranchAddress("true_stop_x", &true_stop_x);
    tree->SetBranchAddress("true_stop_y", &true_stop_y);
    tree->SetBranchAddress("true_stop_z", &true_stop_z);
    tree->SetBranchAddress("stop_process", &stop_process);
    tree->SetBranchAddress("had_inelastic", &had_inelastic);
    tree->SetBranchAddress("had_elastic",   &had_elastic);
    tree->SetBranchAddress("n_prim_daughters", &n_prim_daughters);

    tree->SetBranchAddress("track_id",        &track_id);
    tree->SetBranchAddress("track_parent_id", &track_parent_id);
    tree->SetBranchAddress("track_pdg",       &track_pdg);
    tree->SetBranchAddress("track_process",   &track_process);
    tree->SetBranchAddress("track_ke",        &track_ke);
    tree->SetBranchAddress("track_start_x",   &track_start_x);
    tree->SetBranchAddress("track_start_y",   &track_start_y);
    tree->SetBranchAddress("track_start_z",   &track_start_z);
    tree->SetBranchAddress("track_end_x",     &track_end_x);
    tree->SetBranchAddress("track_end_y",     &track_end_y);
    tree->SetBranchAddress("track_end_z",     &track_end_z);
    tree->SetBranchAddress("track_dir_x",     &track_dir_x);
    tree->SetBranchAddress("track_dir_y",     &track_dir_y);
    tree->SetBranchAddress("track_dir_z",     &track_dir_z);
    tree->SetBranchAddress("track_nhits",     &track_nhits);

    tree->SetBranchAddress("cher_end_len", &cher_end_len);
    tree->SetBranchAddress("cher_end_x",   &cher_end_x);
    tree->SetBranchAddress("cher_end_y",   &cher_end_y);
    tree->SetBranchAddress("cher_end_z",   &cher_end_z);
    tree->SetBranchAddress("cher_n_prim",  &cher_n_prim);
    tree->SetBranchAddress("true_stopvol", &true_stopvol);
    tree->SetBranchAddress("true_exit_ke", &true_exit_ke);
    tree->SetBranchAddress("n_elastic",    &n_elastic);
    tree->SetBranchAddress("n_inelastic",  &n_inelastic);

    // ------------------------------------------------------------------
    // 3) Loop over all events (tree entries) and access the data
    // ------------------------------------------------------------------
    for (Long64_t i = 0; i < n_entries; ++i) {

        // Load the branch data for this entry into the variables above
        tree->GetEntry(i);

        std::cout << "\n--- Event " << i << " ---" << std::endl;
        std::cout << "  event_number = " << event_number
                   << ", true_pdg = " << true_pdg
                   << ", true_E = " << true_E << " MeV"
                   << ", n_digihits = " << n_digihits << std::endl;
        std::cout << "  true vertex = ("
                   << true_vtx_x << ", " << true_vtx_y << ", " << true_vtx_z
                   << ") cm"  << std::endl;
        std::cout << "  true stop point = ("
                   << true_stop_x << ", " << true_stop_y << ", " << true_stop_z
                   << ") cm"  << std::endl;
        std::cout  << " number of elastic scatters = " << n_elastic 
                   << " number of inelastic scatters = " << n_inelastic
                   << std::endl;
        std::cout  << " total Cherenkov photons from primary = " << cher_n_prim
                   << " Cherenkov emission path length = " << cher_end_len << " cm"
                   
                   << std::endl;
        std::cout  << std::endl;

        // Example: loop over the digitized hits of this event
        // (hit_pmt_charges, hit_pmt_calibrated_times, hit_x, ... all have the
        // same length == n_digihits)
        if (hit_pmt_charges) {
            for (size_t h = 0; h < hit_pmt_charges->size(); ++h) {
                // Uncomment to print every single hit (can be verbose!)
                // std::cout << "    hit " << h
                //           << ": q=" << hit_pmt_charges->at(h)
                //           << " t=" << hit_pmt_calibrated_times->at(h)
                //           << " pos=(" << hit_x->at(h) << ", "
                //           << hit_y->at(h) << ", " << hit_z->at(h) << ")"
                //           << " tube=" << hit_tube->at(h)
                //           << " track_id=" << hit_track_id->at(h) << std::endl;
            }
        }

        // Example: loop over every track in this event (primary + all its
        // saved descendants at any depth)
        if (track_id) {
            for (size_t s = 0; s < track_id->size(); ++s) {
                std::cout << "    track " << s
                          << ": id=" << track_id->at(s)
                          << " parent_id=" << track_parent_id->at(s)
                          << " pdg=" << track_pdg->at(s)
                          << " ke=" << track_ke->at(s) << " MeV"
                          << " start=(" << track_start_x->at(s) << ", "
                          << track_start_y->at(s) << ", " << track_start_z->at(s) << ") cm"
                          << " end=(" << track_end_x->at(s) << ", "
                          << track_end_y->at(s) << ", " << track_end_z->at(s) << ") cm"
                          << " nhits=" << track_nhits->at(s)
                          << " process=" << track_process->at(s) << std::endl;
            }
        }
    }

    // ------------------------------------------------------------------
    // 4) Clean up
    // ------------------------------------------------------------------
    file->Close();
}
