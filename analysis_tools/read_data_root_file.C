// read_data_root_file.C
//
// Basic template macro to read the merged WCTE real-data production file:
//   /eos/user/a/acraplet/WCSim/Analysis/data/WCTE_merged_production_R1670.root
//
// Unlike the simulation flat-ntuple (see read_flatten_root_file.C), this is a
// real DAQ production file and contains SEVERAL trees:
//
//   WCTEReadoutWindows        - the main tree: one entry per readout window
//                                (i.e. one triggered "event"). Contains the
//                                digitized WCTE PMT hits (time+charge) plus
//                                all the beam-line instrumentation (VME TDC/
//                                ADC) info used for particle ID and timing.
//   vme_analysis_scalar_results - single-entry tree of run-averaged PID cuts,
//                                TOF/momentum summary statistics.
//   vme_analysis_run_info     - single-entry tree of run-level metadata
//                                (beam momentum, collimator/jaw positions...).
//   DataQualityMetrics        - single-entry tree of data-quality counters.
//   Configuration              - single-entry tree of PMT masks / good-channel
//                                lists used to build this production file.
//
// Run with:
//   root -l read_data_root_file.C
// or from inside a ROOT session:
//   .x read_data_root_file.C

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <vector>
#include <string>
#include <iostream>

void read_data_root_file()
{
    // ------------------------------------------------------------------
    // 1) Open the input file
    // ------------------------------------------------------------------
    const char* file_path =
        "/eos/user/a/acraplet/WCSim/Analysis/data/WCTE_merged_production_R1670.root";

    TFile* file = TFile::Open(file_path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: could not open file " << file_path << std::endl;
        return;
    }
    std::cout << "Opened " << file_path << std::endl;

    // ------------------------------------------------------------------
    // 2) Single-entry "summary" trees: run configuration, PID cuts, run
    //    metadata, data-quality counters. Each of these trees has exactly
    //    ONE entry (they describe the whole run/file, not per-event data).
    //    Rather than declaring every single branch by hand (there are
    //    ~90+ momentum/TOF summary branches in vme_analysis_scalar_results
    //    alone), TTree::Show(entry) is a quick way to dump every branch
    //    value for a given entry - handy for these small "header-like" trees.
    // ------------------------------------------------------------------
    const char* summary_tree_names[] = {
        "Configuration", "vme_analysis_run_info",
        "vme_analysis_scalar_results", "DataQualityMetrics"
    };
    for (const char* name : summary_tree_names) {
        TTree* t = (TTree*)file->Get(name);
        if (!t) {
            std::cerr << "WARNING: tree '" << name << "' not found" << std::endl;
            continue;
        }
        std::cout << "\n=== " << name << " (single-entry summary tree) ===" << std::endl;
        t->GetEntry(0);
        t->Show(0); // dumps every branch name + value for entry 0
    }

    // ------------------------------------------------------------------
    // 3) Main tree: WCTEReadoutWindows
    //    One entry = one readout window (triggered event). Per-hit
    //    quantities are std::vector<...> branches (one WCTE PMT hit each),
    //    while per-window quantities (run/event numbers, VME TDC times for
    //    the beam-line TOF/ACT counters, etc.) are scalar branches.
    // ------------------------------------------------------------------
    TTree* tree = (TTree*)file->Get("WCTEReadoutWindows");
    if (!tree) {
        std::cerr << "ERROR: could not find TTree 'WCTEReadoutWindows'" << std::endl;
        file->Close();
        return;
    }

    Long64_t n_entries = tree->GetEntries();
    std::cout << "\nTree 'WCTEReadoutWindows' has " << n_entries
               << " entries (readout windows)" << std::endl;

    // --- Event/window identification ---
    Double_t  window_time    = 0; // absolute window start time
    Long_t    start_counter  = 0; // hardware clock counter at window start
    Int_t     run_id         = 0;
    Int_t     sub_run_id     = 0;
    Int_t     spill_counter  = 0; // beam spill number
    Int_t     event_number   = 0;
    Int_t     readout_number = 0;
    Int_t     window_data_quality_mask = 0; // bitmask flagging data-quality issues in this window

    // --- WCTE PMT hits (one entry per digitized hit in this window) ---
    // hit_mpmt_slot_ids + hit_pmt_position_ids identify WHICH PMT was hit
    // (hardware IDs), analogous to hit_mpmt/hit_pmt_in_mpmt in the sim tree
    // - but note there is NO position (x,y,z) branch here: the data has no
    // built-in geometry, positions must come from a separate mapping/DB.
    std::vector<double>* hit_pmt_calibrated_times = nullptr; // calibrated hit time (ns)
    std::vector<double>* hit_pmt_charges          = nullptr; // calibrated hit charge (p.e. or ADC-equivalent)
    std::vector<int>*    hit_mpmt_slot_ids        = nullptr; // mPMT module slot ID
    std::vector<int>*    hit_pmt_position_ids     = nullptr; // PMT position index within its mPMT
    std::vector<int>*    hit_pmt_readout_mask     = nullptr; // per-hit readout/quality mask

    // --- Beam-line "trigger board" hits/waveforms (separate hardware from
    //     the WCTE PMTs: digitizes the TOF/ACT/tagger beam instrumentation) ---
    std::vector<int>*    trigger_board_hit_card_ids           = nullptr;
    std::vector<int>*    trigger_board_hit_channel_ids        = nullptr;
    std::vector<float>*  trigger_board_hit_charges            = nullptr;
    std::vector<double>* trigger_board_hit_times              = nullptr;
    std::vector<int>*    trigger_board_waveform_card_ids      = nullptr;
    std::vector<int>*    trigger_board_waveform_channel_ids   = nullptr;
    std::vector<double>* trigger_board_waveform_times         = nullptr;
    std::vector<std::vector<double>>* trigger_board_waveforms = nullptr; // full digitized waveform samples (large!)

    // --- VME beam-instrumentation scalars: used offline for particle ID ---
    // T0/T1/T4/T5 are time-of-flight (TOF) scintillator counters along the
    // beamline; ACT0-5 are threshold Cherenkov counters (left/right PMT);
    // "mu_tag" is a downstream muon-tagging counter.
    Double_t vme_t0_time = 0, vme_t1_time = 0, vme_t4_time = 0;
    Double_t vme_t0_time_second_hit = 0, vme_t1_time_second_hit = 0, vme_t4_time_second_hit = 0;
    Double_t vme_t5_time = 0;
    Double_t vme_t4_l_time = 0, vme_t4_r_time = 0;
    Double_t vme_t4_l_second_hit = 0, vme_t4_r_second_hit = 0;
    Double_t vme_mu_tag_l_time = 0, vme_mu_tag_r_time = 0;
    Double_t vme_mu_tag_l_charge = 0, vme_mu_tag_r_charge = 0, vme_mu_tag_total = 0;
    Double_t vme_ref0_time = 0, vme_ref1_time = 0;

    // T0/T1 each have 4 readout channels; T4 has 2 - stored per-channel time & charge
    Double_t vme_time_t0[4]   = {0}; // vme_time_t0_0..3
    Double_t vme_time_t1[4]   = {0}; // vme_time_t1_0..3
    Double_t vme_time_t4[2]   = {0}; // vme_time_t4_0..1
    Double_t vme_charge_t0[4] = {0}; // vme_charge_t0_0..3
    Double_t vme_charge_t1[4] = {0}; // vme_charge_t1_0..3
    Double_t vme_charge_t4[2] = {0}; // vme_charge_t4_0..1

    // ACT0-5 threshold Cherenkov counters, each read out on left (l) and
    // right (r) PMTs, each with a charge and a time
    Double_t vme_act_l_charge[6] = {0}; // vme_act{0..5}_l_charge
    Double_t vme_act_r_charge[6] = {0}; // vme_act{0..5}_r_charge
    Double_t vme_act_l_time[6]   = {0}; // vme_act{0..5}_l_time
    Double_t vme_act_r_time[6]   = {0}; // vme_act{0..5}_r_time

    // Precomputed time-of-flight differences between counter pairs, and
    // slew/offset-corrected versions - these are what's actually used for PID
    Double_t vme_tof_t0t1 = 0, vme_tof_t0t4 = 0, vme_tof_t4t1 = 0;
    Double_t vme_tof_t0t5 = 0, vme_tof_t1t5 = 0, vme_tof_t4t5 = 0;
    Double_t vme_tof_corr = 0, vme_tof_t0t4_corr = 0;

    Double_t vme_act_eveto  = 0; // "electron veto" ACT sum/threshold quantity
    Double_t vme_act_tagger = 0; // tagger-side ACT quantity

    Int_t vme_event_id             = 0;
    Int_t vme_spill_number         = 0;
    Int_t vme_evt_quality_bitmask  = 0; // VME-side data-quality flags
    Int_t vme_digi_issues_bitmask  = 0; // digitizer-issue flags

    // --- "T5" tagging system info (separate small trigger/tag detector) ---
    Int_t T5_event_nr    = 0;
    Int_t T5_hit_bitmask = 0;

    // Attach each variable to its branch
    tree->SetBranchAddress("window_time", &window_time);
    tree->SetBranchAddress("start_counter", &start_counter);
    tree->SetBranchAddress("run_id", &run_id);
    tree->SetBranchAddress("sub_run_id", &sub_run_id);
    tree->SetBranchAddress("spill_counter", &spill_counter);
    tree->SetBranchAddress("event_number", &event_number);
    tree->SetBranchAddress("readout_number", &readout_number);
    tree->SetBranchAddress("window_data_quality_mask", &window_data_quality_mask);

    tree->SetBranchAddress("hit_pmt_calibrated_times", &hit_pmt_calibrated_times);
    tree->SetBranchAddress("hit_pmt_charges", &hit_pmt_charges);
    tree->SetBranchAddress("hit_mpmt_slot_ids", &hit_mpmt_slot_ids);
    tree->SetBranchAddress("hit_pmt_position_ids", &hit_pmt_position_ids);
    tree->SetBranchAddress("hit_pmt_readout_mask", &hit_pmt_readout_mask);

    tree->SetBranchAddress("trigger_board_hit_card_ids", &trigger_board_hit_card_ids);
    tree->SetBranchAddress("trigger_board_hit_channel_ids", &trigger_board_hit_channel_ids);
    tree->SetBranchAddress("trigger_board_hit_charges", &trigger_board_hit_charges);
    tree->SetBranchAddress("trigger_board_hit_times", &trigger_board_hit_times);
    tree->SetBranchAddress("trigger_board_waveform_card_ids", &trigger_board_waveform_card_ids);
    tree->SetBranchAddress("trigger_board_waveform_channel_ids", &trigger_board_waveform_channel_ids);
    tree->SetBranchAddress("trigger_board_waveform_times", &trigger_board_waveform_times);
    tree->SetBranchAddress("trigger_board_waveforms", &trigger_board_waveforms);

    tree->SetBranchAddress("vme_t0_time", &vme_t0_time);
    tree->SetBranchAddress("vme_t1_time", &vme_t1_time);
    tree->SetBranchAddress("vme_t4_time", &vme_t4_time);
    tree->SetBranchAddress("vme_t0_time_second_hit", &vme_t0_time_second_hit);
    tree->SetBranchAddress("vme_t1_time_second_hit", &vme_t1_time_second_hit);
    tree->SetBranchAddress("vme_t4_time_second_hit", &vme_t4_time_second_hit);
    tree->SetBranchAddress("vme_t5_time", &vme_t5_time);
    tree->SetBranchAddress("vme_t4_l_time", &vme_t4_l_time);
    tree->SetBranchAddress("vme_t4_r_time", &vme_t4_r_time);
    tree->SetBranchAddress("vme_t4_l_second_hit", &vme_t4_l_second_hit);
    tree->SetBranchAddress("vme_t4_r_second_hit", &vme_t4_r_second_hit);
    tree->SetBranchAddress("vme_mu_tag_l_time", &vme_mu_tag_l_time);
    tree->SetBranchAddress("vme_mu_tag_r_time", &vme_mu_tag_r_time);
    tree->SetBranchAddress("vme_mu_tag_l_charge", &vme_mu_tag_l_charge);
    tree->SetBranchAddress("vme_mu_tag_r_charge", &vme_mu_tag_r_charge);
    tree->SetBranchAddress("vme_mu_tag_total", &vme_mu_tag_total);
    tree->SetBranchAddress("vme_ref0_time", &vme_ref0_time);
    tree->SetBranchAddress("vme_ref1_time", &vme_ref1_time);

    // T0/T1: 4 channels each; T4: 2 channels - filled via a small loop
    for (int i = 0; i < 4; i++) {
        tree->SetBranchAddress(Form("vme_time_t0_%d", i), &vme_time_t0[i]);
        tree->SetBranchAddress(Form("vme_time_t1_%d", i), &vme_time_t1[i]);
        tree->SetBranchAddress(Form("vme_charge_t0_%d", i), &vme_charge_t0[i]);
        tree->SetBranchAddress(Form("vme_charge_t1_%d", i), &vme_charge_t1[i]);
    }
    for (int i = 0; i < 2; i++) {
        tree->SetBranchAddress(Form("vme_time_t4_%d", i), &vme_time_t4[i]);
        tree->SetBranchAddress(Form("vme_charge_t4_%d", i), &vme_charge_t4[i]);
    }

    // ACT0-5, left/right charge & time - filled via a small loop
    for (int i = 0; i < 6; i++) {
        tree->SetBranchAddress(Form("vme_act%d_l_charge", i), &vme_act_l_charge[i]);
        tree->SetBranchAddress(Form("vme_act%d_r_charge", i), &vme_act_r_charge[i]);
        tree->SetBranchAddress(Form("vme_act%d_l_time", i), &vme_act_l_time[i]);
        tree->SetBranchAddress(Form("vme_act%d_r_time", i), &vme_act_r_time[i]);
    }

    tree->SetBranchAddress("vme_tof_t0t1", &vme_tof_t0t1);
    tree->SetBranchAddress("vme_tof_t0t4", &vme_tof_t0t4);
    tree->SetBranchAddress("vme_tof_t4t1", &vme_tof_t4t1);
    tree->SetBranchAddress("vme_tof_t0t5", &vme_tof_t0t5);
    tree->SetBranchAddress("vme_tof_t1t5", &vme_tof_t1t5);
    tree->SetBranchAddress("vme_tof_t4t5", &vme_tof_t4t5);
    tree->SetBranchAddress("vme_tof_corr", &vme_tof_corr);
    tree->SetBranchAddress("vme_tof_t0t4_corr", &vme_tof_t0t4_corr);
    tree->SetBranchAddress("vme_act_eveto", &vme_act_eveto);
    tree->SetBranchAddress("vme_act_tagger", &vme_act_tagger);

    tree->SetBranchAddress("vme_event_id", &vme_event_id);
    tree->SetBranchAddress("vme_spill_number", &vme_spill_number);
    tree->SetBranchAddress("vme_evt_quality_bitmask", &vme_evt_quality_bitmask);
    tree->SetBranchAddress("vme_digi_issues_bitmask", &vme_digi_issues_bitmask);

    tree->SetBranchAddress("T5_event_nr", &T5_event_nr);
    tree->SetBranchAddress("T5_hit_bitmask", &T5_hit_bitmask);

    // ------------------------------------------------------------------
    // 4) Loop over readout windows and access the data
    //
    //    NOTE: this file has ~89k windows and the waveform branch alone is
    //    ~1 GB, so only loop over a handful of entries here as a template -
    //    remove the n_to_print cap for a real analysis.
    // ------------------------------------------------------------------
    Long64_t n_to_print = std::min<Long64_t>(n_entries, 3);
    for (Long64_t i = 0; i < n_to_print; ++i) {

        tree->GetEntry(i);

        std::cout << "\n--- Window " << i << " ---" << std::endl;
        std::cout << "  run_id=" << run_id << " event_number=" << event_number
                   << " spill_counter=" << spill_counter
                   << " n_pmt_hits=" << (hit_pmt_charges ? hit_pmt_charges->size() : 0)
                   << std::endl;
        std::cout << "  window_data_quality_mask=" << window_data_quality_mask << std::endl;
        std::cout << "  TOF t0-t1=" << vme_tof_t0t1 << " ns, t0-t4=" << vme_tof_t0t4 << " ns"
                   << " (corrected: " << vme_tof_corr << " ns)" << std::endl;

        // Example: loop over the WCTE PMT hits of this window
        if (hit_pmt_charges) {
            for (size_t h = 0; h < hit_pmt_charges->size(); ++h) {
                // Uncomment to print every single hit (can be verbose!)
                // std::cout << "    hit " << h
                //           << ": q=" << hit_pmt_charges->at(h)
                //           << " t=" << hit_pmt_calibrated_times->at(h)
                //           << " mpmt_slot=" << hit_mpmt_slot_ids->at(h)
                //           << " pmt_pos=" << hit_pmt_position_ids->at(h) << std::endl;
            }
        }
    }

    // ------------------------------------------------------------------
    // 5) Clean up
    // ------------------------------------------------------------------
    file->Close();
}
