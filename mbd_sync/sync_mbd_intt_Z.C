#include <TFile.h>

#include <TTree.h>

#include <TDirectory.h>

// #include "InttEvent.cc"

#include "mbdtree.C"

#include <iostream>

//R__LOAD_LIBRARY(libInttEvent.so)

void sync_mbd_intt_Z() {

    TDirectory * gDir = gDirectory;

    TFile * f_mbd = TFile::Open("/sphenix/tg/tg01/commissioning/INTT/subsystems/MBD/auau2023_v0/beam_seb18-00020869-0000_mbd.root");
    gDirectory = gDir;
    TTree * t_mbd = (TTree * ) f_mbd -> Get("t");
    cout << " " << t_mbd << endl;
    if (!t_mbd) return;
    mbdtree mbdt(t_mbd);

    string folder_direction = "/sphenix/user/ChengWei/INTT/INTT_commissioning/ZeroField/20869/folder_beam_inttall-00020869-0000_event_base_ana_cluster_full_survey_3.32_excludeR20000_200kEvent_3HotCut_advanced";
    TFile * f_intt = TFile::Open(Form("%s/INTT_zvtx.root", folder_direction.c_str()));
    gDirectory = gDir;
    TTree * t_intt = (TTree * ) f_intt -> Get("tree_Z");
    if (!t_intt) return;

    int intt_eID, intt_N_cluster_outer, intt_N_cluster_inner, intt_N_good;
    double intt_zvtx, intt_zvtxE, intt_rangeL, intt_rangeR, intt_width_density;
    Long64_t intt_bco_full; 

    t_intt -> SetBranchAddress("eID", & intt_eID);
    t_intt -> SetBranchAddress("bco_full", & intt_bco_full);
    t_intt -> SetBranchAddress("nclu_inner", & intt_N_cluster_inner);
    t_intt -> SetBranchAddress("nclu_outer", & intt_N_cluster_outer);
    t_intt -> SetBranchAddress("zvtx", & intt_zvtx);
    t_intt -> SetBranchAddress("zvtxE", & intt_zvtxE);
    t_intt -> SetBranchAddress("rangeL", & intt_rangeL);
    t_intt -> SetBranchAddress("rangeR", & intt_rangeR);
    t_intt -> SetBranchAddress("N_good", & intt_N_good);
    t_intt -> SetBranchAddress("Width_density", & intt_width_density);

    cout << t_mbd -> GetEntries() << " " << t_intt -> GetEntries() << endl;

    TFile * out_file = new TFile(Form("%s/INTT_MBD_zvtx.root",folder_direction.c_str()),"RECREATE");

    int out_eID, out_N_cluster_outer, out_N_cluster_inner, out_N_good;
    double out_zvtx, out_zvtxE, out_rangeL, out_rangeR, out_width_density;
    Long64_t out_bco_full; 
    double out_mbd_bz, out_mbd_bqs, out_mbd_bqn;

    TTree * tree_out =  new TTree ("tree_Z", "INTT Z info.");
    tree_out -> Branch("eID",&out_eID);
    tree_out -> Branch("bco_full",&out_bco_full);
    tree_out -> Branch("nclu_inner",&out_N_cluster_inner);
    tree_out -> Branch("nclu_outer",&out_N_cluster_outer);
    tree_out -> Branch("zvtx",&out_zvtx);
    tree_out -> Branch("zvtxE",&out_zvtxE);
    tree_out -> Branch("rangeL",&out_rangeL);
    tree_out -> Branch("rangeR",&out_rangeR);
    tree_out -> Branch("N_good",&out_N_good);
    tree_out -> Branch("Width_density",&out_width_density);
    tree_out -> Branch("mbd_bz",&out_mbd_bz); // note : mbd branch
    tree_out -> Branch("mbd_bqs",&out_mbd_bqs);
    tree_out -> Branch("mbd_bqn",&out_mbd_bqn);


    gDirectory = gDir;

    TH2F * h_qmbd_nintt = new TH2F("h_qmbd_nintt", "BbcQ vs Intt N", 200, 0, 9000, 200, 0, 4000);
    TH2F * intt_mbd_bco = new TH2F("intt_mbd_bco", "INTT - MBD", 100, 0, 50000, 100, -10, 100000);

    int prev_mbdclk = 0;
    ULong64_t prev_bco = 0;

    bool found_firstevt = false;
    int mbd_evt_offset = 0;
    int intt_evt_offset = 0;
    for (int i = 0; i < t_intt -> GetEntries(); i++) {
        mbdt.LoadTree(i + mbd_evt_offset);
        mbdt.GetEntry(i + mbd_evt_offset);
        t_intt -> GetEntry(i + intt_evt_offset);

        float bbcq = mbdt.bqn + mbdt.bqs;

        unsigned short mbdclk = mbdt.femclk;
        ULong64_t bco = intt_bco_full;
        ULong64_t bco16 = bco & 0xFFFF;

        int mbd_prvdif = (mbdclk - prev_mbdclk) & 0xFFFF;
        ULong64_t intt_prvdif = bco - prev_bco;

        prev_mbdclk = mbdclk;
        prev_bco = bco;

        out_mbd_bqn = mbdt.bqn;
        out_mbd_bqs = mbdt.bqs;
        out_mbd_bz  = mbdt.bz;

        out_eID = intt_eID;
        out_bco_full = intt_bco_full;
        out_N_cluster_inner = intt_N_cluster_inner;
        out_N_cluster_outer = intt_N_cluster_outer;
        out_zvtx = intt_zvtx;
        out_zvtxE = intt_zvtxE;
        out_rangeL = intt_rangeL;
        out_rangeR = intt_rangeR;
        out_N_good = intt_N_good;
        out_width_density = intt_width_density;

        if (intt_N_cluster_inner != -1 && intt_N_cluster_outer != -1){
            if (mbdt.bqs > 100 && mbdt.bqn > 100){
                h_qmbd_nintt -> Fill(intt_N_cluster_inner + intt_N_cluster_outer, bbcq);
            }
            
            intt_mbd_bco -> Fill(i, (mbdclk - bco16) & 0xFFFF);

            tree_out -> Fill();
        }


        

        if ((i % 1000) == 0) {
            cout << i << " " << hex << setw(6) << mbdclk << " " << setw(6) << bco16 << " (mbd-intt)" << setw(6) << ((mbdclk - bco16) & 0xFFFF) <<
                "      (femclk-prev)" << setw(6) << mbd_prvdif << " (bco-prev)" << setw(6) << intt_prvdif << dec << endl;
        }


        t_intt -> GetEntry(i + 1 + intt_evt_offset);
        ULong64_t next_bco16 = (intt_bco_full) & 0xFFFF;
        mbdt.LoadTree(i + 1);
        mbdt.GetEntry(i + 1);
        unsigned short next_mbdclk = mbdt.femclk;
        if (((next_mbdclk - next_bco16) & 0xFFFF) != ((mbdclk - bco16) & 0xFFFF)) intt_evt_offset += 1;
    }

    out_file -> cd();
    tree_out->SetDirectory(out_file);
    tree_out -> Write("", TObject::kOverwrite);
    out_file -> Close();
    
    TFile * froot = new TFile(Form("%s/sync.root",folder_direction.c_str()), "recreate");
    h_qmbd_nintt -> Write();
    intt_mbd_bco -> Write();
    froot -> Close();

    
}