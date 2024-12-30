#ifndef SEMICHIPCLUSTERING_H
#define SEMICHIPCLUSTERING_H

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <cctype> // For isdigit

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <THStack.h>
#include <TCanvas.h> // note : for the combined case
#include <TGraph.h>  // note : for the combined case

#include <TKey.h>
#include <TRandom.h> // note : for the offset
#include <TRandom3.h> // note : for the offset

#include <TColor.h>

#include <TObjArray.h>

#include "../Constants.h"

class SemiChipClustering{
    public:
        SemiChipClustering(
            int process_id_in, // note : 1 or 2
            int runnumber_in, // note : still, (54280 or -1)
            int run_nEvents_in,
            std::string input_directory_in,
            std::string input_file_name_in,
            std::string output_directory_in,

            std::string output_file_name_suffix_in,

            bool BcoFullDiffCut_in,
            bool INTT_vtxZ_QA_in,
            bool ApplyHitQA_in = true,
            bool clone_hit_remove_BCO_tag_in = true,
            std::pair<bool, int> cut_HitBcoDiff_in = {true, 55},
            std::vector<int> adc_conversion_vec_in = {35, 45, 60, 90, 120, 150, 180, 210}
        );

        std::string GetOutputFileName() {return output_filename;}
        void MainProcess();
        void EndRun();

    protected:
        struct inttHitstr{
            int hit_server;
            int hit_felix_ch;
            int hit_chip;
            int hit_channel;
            int hit_adc;
            int hit_bco;
            int hit_bco_diff;
        };

        // note : ------------------------- for constructor --------------------------------
        int process_id;
        int runnumber;
        int run_nEvents;
        std::string input_directory;
        std::string input_file_name;
        std::string output_directory;
        std::string output_file_name_suffix;

        bool BcoFullDiffCut;
        bool INTT_vtxZ_QA;
        bool ApplyHitQA;
        bool clone_hit_remove_BCO_tag;
        std::pair<bool, int> cut_HitBcoDiff;
        std::vector<int> adc_conversion_vec;

        // note : ------------------------- for input root file --------------------------------
        TFile * file_in;
        TTree * tree_in;
        std::string tree_name = "EventTree";
        void PrepareInputRootFile();
        std::map<std::string, int> GetInputTreeBranchesMap(TTree * m_tree_in);

        // note : MBD & centrality relevant
        float MBD_z_vtx;
        bool is_min_bias;
        float MBD_centrality;
        float MBD_south_charge_sum;
        float MBD_north_charge_sum;
        float MBD_charge_sum;
        float MBD_charge_asymm;
        int InttBcoFullDiff_next;

        int NClus;

        // note : trigger tag
        int MBDNSg2;
        int MBDNSg2_vtxZ10cm;
        int MBDNSg2_vtxZ30cm;
        int MBDNSg2_vtxZ60cm;

        double INTTvtxZ;
        double INTTvtxZError;
        double NgroupTrapezoidal;
        double NgroupCoarse;
        double TrapezoidalFitWidth;
        double TrapezoidalFWHM;

        std::vector<unsigned long> *InttRawHit_bco;
        std::vector<unsigned int> *InttRawHit_packetid;
        // std::vector<unsigned int> *InttRawHit_word;
        std::vector<unsigned short> *InttRawHit_fee;
        std::vector<unsigned short> *InttRawHit_channel_id;
        std::vector<unsigned short> *InttRawHit_chip_id;
        std::vector<unsigned short> *InttRawHit_adc;
        std::vector<unsigned short> *InttRawHit_FPHX_BCO;
        // std::vector<unsigned short> *InttRawHit_full_FPHX;
        // std::vector<unsigned short> *InttRawHit_full_ROC;
        // std::vector<unsigned short> *InttRawHit_amplitude;

        // note : ------------------------- for analysis --------------------------------
        struct chip_clu_info{
            std::vector<double> adc_vec;
            std::vector<int> size_vec;
            std::vector<double> edge_l_vec;
            std::vector<double> edge_r_vec;
            std::vector<double> pos_vec;
            int largest_size;
        };
        long long evt_INTT_bco_full;
        std::map<std::string, inttHitstr> evt_inttHits_map;
        TH1D * h1D_ClusSizeCount;

        void EvtCleanUp();
        
        bool DoSemiChipCluster(
            int eID_count,
            int NClus_in,
            std::vector<unsigned long> bco_vec_in,
            std::vector<unsigned int> packetid_vec_in,
            std::vector<unsigned short> fee_vec_in,
            std::vector<unsigned short> channel_id_vec_in,
            std::vector<unsigned short> chip_id_vec_in,
            std::vector<unsigned short> adc_vec_in,
            std::vector<unsigned short> FPHX_BCO_vec_in
        );
        SemiChipClustering::chip_clu_info find_Ngroup(TH1D * hist_in);
        double CoM_cluster_pos(TH1D * hist_in, double edge_l, double edge_r);


        // note : ------------------------- for output root file --------------------------------
        TFile * file_out;
        virtual void PrepareOutPutFileName();
        virtual void PrepareOutPutRootFile();
        std::string output_filename;

        // note : ------------------------- for output hist --------------------------------
        void PrepareHists();

        std::map<std::string, TH1D*> h1D_chip_hit_container_map;

        TH1D * h1D_NSize2Clus;
        TH1D * h1D_NSize1Clus;
        TH1D * h1D_ClusSizeCount_all;
        
        TH2D * h2D_NClus_NSize2Clus;
        TH2D * h2D_NClus_NSize1Clus;
        TH2D * h2D_nInttRawHit_NSize2Clus;
        TH2D * h2D_nInttRawHit_NSize1Clus;

        TH2D * h2D_nChipHit_NSize2Clus;
        TH2D * h2D_nChipHit_NSize1Clus;

        // note : -------------------------------- The constant values ------------------------------
        const int nFelix = 8;
        const int Felix_offset = 3001;
        const int nFelix_channel = 14; 
        const int nChip = 26;
        const int nHitBco = 128;

        const int inner_barrel_0 = 3;
        const int inner_barrel_1 = 4;
        const int outer_barrel_0 = 5;
        const int outer_barrel_1 = 6;

        const int south_typeB_ZID = 1;
        const int south_typeA_ZID = 0;
        const int north_typeA_ZID = 2;
        const int north_typeB_ZID = 3;

        std::vector<int> typeA_sensorZID = {0,2}; // note : sensor Z ID for type A // note -> 1, 0, 2, 3 

        std::pair<double, double> cut_GlobalMBDvtxZ    = Constants::cut_GlobalMBDvtxZ;
        std::pair<double, double> cut_vtxZDiff = Constants::cut_vtxZDiff;
        std::pair<double, double> cut_TrapezoidalFitWidth = Constants::cut_TrapezoidalFitWidth;
        std::pair<double, double> cut_TrapezoidalFWHM = Constants::cut_TrapezoidalFWHM;
        std::pair<double, double> cut_INTTvtxZError = Constants::cut_INTTvtxZError;

        int cut_InttBcoFullDIff_next = Constants::cut_InttBcoFullDIff_next;

        std::vector<double> centrality_edges = Constants::centrality_edges;
        int nCentrality_bin;

        int B0L0_index = Constants::B0L0_index;
        int B0L1_index = Constants::B0L1_index;
        int B1L0_index = Constants::B1L0_index;
        int B1L1_index = Constants::B1L1_index;
        int nLadder_inner = Constants::nLadder_inner;
        int nLadder_outer = Constants::nLadder_outer;

};

#endif