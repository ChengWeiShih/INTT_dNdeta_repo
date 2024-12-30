#include "RestDist.h"

RestDist::RestDist(
    int process_id_in,
    int runnumber_in,
    int run_nEvents_in,
    std::string input_directory_in,
    std::string input_file_name_in,
    std::string output_directory_in,

    std::string output_file_name_suffix_in,

    bool Apply_cut_in,
    bool ApplyVtxZReWeighting_in,
    std::pair<bool, int> ApplyEvtBcoFullDiffCut_in,
    std::pair<bool, std::pair<double,double>> RequireVtxZRange_in,
    std::pair<bool, std::pair<double,double>> isClusQA_in // note : true/false, {ADC, phi_size}
):
    process_id(process_id_in),
    runnumber(runnumber_in),
    run_nEvents(run_nEvents_in),
    input_directory(input_directory_in),
    input_file_name(input_file_name_in),
    output_directory(output_directory_in),
    output_file_name_suffix(output_file_name_suffix_in),
    Apply_cut(Apply_cut_in),
    ApplyVtxZReWeighting(ApplyVtxZReWeighting_in),
    ApplyEvtBcoFullDiffCut(ApplyEvtBcoFullDiffCut_in),
    RequireVtxZRange(RequireVtxZRange_in),
    isClusQA(isClusQA_in)
{
    PrepareOutPutFileName();
    PrepareInputFile();

    run_nEvents = (run_nEvents == -1) ? tree_in->GetEntries() : run_nEvents;
    run_nEvents = (run_nEvents > tree_in->GetEntries()) ? tree_in->GetEntries() : run_nEvents;
    nCentrality_bin = centrality_edges.size() - 1;

    PrepareOutputFile();
    PrepareHist();

}

void RestDist::PrepareOutPutFileName()
{
    std::string job_index = std::to_string( process_id );
    int job_index_len = 5;
    job_index.insert(0, job_index_len - job_index.size(), '0');

    std::string runnumber_str = std::to_string( runnumber );
    if (runnumber != -1){
        int runnumber_str_len = 8;
        runnumber_str.insert(0, runnumber_str_len - runnumber_str.size(), '0');
    }

    if (output_file_name_suffix.size() > 0 && output_file_name_suffix[0] != '_') {
        output_file_name_suffix = "_" + output_file_name_suffix;
    }

    output_filename = "RestDist";
    output_filename = (runnumber != -1) ? "Data_" + output_filename : "MC_" + output_filename;

    output_filename += (Apply_cut) ? "_vtxZQA" : "_NovtxZQA";
    output_filename += (ApplyVtxZReWeighting) ? "_VtxZReWeighting" : "";
    output_filename += (runnumber != -1 && ApplyEvtBcoFullDiffCut.first) ? Form("_EvtBcoFullDiffCut%d", ApplyEvtBcoFullDiffCut.second) : "";

    if (RequireVtxZRange.first){
        std::string vtxZRange;

        double vtxZ_range_l = RequireVtxZRange.second.first;
        double vtxZ_range_r = RequireVtxZRange.second.second;

        vtxZRange += "_vtxZRange";
        vtxZRange += (vtxZ_range_l < 0) ? Form("M%.0fp%d",fabs(vtxZ_range_l),int(fabs(vtxZ_range_l)*10)%10) : Form("%.0fp%d",vtxZ_range_l, int(vtxZ_range_l*10)%10);
        vtxZRange += "to";
        vtxZRange += (vtxZ_range_r < 0) ? Form("M%.0fp%d",fabs(vtxZ_range_r),int(fabs(vtxZ_range_r)*10)%10) : Form("%.0fp%d",vtxZ_range_r, int(vtxZ_range_r*10)%10);

        output_filename += vtxZRange;
    }

    output_filename += (isClusQA.first) ? Form("_ClusQAAdc%.0fPhiSize%.0f",isClusQA.second.first,isClusQA.second.second) : "";

    output_filename += output_file_name_suffix;
    output_filename += (runnumber != -1) ? Form("_%s_%s.root",runnumber_str.c_str(),job_index.c_str()) : Form("_%s.root",job_index.c_str());
}

std::map<std::string, int> RestDist::GetInputTreeBranches(TTree * m_tree_in)
{
    std::map<std::string, int> branch_map;
    TObjArray * branch_list = m_tree_in -> GetListOfBranches();
    for (int i = 0; i < branch_list -> GetEntries(); i++)
    {
        TBranch * branch = dynamic_cast<TBranch*>(branch_list->At(i));
        branch_map[branch -> GetName()] = 1;
    }
    return branch_map;
}

void RestDist::PrepareInputFile()
{
    file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_file_name.c_str()));
    tree_in = (TTree*)file_in->Get(tree_name.c_str());

    std::map<std::string, int> branch_map = GetInputTreeBranches(tree_in);

    tree_in -> SetBranchStatus("*", 0);

    tree_in -> SetBranchStatus("is_min_bias", 1);
    tree_in -> SetBranchStatus("MBD_centrality", 1);
    tree_in -> SetBranchStatus("MBD_z_vtx", 1);
    tree_in -> SetBranchStatus("MBD_charge_sum", 1);
    tree_in -> SetBranchStatus("MBD_charge_asymm", 1);
    tree_in -> SetBranchStatus("MBD_south_charge_sum", 1);
    tree_in -> SetBranchStatus("MBD_north_charge_sum", 1);

    tree_in -> SetBranchStatus("INTTvtxZ", 1);
    tree_in -> SetBranchStatus("INTTvtxZError", 1);
    tree_in -> SetBranchStatus("NgroupTrapezoidal", 1);
    tree_in -> SetBranchStatus("NgroupCoarse", 1);
    tree_in -> SetBranchStatus("TrapezoidalFitWidth", 1);
    tree_in -> SetBranchStatus("TrapezoidalFWHM", 1);

    tree_in -> SetBranchStatus("NClus", 1);
    tree_in -> SetBranchStatus("NClus_Layer1", 1);
    
    // note : for data
    if (branch_map.find("MBDNSg2") != branch_map.end()) {
        tree_in -> SetBranchStatus("MBDNSg2", 1);
        m_withTrig = true;
    }
    if (branch_map.find("MBDNSg2_vtxZ10cm") != branch_map.end()) {tree_in -> SetBranchStatus("MBDNSg2_vtxZ10cm", 1);}
    if (branch_map.find("MBDNSg2_vtxZ30cm") != branch_map.end()) {tree_in -> SetBranchStatus("MBDNSg2_vtxZ30cm", 1);}
    
    if (branch_map.find("InttBcoFullDiff_next") != branch_map.end()) {tree_in -> SetBranchStatus("InttBcoFullDiff_next", 1); }
    
    // note :for MC
    if (branch_map.find("NTruthVtx") != branch_map.end()) {tree_in -> SetBranchStatus("NTruthVtx", 1);}
    if (branch_map.find("TruthPV_trig_z") != branch_map.end()) {tree_in -> SetBranchStatus("TruthPV_trig_z", 1);}


    tree_in -> SetBranchStatus("ClusX", 1);
    tree_in -> SetBranchStatus("ClusY", 1);
    tree_in -> SetBranchStatus("ClusZ", 1);

    tree_in -> SetBranchStatus("ClusAdc", 1);
    tree_in -> SetBranchStatus("ClusPhiSize", 1);
    tree_in -> SetBranchStatus("ClusZSize", 1);

    tree_in -> SetBranchStatus("ClusLayer", 1);
    tree_in -> SetBranchStatus("ClusLadderZId", 1);
    tree_in -> SetBranchStatus("ClusLadderPhiId", 1);

    tree_in -> SetBranchStatus("ClusEta_INTTz", 1);
    tree_in -> SetBranchStatus("ClusEta_MBDz", 1);
    tree_in -> SetBranchStatus("ClusPhi_AvgPV", 1);

    // note : for MC
    if (branch_map.find("ClusEta_TrueXYZ") != branch_map.end()) {tree_in -> SetBranchStatus("ClusEta_TrueXYZ", 1);}
    if (branch_map.find("ClusPhi_TrueXY") != branch_map.end()) {tree_in -> SetBranchStatus("ClusPhi_TrueXY", 1);}

    // Division : ---SetBranchAddress-----------------------------------------------------------------------------------------------


    ClusX = 0;
    ClusY = 0;
    ClusZ = 0;
    ClusAdc = 0;
    ClusPhiSize = 0;
    ClusZSize = 0;
    ClusLayer = 0;
    ClusLadderZId = 0;
    ClusLadderPhiId = 0;
    ClusEta_INTTz = 0;
    ClusEta_MBDz = 0;
    ClusPhi_AvgPV = 0;
    ClusEta_TrueXYZ = 0;
    ClusPhi_TrueXY = 0;

    tree_in -> SetBranchAddress("is_min_bias", &is_min_bias);
    tree_in -> SetBranchAddress("MBD_centrality", &MBD_centrality);
    tree_in -> SetBranchAddress("MBD_z_vtx", &MBD_z_vtx);
    tree_in -> SetBranchAddress("MBD_charge_sum", &MBD_charge_sum);
    tree_in -> SetBranchAddress("MBD_charge_asymm", &MBD_charge_asymm);
    tree_in -> SetBranchAddress("MBD_south_charge_sum", &MBD_south_charge_sum);
    tree_in -> SetBranchAddress("MBD_north_charge_sum", &MBD_north_charge_sum);

    tree_in -> SetBranchAddress("INTTvtxZ", &INTTvtxZ);
    tree_in -> SetBranchAddress("INTTvtxZError", &INTTvtxZError);
    tree_in -> SetBranchAddress("NgroupTrapezoidal", &NgroupTrapezoidal);
    tree_in -> SetBranchAddress("NgroupCoarse", &NgroupCoarse);
    tree_in -> SetBranchAddress("TrapezoidalFitWidth", &TrapezoidalFitWidth);
    tree_in -> SetBranchAddress("TrapezoidalFWHM", &TrapezoidalFWHM);

    tree_in -> SetBranchAddress("NClus", &NClus);
    tree_in -> SetBranchAddress("NClus_Layer1", &NClus_Layer1);

    tree_in -> SetBranchAddress("ClusX", &ClusX);
    tree_in -> SetBranchAddress("ClusY", &ClusY);
    tree_in -> SetBranchAddress("ClusZ", &ClusZ);

    tree_in -> SetBranchAddress("ClusAdc", &ClusAdc);
    tree_in -> SetBranchAddress("ClusPhiSize", &ClusPhiSize);
    tree_in -> SetBranchAddress("ClusZSize", &ClusZSize);
    
    tree_in -> SetBranchAddress("ClusLayer", &ClusLayer);
    tree_in -> SetBranchAddress("ClusLadderZId", &ClusLadderZId);
    tree_in -> SetBranchAddress("ClusLadderPhiId", &ClusLadderPhiId);

    tree_in -> SetBranchAddress("ClusEta_INTTz", &ClusEta_INTTz);
    tree_in -> SetBranchAddress("ClusEta_MBDz", &ClusEta_MBDz);
    tree_in -> SetBranchAddress("ClusPhi_AvgPV", &ClusPhi_AvgPV);

    // note : MC
    if (branch_map.find("ClusEta_TrueXYZ") != branch_map.end()) {tree_in -> SetBranchAddress("ClusEta_TrueXYZ", &ClusEta_TrueXYZ);}
    if (branch_map.find("ClusPhi_TrueXY") != branch_map.end()) {tree_in -> SetBranchAddress("ClusPhi_TrueXY", &ClusPhi_TrueXY);}

    // note : for data
    if (branch_map.find("MBDNSg2") != branch_map.end()) {tree_in -> SetBranchAddress("MBDNSg2", &MBDNSg2);}
    if (branch_map.find("MBDNSg2_vtxZ10cm") != branch_map.end()) {tree_in -> SetBranchAddress("MBDNSg2_vtxZ10cm", &MBDNSg2_vtxZ10cm);}
    if (branch_map.find("MBDNSg2_vtxZ30cm") != branch_map.end()) {tree_in -> SetBranchAddress("MBDNSg2_vtxZ30cm", &MBDNSg2_vtxZ30cm);}

    if (branch_map.find("InttBcoFullDiff_next") != branch_map.end()) {tree_in -> SetBranchAddress("InttBcoFullDiff_next", &InttBcoFullDiff_next); }

    // note : for MC
    if (branch_map.find("NTruthVtx") != branch_map.end()) {tree_in -> SetBranchAddress("NTruthVtx", &NTruthVtx);}
    if (branch_map.find("TruthPV_trig_z") != branch_map.end()) {tree_in -> SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);}
}

void RestDist::PrepareOutputFile()
{
    file_out = new TFile(Form("%s/%s", output_directory.c_str(), output_filename.c_str()), "RECREATE");
}

void RestDist::PrepareHist()
{
    h1D_centrality_bin = new TH1D("h1D_centrality_bin","h1D_centrality_bin;Centrality [%];Entries",nCentrality_bin,&centrality_edges[0]); // note : the 0-5%
    
    h1D_map.clear();
    
    h1D_map.insert(
        std::make_pair(
            "h1D_PairDeltaEta",
            new TH1D("h1D_PairDeltaEta","h1D_PairDeltaEta;Pair #Delta#eta;Entries",nDeltaEtaBin,DeltaEtaEdge_min,DeltaEtaEdge_max)
        )
    );

    h1D_map.insert(
        std::make_pair(
            "h1D_PairDeltaPhi",
            new TH1D("h1D_PairDeltaPhi","h1D_PairDeltaPhi;Pair #Delta#phi;Entries",nDeltaPhiBin,DeltaPhiEdge_min,DeltaPhiEdge_max)
        )
    );

    h1D_map.insert(std::make_pair("h1D_ClusPhi", new TH1D("h1D_ClusPhi", "h1D_ClusPhi;Cluster #phi (w.r.t avg. PV);Entries (/0.05)",140,-3.5,3.5)));

    h1D_map.insert(std::make_pair("h1D_Cluster_phi_size", new TH1D("h1D_Cluster_phi_size", "h1D_Cluster_phi_size;Cluster #phi size;Entries (/1)",128,0,128)));
    h1D_map.insert(std::make_pair("h1D_Cluster_adc_size1", new TH1D("h1D_Cluster_adc_size1", "h1D_Cluster_adc_size1;Cluster Adc (PhiSize=1);Entries (/7)",25,0,250)));
    h1D_map.insert(std::make_pair("h1D_Cluster_adc", new TH1D("h1D_Cluster_adc", "h1D_Cluster_adc;Cluster Adc;Entries (/7)",100,0,700)));
    h1D_map.insert(std::make_pair("h1D_Cluster_Z", new TH1D("h1D_Cluster_Z", "h1D_Cluster_Z;Cluster Z [cm];Entries (/0.6)",100,-30,30)));
    h1D_map.insert(std::make_pair("h1D_Cluster_Z_typeA", new TH1D("h1D_Cluster_Z_typeA", "h1D_Cluster_Z_typeA;Cluster Z (Type A) [cm];Entries (/0.6)",100,-30,30)));
    h1D_map.insert(std::make_pair("h1D_Cluster_Z_typeB", new TH1D("h1D_Cluster_Z_typeB", "h1D_Cluster_Z_typeB;Cluster Z (Type B) [cm];Entries (/0.6)",100,-30,30)));


    h1D_map.insert(std::make_pair("h1D_ClusEta_INTTz", new TH1D("h1D_ClusEta_INTTz", "h1D_ClusEta_INTTz;Cluster #eta (INTTz);Entries (/0.1)",60,-3,3)));
    h1D_map.insert(std::make_pair("h1D_ClusEta_INTTz_highNClus", new TH1D("h1D_ClusEta_INTTz_highNClus", "h1D_ClusEta_INTTz_highNClus;Cluster #eta (INTTz);Entries (/0.1)",60,-3,3)));
    
    h1D_map.insert(std::make_pair("h1D_ClusEta_INTTz_typeA", new TH1D("h1D_ClusEta_INTTz_typeA", "h1D_ClusEta_INTTz_typeA;Cluster #eta (INTTz, typeA);Entries (/0.1)",60,-3,3)));
    h1D_map.insert(std::make_pair("h1D_ClusEta_INTTz_typeA_highNClus", new TH1D("h1D_ClusEta_INTTz_typeA_highNClus", "h1D_ClusEta_INTTz_typeA_highNClus;Cluster #eta (INTTz, typeA);Entries (/0.1)",60,-3,3)));
    
    h1D_map.insert(std::make_pair("h1D_ClusEta_INTTz_typeB", new TH1D("h1D_ClusEta_INTTz_typeB", "h1D_ClusEta_INTTz_typeB;Cluster #eta (INTTz, typeB);Entries (/0.1)",60,-3,3)));
    h1D_map.insert(std::make_pair("h1D_ClusEta_INTTz_typeB_highNClus", new TH1D("h1D_ClusEta_INTTz_typeB_highNClus", "h1D_ClusEta_INTTz_typeB_highNClus;Cluster #eta (INTTz, typeB);Entries (/0.1)",60,-3,3)));
    
    h1D_map.insert(std::make_pair("h1D_ClusEta_MBDz", new TH1D("h1D_ClusEta_MBDz", "h1D_ClusEta_MBDz;Cluster #eta (MBDz);Entries (/0.1)",60,-3,3)));
    h1D_map.insert(std::make_pair("h1D_ClusEta_MBDz_highNClus", new TH1D("h1D_ClusEta_MBDz_highNClus", "h1D_ClusEta_MBDz_highNClus;Cluster #eta (MBDz);Entries (/0.1)",60,-3,3)));

    h1D_map.insert(std::make_pair("h1D_ClusEta_MBDz_typeA", new TH1D("h1D_ClusEta_MBDz_typeA", "h1D_ClusEta_MBDz_typeA;Cluster #eta (MBDz, typeA);Entries (/0.1)",60,-3,3)));
    h1D_map.insert(std::make_pair("h1D_ClusEta_MBDz_typeA_highNClus", new TH1D("h1D_ClusEta_MBDz_typeA_highNClus", "h1D_ClusEta_MBDz_typeA_highNClus;Cluster #eta (MBDz, typeA);Entries (/0.1)",60,-3,3)));
    
    h1D_map.insert(std::make_pair("h1D_ClusEta_MBDz_typeB", new TH1D("h1D_ClusEta_MBDz_typeB", "h1D_ClusEta_MBDz_typeB;Cluster #eta (MBDz, typeB);Entries (/0.1)",60,-3,3)));
    h1D_map.insert(std::make_pair("h1D_ClusEta_MBDz_typeB_highNClus", new TH1D("h1D_ClusEta_MBDz_typeB_highNClus", "h1D_ClusEta_MBDz_typeB_highNClus;Cluster #eta (MBDz, typeB);Entries (/0.1)",60,-3,3)));

    if (runnumber == -1){
        h1D_map.insert(std::make_pair("h1D_ClusPhi_TrueXY", new TH1D("h1D_ClusPhi_TrueXY", "h1D_ClusPhi_TrueXY;Cluster #phi (w.r.t TrueXY);Entries (/0.05)",140,-3.5,3.5)));

        h1D_map.insert(std::make_pair("h1D_ClusEta_TrueXYZ", new TH1D("h1D_ClusEta_TrueXYZ", "h1D_ClusEta_TrueXYZ;Cluster #eta (TrueXYZ);Entries (/0.1)",60,-3,3)));
        h1D_map.insert(std::make_pair("h1D_ClusEta_TrueXYZ_highNClus", new TH1D("h1D_ClusEta_TrueXYZ_highNClus", "h1D_ClusEta_TrueXYZ_highNClus;Cluster #eta (TrueXYZ);Entries (/0.1)",60,-3,3)));

        h1D_map.insert(std::make_pair("h1D_ClusEta_TrueXYZ_typeA", new TH1D("h1D_ClusEta_TrueXYZ_typeA", "h1D_ClusEta_TrueXYZ_typeA;Cluster #eta (TrueXYZ, typeA);Entries (/0.1)",60,-3,3)));
        h1D_map.insert(std::make_pair("h1D_ClusEta_TrueXYZ_typeA_highNClus", new TH1D("h1D_ClusEta_TrueXYZ_typeA_highNClus", "h1D_ClusEta_TrueXYZ_typeA_highNClus;Cluster #eta (TrueXYZ, typeA);Entries (/0.1)",60,-3,3)));
        
        h1D_map.insert(std::make_pair("h1D_ClusEta_TrueXYZ_typeB", new TH1D("h1D_ClusEta_TrueXYZ_typeB", "h1D_ClusEta_TrueXYZ_typeB;Cluster #eta (TrueXYZ, typeB);Entries (/0.1)",60,-3,3)));
        h1D_map.insert(std::make_pair("h1D_ClusEta_TrueXYZ_typeB_highNClus", new TH1D("h1D_ClusEta_TrueXYZ_typeB_highNClus", "h1D_ClusEta_TrueXYZ_typeB_highNClus;Cluster #eta (TrueXYZ, typeB);Entries (/0.1)",60,-3,3)));

    }

    
    h1D_map.insert(std::make_pair("h1D_NClus", new TH1D("h1D_NClus", "h1D_NClus;NClus;Entries (/80)",100,0,10000)));
    h1D_map.insert(std::make_pair("h1D_NClus_inner", new TH1D("h1D_NClus_inner", "h1D_NClus_inner;NClus (inner);Entries (/45)",100,0,4500)));
    h1D_map.insert(std::make_pair("h1D_NClus_outer", new TH1D("h1D_NClus_outer", "h1D_NClus_outer;NClus (outer);Entries (/45)",100,0,4500)));
    
    h1D_map.insert(std::make_pair("h1D_inner_typeA_NClus", new TH1D("h1D_inner_typeA_NClus", "h1D_inner_typeA_NClus;NClus (inner, type A);Entries (/32)",100,0,3200)));
    h1D_map.insert(std::make_pair("h1D_outer_typeA_NClus", new TH1D("h1D_outer_typeA_NClus", "h1D_outer_typeA_NClus;NClus (outer, type A);Entries (/32)",100,0,3200)));
    h1D_map.insert(std::make_pair("h1D_inner_typeB_NClus", new TH1D("h1D_inner_typeB_NClus", "h1D_inner_typeB_NClus;NClus (inner, type B);Entries (/32)",100,0,3200)));
    h1D_map.insert(std::make_pair("h1D_outer_typeB_NClus", new TH1D("h1D_outer_typeB_NClus", "h1D_outer_typeB_NClus;NClus (outer, type B);Entries (/32)",100,0,3200)));

    h1D_map.insert(
        std::make_pair(
            "h1D_confirm_INTTvtxZ_Inclusive70",
            new TH1D("h1D_confirm_INTTvtxZ_Inclusive70","h1D_confirm_INTTvtxZ_Inclusive70;INTT vtx Z [cm];Entries (/0.6)",60,-60,60)
        )
    );

    for (auto &hist : h1D_map)
    {
        std::string YTitle = hist.second -> GetYaxis() -> GetTitle();
        std::string XTitle = hist.second -> GetXaxis() -> GetTitle();

        std::string YTitle_post;

        if (YTitle.find("Entries") != std::string::npos) // note : found the (Entries)
        {
            YTitle_post = Form("Entries  (/%.2f)", hist.second -> GetBinWidth(1));
            hist.second -> GetYaxis() -> SetTitle(YTitle_post.c_str());
        }
    }

    // Division : ---h2D-----------------------------------------------------------------------------------------------
    h2D_map.clear();

    h2D_map.insert(
        std::make_pair(
            "h2D_confirm_InttNClus_MbdChargeSum",
            new TH2D("h2D_confirm_InttNClus_MbdChargeSum","h2D_confirm_InttNClus_MbdChargeSum;INTT NClus;MBD charge sum", 200, 0, 10000, 200, 0, 2500)
        )
    );

    h2D_map.insert(
        std::make_pair(
            "h2D_confirm_InttInnerOuterNClus",
            new TH2D("h2D_confirm_InttInnerOuterNClus","h2D_confirm_InttInnerOuterNClus;NClus (inner);NClus (outer)", 200, 0, 5500, 200, 0, 5500)
        )
    );

    h2D_map.insert(std::make_pair("h2D_ClusPhiSize_ClusAdc", new TH2D("h2D_ClusPhiSize_ClusAdc","h2D_ClusPhiSize_ClusAdc;Cluster #phi size;Cluster Adc",128,0,128,100,0,10000)));
    h2D_map.insert(std::make_pair("h2D_ClusXY", new TH2D("h2D_ClusXY","h2D_ClusXY;X [cm];Y [cm]",300,-15,15,300,-15,15)));

    // note : INTTz
    h2D_map.insert(std::make_pair("h2D_INTT_ClusEta_ClusAdc", new TH2D("h2D_INTT_ClusEta_ClusAdc","h2D_INTT_ClusEta_ClusAdc;Cluster #eta (INTTz);Cluster Adc",60,-3,3,50,0,500)));
    h2D_map.insert(std::make_pair("h2D_INTT_ClusEta_ClusAdc_highNClus", new TH2D("h2D_INTT_ClusEta_ClusAdc_highNClus","h2D_INTT_ClusEta_ClusAdc_highNClus;Cluster #eta (INTTz);Cluster Adc",60,-3,3,50,0,500)));

    // h2D_map.insert(std::make_pair("h2D_INTT_MBD_multiplicity", new TH2D("h2D_INTT_MBD_multiplicity","h2D_INTT_MBD_multiplicity;INTT NClus;MBD charge sum",200,0,10000,200,0,5500)));
    // h2D_map.insert(std::make_pair("h2D_inner_outer_INTT_multiplicity", new TH2D("h2D_inner_outer_INTT_multiplicity","h2D_inner_outer_INTT_multiplicity;NClus (inner);NClus (outer)",200,0,4000,200,0,4000)));
    h2D_map.insert(std::make_pair("h2D_ClusZ_LadderZID", new TH2D("h2D_ClusZ_LadderZID","h2D_ClusZ_LadderZID;Cluster Z [cm]; Clus Sensor ZId",100,-30,30,5,0,5)));

    h2D_map.insert(std::make_pair("h2D_sensor_HitMap", new TH2D("h2D_sensor_HitMap","h2D_sensor_HitMap;LayerID * 10 + sensor ZId;Ladder Phi Id",40,30,70,18,0,18)));
    h2D_map.insert(std::make_pair("h2D_sensor_AdcMap", new TH2D("h2D_sensor_AdcMap","h2D_sensor_AdcMap;LayerID * 10 + sensor ZId;Ladder Phi Id",40,30,70,18,0,18)));
    h2D_map.insert(std::make_pair("h2D_sensor_AvgAdcMap", new TH2D("h2D_sensor_AvgAdcMap","h2D_sensor_AvgAdcMap;LayerID * 10 + sensor ZId;Ladder Phi Id",40,30,70,18,0,18)));

    h2D_map.insert(std::make_pair("h2D_ClusXY_LargePhiSize", new TH2D("h2D_ClusXY_LargePhiSize","h2D_ClusXY_LargePhiSize;X [cm];Y [cm]",300,-15,15,300,-15,15)));    
    h2D_map.insert(std::make_pair("h2D_ClusXY_PhiSize43_46", new TH2D("h2D_ClusXY_PhiSize43_46","h2D_ClusXY_PhiSize43_46;X [cm];Y [cm]",300,-15,15,300,-15,15)));
    h2D_map.insert(std::make_pair("h2D_ClusXY_LargerPhiSize49", new TH2D("h2D_ClusXY_LargerPhiSize49","h2D_ClusXY_LargerPhiSize49;X [cm];Y [cm]",300,-15,15,300,-15,15)));    

    h2D_map.insert(std::make_pair("h2D_MBDcharge_south_north",   new TH2D("h2D_MBDcharge_south_north","h2D_MBDcharge_south_north;MBD charge sum (north);MBD charge sum (south)",100,0,1500,100,0,1500)));
    h2D_map.insert(std::make_pair("h2D_MBDvtxZ_MBDchargeSum",   new TH2D("h2D_MBDvtxZ_MBDchargeSum","h2D_MBDvtxZ_MBDchargeSum;MBD vtxZ [cm];MBD charge sum (all)",120,-60,60,100,0,3000)));
    h2D_map.insert(std::make_pair("h2D_MBDvtxZ_MBDchargeSouth", new TH2D("h2D_MBDvtxZ_MBDchargeSouth","h2D_MBDvtxZ_MBDchargeSouth;MBD vtxZ [cm];MBD charge sum (south)",120,-60,60,100,0,1500)));
    h2D_map.insert(std::make_pair("h2D_MBDvtxZ_MBDchargeNorth", new TH2D("h2D_MBDvtxZ_MBDchargeNorth","h2D_MBDvtxZ_MBDchargeNorth;MBD vtxZ [cm];MBD charge sum (north)",120,-60,60,100,0,1500)));
    h2D_map.insert(std::make_pair("h2D_MBD_vtxZ_ChargeAssy", new TH2D("h2D_MBD_vtxZ_ChargeAssy","h2D_MBD_vtxZ_ChargeAssy;MBD vtxZ [cm];MBD charge asymmetry",120,-60,60,100,-1.5,1.5)));

    h2D_map.insert(std::make_pair("h2D_NClus_ClusPhiSize", new TH2D("h2D_NClus_ClusPhiSize","h2D_NClus_ClusPhiSize;NClus;Cluster #phi size",100,0,10000,120,0,120)));
}

void RestDist::EvtCleanUp()
{
    evt_sPH_inner_nocolumn_vec.clear();
    evt_sPH_outer_nocolumn_vec.clear();

    NClus_inner_typeA = 0;
    NClus_inner_typeB = 0;
    NClus_outer_typeA = 0;
    NClus_outer_typeB = 0;
}

void RestDist::PrepareClusterVec()
{
    for (int clu_i = 0; clu_i < ClusX -> size(); clu_i++)
    {
        RestDist::clu_info this_clu;

        this_clu.adc = ClusAdc -> at(clu_i);
        this_clu.phi_size = ClusPhiSize -> at(clu_i);
        this_clu.sensorZID = ClusLadderZId -> at(clu_i);
        this_clu.ladderPhiID = ClusLadderPhiId -> at(clu_i);
        this_clu.layerID = ClusLayer -> at(clu_i);

        this_clu.index = clu_i;

        this_clu.x = ClusX -> at(clu_i);
        this_clu.y = ClusY -> at(clu_i);
        this_clu.z = ClusZ -> at(clu_i);

        this_clu.eta_INTTz = ClusEta_INTTz -> at(clu_i);
        this_clu.eta_MBDz  = ClusEta_MBDz  -> at(clu_i);
        this_clu.phi_avgXY = ClusPhi_AvgPV -> at(clu_i);

        if (runnumber == -1){
            this_clu.eta_TrueXYZ = ClusEta_TrueXYZ -> at(clu_i);
            this_clu.phi_TrueXY  = ClusPhi_TrueXY  -> at(clu_i);
        }

        if (isClusQA.first && this_clu.adc <= isClusQA.second.first) {continue;} // note : adc
        if (isClusQA.first && this_clu.phi_size > isClusQA.second.second) {continue;} // note : phi size


        if (this_clu.layerID == 3 || this_clu.layerID == 4){ // note : inner
            if (this_clu.sensorZID == typeA_sensorZID[0] || this_clu.sensorZID == typeA_sensorZID[1]){ // note : type A
                NClus_inner_typeA++;
            }
            else { // note : type B
                NClus_inner_typeB++;
            }
        }
        else { // note : outer 
            if (this_clu.sensorZID == typeA_sensorZID[0] || this_clu.sensorZID == typeA_sensorZID[1]){ // note : type A
                NClus_outer_typeA++;
            }
            else { // note : type B
                NClus_outer_typeB++;
            }
        }

        std::vector<RestDist::clu_info>* p_evt_sPH_nocolumn_vec =
        (this_clu.layerID == 3 || this_clu.layerID == 4) ? (&evt_sPH_inner_nocolumn_vec) : (&evt_sPH_outer_nocolumn_vec);
        p_evt_sPH_nocolumn_vec -> push_back(this_clu);
    }
}

void RestDist::PrepareEvent()
{
    for (int i = 0; i < run_nEvents; i++)
    {
        tree_in -> GetEntry(i);

        EvtCleanUp();

        if (i % 10 == 0) {std::cout << "Processing event " << i << std::endl;}

        // =======================================================================================================================================================
        // note : hard cut

        // note : for data
        if (runnumber != -1 && ApplyEvtBcoFullDiffCut.first && InttBcoFullDiff_next <= ApplyEvtBcoFullDiffCut.second) {continue;}
        if (runnumber != -1 && MBDNSg2 != 1) {continue;} // todo: assume MC no trigger

        // note : for MC
        if (runnumber == -1 && NTruthVtx != 1) {continue;}

        // note : both data and MC
        if (MBD_z_vtx != MBD_z_vtx) {continue;}
        if (MBD_centrality != MBD_centrality) {continue;}
        if (MBD_centrality < 0 || MBD_centrality > 1) {continue;}
        if (INTTvtxZ != INTTvtxZ) {continue;}
        if (MBD_z_vtx < cut_GlobalMBDvtxZ.first || MBD_z_vtx > cut_GlobalMBDvtxZ.second) {continue;} // todo: the hard cut 60 cm 

        // =======================================================================================================================================================
        // note : optional cut
        if (Apply_cut && (MBD_z_vtx - INTTvtxZ < cut_vtxZDiff.first || MBD_z_vtx - INTTvtxZ > cut_vtxZDiff.second) ) {continue;}
        if (Apply_cut && (TrapezoidalFitWidth < cut_TrapezoidalFitWidth.first || TrapezoidalFitWidth > cut_TrapezoidalFitWidth.second)){continue;}
        if (Apply_cut && (TrapezoidalFWHM < cut_TrapezoidalFWHM.first || TrapezoidalFWHM > cut_TrapezoidalFWHM.second)){continue;}
        if (Apply_cut && (INTTvtxZError < cut_INTTvtxZError.first || INTTvtxZError > cut_INTTvtxZError.second)){continue;}

        // =======================================================================================================================================================
        if (ApplyVtxZReWeighting && runnumber != -1){
            std::cout<<"Should not have the vtxZ weighting from the data"<<std::endl;
            exit(1);
        }

        double INTTvtxZWeighting;
        if (ApplyVtxZReWeighting && h1D_INTT_vtxZ_reweighting != nullptr){
            INTTvtxZWeighting = h1D_INTT_vtxZ_reweighting -> GetBinContent(h1D_INTT_vtxZ_reweighting -> FindBin(INTTvtxZ));
        }
        else if (ApplyVtxZReWeighting && h1D_INTT_vtxZ_reweighting == nullptr){
            std::cout << "ApplyVtxZReWeighting is true, but h1D_INTT_vtxZ_reweighting is nullptr" << std::endl;
            exit(1);
        }
        else {
            INTTvtxZWeighting = 1.0;
        }
        
        int Mbin = h1D_centrality_bin -> Fill(MBD_centrality, INTTvtxZWeighting);
        Mbin = (Mbin == -1) ? -1 : Mbin - 1;
        if (Mbin == -1) {
            std::cout << "Mbin == -1, MBD_centrality = " << MBD_centrality << std::endl;
            continue;
        }

        // =======================================================================================================================================================
        // todo: Only use the INTT for the analysis vtxZ range
        if (RequireVtxZRange.first && (INTTvtxZ < RequireVtxZRange.second.first || INTTvtxZ > RequireVtxZRange.second.second)) {continue;}
        // if (RequireVtxZRange.first && (MBD_z_vtx < RequireVtxZRange.second.first || MBD_z_vtx > RequireVtxZRange.second.second)) {continue;}
         
        PrepareClusterVec();

        if (Mbin <= Semi_inclusive_bin){
            h1D_map["h1D_confirm_INTTvtxZ_Inclusive70"] -> Fill(INTTvtxZ, INTTvtxZWeighting);
        }

        int total_NClus = NClus_inner_typeA + NClus_inner_typeB + NClus_outer_typeA + NClus_outer_typeB;

        h1D_map["h1D_NClus"] -> Fill(total_NClus, INTTvtxZWeighting);
        h1D_map["h1D_NClus_inner"] -> Fill(NClus_inner_typeA + NClus_inner_typeB, INTTvtxZWeighting);
        h1D_map["h1D_NClus_outer"] -> Fill(NClus_outer_typeA + NClus_outer_typeB, INTTvtxZWeighting);
        h1D_map["h1D_inner_typeA_NClus"] -> Fill(NClus_inner_typeA, INTTvtxZWeighting);
        h1D_map["h1D_inner_typeB_NClus"] -> Fill(NClus_inner_typeB, INTTvtxZWeighting);
        h1D_map["h1D_outer_typeA_NClus"] -> Fill(NClus_outer_typeA, INTTvtxZWeighting);
        h1D_map["h1D_outer_typeB_NClus"] -> Fill(NClus_outer_typeB, INTTvtxZWeighting);

        // note : h2D
        h2D_map["h2D_confirm_InttNClus_MbdChargeSum"] -> Fill(total_NClus, MBD_charge_sum, INTTvtxZWeighting);
        h2D_map["h2D_confirm_InttInnerOuterNClus"] -> Fill(NClus_inner_typeA + NClus_inner_typeB, NClus_outer_typeA + NClus_outer_typeB, INTTvtxZWeighting);
        h2D_map["h2D_MBD_vtxZ_ChargeAssy"] -> Fill(MBD_z_vtx, MBD_charge_asymm, INTTvtxZWeighting);
        h2D_map["h2D_MBDcharge_south_north"] -> Fill(MBD_south_charge_sum, MBD_north_charge_sum, INTTvtxZWeighting);
        h2D_map["h2D_MBDvtxZ_MBDchargeSum"] -> Fill(MBD_z_vtx, MBD_charge_sum, INTTvtxZWeighting);
        h2D_map["h2D_MBDvtxZ_MBDchargeSouth"] -> Fill(MBD_z_vtx, MBD_south_charge_sum, INTTvtxZWeighting);
        h2D_map["h2D_MBDvtxZ_MBDchargeNorth"] -> Fill(MBD_z_vtx, MBD_north_charge_sum, INTTvtxZWeighting);

        for (int layer_index = 0; layer_index < 2; layer_index++){ // note : inner and outer barrels
            std::vector<RestDist::clu_info>* p_evt_sPH_nocolumn_vec = (layer_index == 0) ? (&evt_sPH_inner_nocolumn_vec) : (&evt_sPH_outer_nocolumn_vec);

            for (int clu_i = 0; clu_i < p_evt_sPH_nocolumn_vec->size(); clu_i++){
                
                RestDist::clu_info this_clu = p_evt_sPH_nocolumn_vec -> at(clu_i);

                h1D_map["h1D_Cluster_Z"] -> Fill(this_clu.z, INTTvtxZWeighting);
                h1D_map["h1D_ClusEta_MBDz"] -> Fill(this_clu.eta_MBDz, INTTvtxZWeighting);
                h1D_map["h1D_ClusPhi"] -> Fill(this_clu.phi_avgXY, INTTvtxZWeighting);
                h1D_map["h1D_ClusEta_INTTz"] -> Fill(this_clu.eta_INTTz, INTTvtxZWeighting);
                h1D_map["h1D_Cluster_phi_size"] -> Fill(this_clu.phi_size, INTTvtxZWeighting);
                h1D_map["h1D_Cluster_adc"] -> Fill(this_clu.adc, INTTvtxZWeighting);
                if (this_clu.phi_size == 1) {h1D_map["h1D_Cluster_adc_size1"] -> Fill(this_clu.adc, INTTvtxZWeighting);}
                if (runnumber == -1) {h1D_map["h1D_ClusPhi_TrueXY"] -> Fill(this_clu.phi_TrueXY), INTTvtxZWeighting;}
                if (runnumber == -1) {h1D_map["h1D_ClusEta_TrueXYZ"] -> Fill(this_clu.eta_TrueXYZ), INTTvtxZWeighting;}

                

                if (total_NClus > HighNClus){
                    h1D_map["h1D_ClusEta_MBDz_highNClus"] -> Fill(this_clu.eta_MBDz, INTTvtxZWeighting);
                    h1D_map["h1D_ClusEta_INTTz_highNClus"] -> Fill(this_clu.eta_INTTz, INTTvtxZWeighting);

                    if (runnumber == -1) { h1D_map["h1D_ClusEta_TrueXYZ_highNClus"] -> Fill(this_clu.eta_TrueXYZ), INTTvtxZWeighting;}
                }

                // note : h2D
                h2D_map["h2D_ClusPhiSize_ClusAdc"] -> Fill(this_clu.phi_size, this_clu.adc, INTTvtxZWeighting);
                h2D_map["h2D_ClusXY"] -> Fill(this_clu.x, this_clu.y, INTTvtxZWeighting);
                h2D_map["h2D_INTT_ClusEta_ClusAdc"] -> Fill(this_clu.eta_INTTz, this_clu.adc, INTTvtxZWeighting);
                h2D_map["h2D_ClusZ_LadderZID"] -> Fill(this_clu.z, this_clu.sensorZID, INTTvtxZWeighting);
                h2D_map["h2D_sensor_HitMap"] -> Fill(this_clu.layerID * 10 + this_clu.sensorZID, this_clu.ladderPhiID, INTTvtxZWeighting);
                h2D_map["h2D_sensor_AdcMap"] -> Fill(this_clu.layerID * 10 + this_clu.sensorZID, this_clu.ladderPhiID, this_clu.adc * INTTvtxZWeighting);

                if (this_clu.phi_size > 30) {h2D_map["h2D_ClusXY_LargePhiSize"]->Fill(this_clu.x, this_clu.y, INTTvtxZWeighting);}
                if (this_clu.phi_size == 43 || this_clu.phi_size == 46){h2D_map["h2D_ClusXY_PhiSize43_46"] -> Fill(this_clu.x, this_clu.y, INTTvtxZWeighting);}
                if (this_clu.phi_size >= 49){h2D_map["h2D_ClusXY_LargerPhiSize49"]->Fill(this_clu.x, this_clu.y, INTTvtxZWeighting);}
                h2D_map["h2D_NClus_ClusPhiSize"] -> Fill(total_NClus, this_clu.phi_size, INTTvtxZWeighting);

                if (total_NClus > HighNClus){
                    h2D_map["h2D_INTT_ClusEta_ClusAdc_highNClus"] -> Fill(this_clu.eta_INTTz, this_clu.adc, INTTvtxZWeighting);
                }

                // Division : ---type--------------------------------------------------------------------------------------------------------------------------
                if (this_clu.sensorZID == typeA_sensorZID[0] || this_clu.sensorZID == typeA_sensorZID[1]){ // note : type A
                    h1D_map["h1D_Cluster_Z_typeA"] -> Fill(this_clu.z, INTTvtxZWeighting);
                    h1D_map["h1D_ClusEta_MBDz_typeA"] -> Fill(this_clu.eta_MBDz, INTTvtxZWeighting);
                    h1D_map["h1D_ClusEta_INTTz_typeA"] -> Fill(this_clu.eta_INTTz, INTTvtxZWeighting);

                    if (runnumber == -1) {h1D_map["h1D_ClusEta_TrueXYZ_typeA"] -> Fill(this_clu.eta_TrueXYZ), INTTvtxZWeighting;}

                    if (total_NClus > HighNClus){
                        h1D_map["h1D_ClusEta_MBDz_typeA_highNClus"] -> Fill(this_clu.eta_MBDz, INTTvtxZWeighting);
                        h1D_map["h1D_ClusEta_INTTz_typeA_highNClus"] -> Fill(this_clu.eta_INTTz, INTTvtxZWeighting);

                        if (runnumber == -1) {h1D_map["h1D_ClusEta_TrueXYZ_typeA_highNClus"] -> Fill(this_clu.eta_TrueXYZ, INTTvtxZWeighting);}
                    }
                }
                else { // note : type B
                    h1D_map["h1D_Cluster_Z_typeB"] -> Fill(this_clu.z, INTTvtxZWeighting);
                    h1D_map["h1D_ClusEta_MBDz_typeB"] -> Fill(this_clu.eta_MBDz, INTTvtxZWeighting);
                    h1D_map["h1D_ClusEta_INTTz_typeB"] -> Fill(this_clu.eta_INTTz, INTTvtxZWeighting);

                    if (runnumber == -1) {h1D_map["h1D_ClusEta_TrueXYZ_typeB"] -> Fill(this_clu.eta_TrueXYZ), INTTvtxZWeighting;}

                    if (total_NClus > HighNClus){
                        h1D_map["h1D_ClusEta_MBDz_typeB_highNClus"] -> Fill(this_clu.eta_MBDz, INTTvtxZWeighting);
                        h1D_map["h1D_ClusEta_INTTz_typeB_highNClus"] -> Fill(this_clu.eta_INTTz, INTTvtxZWeighting);

                        if (runnumber == -1) {h1D_map["h1D_ClusEta_TrueXYZ_typeB_highNClus"] -> Fill(this_clu.eta_TrueXYZ, INTTvtxZWeighting);}
                    }
                }

            }
        }


        // Division : --- For pairs ---------------------------------------------------------------------------------------------------------------------
        inner_clu_phi_map.clear();
        outer_clu_phi_map.clear();
        inner_clu_phi_map = std::vector<std::vector<std::pair<bool,RestDist::clu_info>>>(360);
        outer_clu_phi_map = std::vector<std::vector<std::pair<bool,RestDist::clu_info>>>(360);

        for (int inner_i = 0; inner_i < int(evt_sPH_inner_nocolumn_vec.size()); inner_i++) {
            RestDist::clu_info this_clu = evt_sPH_inner_nocolumn_vec[inner_i];
            double Clus_InnerPhi_Offset = (this_clu.phi_avgXY < 0) ? this_clu.phi_avgXY * (180./TMath::Pi()) + 360 : this_clu.phi_avgXY * (180./TMath::Pi());
            inner_clu_phi_map[ int(Clus_InnerPhi_Offset) ].push_back({false,this_clu});
        }
        for (int outer_i = 0; outer_i < int(evt_sPH_outer_nocolumn_vec.size()); outer_i++) {
            RestDist::clu_info this_clu = evt_sPH_outer_nocolumn_vec[outer_i];
            double Clus_OuterPhi_Offset = (this_clu.phi_avgXY < 0) ? this_clu.phi_avgXY * (180./TMath::Pi()) + 360 : this_clu.phi_avgXY * (180./TMath::Pi());
            outer_clu_phi_map[ int(Clus_OuterPhi_Offset) ].push_back({false,this_clu});            
        }

        for (int inner_phi_i = 0; inner_phi_i < 360; inner_phi_i++) // note : each phi cell (1 degree)
        {
            // note : N cluster in this phi cell
            for (int inner_phi_clu_i = 0; inner_phi_clu_i < int(inner_clu_phi_map[inner_phi_i].size()); inner_phi_clu_i++)
            {
                RestDist::clu_info inner_clu = inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].second;

                // todo: change the outer phi scan range
                // note : the outer phi index, -1, 0, 1
                // note : the outer phi index, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 for the scan test
                for (int scan_i = -5; scan_i < 6; scan_i++)
                {
                    int true_scan_i = ((inner_phi_i + scan_i) < 0) ? 360 + (inner_phi_i + scan_i) : ((inner_phi_i + scan_i) > 359) ? (inner_phi_i + scan_i)-360 : inner_phi_i + scan_i;

                    // note : N clusters in that outer phi cell
                    for (int outer_phi_clu_i = 0; outer_phi_clu_i < int(outer_clu_phi_map[true_scan_i].size()); outer_phi_clu_i++)
                    {                        
                        RestDist::clu_info outer_clu = outer_clu_phi_map[true_scan_i][outer_phi_clu_i].second;

                        double deltaEta = inner_clu.eta_INTTz - outer_clu.eta_INTTz;
                        double deltaPhi = inner_clu.phi_avgXY - outer_clu.phi_avgXY;

                        h1D_map["h1D_PairDeltaEta"] -> Fill(deltaEta, INTTvtxZWeighting);
                        h1D_map["h1D_PairDeltaPhi"] -> Fill(deltaPhi, INTTvtxZWeighting);

                    }
                }
            }
        }         


        // for (int inner_i = 0; inner_i < evt_sPH_inner_nocolumn_vec.size(); inner_i++){
        //     for (int outer_i = 0; outer_i < evt_sPH_outer_nocolumn_vec.size(); outer_i++){

        //         RestDist::clu_info inner_clu = evt_sPH_inner_nocolumn_vec.at(inner_i);
        //         RestDist::clu_info outer_clu = evt_sPH_outer_nocolumn_vec.at(outer_i);

        //         double deltaEta = inner_clu.eta_INTTz - outer_clu.eta_INTTz;
        //         double deltaPhi = inner_clu.phi_avgXY - outer_clu.phi_avgXY;

        //         h1D_map["h1D_PairDeltaEta"] -> Fill(deltaEta, INTTvtxZWeighting);
        //         h1D_map["h1D_PairDeltaPhi"] -> Fill(deltaPhi, INTTvtxZWeighting);
        //     }
        // }
        
    }// note : end of event loop

}

void RestDist::EndRun()
{
    file_out -> cd();
    h1D_centrality_bin -> Write();

    for (auto &hist : h1D_map)
    {
        hist.second -> Write();
    }

    for (auto &hist : h2D_map)
    {
        hist.second -> Write();
    }

    h2D_map["h2D_sensor_AdcMap"] -> Sumw2(true);
    h2D_map["h2D_sensor_HitMap"] -> Sumw2(true);
    h2D_map["h2D_sensor_AvgAdcMap"] -> Divide(h2D_map["h2D_sensor_AdcMap"], h2D_map["h2D_sensor_HitMap"], 1, 1);
    h2D_map["h2D_sensor_AvgAdcMap"] -> Write();

    file_out -> Close();
}