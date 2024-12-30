#include "SemiChipClustering.h"

SemiChipClustering::SemiChipClustering(
    int process_id_in, // note : 1 or 2
    int runnumber_in, // note : still, (54280 or -1)
    int run_nEvents_in,
    std::string input_directory_in,
    std::string input_file_name_in,
    std::string output_directory_in,

    std::string output_file_name_suffix_in,

    bool BcoFullDiffCut_in,
    bool INTT_vtxZ_QA_in,
    bool ApplyHitQA_in,
    bool clone_hit_remove_BCO_tag_in,
    std::pair<bool, int> cut_HitBcoDiff_in,
    std::vector<int> adc_conversion_vec_in
):
    process_id(process_id_in),
    runnumber(runnumber_in),
    run_nEvents(run_nEvents_in),
    input_directory(input_directory_in),
    input_file_name(input_file_name_in),
    output_directory(output_directory_in),
    output_file_name_suffix(output_file_name_suffix_in),
    BcoFullDiffCut(BcoFullDiffCut_in),
    INTT_vtxZ_QA(INTT_vtxZ_QA_in),
    ApplyHitQA(ApplyHitQA_in),
    clone_hit_remove_BCO_tag(clone_hit_remove_BCO_tag_in),
    cut_HitBcoDiff(cut_HitBcoDiff_in),
    adc_conversion_vec(adc_conversion_vec_in)
{
    PrepareInputRootFile();

    run_nEvents = (run_nEvents == -1) ? tree_in->GetEntries() : run_nEvents;
    run_nEvents = (run_nEvents > tree_in->GetEntries()) ? tree_in->GetEntries() : run_nEvents;

    PrepareOutPutFileName();
    PrepareOutPutRootFile();
}

std::map<std::string, int> SemiChipClustering::GetInputTreeBranchesMap(TTree * m_tree_in)
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

void SemiChipClustering::PrepareInputRootFile(){
    file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_file_name.c_str()));
    if (!file_in || file_in -> IsZombie() || file_in == nullptr) {
        std::cout << "Error: cannot open file: " << input_file_name << std::endl;
        exit(1);
    }

    tree_in = (TTree*)file_in -> Get(tree_name.c_str());
    
    std::map<std::string, int> branch_map = GetInputTreeBranchesMap(tree_in);

    tree_in -> SetBranchStatus("*", 0);

    tree_in -> SetBranchStatus("MBD_z_vtx", 1);
    tree_in -> SetBranchStatus("is_min_bias", 1);
    tree_in -> SetBranchStatus("MBD_centrality", 1);
    tree_in -> SetBranchStatus("MBD_south_charge_sum", 1);
    tree_in -> SetBranchStatus("MBD_north_charge_sum", 1);
    tree_in -> SetBranchStatus("MBD_charge_sum", 1);
    tree_in -> SetBranchStatus("MBD_charge_asymm", 1);
    tree_in -> SetBranchStatus("InttBcoFullDiff_next", 1);

    tree_in -> SetBranchStatus("MBDNSg2",1);
    tree_in -> SetBranchStatus("MBDNSg2_vtxZ10cm",1);
    tree_in -> SetBranchStatus("MBDNSg2_vtxZ30cm",1);
    tree_in -> SetBranchStatus("MBDNSg2_vtxZ60cm",1);

    tree_in -> SetBranchStatus("InttRawHit_bco", 1);
    tree_in -> SetBranchStatus("InttRawHit_packetid", 1);
    tree_in -> SetBranchStatus("InttRawHit_fee", 1);
    tree_in -> SetBranchStatus("InttRawHit_channel_id", 1);
    tree_in -> SetBranchStatus("InttRawHit_chip_id", 1);
    tree_in -> SetBranchStatus("InttRawHit_adc", 1);
    tree_in -> SetBranchStatus("InttRawHit_FPHX_BCO", 1);

    // note : INTT vertex Z
    tree_in -> SetBranchStatus("INTTvtxZ", 1);
    tree_in -> SetBranchStatus("INTTvtxZError", 1);
    tree_in -> SetBranchStatus("NgroupTrapezoidal", 1);
    tree_in -> SetBranchStatus("NgroupCoarse", 1);
    tree_in -> SetBranchStatus("TrapezoidalFitWidth", 1);
    tree_in -> SetBranchStatus("TrapezoidalFWHM", 1);

    tree_in -> SetBranchStatus("NClus", 1);


    InttRawHit_bco = 0;
    InttRawHit_packetid = 0;
    InttRawHit_fee = 0;
    InttRawHit_channel_id = 0;
    InttRawHit_chip_id = 0;
    InttRawHit_adc = 0;
    InttRawHit_FPHX_BCO = 0;

    tree_in -> SetBranchAddress("MBD_z_vtx", &MBD_z_vtx);
    tree_in -> SetBranchAddress("is_min_bias", &is_min_bias);
    tree_in -> SetBranchAddress("MBD_centrality", &MBD_centrality);
    tree_in -> SetBranchAddress("MBD_south_charge_sum", &MBD_south_charge_sum);
    tree_in -> SetBranchAddress("MBD_north_charge_sum", &MBD_north_charge_sum);
    tree_in -> SetBranchAddress("MBD_charge_sum", &MBD_charge_sum);
    tree_in -> SetBranchAddress("MBD_charge_asymm", &MBD_charge_asymm);
    tree_in -> SetBranchAddress("InttBcoFullDiff_next", &InttBcoFullDiff_next);

    tree_in -> SetBranchAddress("MBDNSg2",&MBDNSg2);
    tree_in -> SetBranchAddress("MBDNSg2_vtxZ10cm",&MBDNSg2_vtxZ10cm);
    tree_in -> SetBranchAddress("MBDNSg2_vtxZ30cm",&MBDNSg2_vtxZ30cm);
    tree_in -> SetBranchAddress("MBDNSg2_vtxZ60cm",&MBDNSg2_vtxZ60cm);

    tree_in -> SetBranchAddress("InttRawHit_bco", &InttRawHit_bco);
    tree_in -> SetBranchAddress("InttRawHit_packetid", &InttRawHit_packetid);
    tree_in -> SetBranchAddress("InttRawHit_fee", &InttRawHit_fee);
    tree_in -> SetBranchAddress("InttRawHit_channel_id", &InttRawHit_channel_id);
    tree_in -> SetBranchAddress("InttRawHit_chip_id", &InttRawHit_chip_id);
    tree_in -> SetBranchAddress("InttRawHit_adc", &InttRawHit_adc);
    tree_in -> SetBranchAddress("InttRawHit_FPHX_BCO", &InttRawHit_FPHX_BCO);

    // note : INTT vertex Z
    tree_in -> SetBranchAddress("INTTvtxZ", &INTTvtxZ);
    tree_in -> SetBranchAddress("INTTvtxZError", &INTTvtxZError);
    tree_in -> SetBranchAddress("NgroupTrapezoidal", &NgroupTrapezoidal);
    tree_in -> SetBranchAddress("NgroupCoarse", &NgroupCoarse);
    tree_in -> SetBranchAddress("TrapezoidalFitWidth", &TrapezoidalFitWidth);
    tree_in -> SetBranchAddress("TrapezoidalFWHM", &TrapezoidalFWHM);

    tree_in -> SetBranchAddress("NClus", &NClus);


}

void SemiChipClustering::PrepareOutPutFileName()
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

    output_filename = "SemiChipClus";
    output_filename = (runnumber != -1) ? "Data_" + output_filename : "MC_" + output_filename;
    
    output_filename += (BcoFullDiffCut && runnumber != -1) ? "_BcoFullDiffCut" : "";
    output_filename += (INTT_vtxZ_QA) ? "_VtxZQA" : "";

    output_filename += (ApplyHitQA) ? "_HitQA" : "";
    output_filename += (cut_HitBcoDiff.first) ? Form("_HitBcoDiff%d",cut_HitBcoDiff.second) : "";
    output_filename += (clone_hit_remove_BCO_tag) ? "_CloneHitRm" : "";

    output_filename += output_file_name_suffix;
    output_filename += (runnumber != -1) ? Form("_%s_%s.root",runnumber_str.c_str(),job_index.c_str()) : Form("_%s.root",job_index.c_str());
}

void SemiChipClustering::PrepareOutPutRootFile()
{
    file_out = new TFile(Form("%s/%s", output_directory.c_str(), output_filename.c_str()), "RECREATE");
}

void SemiChipClustering::PrepareHists()
{
    h1D_NSize2Clus = new TH1D("h1D_NSize2Clus","h1D_NSize2Clus;NSize2Clus;Entries", 30, 0, 30);
    h1D_NSize1Clus = new TH1D("h1D_NSize1Clus","h1D_NSize1Clus;NSize1Clus;Entries", 30, 0, 30);
    
    h2D_NClus_NSize2Clus = new TH2D("h2D_NClus_NSize2Clus","h2D_NClus_NSize2Clus;NClus;NSize2Clus", 200, 0, 10000, 30, 0, 30);
    h2D_NClus_NSize1Clus = new TH2D("h2D_NClus_NSize1Clus","h2D_NClus_NSize1Clus;NClus;NSize1Clus", 200, 0, 10000, 30, 0, 30);
    h2D_nInttRawHit_NSize2Clus = new TH2D("h2D_nInttRawHit_NSize2Clus","h2D_nInttRawHit_NSize2Clus;NInttRawHits;NSize2Clus",200, 0, 25000, 30, 0, 30);
    h2D_nInttRawHit_NSize1Clus = new TH2D("h2D_nInttRawHit_NSize1Clus","h2D_nInttRawHit_NSize1Clus;NInttRawHits;NSize1Clus",200, 0, 25000, 30, 0, 30);

    h2D_nChipHit_NSize2Clus = new TH2D("h2D_nChipHit_NSize2Clus","h2D_nChipHit_NSize2Clus;NChipHits;NSize2Clus",128, 0, 128, 30, 0, 30);
    h2D_nChipHit_NSize1Clus = new TH2D("h2D_nChipHit_NSize1Clus","h2D_nChipHit_NSize1Clus;NChipHits;NSize1Clus",128, 0, 128, 30, 0, 30);

    h1D_ClusSizeCount = new TH1D("h1D_ClusSizeCount","h1D_ClusSizeCount;ClusSize;Entries", 130,0,130);
    h1D_ClusSizeCount_all = new TH1D("h1D_ClusSizeCount_all","h1D_ClusSizeCount_all;ClusSize;Entries", 130,0,130);

    h1D_chip_hit_container_map.clear();
    for (int server_i = 0; server_i < nFelix; server_i++){
        for (int Fch_i = 0; Fch_i < nFelix_channel; Fch_i++){
            for (int chip_i = 0; chip_i < nChip; chip_i++){
                
                h1D_chip_hit_container_map[Form("server%d_Fch%d_chip%d",server_i,Fch_i,chip_i)] = 
                new TH1D(Form("server%d_Fch%d_chip%d",server_i,Fch_i,chip_i), Form("server%d_Fch%d_chip%d",server_i,Fch_i,chip_i), 128, 0, 128);

            }
        }
    }
}

void SemiChipClustering::EvtCleanUp()
{
    evt_inttHits_map.clear();

    for ( auto &pair : h1D_chip_hit_container_map)
    {
        pair.second -> Reset("ICESM");
    }

    h1D_ClusSizeCount -> Reset("ICESM");
}

bool SemiChipClustering::DoSemiChipCluster(
    int eID_count,
    int NClus_in,
    std::vector<unsigned long> bco_vec_in,
    std::vector<unsigned int> packetid_vec_in,
    std::vector<unsigned short> fee_vec_in,
    std::vector<unsigned short> channel_id_vec_in,
    std::vector<unsigned short> chip_id_vec_in,
    std::vector<unsigned short> adc_vec_in,
    std::vector<unsigned short> FPHX_BCO_vec_in
)
{
    EvtCleanUp();

    for (int hit_i = 0; hit_i < bco_vec_in.size(); hit_i++)
    {
        if (hit_i == 0) {evt_INTT_bco_full = bco_vec_in[hit_i];}

        int bco_full = bco_vec_in[hit_i];

        int server = packetid_vec_in[hit_i] - Felix_offset; // note : the felix server ID
        int felix_ch = fee_vec_in[hit_i]; // note : the felix channel ID 0 - 13
        int chip = (chip_id_vec_in[hit_i] - 1) % 26; // note : chip ID 0 - 25
        int channel = channel_id_vec_in[hit_i]; // note : channel ID 0 - 127
        int adc = adc_vec_in[hit_i]; // note : adc value
        int bco = FPHX_BCO_vec_in[hit_i];
        int bco_diff = (bco - (bco_full & 0x7fU) + 128) % 128;


        if (eID_count % 1000 == 0 && hit_i%25 == 0) {
            std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - server: "<< server << " felix_ch: "<< felix_ch << " chip: "<< chip << " channel: "<< channel << " adc: "<< adc << " bco: "<< bco << " bco_diff: "<< bco_diff << std::endl;
        }

        if (ApplyHitQA && (server < 0 || server > 7)) // note : server ID 0 - 7
        {
            std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - server out of range, this server is : "<< server << std::endl;
            continue;
        }

        if (ApplyHitQA && (felix_ch < 0 || felix_ch > 13)) // note : felix channel ID 0 - 13
        {
            std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - felix channel out of range, this channel is : "<< felix_ch << std::endl;
            continue;
        }

        if (ApplyHitQA && (chip < 0 || chip > 25)) // note : chip ID 0 - 25
        {
            std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - chip out of range, this chip is : "<< chip << std::endl;
            continue;
        }

        if (ApplyHitQA && (channel < 0 || channel > 127)) // note : channel ID 0 - 127
        {
            std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - channel out of range, this channel is : "<< channel << std::endl;
            continue;
        }

        if (ApplyHitQA && (adc < 0 || adc > 7)) // note : adc ID 0 - 7
        {
            std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - adc out of range, this adc is : "<< adc << std::endl;
            continue;
        }
        
        if (ApplyHitQA && (bco < 0 || bco > 127)) // note : bco ID 0 - 127
        {
            std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - bco out of range, this bco is : "<< bco << std::endl;
            continue;
        }

        if (cut_HitBcoDiff.first && (bco_diff < cut_HitBcoDiff.second - 1 || bco_diff > cut_HitBcoDiff.second + 1)) // note : bco_diff cut
        {
            // std::cout << "eID: "<< eID_count <<" INTTBcoResolution::PrepareINTT - bco_diff out of range, this bco_diff is : "<< bco_diff << std::endl;
            continue;
        }

        SemiChipClustering::inttHitstr this_hit;
        
        this_hit.hit_server = server;
        this_hit.hit_felix_ch = felix_ch;
        this_hit.hit_chip = chip;
        this_hit.hit_channel = channel;
        this_hit.hit_adc = adc_conversion_vec[adc];
        this_hit.hit_bco = bco;
        this_hit.hit_bco_diff = bco_diff;

        if (clone_hit_remove_BCO_tag) // note : can only be one of each
        {
            if (evt_inttHits_map.find(Form("%d_%d_%d_%d_%d",server,felix_ch,chip,channel,bco)) == evt_inttHits_map.end()) // note : if it's not found, we just add it
            {
                evt_inttHits_map[Form("%d_%d_%d_%d_%d",server,felix_ch,chip,channel,bco)] = this_hit;
            }
        }
        else // note : if we don't want to remove the clone hits
        {
            evt_inttHits_map[Form("%d",hit_i)] = this_hit; // note : only index i to make the key unique
        }

    } // note : end of loop over hits

    for (auto pair : evt_inttHits_map)
    {
        SemiChipClustering::inttHitstr this_hit = pair.second;
        int server = this_hit.hit_server;
        int felix_ch = this_hit.hit_felix_ch;
        int chip = this_hit.hit_chip;
        int channel = this_hit.hit_channel;
        int adc = this_hit.hit_adc;
        int bco = this_hit.hit_bco;
        int bco_diff = this_hit.hit_bco_diff;

        h1D_chip_hit_container_map[Form("server%d_Fch%d_chip%d",server,felix_ch,chip)] -> Fill(channel, adc);
    }

    for (auto pair : h1D_chip_hit_container_map)
    {
        if (pair.second -> GetEntries() == 0) {continue;}
        
        SemiChipClustering::chip_clu_info clu_info = find_Ngroup(pair.second);
        h1D_ClusSizeCount -> Reset("ICESM");

        for (auto size : clu_info.size_vec)
        {
            h1D_ClusSizeCount -> Fill(size);
            h1D_ClusSizeCount_all -> Fill(size);
        }


        h1D_NSize2Clus -> Fill(h1D_ClusSizeCount -> GetBinContent(3));
        h1D_NSize1Clus -> Fill(h1D_ClusSizeCount -> GetBinContent(2));
        
        h2D_NClus_NSize2Clus -> Fill(NClus_in, h1D_ClusSizeCount -> GetBinContent(3));
        h2D_NClus_NSize1Clus -> Fill(NClus_in, h1D_ClusSizeCount -> GetBinContent(2));
        
        h2D_nInttRawHit_NSize2Clus -> Fill(evt_inttHits_map.size(), h1D_ClusSizeCount -> GetBinContent(3));
        h2D_nInttRawHit_NSize1Clus -> Fill(evt_inttHits_map.size(), h1D_ClusSizeCount -> GetBinContent(2));
        
        h2D_nChipHit_NSize2Clus -> Fill(pair.second -> GetEntries(), h1D_ClusSizeCount -> GetBinContent(3));
        h2D_nChipHit_NSize1Clus -> Fill(pair.second -> GetEntries(), h1D_ClusSizeCount -> GetBinContent(2));
    
    }

    return true;
}

SemiChipClustering::chip_clu_info SemiChipClustering::find_Ngroup(TH1D * hist_in)
{
    double Highest_bin_Content __attribute__((unused)) = hist_in -> GetBinContent(hist_in -> GetMaximumBin());
    double Highest_bin_Center   = hist_in -> GetBinCenter(hist_in -> GetMaximumBin());

    int group_Nbin = 0;
    int peak_group_ID = -9999;
    double group_entry = 0;
    double peak_group_ratio __attribute__((unused));
    std::vector<int> group_Nbin_vec; group_Nbin_vec.clear();
    std::vector<double> group_entry_vec; group_entry_vec.clear();
    std::vector<double> group_widthL_vec; group_widthL_vec.clear();
    std::vector<double> group_widthR_vec; group_widthR_vec.clear();
    std::vector<double> clu_pos_vec; clu_pos_vec.clear();

    for (int i = 0; i < hist_in -> GetNbinsX(); i++){
        // todo : the background rejection is here : Highest_bin_Content/2. for the time being
        // double bin_content = ( hist_in -> GetBinContent(i+1) <= Highest_bin_Content/2.) ? 0. : ( hist_in -> GetBinContent(i+1) - Highest_bin_Content/2. );

        double bin_content = hist_in -> GetBinContent(i+1);

        if (bin_content != 0){
            
            if (group_Nbin == 0) {
                group_widthL_vec.push_back(hist_in -> GetBinCenter(i+1) - (hist_in -> GetBinWidth(i+1)/2.));
            }

            group_Nbin += 1; 
            group_entry += bin_content;
        }
        else if (bin_content == 0 && group_Nbin != 0){
            group_widthR_vec.push_back(hist_in -> GetBinCenter(i+1) - (hist_in -> GetBinWidth(i+1)/2.));
            group_Nbin_vec.push_back(group_Nbin);
            group_entry_vec.push_back(group_entry);
            group_Nbin = 0;
            group_entry = 0;
        }
    }
    if (group_Nbin != 0) {
        group_Nbin_vec.push_back(group_Nbin);
        group_entry_vec.push_back(group_entry);
        group_widthR_vec.push_back(hist_in -> GetXaxis()->GetXmax());
    } // note : the last group at the edge

    // note : find the peak group
    for (int i = 0; i < int(group_Nbin_vec.size()); i++){
        if (group_widthL_vec[i] < Highest_bin_Center && Highest_bin_Center < group_widthR_vec[i]){
            peak_group_ID = i;
            break;
        }
    }
    
    // note : On Nov 6, 2024, we no longer need to calculate the ratio of the peak group
    // double denominator = (accumulate( group_entry_vec.begin(), group_entry_vec.end(), 0.0 ));
    // denominator = (denominator <= 0) ? 1. : denominator;
    // peak_group_ratio = group_entry_vec[peak_group_ID] / denominator;
    peak_group_ratio = -999;

    double peak_group_left __attribute__((unused))  = (double(group_Nbin_vec.size()) == 0) ? -999 : group_widthL_vec[peak_group_ID];
    double peak_group_right __attribute__((unused)) = (double(group_Nbin_vec.size()) == 0) ? 999  : group_widthR_vec[peak_group_ID];

    for (int i = 0; i < int(group_Nbin_vec.size()); i++){
        clu_pos_vec.push_back( CoM_cluster_pos(hist_in, group_widthL_vec[i], group_widthR_vec[i]) );
    }

    SemiChipClustering::chip_clu_info out_str;

    out_str.adc_vec = group_entry_vec;
    out_str.size_vec = group_Nbin_vec;
    out_str.edge_l_vec = group_widthL_vec;
    out_str.edge_r_vec = group_widthR_vec;
    out_str.pos_vec = clu_pos_vec;
    out_str.largest_size = *(std::max_element(group_Nbin_vec.begin(), group_Nbin_vec.end()));

    // for (int i = 0; i < group_Nbin_vec.size(); i++)
    // {
    //     std::cout<<" "<<std::endl;
    //     std::cout<<"group width : "<<group_Nbin_vec[i]<<std::endl;
    //     std::cout<<"group adc : "<<group_entry_vec[i]<<std::endl;
    //     std::cout<<group_widthL_vec[i]<<" "<<group_widthR_vec[i]<<std::endl;
    // }

    // cout<<" "<<endl;
    // cout<<"N group : "<<group_Nbin_vec.size()<<endl;
    // cout<<"Peak group ID : "<<peak_group_ID<<endl;
    // cout<<"peak group width : "<<group_widthL_vec[peak_group_ID]<<" "<<group_widthR_vec[peak_group_ID]<<endl;
    // cout<<"ratio : "<<peak_group_ratio<<endl;
    
    // note : {N_group, ratio (if two), peak widthL, peak widthR}
    return out_str;
}

double SemiChipClustering::CoM_cluster_pos(TH1D * hist_in, double edge_l, double edge_r)
{
    double sum = 0;
    double sum_weight = 0;

    for (int i = 0; i < hist_in -> GetNbinsX(); i++){
        double bin_content = hist_in -> GetBinContent(i+1);
        double bin_center = hist_in -> GetBinCenter(i+1);

        if (edge_l < bin_center && bin_center < edge_r){
            sum += bin_content * bin_center;
            sum_weight += bin_content;
        }
    }

    return (sum_weight == 0) ? -999 : (sum / sum_weight) - hist_in -> GetBinWidth(1)/2.;
}

void SemiChipClustering::MainProcess()
{
    for (int i = 0; i < run_nEvents; i++)
    {
        tree_in -> GetEntry(i);
        EvtCleanUp();

        if (i % 10 == 0) {std::cout << "Processing event " << i<<", NClus : "<< NClus << std::endl;}

        // =======================================================================================================================================================
        // note : hard cut

        // note : for data
        if (runnumber != -1 && BcoFullDiffCut && InttBcoFullDiff_next <= cut_InttBcoFullDIff_next) {continue;}
        if (runnumber != -1 && MBDNSg2 != 1) {continue;} // todo: assume MC no trigger

        // note : for MC
        // if (runnumber == -1 && NTruthVtx != 1) {continue;}

        // note : both data and MC
        if (is_min_bias != 1) {continue;}
        if (MBD_z_vtx != MBD_z_vtx) {continue;}
        if (MBD_centrality != MBD_centrality) {continue;}
        if (MBD_centrality < 0 || MBD_centrality > 1) {continue;}
        if (INTTvtxZ != INTTvtxZ) {continue;}
        if (MBD_z_vtx < cut_GlobalMBDvtxZ.first || MBD_z_vtx > cut_GlobalMBDvtxZ.second) {continue;} // todo: the hard cut 60 cm 
        // =======================================================================================================================================================

        // =======================================================================================================================================================
        // note : optional cut
        if (INTT_vtxZ_QA && (MBD_z_vtx - INTTvtxZ < cut_vtxZDiff.first || MBD_z_vtx - INTTvtxZ > cut_vtxZDiff.second) ) {continue;}
        if (INTT_vtxZ_QA && (TrapezoidalFitWidth < cut_TrapezoidalFitWidth.first || TrapezoidalFitWidth > cut_TrapezoidalFitWidth.second)){continue;}
        if (INTT_vtxZ_QA && (TrapezoidalFWHM < cut_TrapezoidalFWHM.first || TrapezoidalFWHM > cut_TrapezoidalFWHM.second)){continue;}
        if (INTT_vtxZ_QA && (INTTvtxZError < cut_INTTvtxZError.first || INTTvtxZError > cut_INTTvtxZError.second)){continue;}
        // =======================================================================================================================================================

        DoSemiChipCluster(
            i,
            NClus,
            *InttRawHit_bco,
            *InttRawHit_packetid,
            *InttRawHit_fee,
            *InttRawHit_channel_id,
            *InttRawHit_chip_id,
            *InttRawHit_adc,
            *InttRawHit_FPHX_BCO
        );
        

    }
}

void SemiChipClustering::EndRun()
{
    file_out -> cd();

    h1D_ClusSizeCount_all -> Write();
    h1D_NSize2Clus -> Write();
    h1D_NSize1Clus -> Write();
    h2D_NClus_NSize2Clus -> Write();
    h2D_NClus_NSize1Clus -> Write();
    h2D_nInttRawHit_NSize2Clus -> Write();
    h2D_nInttRawHit_NSize1Clus -> Write();
    h2D_nChipHit_NSize2Clus -> Write();
    h2D_nChipHit_NSize1Clus -> Write();

    file_out -> Close();
}