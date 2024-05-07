#ifndef EtaDistReader_h
#define EtaDistReader_h

#include "INTTEta.h"


// note : the code works with the histograms, so the centrality bins are already defined in the histograms, ranging from 0 to something big number
// todo : the signal_region and centrality_region is hard_written in the INTTEta.h, which is written in the "ana_map_v1.h", so modify there


class EtaDistReader : public INTTEta
{
    public:
        EtaDistReader(string run_type, string out_folder_directory, vector<pair<int,vector<int>>> included_eta_z_map, string input_file_directory, bool centrality_Z_map_bool = 0) : 
        INTTEta(run_type, out_folder_directory), included_eta_z_map(included_eta_z_map), input_file_directory(input_file_directory), centrality_Z_map_bool(centrality_Z_map_bool)
        {
            gErrorIgnoreLevel = kError;

            cout<<"In EtaDistReader, centrality region size : "<<centrality_region.size()<<endl;    
            cout<<"In EtaDistReader, plot_text: "<<plot_text<<endl;
            cout<<"In EtaDistReader, signal_region: "<<signal_region<<endl;


            file_in = TFile::Open(input_file_directory.c_str());

            N_centrality_bin = centrality_region.size();
            eta_correction = included_eta_z_map[0].first - 1; // todo : the eta bin starts from 15, so the correction is 14
            
            SignalNTrack_eta_z_Multi_2D.clear();
            SignalNTrack_eta_z_Single_2D.clear();
            TrueNtrack_eta_z_MC_2D.clear();
            dNdeta_1D_MC.clear();
            dNdeta_1D_reco_single.clear();
            dNdeta_1D_reco_multi.clear();
            DeltaPhi_Multi_1D.clear();
            DeltaPhi_Multi_Stack.clear();
            SignalNTrack_Single.clear(); 
            SignalNTrack_Multi.clear();
            N_event_counting = vector<vector<int>>(N_centrality_bin, vector<int>(included_eta_z_map.size(),0));
            N_event_counting_MC = vector<vector<int>>(N_centrality_bin, vector<int>(included_eta_z_map.size(),0));

            DeltaPhi_Multi_Stack_hist_out.clear();
            
            solid_line = new TLine();
            solid_line -> SetLineWidth(1);
            solid_line -> SetLineColor(2);
            solid_line -> SetLineStyle(1);

            // file_out = new TFile(Form("%s/alpha_correction_map.root",out_folder_directory.c_str()),"RECREATE");

            legend = new TLegend(0.45,0.8,0.70,0.9);
            // legend -> SetMargin(0);
            legend->SetTextSize(0.03);

            ReadFileHist();
            InitHist();
            MainPreparation();
            // FinaldNdEta();


            return;
        };

        vector<TH1F *> GetDeltaPhi_Multi_stack_1D();
        vector<TH1F *> GetdNdeta_1D_MC();
        vector<TH1F *> GetdNdeta_1D_reco_single();
        vector<TH1F *> GetdNdeta_1D_reco_multi();
        vector<string> Get_centrality_region() {return centrality_region;};
        vector<TH2F *> Get_alpha_corr_map() {return {eta_Mbin_correction_loose, eta_Mbin_correction_tight};};
        void FindCoveredRegion();
        void Set_alpha_corr_map(TH2F * hist_in_loose, TH2F * hist_in_tight);
        void FinaldNdEta();
        void EndRun();

        ~EtaDistReader()
        {
            file_in -> Close();
            cout<<"EtaDistReader done, goodbye"<<endl;
        };

    protected:

        bool centrality_Z_map_bool;
        vector<pair<int,vector<int>>> included_eta_z_map;
        string input_file_directory;
        int eta_correction; // todo : the eta bin starts from 15, so the correction is 14
        int N_centrality_bin;

        TFile * file_in;
        string color_code[5] = {"#fedc97", "#b5b682", "#7c9885", "#28666e", "#033f63"};
        // vector<string> centrality_region;

        vector<TH2F *> SignalNTrack_eta_z_Multi_2D;
        vector<TH2F *> SignalNTrack_eta_z_Single_2D;
        map<string,TH1F *> DeltaPhi_Multi_1D;
        vector<TH2F *> TrueNtrack_eta_z_MC_2D;
        TH2F * centrality_Z_map;
        TH2F * centrality_Z_map_MC;
        TH2F * eta_z_ref;
        TH2F * eta_z_ref_used_check;

        TH2F * eta_Mbin_correction_tight;
        TH2F * eta_Mbin_correction_loose;
        TH2F * eta_Mbin_correction_loose_noUPC; // note : no ultra peripheral collision bin
        TH2F * input_eta_Mbin_correction_tight;
        TH2F * input_eta_Mbin_correction_loose;


        vector<TH1F *> dNdeta_1D_MC;
        vector<TH1F *> dNdeta_1D_reco_single;
        vector<TH1F *> dNdeta_1D_reco_multi;

        map<string, THStack *> DeltaPhi_Multi_Stack;
        vector<TH1F *> DeltaPhi_Multi_Stack_hist_out;
        map<string, int> SignalNTrack_Single;
        map<string, int> SignalNTrack_Multi;
        vector<vector<int>> N_event_counting;
        vector<vector<int>> N_event_counting_MC;

        TLegend * legend;
        TLine * solid_line;
        TH1F * temp_hist;

        // TFile * file_out; 

        void ReadFileHist();
        void InitHist();
        void MainPreparation();
        void DrawCoordLine(TH2F * hist_in);
        void DrawEtaCoverBox(TH2F * hist_in);
        

};

void EtaDistReader::DrawCoordLine(TH2F * hist_in)
{
    // note : draw the vertical line, to segment the eta region
    for (int i = 0; i < hist_in -> GetNbinsX(); i++) { 
        coord_line -> DrawLine(
            hist_in->GetXaxis()->GetBinLowEdge(i+1), 
            hist_in->GetYaxis()->GetBinLowEdge(1), 
            hist_in->GetXaxis()->GetBinLowEdge(i+1), 
            hist_in->GetYaxis()->GetBinLowEdge(hist_in -> GetNbinsY()) + hist_in->GetYaxis()->GetBinWidth(hist_in -> GetNbinsY())
        ); 
    }
    coord_line -> DrawLine(
        hist_in->GetXaxis()->GetBinLowEdge(hist_in->GetNbinsX()) + hist_in->GetXaxis()->GetBinWidth(hist_in->GetNbinsX()), 
        hist_in->GetYaxis()->GetBinLowEdge(1), 
        hist_in->GetXaxis()->GetBinLowEdge(hist_in->GetNbinsX()) + hist_in->GetXaxis()->GetBinWidth(hist_in->GetNbinsX()), 
        hist_in->GetYaxis()->GetBinLowEdge(hist_in -> GetNbinsY()) + hist_in->GetYaxis()->GetBinWidth(hist_in -> GetNbinsY())
    );

    // note : draw the horizontal line, to segment the z region
    for (int i = 0; i < hist_in -> GetNbinsY(); i++) { 
        coord_line -> DrawLine(
            hist_in->GetXaxis()->GetBinLowEdge(1), 
            hist_in->GetYaxis()->GetBinLowEdge(i+1), 
            hist_in->GetXaxis()->GetBinLowEdge(hist_in->GetNbinsX()) + hist_in->GetXaxis()->GetBinWidth(hist_in->GetNbinsX()), 
            hist_in->GetYaxis()->GetBinLowEdge(i+1)
        ); 
    }
    coord_line -> DrawLine(
        hist_in->GetXaxis()->GetBinLowEdge(1), 
        hist_in->GetYaxis()->GetBinLowEdge(hist_in -> GetNbinsY()) + hist_in->GetYaxis()->GetBinWidth(hist_in -> GetNbinsY()), 
        hist_in->GetXaxis()->GetBinLowEdge(hist_in->GetNbinsX()) + hist_in->GetXaxis()->GetBinWidth(hist_in->GetNbinsX()), 
        hist_in->GetYaxis()->GetBinLowEdge(hist_in -> GetNbinsY()) + hist_in->GetYaxis()->GetBinWidth(hist_in -> GetNbinsY())
    ); 
}

void EtaDistReader::DrawEtaCoverBox(TH2F * hist_in)
{
    double INTT_layer_R[4] = {71.88, 77.32, 96.8, 102.62}; // note : the radii of the INTT layers

    // note : low  Z edge for left  eta
    // note : high Z edge for right eta
    for (int i = 0; i < hist_in -> GetNbinsY(); i++)
    {
        double bin_low_zvtx  = hist_in->GetYaxis()->GetBinLowEdge(i+1);
        double bin_high_zvtx = hist_in->GetYaxis()->GetBinLowEdge(i+1) + hist_in->GetYaxis()->GetBinWidth(i+1);

        double INTT_eta_acceptance_l = -0.5 * TMath::Log((sqrt(pow(-230.-bin_low_zvtx,2)+pow(INTT_layer_R[0],2))-(-230.-bin_low_zvtx)) / (sqrt(pow(-230.-bin_low_zvtx,2)+pow(INTT_layer_R[0],2))+(-230.-bin_low_zvtx))); // note : left
        double INTT_eta_acceptance_r = -0.5 * TMath::Log((sqrt(pow(230.-bin_high_zvtx,2)+pow(INTT_layer_R[3],2))-(230.-bin_high_zvtx)) / (sqrt(pow(230.-bin_high_zvtx,2)+pow(INTT_layer_R[3],2))+(230.-bin_high_zvtx))); // note : right
        
        solid_line -> DrawLine(INTT_eta_acceptance_l, bin_low_zvtx, INTT_eta_acceptance_l, bin_high_zvtx); // note : vertical
        solid_line -> DrawLine(INTT_eta_acceptance_r, bin_low_zvtx, INTT_eta_acceptance_r, bin_high_zvtx); // note : vertical
        solid_line -> DrawLine(INTT_eta_acceptance_l, bin_low_zvtx, INTT_eta_acceptance_r, bin_low_zvtx); // note : horizontal
        solid_line -> DrawLine(INTT_eta_acceptance_l, bin_high_zvtx, INTT_eta_acceptance_r, bin_high_zvtx); // note : horizontal
        
        cout<<"bin Yaxis ID: "<< i+1 <<" low_edge "<<bin_low_zvtx<<" high edge: "<<bin_high_zvtx<<" eta region: "<<INTT_eta_acceptance_l<<"~"<<INTT_eta_acceptance_r<<endl;   
    }
}

void EtaDistReader::Set_alpha_corr_map(TH2F * hist_in_loose, TH2F * hist_in_tight)
{
    input_eta_Mbin_correction_loose = hist_in_loose;
    input_eta_Mbin_correction_tight = hist_in_tight;

    std::cout<<"in EtaDistReader, input map loose check "<<input_eta_Mbin_correction_loose -> GetBinContent(3,3)<<endl;
    std::cout<<"in EtaDistReader, input map tight check "<<input_eta_Mbin_correction_tight -> GetBinContent(3,3)<<endl;

    return;
}

void EtaDistReader::ReadFileHist()
{

    // note : the histograms of the number of tracklets in the signal region for both methods
    // note : the 1D vector is for different centrality bins
    for (int Mbin = 0; Mbin < N_centrality_bin; Mbin++)
    {
        // note : the entries with the deltaphi within 1 degree, which means the background entries are also included. 
        // note : the loose method
        SignalNTrack_eta_z_Multi_2D.push_back((TH2F *) file_in -> Get(Form("Reco_SignalNTracklet_Multi_MBin%i",Mbin)));
        // note : the tight method
        SignalNTrack_eta_z_Single_2D.push_back((TH2F *) file_in -> Get(Form("Reco_SignalNTracklet_Single_MBin%i",Mbin)));
    }

    // note : the 1D histograms of the DeltaPhi distributions of each centrality-eta-z bin (the multi method, for the background subtraction)
    // note : the format of the key is "MBin_EtaBin_ZBin"
    for (int Mbin = 0; Mbin < N_centrality_bin; Mbin++) // note : centrality bin
    {
        for (auto eta_z : included_eta_z_map) // note : eta bin
        {
            for(int zbin = 0; zbin < eta_z.second.size(); zbin++) // note :  Z bin
            {
                DeltaPhi_Multi_1D[Form("%i_%i_%i",Mbin,eta_z.first,eta_z.second[zbin])] = (TH1F *) file_in -> Get(Form("Reco_DeltaPhi1D_Multi_MBin%i_Eta%i_Z%i",Mbin,eta_z.first,eta_z.second[zbin]));
                DeltaPhi_Multi_1D[Form("%i_%i_%i",Mbin,eta_z.first,eta_z.second[zbin])] -> SetFillColor(TColor::GetColor(color_code[zbin%5].c_str()));
            }
        }
    }

    // note : the histograms of the number of True track from the MC
    // note : the 1D vector is for different centrality bins
    for (int Mbin = 0; Mbin < N_centrality_bin; Mbin++)
    {
        TrueNtrack_eta_z_MC_2D.push_back((TH2F *) file_in -> Get(Form("NTrueTrack_MBin%i",Mbin)));
    }

    // note : the centrality Z map, to keep the N event counting
    string centrality_Z_map_name = (centrality_Z_map_bool == 0) ? "MBin_Z_evt_map" : "Z_MBin_evt_map"; // todo : the correct name should be "MBin_Z_evt_map", correct next time
    centrality_Z_map    = (TH2F *) file_in -> Get(centrality_Z_map_name.c_str());
    centrality_Z_map_MC = (TH2F *) file_in -> Get((centrality_Z_map_name+"_MC").c_str()); 

    // note : the reference map of Eta-Z 2D histogram
    eta_z_ref = (TH2F *) file_in -> Get("Eta_Z_reference");

    eta_z_ref_used_check = (TH2F*) eta_z_ref -> Clone("eta_z_ref_used_check");
    
    cout<<"ReadFileHist done"<<endl;

    return;

}

void EtaDistReader::FindCoveredRegion()
{
    c1 -> cd();
    eta_z_ref -> Draw("colz0");
    DrawCoordLine(eta_z_ref);    
    DrawEtaCoverBox(eta_z_ref);
    c1 -> Print(Form("%s/eta_z_ref.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

    // note : for the used bins check
    c1 -> cd();
    double eta_z_ref_used_check_max_entry = eta_z_ref_used_check->GetMaximum();
    for (auto eta_z : included_eta_z_map){
        for (auto zbin : eta_z.second){
            eta_z_ref_used_check -> SetBinContent(eta_z.first, zbin, eta_z_ref_used_check_max_entry * 1000);        
        }
    }
    eta_z_ref_used_check -> Draw("colz0");
    DrawCoordLine(eta_z_ref_used_check);    
    DrawEtaCoverBox(eta_z_ref_used_check);
    c1 -> Print(Form("%s/eta_z_ref_used_check.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

}

void EtaDistReader::InitHist()
{
    for (int i = 0; i < N_centrality_bin; i++)
    {
        dNdeta_1D_MC.push_back(new TH1F("","",included_eta_z_map.size(), eta_z_ref -> GetXaxis() ->GetBinLowEdge(included_eta_z_map[0].first), eta_z_ref -> GetXaxis() ->GetBinCenter(included_eta_z_map[included_eta_z_map.size() - 1].first) + eta_z_ref -> GetXaxis() ->GetBinWidth(included_eta_z_map[included_eta_z_map.size() - 1].first) / 2.));
        dNdeta_1D_MC[i] -> SetMarkerStyle(20);
        dNdeta_1D_MC[i] -> SetMarkerSize(0.8);
        dNdeta_1D_MC[i] -> SetMarkerColor(TColor::GetColor("#1A3947"));
        dNdeta_1D_MC[i] -> SetLineColor(TColor::GetColor("#1A3947"));
        dNdeta_1D_MC[i] -> SetLineWidth(2);
        dNdeta_1D_MC[i] -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");
        dNdeta_1D_MC[i] -> GetXaxis() -> SetTitle("#eta");
        // dNdeta_1D_MC[i] -> GetYaxis() -> SetRangeUser(0,50);
        dNdeta_1D_MC[i] -> SetTitleSize(0.06, "X");
        dNdeta_1D_MC[i] -> SetTitleSize(0.06, "Y");
        dNdeta_1D_MC[i] -> GetXaxis() -> SetTitleOffset (0.71);
        dNdeta_1D_MC[i] -> GetYaxis() -> SetTitleOffset (1.1);
        dNdeta_1D_MC[i] -> GetXaxis() -> CenterTitle(true);
        dNdeta_1D_MC[i] -> GetYaxis() -> CenterTitle(true);

        dNdeta_1D_reco_single.push_back(new TH1F("","",included_eta_z_map.size(), eta_z_ref -> GetXaxis() ->GetBinLowEdge(included_eta_z_map[0].first), eta_z_ref -> GetXaxis() ->GetBinCenter(included_eta_z_map[included_eta_z_map.size() - 1].first) + eta_z_ref -> GetXaxis() ->GetBinWidth(included_eta_z_map[included_eta_z_map.size() - 1].first) / 2.));
        dNdeta_1D_reco_single[i] -> SetMarkerStyle(20);
        dNdeta_1D_reco_single[i] -> SetMarkerSize(0.8);
        dNdeta_1D_reco_single[i] -> SetMarkerColor(TColor::GetColor("#1A3947"));
        dNdeta_1D_reco_single[i] -> SetLineColor(TColor::GetColor("#1A3947"));
        dNdeta_1D_reco_single[i] -> SetLineWidth(2);
        dNdeta_1D_reco_single[i] -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");
        dNdeta_1D_reco_single[i] -> GetXaxis() -> SetTitle("#eta");
        // dNdeta_1D_reco_single[i] -> GetYaxis() -> SetRangeUser(0,50);
        dNdeta_1D_reco_single[i] -> SetTitleSize(0.06, "X");
        dNdeta_1D_reco_single[i] -> SetTitleSize(0.06, "Y");
        dNdeta_1D_reco_single[i] -> GetXaxis() -> SetTitleOffset (0.71);
        dNdeta_1D_reco_single[i] -> GetYaxis() -> SetTitleOffset (1.1);
        dNdeta_1D_reco_single[i] -> GetXaxis() -> CenterTitle(true);
        dNdeta_1D_reco_single[i] -> GetYaxis() -> CenterTitle(true);

        dNdeta_1D_reco_multi.push_back(new TH1F("","",included_eta_z_map.size(), eta_z_ref -> GetXaxis() ->GetBinLowEdge(included_eta_z_map[0].first), eta_z_ref -> GetXaxis() ->GetBinCenter(included_eta_z_map[included_eta_z_map.size() - 1].first) + eta_z_ref -> GetXaxis() ->GetBinWidth(included_eta_z_map[included_eta_z_map.size() - 1].first) / 2.));
        dNdeta_1D_reco_multi[i] -> SetMarkerStyle(20);
        dNdeta_1D_reco_multi[i] -> SetMarkerSize(0.8);
        dNdeta_1D_reco_multi[i] -> SetMarkerColor(TColor::GetColor("#c48045"));
        dNdeta_1D_reco_multi[i] -> SetLineColor(TColor::GetColor("#c48045"));
        dNdeta_1D_reco_multi[i] -> SetLineWidth(2);
        dNdeta_1D_reco_multi[i] -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");
        dNdeta_1D_reco_multi[i] -> GetXaxis() -> SetTitle("#eta");
        // dNdeta_1D_reco_multi[i] -> GetYaxis() -> SetRangeUser(0,50);
        dNdeta_1D_reco_multi[i] -> SetTitleSize(0.06, "X");
        dNdeta_1D_reco_multi[i] -> SetTitleSize(0.06, "Y");
        dNdeta_1D_reco_multi[i] -> GetXaxis() -> SetTitleOffset (0.71);
        dNdeta_1D_reco_multi[i] -> GetYaxis() -> SetTitleOffset (1.1);
        dNdeta_1D_reco_multi[i] -> GetXaxis() -> CenterTitle(true);
        dNdeta_1D_reco_multi[i] -> GetYaxis() -> CenterTitle(true);
    }

    eta_Mbin_correction_tight = new TH2F(
        "",
        "eta_Mbin_correction_tight;#eta;Mbin;#alpha factor",
        included_eta_z_map.size(), 
        eta_z_ref -> GetXaxis() ->GetBinLowEdge(included_eta_z_map[0].first), 
        eta_z_ref -> GetXaxis() ->GetBinCenter(included_eta_z_map[included_eta_z_map.size() - 1].first) + eta_z_ref -> GetXaxis() ->GetBinWidth(included_eta_z_map[included_eta_z_map.size() - 1].first) / 2.,
        N_centrality_bin,
        0, 
        N_centrality_bin
    );

    eta_Mbin_correction_loose = new TH2F(
        "",
        "eta_Mbin_correction_loose;#eta;Mbin;#alpha factor",
        included_eta_z_map.size(), 
        eta_z_ref -> GetXaxis() ->GetBinLowEdge(included_eta_z_map[0].first), 
        eta_z_ref -> GetXaxis() ->GetBinCenter(included_eta_z_map[included_eta_z_map.size() - 1].first) + eta_z_ref -> GetXaxis() ->GetBinWidth(included_eta_z_map[included_eta_z_map.size() - 1].first) / 2.,
        N_centrality_bin,
        0, 
        N_centrality_bin
    );

    eta_Mbin_correction_loose_noUPC = new TH2F(
        "",
        "eta_Mbin_correction_loose_noUPC;#eta;Mbin;#alpha factor",
        included_eta_z_map.size(), 
        eta_z_ref -> GetXaxis() ->GetBinLowEdge(included_eta_z_map[0].first), 
        eta_z_ref -> GetXaxis() ->GetBinCenter(included_eta_z_map[included_eta_z_map.size() - 1].first) + eta_z_ref -> GetXaxis() ->GetBinWidth(included_eta_z_map[included_eta_z_map.size() - 1].first) / 2.,
        N_centrality_bin,
        0, 
        N_centrality_bin
    );

    input_eta_Mbin_correction_tight = nullptr;
    input_eta_Mbin_correction_loose = nullptr;

    cout<<"InitHist done"<<endl;
}




// note : to count the number of reco tracklets in the signal region with the selected eta bin and z bin, for each centrality bin
// note : for both methods
// note : to increase the statistic, try to stack the DeltaPhi distributions with same eta bin but different z bin
// note : the key is "MBin_EtaBin"
// note : for each eta bin and each centrality bin, count the number of event. (Each eta bin has different selected z region, so the number of event is different for each eta bin)
void EtaDistReader::MainPreparation()
{
    c1 -> cd();
    c1 -> Print(Form("%s/Reco_DeltaPhi_Multi_Stack.pdf(", out_folder_directory.c_str())); 
    for (int Mbin = 0; Mbin < N_centrality_bin; Mbin++) // note : centrality bin
    {
        for (auto eta_z : included_eta_z_map) // note : eta bin
        {
            DeltaPhi_Multi_Stack[Form("%i_%i",Mbin,eta_z.first)] = new THStack(Form("DeltaPhi_Multi_Stack_MBin%i_Eta%i",Mbin,eta_z.first),Form("DeltaPhi_Multi_Stack_MBin%i_Eta%i",Mbin,eta_z.first));
            for(int zbin : eta_z.second) // note :  Z bin
            {
                // note : the key to the map: "MBin_EtaBin"
                DeltaPhi_Multi_Stack[Form("%i_%i",Mbin,eta_z.first)] -> Add(DeltaPhi_Multi_1D[Form("%i_%i_%i",Mbin,eta_z.first,zbin)]);
                SignalNTrack_Single[Form("%i_%i",Mbin,eta_z.first)] += SignalNTrack_eta_z_Single_2D[Mbin] -> GetBinContent(eta_z.first, zbin); // note : the signal counting for the single method
                SignalNTrack_Multi[Form("%i_%i",Mbin,eta_z.first)]  +=  SignalNTrack_eta_z_Multi_2D[Mbin] -> GetBinContent(eta_z.first, zbin); // note : the signal counting for the multi method

                // note : to count the TrueNTtrack in the true level information
                dNdeta_1D_MC[Mbin]->SetBinContent(eta_z.first - eta_correction, TrueNtrack_eta_z_MC_2D[Mbin] -> GetBinContent(eta_z.first, zbin) + dNdeta_1D_MC[Mbin]->GetBinContent(eta_z.first - eta_correction));

                N_event_counting[Mbin][eta_z.first - eta_correction - 1]    += centrality_Z_map -> GetBinContent(Mbin + 1, zbin);
                N_event_counting_MC[Mbin][eta_z.first - eta_correction - 1] += centrality_Z_map_MC -> GetBinContent(Mbin + 1, zbin);
            }

            // note : Background fit on the DeltaPhi distributions, to know the background level
            // note : print the DeltaPhi distributions for each centrality bin adn each eta bin, after stacking
            temp_hist = (TH1F *) DeltaPhi_Multi_Stack[Form("%i_%i",Mbin,eta_z.first)] -> GetStack() -> Last();
            
            DeltaPhi_Multi_Stack_hist_out.push_back(temp_hist);
            DeltaPhi_Multi_Stack_hist_out[DeltaPhi_Multi_Stack_hist_out.size() - 1] -> SetTitle(Form("MBin%i_EtaBin%i",Mbin,eta_z.first));

            double hist_offset = get_dist_offset(temp_hist, 15);

            bkg_fit_pol2 -> SetParameters(hist_offset, 0, -0.2, 0, signal_region);
            bkg_fit_pol2 -> FixParameter(4, signal_region);
            bkg_fit_pol2 -> SetParLimits(2, -100, 0);
            temp_hist -> Fit(bkg_fit_pol2,"NQ");
            // note : extract the background region (which includes the signal region also)
            draw_pol2_line -> SetParameters(bkg_fit_pol2 -> GetParameter(0), bkg_fit_pol2 -> GetParameter(1), bkg_fit_pol2 -> GetParameter(2), bkg_fit_pol2 -> GetParameter(3));

            DeltaPhi_Multi_Stack[Form("%i_%i",Mbin,eta_z.first)] -> SetMinimum(0);
            DeltaPhi_Multi_Stack[Form("%i_%i",Mbin,eta_z.first)] -> SetMaximum( temp_hist -> GetBinContent(temp_hist -> GetMaximumBin()) * 1.5);
            DeltaPhi_Multi_Stack[Form("%i_%i",Mbin,eta_z.first)] -> Draw(); 

            draw_pol2_line -> Draw("lsame");
            
            double pol2_bkg_integral = fabs(draw_pol2_line -> Integral( -1. * signal_region, signal_region )) / temp_hist -> GetBinWidth(1);
            // cout<<i<<" "<<i1<<" pol2_bkg integral: "<<pol2_bkg_integral<<endl;

            // note : for the multi method, count the number of tracklets after the back subtraction
            dNdeta_1D_reco_multi[Mbin]  -> SetBinContent(eta_z.first - eta_correction, SignalNTrack_Multi[Form("%i_%i",Mbin,eta_z.first)] - pol2_bkg_integral);
            dNdeta_1D_reco_single[Mbin] -> SetBinContent(eta_z.first - eta_correction, SignalNTrack_Single[Form("%i_%i",Mbin,eta_z.first)]);

            coord_line -> DrawLine(-1*signal_region, 0, -1 * signal_region, temp_hist -> GetBinContent(temp_hist -> GetMaximumBin()) * 1.5);
            coord_line -> DrawLine(signal_region, 0, signal_region, temp_hist -> GetBinContent(temp_hist -> GetMaximumBin()) * 1.5);
            
            ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
            
            draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s, #eta: %.2f ~ %.2f",centrality_region[Mbin].c_str(), 
                eta_z_ref -> GetXaxis() -> GetBinCenter(eta_z.first) - eta_z_ref -> GetXaxis() -> GetBinWidth(eta_z.first)/2., 
                eta_z_ref -> GetXaxis() -> GetBinCenter(eta_z.first) + eta_z_ref -> GetXaxis() -> GetBinWidth(eta_z.first)/2.)
            );
            
            draw_text -> DrawLatex(0.21, 0.86, Form("MBin: %i, #eta_bin: %i, #Delta#Phi_bin width: %.2f", Mbin, eta_z.first, temp_hist -> GetBinWidth(1)));
            draw_text -> DrawLatex(0.21, 0.82, Form("pol2: %.2f + %.2f(x-%.2f) + %.2f(x-%.2f)^{2}", bkg_fit_pol2 -> GetParameter(0), bkg_fit_pol2 -> GetParameter(1), bkg_fit_pol2 -> GetParameter(3), bkg_fit_pol2 -> GetParameter(2), bkg_fit_pol2 -> GetParameter(3)));
            draw_text -> DrawLatex(0.21, 0.78, Form("Signal size: %i, pol2 bkg size: %.2f", SignalNTrack_Multi[Form("%i_%i",Mbin,eta_z.first)], pol2_bkg_integral));

            c1 -> Print(Form("%s/Reco_DeltaPhi_Multi_Stack.pdf", out_folder_directory.c_str()));
            c1 -> Clear();
        }
    }
    c1 -> Print(Form("%s/Reco_DeltaPhi_Multi_Stack.pdf)", out_folder_directory.c_str()));
    c1 -> Clear();

    cout<<"MainPreparation done"<<endl;
    return;
}

void EtaDistReader::FinaldNdEta()
{
    double MC_hist_counting = 0;
    double data_tight_counting = 0;
    double data_loose_counting = 0;

    c1 -> cd();
    c1 -> Print(Form("%s/dNdeta_combine_final_no_correction.pdf(", out_folder_directory.c_str()));
    for (int Mbin = 0; Mbin < N_centrality_bin; Mbin++)
    {   
        // // note : check the ratio of the bin error to the bin content of the three histograms
        // for (int i = 0; i < included_eta_z_map.size(); i++)
        // {
        //     cout<<" "<<endl;
        //     cout<<"1MBin: "<<Mbin<<" eta bin: "<<i<<" MC: "<<dNdeta_1D_MC[Mbin] -> GetBinError(i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(i+1)<<" reco_single: "<<dNdeta_1D_reco_single[Mbin] -> GetBinError(i+1) / dNdeta_1D_reco_single[Mbin] -> GetBinContent(i+1)<<" reco_multi: "<<dNdeta_1D_reco_multi[Mbin] -> GetBinError(i+1) / dNdeta_1D_reco_multi[Mbin] -> GetBinContent(i+1)<<endl;
        // }

        std::cout<<" "<<std::endl;

        dNdeta_1D_MC[Mbin]          -> Scale(1./double(dNdeta_1D_MC[Mbin]          -> GetBinWidth(1) ));
        dNdeta_1D_reco_single[Mbin] -> Scale(1./double(dNdeta_1D_reco_single[Mbin] -> GetBinWidth(1) ));
        dNdeta_1D_reco_multi[Mbin]  -> Scale(1./double(dNdeta_1D_reco_multi[Mbin]  -> GetBinWidth(1) ));

        // // note : check the ratio of the bin error to the bin content of the three histograms
        // for (int i = 0; i < included_eta_z_map.size(); i++)
        // {
        //     cout<<" "<<endl;
        //     cout<<"2MBin: "<<Mbin<<" eta bin: "<<i<<" MC: "<<dNdeta_1D_MC[Mbin] -> GetBinError(i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(i+1)<<" reco_single: "<<dNdeta_1D_reco_single[Mbin] -> GetBinError(i+1) / dNdeta_1D_reco_single[Mbin] -> GetBinContent(i+1)<<" reco_multi: "<<dNdeta_1D_reco_multi[Mbin] -> GetBinError(i+1) / dNdeta_1D_reco_multi[Mbin] -> GetBinContent(i+1)<<endl;
        // }

        for (int i = 0; i < included_eta_z_map.size(); i++)
        {
            dNdeta_1D_MC[Mbin] -> SetBinContent(i+1, dNdeta_1D_MC[Mbin] -> GetBinContent(i+1) / double(N_event_counting_MC[Mbin][i]));
            dNdeta_1D_MC[Mbin] -> SetBinError(i+1, dNdeta_1D_MC[Mbin] -> GetBinError(i+1) / double(N_event_counting_MC[Mbin][i]));

            dNdeta_1D_reco_single[Mbin] -> SetBinContent(i+1, dNdeta_1D_reco_single[Mbin] -> GetBinContent(i+1) / double(N_event_counting[Mbin][i]));
            dNdeta_1D_reco_single[Mbin] -> SetBinError(i+1, dNdeta_1D_reco_single[Mbin] -> GetBinError(i+1) / double(N_event_counting[Mbin][i]));

            dNdeta_1D_reco_multi[Mbin]  -> SetBinContent(i+1, dNdeta_1D_reco_multi[Mbin] -> GetBinContent(i+1) / double(N_event_counting[Mbin][i]));
            dNdeta_1D_reco_multi[Mbin]  -> SetBinError(i+1, dNdeta_1D_reco_multi[Mbin] -> GetBinError(i+1) / double(N_event_counting[Mbin][i]));
        
            if (input_eta_Mbin_correction_loose != nullptr)
            {
                dNdeta_1D_reco_multi[Mbin] -> SetBinContent(i+1, dNdeta_1D_reco_multi[Mbin] -> GetBinContent(i+1) / input_eta_Mbin_correction_loose -> GetBinContent(i+1, Mbin+1));
                dNdeta_1D_reco_multi[Mbin] -> SetBinError(i+1, dNdeta_1D_reco_multi[Mbin] -> GetBinError(i+1) /     input_eta_Mbin_correction_loose -> GetBinContent(i+1, Mbin+1));
            }

            if (input_eta_Mbin_correction_tight != nullptr)
            {
                dNdeta_1D_reco_single[Mbin] -> SetBinContent(i+1, dNdeta_1D_reco_single[Mbin] -> GetBinContent(i+1) / input_eta_Mbin_correction_tight -> GetBinContent(i+1, Mbin+1));
                dNdeta_1D_reco_single[Mbin] -> SetBinError(i+1,   dNdeta_1D_reco_single[Mbin] -> GetBinError(i+1)   / input_eta_Mbin_correction_tight -> GetBinContent(i+1, Mbin+1));
            }

        }
        
        cout<<"----------- for the case of method tight ----------------" <<endl;
        // note : to print the ratio between reco track and MC track
        std::cout << "centrality bin : "<<Mbin<<", ";
        for (int bin_i = 0; bin_i < dNdeta_1D_MC[Mbin] -> GetNbinsX(); bin_i++) {
            MC_hist_counting += dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1);
            data_tight_counting += dNdeta_1D_reco_single[Mbin] -> GetBinContent(bin_i+1);
            std::cout <<"~~"<<dNdeta_1D_reco_single[Mbin] -> GetBinContent(bin_i+1) <<", "<< dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1)<< ", " << Form("%.3f",dNdeta_1D_reco_single[Mbin] -> GetBinContent(bin_i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1)) <<"~~, ";

            eta_Mbin_correction_tight -> SetBinContent(bin_i + 1, Mbin+1, dNdeta_1D_reco_single[Mbin] -> GetBinContent(bin_i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1));
        }
        std::cout << std::endl;

        cout<<"----------- for the case of method inclusive ----------------" <<endl;
        // note : to print the ratio between reco track and MC track
        std::cout << "centrality bin : "<<Mbin<<", ";
        for (int bin_i = 0; bin_i < dNdeta_1D_MC[Mbin] -> GetNbinsX(); bin_i++) {
            data_loose_counting += dNdeta_1D_reco_multi[Mbin] -> GetBinContent(bin_i+1);
            std::cout <<"~~"<<dNdeta_1D_reco_multi[Mbin] -> GetBinContent(bin_i+1) <<", "<< dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1)<< ", " << Form( "%.3f", dNdeta_1D_reco_multi[Mbin] -> GetBinContent(bin_i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1)) <<"~~, ";

            eta_Mbin_correction_loose       -> SetBinContent(bin_i + 1, Mbin + 1, dNdeta_1D_reco_multi[Mbin] -> GetBinContent(bin_i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1));
            eta_Mbin_correction_loose_noUPC -> SetBinContent(bin_i + 1, Mbin + 1, dNdeta_1D_reco_multi[Mbin] -> GetBinContent(bin_i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(bin_i+1));
        }
        std::cout << std::endl;

        // note : check the bin content of the three histograms
        // note : and check the bin error of the three histograms
        // for (int i = 0; i < included_eta_z_map.size(); i++)
        // {
        //     cout<<" "<<endl;
        //     cout<<"3MBin: "<<Mbin<<" eta bin: "<<i<<" MC: "<<dNdeta_1D_MC[Mbin] -> GetBinContent(i+1)<<" reco_single: "<<dNdeta_1D_reco_single[Mbin] -> GetBinContent(i+1)<<" reco_multi: "<<dNdeta_1D_reco_multi[Mbin] -> GetBinContent(i+1)<<endl;
        //     cout<<"3MBin: "<<Mbin<<" eta bin: "<<i<<" MC: "<<dNdeta_1D_MC[Mbin] -> GetBinError(i+1)<<" reco_single: "<<dNdeta_1D_reco_single[Mbin] -> GetBinError(i+1)<<" reco_multi: "<<dNdeta_1D_reco_multi[Mbin] -> GetBinError(i+1)<<endl;
        // }

        // // note : check the ratio of the bin error to the bin content of the three histograms
        // for (int i = 0; i < included_eta_z_map.size(); i++)
        // {
        //     cout<<" "<<endl;
        //     cout<<"3MBin: "<<Mbin<<" eta bin: "<<i<<" MC: "<<dNdeta_1D_MC[Mbin] -> GetBinError(i+1) / dNdeta_1D_MC[Mbin] -> GetBinContent(i+1)<<" reco_single: "<<dNdeta_1D_reco_single[Mbin] -> GetBinError(i+1) / dNdeta_1D_reco_single[Mbin] -> GetBinContent(i+1)<<" reco_multi: "<<dNdeta_1D_reco_multi[Mbin] -> GetBinError(i+1) / dNdeta_1D_reco_multi[Mbin] -> GetBinContent(i+1)<<endl;
        // }

        if (Mbin == 0)
        {
            legend -> AddEntry(dNdeta_1D_MC[Mbin], "MC","f");
            legend -> AddEntry(dNdeta_1D_reco_single[Mbin], "Reco. Method tight","lep");
            legend -> AddEntry(dNdeta_1D_reco_multi[Mbin], "Reco. Method loose","lep");
        }


        dNdeta_1D_MC[Mbin] -> GetYaxis() -> SetRangeUser(0, dNdeta_1D_MC[Mbin] -> GetMaximum() * 1.5);

        dNdeta_1D_MC[Mbin] -> Draw("hist");
        dNdeta_1D_reco_single[Mbin] -> Draw("p same");
        dNdeta_1D_reco_multi[Mbin]  -> Draw("p same");

        legend -> Draw("same");

        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[Mbin].c_str()));
       
        c1 -> Print(Form("%s/dNdeta_combine_final_no_correction.pdf", out_folder_directory.c_str()));
        c1 -> Clear();

        cout<<"Mbin : "<<Mbin<<" MC_hist_counting: "<<MC_hist_counting<<" data_tight_counting: "<<data_tight_counting<<" data_loose_counting: "<<data_loose_counting<<endl;

        MC_hist_counting = 0;
        data_tight_counting = 0;
        data_loose_counting = 0;
    }
    c1 -> Print(Form("%s/dNdeta_combine_final_no_correction.pdf)", out_folder_directory.c_str()));
    c1 -> Clear();

    gStyle->SetPaintTextFormat("1.3f");

    c1 -> cd();
    eta_Mbin_correction_loose -> Draw("colz0");
    eta_Mbin_correction_loose -> SetMarkerSize(0.7);
    eta_Mbin_correction_loose -> Draw("HIST TEXT45 SAME");
    c1 -> Print(Form("%s/eta_Mbin_correction_loose.pdf", out_folder_directory.c_str()));
    c1 -> Clear();

    c1 -> cd();
    for (int i = 0; i < eta_Mbin_correction_loose_noUPC -> GetNbinsX(); i++) {eta_Mbin_correction_loose_noUPC -> SetBinContent(i+1, eta_Mbin_correction_loose_noUPC -> GetNbinsY()-1, 0);}
    eta_Mbin_correction_loose_noUPC -> Draw("colz0");
    eta_Mbin_correction_loose_noUPC -> SetMarkerSize(0.7);
    eta_Mbin_correction_loose_noUPC -> Draw("HIST TEXT45 SAME");
    c1 -> Print(Form("%s/eta_Mbin_correction_loose_noUPC.pdf", out_folder_directory.c_str()));
    c1 -> Clear();

    c1 -> cd();
    eta_Mbin_correction_tight -> Draw("colz0");
    eta_Mbin_correction_tight -> SetMarkerSize(0.7);
    eta_Mbin_correction_tight -> Draw("HIST TEXT45 SAME");
    c1 -> Print(Form("%s/eta_Mbin_correction_tight.pdf", out_folder_directory.c_str()));
    c1 -> Clear();

    cout<<"FinaldNdEta done"<<endl;
    return;
}

void EtaDistReader::EndRun()
{
    TFile * file_out = new TFile(Form("%s/alpha_correction_map.root",out_folder_directory.c_str()),"RECREATE");
    // file_out -> cd();
    eta_Mbin_correction_loose -> Write("eta_Mbin_correction_loose");
    eta_Mbin_correction_loose_noUPC -> Write("eta_Mbin_correction_loose_noUPC");
    eta_Mbin_correction_tight -> Write("eta_Mbin_correction_tight");

    file_out -> Close();
}

vector<TH1F *> EtaDistReader::GetDeltaPhi_Multi_stack_1D()
{
    return DeltaPhi_Multi_Stack_hist_out;
}
vector<TH1F *> EtaDistReader::GetdNdeta_1D_MC()
{
    return dNdeta_1D_MC;    
}
vector<TH1F *> EtaDistReader::GetdNdeta_1D_reco_single()
{
    return dNdeta_1D_reco_single;    
}
vector<TH1F *> EtaDistReader::GetdNdeta_1D_reco_multi()
{
    return dNdeta_1D_reco_multi;    
}

#endif