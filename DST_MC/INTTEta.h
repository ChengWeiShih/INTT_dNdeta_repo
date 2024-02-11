#ifndef INTTEta_h
#define INTTEta_h

#include "INTTXYvtxEvt.h"

class INTTEta : public INTTXYvtxEvt
{
    public:
        INTTEta(string run_type, string out_folder_directory, pair<double,double> beam_origin, int geo_mode_id, double phi_diff_cut = 0.11, pair<double, double> DCA_cut = {-1,1}, int N_clu_cutl = 20, int N_clu_cut = 10000, bool draw_event_display = true, double peek = 3.32405, double angle_diff_new_l = 0.0, double angle_diff_new_r = 3, bool print_message_opt = true):
        INTTXYvtxEvt(run_type, out_folder_directory, beam_origin, geo_mode_id, phi_diff_cut, DCA_cut, N_clu_cutl, N_clu_cut, draw_event_display, peek, angle_diff_new_l, angle_diff_new_r, print_message_opt, 2.5, 9)
        {
            track_cluster_ratio_1D.clear();
            dNdeta_1D.clear();
            track_eta_phi_2D.clear();
            track_eta_z_2D.clear();
            dNdeta_1D_MC.clear();
            final_dNdeta_1D.clear();
            final_dNdeta_1D_MC.clear();
            track_eta_z_2D_MC.clear();
            track_delta_phi_1D.clear();
            track_DCA_distance.clear();
            final_track_delta_phi_1D.clear();
            // for (int i = 0; i < 360; i++)
            // {
            //     inner_clu_phi_map[i].clear();
            //     outer_clu_phi_map[i].clear();
            // }

            std::fill(std::begin(inner_clu_phi_map), std::end(inner_clu_phi_map), std::vector<pair<bool,clu_info>>());
            std::fill(std::begin(outer_clu_phi_map), std::end(outer_clu_phi_map), std::vector<pair<bool,clu_info>>());

            // note : 10 is for the centrality bin
            N_GoodEvent_vec = vector<long long>(10,0); // todo : if the centrality bin is changed, the vector should be updated 

            InitHist();
            InitGraph();

            out_folder_directory_evt = out_folder_directory + "/evt_display";
            system(Form("mkdir %s", out_folder_directory_evt.c_str()));

            loop_NGoodPair = 0;
            evt_NTrack = 0;
            evt_NTrack_MC = 0;
            N_GoodEvent = 0;
            return;
        };
        
        void ProcessEvt(int event_i, vector<clu_info> temp_sPH_inner_nocolumn_vec, vector<clu_info> temp_sPH_outer_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_rz_vec, int NvtxMC, vector<double> TrigvtxMC, bool PhiCheckTag, Long64_t bco_full, pair<double,double> evt_z, int centrality_bin, vector<vector<float>> true_track_info); 
        // void ProcessEvtMC(int event_i, vector<float> true_track_info, vector<double> TrigvtxMC);
        void ClearEvt() override;
        void PrintPlots() override;
        void EndRun() override;

    protected:
        // note : {5, 15, 25, 35, 45, 55, 65, 75, 85, 95, inclusive}
        vector<TH1F *> track_cluster_ratio_1D;
        vector<TH1F *> dNdeta_1D;
        vector<TH2F *> track_eta_phi_2D; 
        vector<TH2F *> track_eta_z_2D; 

        vector<TH1F *> dNdeta_1D_MC;
        vector<TH2F *> track_eta_z_2D_MC; 

        vector<TH1F *> track_delta_eta_1D;
        vector<TH1F *> track_delta_eta_1D_post; 
        vector<TH1F *> track_delta_phi_1D;
        vector<TH1F *> track_DCA_distance;
        vector<TH2F *> track_phi_DCA_2D;

        vector<TH2F *> track_DeltaPhi_eta_2D;
        vector<TH2F *> track_DeltaPhi_DeltaEta_2D;

        TH2F * track_correlation_2D;
        vector<TH1F *> track_ratio_1D; // note : N reco track /N true track

        TH2F * reco_eta_correlation_2D;
        TH2F * reco_eta_diff_reco3P_2D;
        TH1F * reco_eta_diff_1D;

        vector<vector<TH1F *>> final_track_delta_phi_1D; // note : centrality_bin and different eta
        vector<TH1F *> final_dNdeta_1D;
        vector<TH1F *> final_dNdeta_1D_MC;

        TGraph * evt_reco_track_gr_All; // note : all the tracklets in one event, no selection
        TGraph * evt_reco_track_gr_PhiLoose; // note : the tracklets that pass the loose phi selection
        TGraph * evt_reco_track_gr_Z; // note : the tracklets that pass the z selection and the loose phi selection
        TGraph * evt_reco_track_gr_ZDCA; // note : the tracklets that pass the z and phi selection, and the loose phi selection
        TGraph * evt_reco_track_gr_ZDCAPhi; // note : the tracklets that pass the z, phi and DCA selection, and the loose phi selection
        TGraph * evt_true_track_gr; // note : the true tracklets in one event, no selection

        vector<pair<bool,clu_info>> inner_clu_phi_map[360]; // note: phi
        vector<pair<bool,clu_info>> outer_clu_phi_map[360]; // note: phi

        string out_folder_directory_evt;

        int draw_evt_cut = 15; // todo : change the cut value (number of true track)
        int print_plot_ratio = 3; // note : print the event display every 3 events
        int inner_NClu;
        int outer_NClu;
        int loop_NGoodPair;
        int evt_NTrack;
        int evt_NTrack_MC;
        long long N_GoodEvent;
        // double INTT_layer_R = 71.88; // note: the innermost, B0L0, layer radius
        double INTT_layer_R[4] = {71.88, 77.32, 96.8, 102.62}; // note : the radii of the INTT layers
        
        vector<long long> N_GoodEvent_vec;
        pair<double,double> Get_eta_pair;
        vector<vector<double>> final_eta_entry;

        TH2F * track_cluster_ratio_multiplicity_2D; // note : x : NClus / NTracks, y : ratio
        TGraphErrors * track_gr;

        // todo : if the centrality bin is changed, the following map and vector should be updated
        map<int,int> centrality_map = {
            {5, 0},
            {15, 1},
            {25, 2},
            {35, 3},
            {45, 4},
            {55, 5},
            {65, 6},
            {75, 7},
            {85, 8},
            {95, 9}
        };

        vector<string> centrality_region = {
            "0-5%",
            "5-15%",
            "15-25%",
            "25-35%",
            "35-45%",
            "45-55%",
            "55-65%",
            "65-75%",
            "75-85%",
            "85-95%",
            "0-95%"
        };

        TH1F * eta_region_hist;
        vector<double> eta_region = { // todo: if the eta region is changed, the following vector should be updated
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1,
            1.2,
            1.4,
            1.6,
            1.8,
            2,
            2.2
        };
        
        void InitHist();
        void InitGraph();

        double grEY_stddev(TGraphErrors * input_grr);
        pair<double, double> mirrorPolynomial(double a, double b);
        pair<double,double> Get_eta(vector<double>p0, vector<double>p1, vector<double>p2);
        double convertTo360(double radian);
        void print_evt_plot(int event_i, int NTrueTrack, int innerNClu, int outerNClu);
        double get_clu_eta(vector<double> vertex, vector<double> clu_pos);
        double get_dist_offset(TH1F * hist_in, int check_N_bin);
};

void INTTEta::InitHist()
{
    // double Eta_NBin = 145;
    // double Eta_Min = -2.9;
    // double Eta_Max = 2.9;

    double Eta_NBin = 61;
    double Eta_Min = -3.05;
    double Eta_Max = 3.05;

    final_track_delta_phi_1D = vector<vector<TH1F *>>(centrality_region.size());
    final_eta_entry = vector<vector<double>>(centrality_region.size());

    for (int i = 0; i < centrality_region.size(); i++) 
    {
        track_cluster_ratio_1D.push_back(new TH1F("","track_cluster_ratio_1D",200,0,15));
        track_cluster_ratio_1D[i] -> GetXaxis() -> SetTitle("NClus / NTracks");
        track_cluster_ratio_1D[i] -> GetYaxis() -> SetTitle("Entry");
        track_cluster_ratio_1D[i] -> GetXaxis() -> SetNdivisions(505);

        dNdeta_1D.push_back(new TH1F("","",Eta_NBin, Eta_Min, Eta_Max));
        dNdeta_1D[i] -> SetMarkerStyle(20);
        dNdeta_1D[i] -> SetMarkerSize(0.8);
        dNdeta_1D[i] -> SetMarkerColor(TColor::GetColor("#1A3947"));
        dNdeta_1D[i] -> SetLineColor(TColor::GetColor("#1A3947"));
        dNdeta_1D[i] -> SetLineWidth(2);
        dNdeta_1D[i] -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");
        dNdeta_1D[i] -> GetXaxis() -> SetTitle("#eta");
        dNdeta_1D[i] -> GetYaxis() -> SetRangeUser(0,50);
        dNdeta_1D[i] -> SetTitleSize(0.06, "X");
        dNdeta_1D[i] -> SetTitleSize(0.06, "Y");
        dNdeta_1D[i] -> GetXaxis() -> SetTitleOffset (0.71);
        dNdeta_1D[i] -> GetYaxis() -> SetTitleOffset (1.1);
        dNdeta_1D[i] -> GetXaxis() -> CenterTitle(true);
        dNdeta_1D[i] -> GetYaxis() -> CenterTitle(true);

        dNdeta_1D_MC.push_back(new TH1F("","",Eta_NBin, Eta_Min, Eta_Max));
        dNdeta_1D_MC[i] -> SetMarkerStyle(20);
        dNdeta_1D_MC[i] -> SetMarkerSize(0.8);
        dNdeta_1D_MC[i] -> SetMarkerColor(TColor::GetColor("#1A3947"));
        dNdeta_1D_MC[i] -> SetLineColor(TColor::GetColor("#1A3947"));
        dNdeta_1D_MC[i] -> SetLineWidth(2);
        dNdeta_1D_MC[i] -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");
        dNdeta_1D_MC[i] -> GetXaxis() -> SetTitle("#eta");
        dNdeta_1D_MC[i] -> GetYaxis() -> SetRangeUser(0,50);
        dNdeta_1D_MC[i] -> SetTitleSize(0.06, "X");
        dNdeta_1D_MC[i] -> SetTitleSize(0.06, "Y");
        dNdeta_1D_MC[i] -> GetXaxis() -> SetTitleOffset (0.71);
        dNdeta_1D_MC[i] -> GetYaxis() -> SetTitleOffset (1.1);
        dNdeta_1D_MC[i] -> GetXaxis() -> CenterTitle(true);
        dNdeta_1D_MC[i] -> GetYaxis() -> CenterTitle(true);

        final_dNdeta_1D.push_back(new TH1F("","",eta_region.size() - 1, &eta_region[0]));
        final_dNdeta_1D[i] -> SetMarkerStyle(20);
        final_dNdeta_1D[i] -> SetMarkerSize(0.8);
        final_dNdeta_1D[i] -> SetMarkerColor(TColor::GetColor("#1A3947"));
        final_dNdeta_1D[i] -> SetLineColor(TColor::GetColor("#1A3947"));
        final_dNdeta_1D[i] -> SetLineWidth(2);
        final_dNdeta_1D[i] -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");
        final_dNdeta_1D[i] -> GetXaxis() -> SetTitle("#eta");
        final_dNdeta_1D[i] -> GetYaxis() -> SetRangeUser(0,50);
        final_dNdeta_1D[i] -> SetTitleSize(0.06, "X");
        final_dNdeta_1D[i] -> SetTitleSize(0.06, "Y");
        final_dNdeta_1D[i] -> GetXaxis() -> SetTitleOffset (0.71);
        final_dNdeta_1D[i] -> GetYaxis() -> SetTitleOffset (1.1);
        final_dNdeta_1D[i] -> GetXaxis() -> CenterTitle(true);
        final_dNdeta_1D[i] -> GetYaxis() -> CenterTitle(true);

        final_dNdeta_1D_MC.push_back(new TH1F("","",eta_region.size() - 1, &eta_region[0]));
        final_dNdeta_1D_MC[i] -> SetMarkerStyle(20);
        final_dNdeta_1D_MC[i] -> SetMarkerSize(0.8);
        final_dNdeta_1D_MC[i] -> SetMarkerColor(TColor::GetColor("#1A3947"));
        final_dNdeta_1D_MC[i] -> SetLineColor(TColor::GetColor("#1A3947"));
        final_dNdeta_1D_MC[i] -> SetLineWidth(2);
        final_dNdeta_1D_MC[i] -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");
        final_dNdeta_1D_MC[i] -> GetXaxis() -> SetTitle("#eta");
        final_dNdeta_1D_MC[i] -> GetYaxis() -> SetRangeUser(0,50);
        final_dNdeta_1D_MC[i] -> SetTitleSize(0.06, "X");
        final_dNdeta_1D_MC[i] -> SetTitleSize(0.06, "Y");
        final_dNdeta_1D_MC[i] -> GetXaxis() -> SetTitleOffset (0.71);
        final_dNdeta_1D_MC[i] -> GetYaxis() -> SetTitleOffset (1.1);
        final_dNdeta_1D_MC[i] -> GetXaxis() -> CenterTitle(true);
        final_dNdeta_1D_MC[i] -> GetYaxis() -> CenterTitle(true);

        track_eta_phi_2D.push_back(new TH2F("","track_eta_phi_2D", 200, 0, 360, Eta_NBin, Eta_Min, Eta_Max));
        track_eta_phi_2D[i] -> GetXaxis() -> SetTitle("tracklet #phi");
        track_eta_phi_2D[i] -> GetYaxis() -> SetTitle("tracklet #eta");
        track_eta_phi_2D[i] -> GetXaxis() -> SetNdivisions(505);

        track_eta_z_2D.push_back(new TH2F("","track_eta_z_2D", 300, -450, 0, Eta_NBin, Eta_Min, Eta_Max));
        track_eta_z_2D[i] -> GetXaxis() -> SetTitle("z vertex [mm]");
        track_eta_z_2D[i] -> GetYaxis() -> SetTitle("tracklet #eta");
        track_eta_z_2D[i] -> GetXaxis() -> SetNdivisions(505);

        track_eta_z_2D_MC.push_back(new TH2F("","track_eta_z_2D", 300, -450, 0, Eta_NBin, Eta_Min, Eta_Max));
        track_eta_z_2D_MC[i] -> GetXaxis() -> SetTitle("z vertex [mm]");
        track_eta_z_2D_MC[i] -> GetYaxis() -> SetTitle("tracklet #eta");
        track_eta_z_2D_MC[i] -> GetXaxis() -> SetNdivisions(505);

        track_delta_eta_1D.push_back(new TH1F("","track_delta_eta_1D", 200, -1.5, 1.5));
        track_delta_eta_1D[i] -> GetXaxis() -> SetTitle("Clus #Delta#eta");
        track_delta_eta_1D[i] -> GetYaxis() -> SetTitle("Entry");
        track_delta_eta_1D[i] -> GetXaxis() -> SetNdivisions(505);

        track_delta_eta_1D_post.push_back(new TH1F("","track_delta_eta_1D_ potrack_delta_eta_1D_post", 200, -1.5, 1.5));
        track_delta_eta_1D_post[i] -> GetXaxis() -> SetTitle("Clus #Delta#eta");
        track_delta_eta_1D_post[i] -> GetYaxis() -> SetTitle("Entry");
        track_delta_eta_1D_post[i] -> GetXaxis() -> SetNdivisions(505);

        track_delta_phi_1D.push_back(new TH1F("","track_delta_phi_1D", 200, -1.5, 1.5));
        track_delta_phi_1D[i] -> GetXaxis() -> SetTitle("Clus #Delta#phi [degree]");
        track_delta_phi_1D[i] -> GetYaxis() -> SetTitle("Entry");
        track_delta_phi_1D[i] -> GetXaxis() -> SetNdivisions(505);

        track_DCA_distance.push_back(new TH1F("","track_DCA_distance", 200, -5, 5));
        track_DCA_distance[i] -> GetXaxis() -> SetTitle("DCA distance [mm]");
        track_DCA_distance[i] -> GetYaxis() -> SetTitle("Entry");
        track_DCA_distance[i] -> GetXaxis() -> SetNdivisions(505);

        track_phi_DCA_2D.push_back(new TH2F("","track_phi_DCA_2D", 200, -1.5, 1.5, 200, -3, 3));
        track_phi_DCA_2D[i] -> GetXaxis() -> SetTitle("Clus #Delta#phi [degree]");
        track_phi_DCA_2D[i] -> GetYaxis() -> SetTitle("DCA distance [mm]");
        track_phi_DCA_2D[i] -> GetXaxis() -> SetNdivisions(505);

        track_ratio_1D.push_back(new TH1F("","track_ratio_1D", 200, 0, 1.5));
        track_ratio_1D[i] -> GetXaxis() -> SetTitle("N reco track / N true track");
        track_ratio_1D[i] -> GetYaxis() -> SetTitle("Entry");
        track_ratio_1D[i] -> GetXaxis() -> SetNdivisions(505);

        track_DeltaPhi_eta_2D.push_back(new TH2F("","track_DeltaPhi_eta_2D", 200, -1.5, 1.5, 200, -4, 4));
        track_DeltaPhi_eta_2D[i] -> GetXaxis() -> SetTitle("Clus #Delta#phi [degree]");
        track_DeltaPhi_eta_2D[i] -> GetYaxis() -> SetTitle("Tracklet #eta");
        track_DeltaPhi_eta_2D[i] -> GetXaxis() -> SetNdivisions(505);

        track_DeltaPhi_DeltaEta_2D.push_back(new TH2F("","track_DeltaPhi_DeltaEta_2D", 200, -1.5, 1.5, 200, -1, 1));
        track_DeltaPhi_DeltaEta_2D[i] -> GetXaxis() -> SetTitle("Clus #Delta#phi [degree]");
        track_DeltaPhi_DeltaEta_2D[i] -> GetYaxis() -> SetTitle("#Delta#eta");
        track_DeltaPhi_DeltaEta_2D[i] -> GetXaxis() -> SetNdivisions(505);

        final_track_delta_phi_1D[i].clear();
        for (int eta_i = 0; eta_i < eta_region.size() - 1; eta_i++)
        {
            final_track_delta_phi_1D[i].push_back(new TH1F("","final_track_delta_phi_1D", 75, -1.5, 1.5));
            final_track_delta_phi_1D[i][eta_i] -> GetXaxis() -> SetTitle("Clus #Delta#phi [degree]");
            final_track_delta_phi_1D[i][eta_i] -> GetYaxis() -> SetTitle("Entry");
            final_track_delta_phi_1D[i][eta_i] -> GetXaxis() -> SetNdivisions(505);
        }
    }
    
    track_cluster_ratio_multiplicity_2D = new TH2F("","track_cluster_ratio_multiplicity_2D", 200, 0, 9000, 200, 0, 15);
    track_cluster_ratio_multiplicity_2D -> GetXaxis() -> SetTitle("NClus");
    track_cluster_ratio_multiplicity_2D -> GetYaxis() -> SetTitle("NClus / NTracks");
    track_cluster_ratio_multiplicity_2D -> GetXaxis() -> SetNdivisions(505);

    track_correlation_2D = new TH2F("","track_correlation_2D", 200, 0, 3000, 200, 0, 3000);
    track_correlation_2D -> GetXaxis() -> SetTitle("N true track");
    track_correlation_2D -> GetYaxis() -> SetTitle("N reco track");
    track_correlation_2D -> GetXaxis() -> SetNdivisions(505);

    reco_eta_correlation_2D = new TH2F("","reco_eta_correlation_2D", 200, -3, 3, 200, -3, 3);
    reco_eta_correlation_2D -> GetXaxis() -> SetTitle("Reco 3P #eta");
    reco_eta_correlation_2D -> GetYaxis() -> SetTitle("Reco 2P #eta");
    reco_eta_correlation_2D -> GetXaxis() -> SetNdivisions(505);

    reco_eta_diff_reco3P_2D = new TH2F("","reco_eta_diff_reco3P_2D", 200, -3, 3, 200, -0.35, 0.35);
    reco_eta_diff_reco3P_2D -> GetXaxis() -> SetTitle("Reco 3P #eta");
    reco_eta_diff_reco3P_2D -> GetYaxis() -> SetTitle("Reco #eta diff");
    reco_eta_diff_reco3P_2D -> GetXaxis() -> SetNdivisions(505);

    reco_eta_diff_1D = new TH1F("","reco_eta_diff_1D", 200, -0.35, 0.35);
    reco_eta_diff_1D -> GetXaxis() -> SetTitle("Reco #eta diff");
    reco_eta_diff_1D -> GetYaxis() -> SetTitle("Entry");
    reco_eta_diff_1D -> GetXaxis() -> SetNdivisions(505);

    eta_region_hist = new TH1F("","", eta_region.size() - 1, &eta_region[0]);

    cout<<" running ? in INTTEta, InitHist"<<endl;

}

void INTTEta::InitGraph()
{
    track_gr = new TGraphErrors();

    evt_reco_track_gr_All = new TGraph();
    evt_reco_track_gr_All -> SetMarkerStyle(5);
    evt_reco_track_gr_All -> SetMarkerSize(1);
    evt_reco_track_gr_All -> SetMarkerColor(4);

    evt_reco_track_gr_PhiLoose = new TGraph();
    evt_reco_track_gr_PhiLoose -> SetMarkerStyle(5);
    evt_reco_track_gr_PhiLoose -> SetMarkerSize(1);
    evt_reco_track_gr_PhiLoose -> SetMarkerColor(4);

    evt_reco_track_gr_Z = new TGraph();
    evt_reco_track_gr_Z -> SetMarkerStyle(5);
    evt_reco_track_gr_Z -> SetMarkerSize(1);
    evt_reco_track_gr_Z -> SetMarkerColor(4);
    
    evt_reco_track_gr_ZDCA = new TGraph();
    evt_reco_track_gr_ZDCA -> SetMarkerStyle(5);
    evt_reco_track_gr_ZDCA -> SetMarkerSize(1);
    evt_reco_track_gr_ZDCA -> SetMarkerColor(4);
    
    evt_reco_track_gr_ZDCAPhi = new TGraph();
    evt_reco_track_gr_ZDCAPhi -> SetMarkerStyle(5);
    evt_reco_track_gr_ZDCAPhi -> SetMarkerSize(1);
    evt_reco_track_gr_ZDCAPhi -> SetMarkerColor(4);

    evt_true_track_gr = new TGraph();
    evt_true_track_gr -> SetMarkerStyle(20);
    evt_true_track_gr -> SetMarkerSize(1);
    // evt_true_track_gr -> SetMarkerColor(2);
    evt_true_track_gr -> SetMarkerColorAlpha(2,0.5);


    cout<<" running ? "<<endl;
    return;
}

void INTTEta::print_evt_plot(int event_i, int NTrueTrack, int innerNClu, int outerNClu)
{   
    c1 -> Clear();
    c1 -> Print( Form("%s/evt_track_eID_%i.pdf(", out_folder_directory_evt.c_str(), event_i) );
    c1 -> cd();

    evt_reco_track_gr_All -> GetXaxis() -> SetTitle("#eta");
    evt_reco_track_gr_All -> GetYaxis() -> SetTitle("#phi");
    evt_reco_track_gr_All -> GetXaxis() -> SetLimits(-3.5, 3.5);
    evt_reco_track_gr_All -> GetYaxis() -> SetRangeUser(-10, 420);
    evt_reco_track_gr_All -> GetXaxis() -> SetNdivisions(505);
    evt_reco_track_gr_All -> Draw("ap");
    evt_true_track_gr -> Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    draw_text -> DrawLatex(0.21, 0.90, Form("NTrueTrack: %i, innerNClu: %i, outerNClu: %i", NTrueTrack, innerNClu, outerNClu));
    draw_text -> DrawLatex(0.21, 0.85, Form("NReco_tracklet: %i", evt_reco_track_gr_All->GetN()));
    coord_line -> DrawLine(-3.5, 0, 3.5, 0);
    coord_line -> DrawLine(-3.5, 360, 3.5, 360);
    c1 -> Print( Form("%s/evt_track_eID_%i.pdf", out_folder_directory_evt.c_str(), event_i) );
    c1 -> Clear();

    evt_reco_track_gr_PhiLoose -> GetXaxis() -> SetTitle("#eta");
    evt_reco_track_gr_PhiLoose -> GetYaxis() -> SetTitle("#phi");
    evt_reco_track_gr_PhiLoose -> GetXaxis() -> SetLimits(-3.5, 3.5);
    evt_reco_track_gr_PhiLoose -> GetYaxis() -> SetRangeUser(-10, 420);
    evt_reco_track_gr_PhiLoose -> GetXaxis() -> SetNdivisions(505);
    evt_reco_track_gr_PhiLoose -> Draw("ap");
    evt_true_track_gr -> Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    draw_text -> DrawLatex(0.21, 0.90, Form("NTrueTrack: %i, innerNClu: %i, outerNClu: %i", NTrueTrack, innerNClu, outerNClu));
    draw_text -> DrawLatex(0.21, 0.85, Form("NReco_tracklet: %i", evt_reco_track_gr_PhiLoose->GetN()));
    coord_line -> DrawLine(-3.5, 0, 3.5, 0);
    coord_line -> DrawLine(-3.5, 360, 3.5, 360);
    c1 -> Print( Form("%s/evt_track_eID_%i.pdf", out_folder_directory_evt.c_str(), event_i) );
    c1 -> Clear();

    evt_reco_track_gr_Z -> GetXaxis() -> SetTitle("#eta");
    evt_reco_track_gr_Z -> GetYaxis() -> SetTitle("#phi");
    evt_reco_track_gr_Z -> GetXaxis() -> SetLimits(-3.5, 3.5);
    evt_reco_track_gr_Z -> GetYaxis() -> SetRangeUser(-10, 420);
    evt_reco_track_gr_Z -> GetXaxis() -> SetNdivisions(505);
    evt_reco_track_gr_Z -> Draw("ap");
    evt_true_track_gr -> Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    draw_text -> DrawLatex(0.21, 0.90, Form("NTrueTrack: %i, innerNClu: %i, outerNClu: %i", NTrueTrack, innerNClu, outerNClu));
    draw_text -> DrawLatex(0.21, 0.85, Form("NReco_tracklet: %i", evt_reco_track_gr_Z->GetN()));
    coord_line -> DrawLine(-3.5, 0, 3.5, 0);
    coord_line -> DrawLine(-3.5, 360, 3.5, 360);
    c1 -> Print( Form("%s/evt_track_eID_%i.pdf", out_folder_directory_evt.c_str(), event_i) );
    c1 -> Clear();

    evt_reco_track_gr_ZDCA -> GetXaxis() -> SetTitle("#eta");
    evt_reco_track_gr_ZDCA -> GetYaxis() -> SetTitle("#phi");
    evt_reco_track_gr_ZDCA -> GetXaxis() -> SetLimits(-3.5, 3.5);
    evt_reco_track_gr_ZDCA -> GetYaxis() -> SetRangeUser(-10, 420);
    evt_reco_track_gr_ZDCA -> GetXaxis() -> SetNdivisions(505);
    evt_reco_track_gr_ZDCA -> Draw("ap");
    evt_true_track_gr -> Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    draw_text -> DrawLatex(0.21, 0.90, Form("NTrueTrack: %i, innerNClu: %i, outerNClu: %i", NTrueTrack, innerNClu, outerNClu));
    draw_text -> DrawLatex(0.21, 0.85, Form("NReco_tracklet: %i", evt_reco_track_gr_ZDCA->GetN()));
    coord_line -> DrawLine(-3.5, 0, 3.5, 0);
    coord_line -> DrawLine(-3.5, 360, 3.5, 360);
    c1 -> Print( Form("%s/evt_track_eID_%i.pdf", out_folder_directory_evt.c_str(), event_i) );
    c1 -> Clear();

    evt_reco_track_gr_ZDCAPhi -> GetXaxis() -> SetTitle("#eta");
    evt_reco_track_gr_ZDCAPhi -> GetYaxis() -> SetTitle("#phi");
    evt_reco_track_gr_ZDCAPhi -> GetXaxis() -> SetLimits(-3.5, 3.5);
    evt_reco_track_gr_ZDCAPhi -> GetYaxis() -> SetRangeUser(-10, 420);
    evt_reco_track_gr_ZDCAPhi -> GetXaxis() -> SetNdivisions(505);
    evt_reco_track_gr_ZDCAPhi -> Draw("ap");
    evt_true_track_gr -> Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    draw_text -> DrawLatex(0.21, 0.90, Form("NTrueTrack: %i, innerNClu: %i, outerNClu: %i", NTrueTrack, innerNClu, outerNClu));
    draw_text -> DrawLatex(0.21, 0.85, Form("NReco_tracklet: %i", evt_reco_track_gr_ZDCAPhi->GetN()));
    coord_line -> DrawLine(-3.5, 0, 3.5, 0);
    coord_line -> DrawLine(-3.5, 360, 3.5, 360);
    c1 -> Print( Form("%s/evt_track_eID_%i.pdf)", out_folder_directory_evt.c_str(), event_i) );
    c1 -> Clear();
}

// note : this function is for the event by event vertex calculation
void INTTEta::ProcessEvt(int event_i, vector<clu_info> temp_sPH_inner_nocolumn_vec, vector<clu_info> temp_sPH_outer_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_rz_vec, int NvtxMC, vector<double> TrigvtxMC, bool PhiCheckTag, Long64_t bco_full, pair<double,double> evt_z, int centrality_bin, vector<vector<float>> true_track_info) // note : evt_z : {z, width} 
{   
    return_tag = 0;

    if (event_i%1 == 0) {cout<<"In INTTEta class, running event : "<<event_i<<endl;}

    inner_NClu = temp_sPH_inner_nocolumn_vec.size();
    outer_NClu = temp_sPH_outer_nocolumn_vec.size();
    total_NClus = inner_NClu + outer_NClu;

    // cout<<"test_0"<<endl;
    if (total_NClus < zvtx_cal_require) {return; cout<<"return confirmation"<<endl;}   
    
    if (run_type == "MC" && NvtxMC != 1) { return; cout<<"In INTTEta class, event : "<<event_i<<" Nvtx : "<<NvtxMC<<" Nvtx more than one "<<endl;}
    if (PhiCheckTag == false)            { return; cout<<"In INTTEta class, event : "<<event_i<<" Nvtx : "<<NvtxMC<<" Not full phi has hits "<<endl;}
    
    if (inner_NClu < 10 || outer_NClu < 10 || total_NClus > N_clu_cut || total_NClus < N_clu_cutl)
    {
        return;
        printf("In INTTEta class, event : %i, low clu continue, NClus : %lu \n", event_i, total_NClus); 
    }

    // todo : the z vertex range is here
    if (-220 > evt_z.first + evt_z.second || -180 < evt_z.first - evt_z.second) {return;}
    
    
    N_GoodEvent += 1;
    N_GoodEvent_vec[centrality_map[centrality_bin]] += 1;

    // cout<<"N inner cluster : "<<inner_NClu<<" N outer cluster : "<<outer_NClu<<endl;

    if (run_type == "MC")
    {
        // note : for the true track case ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        double INTT_eta_acceptance_l = -0.5 * TMath::Log((sqrt(pow(-230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))-(-230.-TrigvtxMC[2]*10.)) / (sqrt(pow(-230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))+(-230.-TrigvtxMC[2]*10.))); // note : left
        double INTT_eta_acceptance_r =  -0.5 * TMath::Log((sqrt(pow(230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))-(230.-TrigvtxMC[2]*10.)) / (sqrt(pow(230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))+(230.-TrigvtxMC[2]*10.))); // note : right
        // if (event_i % 100 == 0){cout<<"z : "<<TrigvtxMC[2]*10.<<" eta : "<<INTT_eta_acceptance_l<<" "<<INTT_eta_acceptance_r<<endl;}

        // cout<<"true_track_info : "<<true_track_info.size()<<endl;    
        for (int track_i = 0; track_i < true_track_info.size(); track_i++)
        {
            if (true_track_info[track_i][2] == 111 || true_track_info[track_i][2] == 22 || abs(true_track_info[track_i][2]) == 2112){continue;}

            if (true_track_info[track_i][0] > INTT_eta_acceptance_l && true_track_info[track_i][0] < INTT_eta_acceptance_r)
            {
                dNdeta_1D_MC[centrality_map[centrality_bin]] -> Fill(true_track_info[track_i][0]);
                dNdeta_1D_MC[dNdeta_1D_MC.size() - 1]        -> Fill(true_track_info[track_i][0]);   
                final_dNdeta_1D_MC[centrality_map[centrality_bin]] -> Fill(true_track_info[track_i][0]);
                final_dNdeta_1D_MC[dNdeta_1D_MC.size() - 1]        -> Fill(true_track_info[track_i][0]);
                track_eta_z_2D_MC[centrality_map[centrality_bin]] -> Fill(TrigvtxMC[2]*10., true_track_info[track_i][0]);
                track_eta_z_2D_MC[track_eta_z_2D_MC.size() - 1]   -> Fill(TrigvtxMC[2]*10., true_track_info[track_i][0]);    

                // cout<<"true track eta : "<<true_track_info[track_i][0]<<" phi : "<<convertTo360(true_track_info[track_i][1])<<endl;
                // cout<<"("<<true_track_info[track_i][0]<<", "<<convertTo360(true_track_info[track_i][1])<<"), ";

                evt_true_track_gr -> SetPoint(evt_true_track_gr->GetN(),true_track_info[track_i][0], convertTo360(true_track_info[track_i][1]));

                evt_NTrack_MC += 1;
            }
        }
        // if (evt_NTrack_MC < 10) {cout<<"evt : "<<event_i<<" ---- N reco track : "<<evt_NTrack<<" N true track : "<<evt_NTrack_MC<<" ratio : "<<double(evt_NTrack) / double(evt_NTrack_MC)<<endl;}
    }


    // cout<<"test_1"<<endl;
    for ( int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++ )
    {    
        loop_NGoodPair = 0;
        
        for ( int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++ )
        {
            // cout<<event_i<<" test_5 "<<inner_i<<" "<<outer_i<<endl;
            // note : calculate the cluster phi after the vertex correction which can enhence the purity of the tracklet selection
            Clus_InnerPhi_Offset = (temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second < 0) ? atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi());
            Clus_OuterPhi_Offset = (temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second < 0) ? atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi());
            double delta_phi = Clus_InnerPhi_Offset - Clus_OuterPhi_Offset;
            double track_phi = (Clus_InnerPhi_Offset + Clus_OuterPhi_Offset)/2.;
            double inner_clu_eta = get_clu_eta({beam_origin.first, beam_origin.second, evt_z.first},{temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y, temp_sPH_inner_nocolumn_vec[inner_i].z});
            double outer_clu_eta = get_clu_eta({beam_origin.first, beam_origin.second, evt_z.first},{temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y, temp_sPH_outer_nocolumn_vec[outer_i].z});
            double delta_eta = inner_clu_eta - outer_clu_eta;
            
            // cout<<"test_2"<<endl;
            // note : before calculating the possible z vertex range of one tracklet, the vertex correction was applied
            pair<double,double> z_range_info = Get_possible_zvtx( 
                0., // get_radius(beam_origin.first,beam_origin.second), 
                {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z}, // note : unsign radius
                {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}  // note : unsign radius
            );

            if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0)
            {
                Get_eta_pair = Get_eta(
                    {get_radius(beam_origin.first,beam_origin.second), evt_z.first,evt_z.second},
                    {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z},
                    {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}
                );

                evt_reco_track_gr_All -> SetPoint(evt_reco_track_gr_All->GetN(),Get_eta_pair.second, track_phi);
            }

            
            // if ( Get_eta_pair.second > 1.4 && Get_eta_pair.second < 1.7 && track_phi > 50 && track_phi < 60 ) {
            //     cout<<" "<<endl;
            //     cout<<"reco eta : "<<Get_eta_pair.second<<" reco phi : "<<track_phi<<" delta eta : "<<delta_eta<<endl;
            //     cout<<"all pair info, inner clu pos "<<temp_sPH_inner_nocolumn_vec[inner_i].x<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].y<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].z<<" phi : "<<Clus_InnerPhi_Offset<<endl;
            //     cout<<"all pair info, outer clu pos "<<temp_sPH_outer_nocolumn_vec[outer_i].x<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].y<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].z<<" phi : "<<Clus_OuterPhi_Offset<<endl;
            //     cout<<"z vertex range : "<<z_range_info.first - z_range_info.second<<" "<<z_range_info.first + z_range_info.second<<endl;
            //     cout<<"accept z range : "<<evt_z.first - evt_z.second<<" "<<evt_z.first + evt_z.second<<endl;
                
            // }

            if (fabs(delta_phi) > 5.72) {continue;}
            if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_PhiLoose -> SetPoint(evt_reco_track_gr_PhiLoose->GetN(),Get_eta_pair.second, track_phi);}

            track_delta_eta_1D[centrality_map[centrality_bin]] -> Fill( delta_eta );
            track_delta_eta_1D[track_delta_eta_1D.size() - 1]  -> Fill( delta_eta );

            track_DeltaPhi_DeltaEta_2D[centrality_map[centrality_bin]]        -> Fill(delta_phi, delta_eta);
            track_DeltaPhi_DeltaEta_2D[track_DeltaPhi_DeltaEta_2D.size() - 1] -> Fill(delta_phi, delta_eta);

            // cout<<event_i<<" test_6 "<<inner_i<<" "<<outer_i<<endl;
            // note : this is a cut to constraint on the z vertex, only if the tracklets with the range that covers the z vertex can pass this cut
            if (z_range_info.first - z_range_info.second > evt_z.first + evt_z.second || z_range_info.first + z_range_info.second < evt_z.first - evt_z.second) {continue;}
            // cout<<"range test: "<<z_range_info.first - z_range_info.second<<" "<<z_range_info.first + z_range_info.second<<" vertex: "<<evt_z.first<<" "<<evt_z.second<<endl;
            // cout<<"test_3"<<endl;

            if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_Z -> SetPoint(evt_reco_track_gr_Z->GetN(),Get_eta_pair.second, track_phi);}
            
            // note : the reason to use the DCA calculation with sign is because that the distribution of 1D DCA will be single peak only
            double DCA_sign = calculateAngleBetweenVectors(
                temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y,
                temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y,
                beam_origin.first, beam_origin.second
            );

            // vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
            //     temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y,
            //     temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y,
            //     beam_origin.first, beam_origin.second
            // );
            // if (DCA_info_vec[0] != fabs(DCA_sign) && fabs( DCA_info_vec[0] - fabs(DCA_sign) ) > 0.1 ){
            //     cout<<"different DCA : "<<DCA_info_vec[0]<<" "<<DCA_sign<<" diff : "<<DCA_info_vec[0] - fabs(DCA_sign)<<endl;
            // }

            track_delta_eta_1D_post[centrality_map[centrality_bin]]     -> Fill( delta_eta );
            track_delta_eta_1D_post[track_delta_eta_1D_post.size() - 1] -> Fill( delta_eta );

            track_delta_phi_1D[centrality_map[centrality_bin]] -> Fill( delta_phi );
            track_delta_phi_1D[track_delta_phi_1D.size() - 1]  -> Fill( delta_phi );
            
            track_DCA_distance[centrality_map[centrality_bin]] -> Fill( DCA_sign );
            track_DCA_distance[track_DCA_distance.size() - 1]  -> Fill( DCA_sign ); 

            track_phi_DCA_2D[centrality_map[centrality_bin]] -> Fill( delta_phi, DCA_sign );
            track_phi_DCA_2D[track_phi_DCA_2D.size() - 1]    -> Fill( delta_phi, DCA_sign );

            // cout<<" "<<endl;
            // cout<<"inner_i : "<<inner_i<<" outer_i "<<outer_i<<" true z : "<<TrigvtxMC[2]*10.<<" reco z: "<<evt_z.first<<endl;
            // cout<<"all pair info, inner clu pos "<<temp_sPH_inner_nocolumn_vec[inner_i].x<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].y<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].z<<" phi : "<<Clus_InnerPhi_Offset<<endl;
            // cout<<"all pair info, outer clu pos "<<temp_sPH_outer_nocolumn_vec[outer_i].x<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].y<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].z<<" phi : "<<Clus_OuterPhi_Offset<<endl;
            // if (fabs(delta_phi) < 4) {cout<<"("<<Get_eta_pair.second<<", "<<track_phi<<"), ";}

            // if (DCA_cut.first < DCA_sign && DCA_sign < DCA_cut.second)
            if (true)
            {
                if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_ZDCA -> SetPoint(evt_reco_track_gr_ZDCA->GetN(),Get_eta_pair.second, track_phi);}

                if (fabs(delta_phi) < phi_diff_cut)
                {
                    if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_ZDCAPhi -> SetPoint(evt_reco_track_gr_ZDCAPhi->GetN(),Get_eta_pair.second, track_phi);}    
                    // cout<<"reco eta : "<<Get_eta_pair.second<<" reco phi : "<<track_phi<<endl;

                    // if (fabs(Get_eta_pair.second) < 0.1) { cout<<"eta "<<Get_eta_pair.second<<", three points : "<<get_radius(beam_origin.first,beam_origin.second)<<" "<<evt_z.first<<" "<<evt_z.second<<" "<<
                    // get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second)<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].z<<" "<<
                    // get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second)<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].z<<endl;}

                    Get_eta_pair = Get_eta(
                        {get_radius(beam_origin.first,beam_origin.second), evt_z.first,evt_z.second},
                        {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z},
                        {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}
                    );

                    double eta_bin = eta_region_hist -> Fill(Get_eta_pair.second);
                    // cout<<"test1 "<<eta_bin<<endl;
                    if (eta_bin != -1)
                    {
                        final_track_delta_phi_1D[centrality_map[centrality_bin]][eta_bin - 1]      -> Fill(delta_phi);
                        final_track_delta_phi_1D[final_track_delta_phi_1D.size() - 1][eta_bin - 1] -> Fill(delta_phi);
                    }
                    

                    reco_eta_correlation_2D -> Fill(Get_eta_pair.second, (inner_clu_eta + outer_clu_eta)/2.);
                    reco_eta_diff_reco3P_2D -> Fill(Get_eta_pair.second, Get_eta_pair.second - (inner_clu_eta + outer_clu_eta)/2.);
                    reco_eta_diff_1D        -> Fill(Get_eta_pair.second - (inner_clu_eta + outer_clu_eta)/2.);

                    track_DeltaPhi_eta_2D[centrality_map[centrality_bin]]   -> Fill(delta_phi, Get_eta_pair.second);
                    track_DeltaPhi_eta_2D[track_DeltaPhi_eta_2D.size() - 1] -> Fill(delta_phi, Get_eta_pair.second);

                    // cout<<"test_5"<<endl;
                    dNdeta_1D[centrality_map[centrality_bin]] -> Fill(Get_eta_pair.second);
                    dNdeta_1D[dNdeta_1D.size() - 1]           -> Fill(Get_eta_pair.second);

                    track_eta_phi_2D[centrality_map[centrality_bin]] -> Fill(track_phi, Get_eta_pair.second);
                    track_eta_phi_2D[track_eta_phi_2D.size() - 1]    -> Fill(track_phi, Get_eta_pair.second);
                    
                    track_eta_z_2D[centrality_map[centrality_bin]] -> Fill(evt_z.first, Get_eta_pair.second);
                    track_eta_z_2D[track_eta_z_2D.size() - 1]      -> Fill(evt_z.first, Get_eta_pair.second);
            
                    evt_NTrack += 1;
                    // cout<<event_i<<" test_10 "<<inner_i<<" "<<outer_i<<endl;
                }

            }
            
        }
    }

    // cout<<" "<<endl;
    // cout<<" "<<endl;
    // cout<<"test_8"<<endl;
    if (run_type == "MC")
    {
        track_correlation_2D -> Fill(evt_NTrack_MC, evt_NTrack);
        track_ratio_1D[centrality_map[centrality_bin]] -> Fill( double(evt_NTrack) / double(evt_NTrack_MC) );
        track_ratio_1D[track_ratio_1D.size() - 1]      -> Fill( double(evt_NTrack) / double(evt_NTrack_MC) );
    }
        
    track_cluster_ratio_multiplicity_2D -> Fill( total_NClus, double(total_NClus) / double(evt_NTrack) );
    track_cluster_ratio_1D[centrality_map[centrality_bin]]    -> Fill( double(total_NClus) / double(evt_NTrack) );
    track_cluster_ratio_1D[track_cluster_ratio_1D.size() - 1] -> Fill( double(total_NClus) / double(evt_NTrack) );

    if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0){print_evt_plot(event_i, evt_NTrack_MC, inner_NClu, outer_NClu);}

    return_tag = 1;
}

// // note : this function is for the event by event vertex calculation
// void INTTEta::ProcessEvt(int event_i, vector<clu_info> temp_sPH_inner_nocolumn_vec, vector<clu_info> temp_sPH_outer_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_rz_vec, int NvtxMC, vector<double> TrigvtxMC, bool PhiCheckTag, Long64_t bco_full, pair<double,double> evt_z, int centrality_bin, vector<vector<float>> true_track_info) // note : evt_z : {z, width} 
// {   
//     return_tag = 0;

//     if (event_i%1 == 0) {cout<<"In INTTEta class, running event : "<<event_i<<endl;}

//     inner_NClu = temp_sPH_inner_nocolumn_vec.size();
//     outer_NClu = temp_sPH_outer_nocolumn_vec.size();
//     total_NClus = inner_NClu + outer_NClu;

//     // cout<<"test_0"<<endl;
//     if (total_NClus < zvtx_cal_require) {return; cout<<"return confirmation"<<endl;}   
    
//     if (run_type == "MC" && NvtxMC != 1) { return; cout<<"In INTTEta class, event : "<<event_i<<" Nvtx : "<<NvtxMC<<" Nvtx more than one "<<endl;}
//     if (PhiCheckTag == false)            { return; cout<<"In INTTEta class, event : "<<event_i<<" Nvtx : "<<NvtxMC<<" Not full phi has hits "<<endl;}
    
//     if (inner_NClu < 10 || outer_NClu < 10 || total_NClus > N_clu_cut || total_NClus < N_clu_cutl)
//     {
//         return;
//         printf("In INTTEta class, event : %i, low clu continue, NClus : %lu \n", event_i, total_NClus); 
//     }

//     // todo : the z vertex range is here
//     if (-220 > evt_z.first + evt_z.second || -180 < evt_z.first - evt_z.second) {return;}
    
    
//     N_GoodEvent += 1;
//     N_GoodEvent_vec[centrality_map[centrality_bin]] += 1;

//     // cout<<"N inner cluster : "<<inner_NClu<<" N outer cluster : "<<outer_NClu<<endl;

//     if (run_type == "MC")
//     {
//         // note : for the true track case ----------------------------------------------------------------------------------------------------------------------------------------------------------------
//         double INTT_eta_acceptance_l = -0.5 * TMath::Log((sqrt(pow(-230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))-(-230.-TrigvtxMC[2]*10.)) / (sqrt(pow(-230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))+(-230.-TrigvtxMC[2]*10.))); // note : left
//         double INTT_eta_acceptance_r =  -0.5 * TMath::Log((sqrt(pow(230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))-(230.-TrigvtxMC[2]*10.)) / (sqrt(pow(230.-TrigvtxMC[2]*10.,2)+pow(INTT_layer_R[3],2))+(230.-TrigvtxMC[2]*10.))); // note : right
//         // if (event_i % 100 == 0){cout<<"z : "<<TrigvtxMC[2]*10.<<" eta : "<<INTT_eta_acceptance_l<<" "<<INTT_eta_acceptance_r<<endl;}

//         // cout<<"true_track_info : "<<true_track_info.size()<<endl;    
//         for (int track_i = 0; track_i < true_track_info.size(); track_i++)
//         {
//             if (true_track_info[track_i][2] == 111 || true_track_info[track_i][2] == 22 || abs(true_track_info[track_i][2]) == 2112){continue;}

//             if (true_track_info[track_i][0] > INTT_eta_acceptance_l && true_track_info[track_i][0] < INTT_eta_acceptance_r)
//             {
//                 dNdeta_1D_MC[centrality_map[centrality_bin]] -> Fill(true_track_info[track_i][0]);
//                 dNdeta_1D_MC[dNdeta_1D_MC.size() - 1]        -> Fill(true_track_info[track_i][0]);   
//                 final_dNdeta_1D_MC[centrality_map[centrality_bin]] -> Fill(true_track_info[track_i][0]);
//                 final_dNdeta_1D_MC[dNdeta_1D_MC.size() - 1]        -> Fill(true_track_info[track_i][0]);
//                 track_eta_z_2D_MC[centrality_map[centrality_bin]] -> Fill(TrigvtxMC[2]*10., true_track_info[track_i][0]);
//                 track_eta_z_2D_MC[track_eta_z_2D_MC.size() - 1]   -> Fill(TrigvtxMC[2]*10., true_track_info[track_i][0]);    

//                 // cout<<"true track eta : "<<true_track_info[track_i][0]<<" phi : "<<convertTo360(true_track_info[track_i][1])<<endl;
//                 // cout<<"("<<true_track_info[track_i][0]<<", "<<convertTo360(true_track_info[track_i][1])<<"), ";

//                 evt_true_track_gr -> SetPoint(evt_true_track_gr->GetN(),true_track_info[track_i][0], convertTo360(true_track_info[track_i][1]));

//                 evt_NTrack_MC += 1;
//             }
//         }
//         // if (evt_NTrack_MC < 10) {cout<<"evt : "<<event_i<<" ---- N reco track : "<<evt_NTrack<<" N true track : "<<evt_NTrack_MC<<" ratio : "<<double(evt_NTrack) / double(evt_NTrack_MC)<<endl;}
//     }

//     // note : put the cluster into the phi map, the first bool is for the cluster usage.
//     // note : false means the cluster is not used
//     for (int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++) {
//         Clus_InnerPhi_Offset = (temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second < 0) ? atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi());
//         inner_clu_phi_map[ int(Clus_InnerPhi_Offset) ].push_back({false,temp_sPH_inner_nocolumn_vec[inner_i]});
//     }
//     for (int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++) {
//         Clus_OuterPhi_Offset = (temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second < 0) ? atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi());
//         outer_clu_phi_map[ int(Clus_OuterPhi_Offset) ].push_back({false,temp_sPH_outer_nocolumn_vec[outer_i]});
//     }

//     // note : for two-cluster tracklets only
//     for (int inner_phi_i = 0; inner_phi_i < 360; inner_phi_i++) // note : each phi cell (1 degree)
//     {
//         // note : number of cluster in this phi cell
//         for (int inner_phi_clu_i = 0; inner_phi_clu_i < inner_clu_phi_map[inner_phi_i].size(); inner_phi_clu_i++)
//         {
//             if (inner_clu_phi_map[inner_phi_i][inner_phi_clu_i].first == true) {continue;}

//             if (inner_phi_i == 0 || inner_phi_i == 359)
//             {

//             }
//             else 
//             {
//                 // note : the outer phi index, -1, 0, 1
//                 for (int scan_i = -1; scan_i < 2; scan_i++)
//                 {

//                 }
//             }

             
//         }
//     }

//     // cout<<"test_1"<<endl;
//     for ( int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++ )
//     {    
//         loop_NGoodPair = 0;
        
//         for ( int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++ )
//         {
//             // cout<<event_i<<" test_5 "<<inner_i<<" "<<outer_i<<endl;
//             // note : calculate the cluster phi after the vertex correction which can enhence the purity of the tracklet selection
//             Clus_InnerPhi_Offset = (temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second < 0) ? atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi());
//             Clus_OuterPhi_Offset = (temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second < 0) ? atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi());
//             double delta_phi = Clus_InnerPhi_Offset - Clus_OuterPhi_Offset;
//             double track_phi = (Clus_InnerPhi_Offset + Clus_OuterPhi_Offset)/2.;
//             double inner_clu_eta = get_clu_eta({beam_origin.first, beam_origin.second, evt_z.first},{temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y, temp_sPH_inner_nocolumn_vec[inner_i].z});
//             double outer_clu_eta = get_clu_eta({beam_origin.first, beam_origin.second, evt_z.first},{temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y, temp_sPH_outer_nocolumn_vec[outer_i].z});
//             double delta_eta = inner_clu_eta - outer_clu_eta;
            
//             // cout<<"test_2"<<endl;
//             // note : before calculating the possible z vertex range of one tracklet, the vertex correction was applied
//             pair<double,double> z_range_info = Get_possible_zvtx( 
//                 0., // get_radius(beam_origin.first,beam_origin.second), 
//                 {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z}, // note : unsign radius
//                 {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}  // note : unsign radius
//             );

//             if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0)
//             {
//                 Get_eta_pair = Get_eta(
//                     {get_radius(beam_origin.first,beam_origin.second), evt_z.first,evt_z.second},
//                     {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z},
//                     {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}
//                 );

//                 evt_reco_track_gr_All -> SetPoint(evt_reco_track_gr_All->GetN(),Get_eta_pair.second, track_phi);
//             }

            
//             // if ( Get_eta_pair.second > 1.4 && Get_eta_pair.second < 1.7 && track_phi > 50 && track_phi < 60 ) {
//             //     cout<<" "<<endl;
//             //     cout<<"reco eta : "<<Get_eta_pair.second<<" reco phi : "<<track_phi<<" delta eta : "<<delta_eta<<endl;
//             //     cout<<"all pair info, inner clu pos "<<temp_sPH_inner_nocolumn_vec[inner_i].x<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].y<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].z<<" phi : "<<Clus_InnerPhi_Offset<<endl;
//             //     cout<<"all pair info, outer clu pos "<<temp_sPH_outer_nocolumn_vec[outer_i].x<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].y<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].z<<" phi : "<<Clus_OuterPhi_Offset<<endl;
//             //     cout<<"z vertex range : "<<z_range_info.first - z_range_info.second<<" "<<z_range_info.first + z_range_info.second<<endl;
//             //     cout<<"accept z range : "<<evt_z.first - evt_z.second<<" "<<evt_z.first + evt_z.second<<endl;
                
//             // }

//             if (fabs(delta_phi) > 5.72) {continue;}
//             if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_PhiLoose -> SetPoint(evt_reco_track_gr_PhiLoose->GetN(),Get_eta_pair.second, track_phi);}

//             track_delta_eta_1D[centrality_map[centrality_bin]] -> Fill( delta_eta );
//             track_delta_eta_1D[track_delta_eta_1D.size() - 1]  -> Fill( delta_eta );

//             track_DeltaPhi_DeltaEta_2D[centrality_map[centrality_bin]]        -> Fill(delta_phi, delta_eta);
//             track_DeltaPhi_DeltaEta_2D[track_DeltaPhi_DeltaEta_2D.size() - 1] -> Fill(delta_phi, delta_eta);

//             // cout<<event_i<<" test_6 "<<inner_i<<" "<<outer_i<<endl;
//             // note : this is a cut to constraint on the z vertex, only if the tracklets with the range that covers the z vertex can pass this cut
//             if (z_range_info.first - z_range_info.second > evt_z.first + evt_z.second || z_range_info.first + z_range_info.second < evt_z.first - evt_z.second) {continue;}
//             // cout<<"range test: "<<z_range_info.first - z_range_info.second<<" "<<z_range_info.first + z_range_info.second<<" vertex: "<<evt_z.first<<" "<<evt_z.second<<endl;
//             // cout<<"test_3"<<endl;

//             if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_Z -> SetPoint(evt_reco_track_gr_Z->GetN(),Get_eta_pair.second, track_phi);}
            
//             // note : the reason to use the DCA calculation with sign is because that the distribution of 1D DCA will be single peak only
//             double DCA_sign = calculateAngleBetweenVectors(
//                 temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y,
//                 temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y,
//                 beam_origin.first, beam_origin.second
//             );

//             // vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
//             //     temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y,
//             //     temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y,
//             //     beam_origin.first, beam_origin.second
//             // );
//             // if (DCA_info_vec[0] != fabs(DCA_sign) && fabs( DCA_info_vec[0] - fabs(DCA_sign) ) > 0.1 ){
//             //     cout<<"different DCA : "<<DCA_info_vec[0]<<" "<<DCA_sign<<" diff : "<<DCA_info_vec[0] - fabs(DCA_sign)<<endl;
//             // }

//             track_delta_eta_1D_post[centrality_map[centrality_bin]]     -> Fill( delta_eta );
//             track_delta_eta_1D_post[track_delta_eta_1D_post.size() - 1] -> Fill( delta_eta );

//             track_delta_phi_1D[centrality_map[centrality_bin]] -> Fill( delta_phi );
//             track_delta_phi_1D[track_delta_phi_1D.size() - 1]  -> Fill( delta_phi );
            
//             track_DCA_distance[centrality_map[centrality_bin]] -> Fill( DCA_sign );
//             track_DCA_distance[track_DCA_distance.size() - 1]  -> Fill( DCA_sign ); 

//             track_phi_DCA_2D[centrality_map[centrality_bin]] -> Fill( delta_phi, DCA_sign );
//             track_phi_DCA_2D[track_phi_DCA_2D.size() - 1]    -> Fill( delta_phi, DCA_sign );

//             // cout<<" "<<endl;
//             // cout<<"inner_i : "<<inner_i<<" outer_i "<<outer_i<<" true z : "<<TrigvtxMC[2]*10.<<" reco z: "<<evt_z.first<<endl;
//             // cout<<"all pair info, inner clu pos "<<temp_sPH_inner_nocolumn_vec[inner_i].x<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].y<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].z<<" phi : "<<Clus_InnerPhi_Offset<<endl;
//             // cout<<"all pair info, outer clu pos "<<temp_sPH_outer_nocolumn_vec[outer_i].x<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].y<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].z<<" phi : "<<Clus_OuterPhi_Offset<<endl;
//             // if (fabs(delta_phi) < 4) {cout<<"("<<Get_eta_pair.second<<", "<<track_phi<<"), ";}

//             // if (DCA_cut.first < DCA_sign && DCA_sign < DCA_cut.second)
//             if (true)
//             {
//                 if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_ZDCA -> SetPoint(evt_reco_track_gr_ZDCA->GetN(),Get_eta_pair.second, track_phi);}

//                 if (fabs(delta_phi) < phi_diff_cut)
//                 {
//                     if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0) {evt_reco_track_gr_ZDCAPhi -> SetPoint(evt_reco_track_gr_ZDCAPhi->GetN(),Get_eta_pair.second, track_phi);}    
//                     // cout<<"reco eta : "<<Get_eta_pair.second<<" reco phi : "<<track_phi<<endl;

//                     // if (fabs(Get_eta_pair.second) < 0.1) { cout<<"eta "<<Get_eta_pair.second<<", three points : "<<get_radius(beam_origin.first,beam_origin.second)<<" "<<evt_z.first<<" "<<evt_z.second<<" "<<
//                     // get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second)<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].z<<" "<<
//                     // get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second)<<" "<<temp_sPH_outer_nocolumn_vec[outer_i].z<<endl;}

//                     Get_eta_pair = Get_eta(
//                         {get_radius(beam_origin.first,beam_origin.second), evt_z.first,evt_z.second},
//                         {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z},
//                         {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}
//                     );

//                     double eta_bin = eta_region_hist -> Fill(Get_eta_pair.second);
//                     // cout<<"test1 "<<eta_bin<<endl;
//                     if (eta_bin != -1)
//                     {
//                         final_track_delta_phi_1D[centrality_map[centrality_bin]][eta_bin - 1]      -> Fill(delta_phi);
//                         final_track_delta_phi_1D[final_track_delta_phi_1D.size() - 1][eta_bin - 1] -> Fill(delta_phi);
//                     }
                    

//                     reco_eta_correlation_2D -> Fill(Get_eta_pair.second, (inner_clu_eta + outer_clu_eta)/2.);
//                     reco_eta_diff_reco3P_2D -> Fill(Get_eta_pair.second, Get_eta_pair.second - (inner_clu_eta + outer_clu_eta)/2.);
//                     reco_eta_diff_1D        -> Fill(Get_eta_pair.second - (inner_clu_eta + outer_clu_eta)/2.);

//                     track_DeltaPhi_eta_2D[centrality_map[centrality_bin]]   -> Fill(delta_phi, Get_eta_pair.second);
//                     track_DeltaPhi_eta_2D[track_DeltaPhi_eta_2D.size() - 1] -> Fill(delta_phi, Get_eta_pair.second);

//                     // cout<<"test_5"<<endl;
//                     dNdeta_1D[centrality_map[centrality_bin]] -> Fill(Get_eta_pair.second);
//                     dNdeta_1D[dNdeta_1D.size() - 1]           -> Fill(Get_eta_pair.second);

//                     track_eta_phi_2D[centrality_map[centrality_bin]] -> Fill(track_phi, Get_eta_pair.second);
//                     track_eta_phi_2D[track_eta_phi_2D.size() - 1]    -> Fill(track_phi, Get_eta_pair.second);
                    
//                     track_eta_z_2D[centrality_map[centrality_bin]] -> Fill(evt_z.first, Get_eta_pair.second);
//                     track_eta_z_2D[track_eta_z_2D.size() - 1]      -> Fill(evt_z.first, Get_eta_pair.second);
            
//                     evt_NTrack += 1;
//                     // cout<<event_i<<" test_10 "<<inner_i<<" "<<outer_i<<endl;
//                 }

//             }
            
//         }
//     }

//     // cout<<" "<<endl;
//     // cout<<" "<<endl;
//     // cout<<"test_8"<<endl;
//     if (run_type == "MC")
//     {
//         track_correlation_2D -> Fill(evt_NTrack_MC, evt_NTrack);
//         track_ratio_1D[centrality_map[centrality_bin]] -> Fill( double(evt_NTrack) / double(evt_NTrack_MC) );
//         track_ratio_1D[track_ratio_1D.size() - 1]      -> Fill( double(evt_NTrack) / double(evt_NTrack_MC) );
//     }
        
//     track_cluster_ratio_multiplicity_2D -> Fill( total_NClus, double(total_NClus) / double(evt_NTrack) );
//     track_cluster_ratio_1D[centrality_map[centrality_bin]]    -> Fill( double(total_NClus) / double(evt_NTrack) );
//     track_cluster_ratio_1D[track_cluster_ratio_1D.size() - 1] -> Fill( double(total_NClus) / double(evt_NTrack) );

//     if (run_type == "MC" && evt_NTrack_MC < draw_evt_cut && event_i % print_plot_ratio == 0){print_evt_plot(event_i, evt_NTrack_MC, inner_NClu, outer_NClu);}

//     return_tag = 1;
// }

void INTTEta::ClearEvt()
{
    if (evt_reco_track_gr_All -> GetN() != 0) {evt_reco_track_gr_All -> Set(0);}
    if (evt_reco_track_gr_PhiLoose -> GetN() != 0) {evt_reco_track_gr_PhiLoose -> Set(0);}
    if (evt_reco_track_gr_Z -> GetN() != 0) {evt_reco_track_gr_Z -> Set(0);}
    if (evt_reco_track_gr_ZDCA -> GetN() != 0) {evt_reco_track_gr_ZDCA -> Set(0);}
    if (evt_reco_track_gr_ZDCAPhi -> GetN() != 0) {evt_reco_track_gr_ZDCAPhi -> Set(0);}
    if (evt_true_track_gr -> GetN() != 0) {evt_true_track_gr -> Set(0);}
    if (track_gr -> GetN() != 0) {track_gr -> Set(0);}

    // for (int i = 0; i < 360; i++)
    // {
    //     inner_clu_phi_map[i].clear();
    //     outer_clu_phi_map[i].clear();
    // }

    std::fill(std::begin(inner_clu_phi_map), std::end(inner_clu_phi_map), std::vector<pair<bool,clu_info>>());
    std::fill(std::begin(outer_clu_phi_map), std::end(outer_clu_phi_map), std::vector<pair<bool,clu_info>>());

    return_tag = 0;
    evt_NTrack = 0;
    evt_NTrack_MC = 0;
    return;
}

double INTTEta::get_dist_offset(TH1F * hist_in, int check_N_bin) // note : check_N_bin 1 to N bins of hist
{
    if (check_N_bin < 0 || check_N_bin > hist_in -> GetNbinsX()) {cout<<" wrong check_N_bin "<<endl; exit(1);}
    double total_entry = 0;
    for (int i = 0; i < check_N_bin; i++)
    {
        total_entry += hist_in -> GetBinContent(i+1); // note : 1, 2, 3.....
        total_entry += hist_in -> GetBinContent(hist_in -> GetNbinsX() - i);
    }

    return total_entry/double(2. * check_N_bin);
}

void INTTEta::PrintPlots()
{
    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_cluster_ratio_1D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_cluster_ratio_1D.size(); i++)
    {
        track_cluster_ratio_1D[i] -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/track_cluster_ratio_1D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_cluster_ratio_1D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/dNdeta_1D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < dNdeta_1D.size(); i++)
    {   
        double N_correction_evt = (i == dNdeta_1D_MC.size() - 1) ? N_GoodEvent : N_GoodEvent_vec[i];
        dNdeta_1D_MC[i] -> Scale(1./double(dNdeta_1D_MC[i] -> GetBinWidth(1) ));
        dNdeta_1D_MC[i] -> Scale(1./double(N_correction_evt));
        
        dNdeta_1D[i] -> Scale(1./double(dNdeta_1D[i] -> GetBinWidth(1) ));
        dNdeta_1D[i] -> Scale(1./double(N_correction_evt));
        dNdeta_1D[i] -> GetYaxis() -> SetRangeUser(0,800);
        dNdeta_1D[i] -> Draw("ep");
        dNdeta_1D_MC[i] -> Draw("hist same");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/dNdeta_1D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/dNdeta_1D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_eta_phi_2D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_eta_phi_2D.size(); i++)
    {
        track_eta_phi_2D[i] -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/track_eta_phi_2D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_eta_phi_2D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_eta_z_2D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_eta_z_2D.size(); i++)
    {
        track_eta_z_2D[i] -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/track_eta_z_2D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_eta_z_2D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_delta_phi_1D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_delta_phi_1D.size(); i++)
    {
        track_delta_phi_1D[i] -> SetMinimum(0);
        track_delta_phi_1D[i] -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        coord_line -> DrawLine(-1*phi_diff_cut, 0, -1 * phi_diff_cut, track_delta_phi_1D[i]->GetMaximum() * 1.1);
        coord_line -> DrawLine(phi_diff_cut, 0, phi_diff_cut, track_delta_phi_1D[i]->GetMaximum() * 1.1);
        c1 -> Print(Form("%s/track_delta_phi_1D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_delta_phi_1D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_delta_eta_1D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_delta_eta_1D.size(); i++)
    {
        track_delta_eta_1D[i] -> SetMinimum(0);
        track_delta_eta_1D[i] -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        // coord_line -> DrawLine(-1*phi_diff_cut, 0, -1 * phi_diff_cut, track_delta_eta_1D[i]->GetMaximum() * 1.1);
        // coord_line -> DrawLine(phi_diff_cut, 0, phi_diff_cut, track_delta_eta_1D[i]->GetMaximum() * 1.1);
        c1 -> Print(Form("%s/track_delta_eta_1D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_delta_eta_1D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_delta_eta_1D_post.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_delta_eta_1D_post.size(); i++)
    {
        track_delta_eta_1D_post[i] -> SetMinimum(0);
        track_delta_eta_1D_post[i] -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        // coord_line -> DrawLine(-1*phi_diff_cut, 0, -1 * phi_diff_cut, track_delta_eta_1D_post[i]->GetMaximum() * 1.1);
        // coord_line -> DrawLine(phi_diff_cut, 0, phi_diff_cut, track_delta_eta_1D_post[i]->GetMaximum() * 1.1);
        c1 -> Print(Form("%s/track_delta_eta_1D_post.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_delta_eta_1D_post.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_DCA_distance.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_DCA_distance.size(); i++)
    {
        track_DCA_distance[i] -> SetMinimum(0);
        track_DCA_distance[i] -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        coord_line -> DrawLine(DCA_cut.first, 0, DCA_cut.first, track_DCA_distance[i]->GetMaximum() * 1.05);
        coord_line -> DrawLine(DCA_cut.second, 0, DCA_cut.second, track_DCA_distance[i]->GetMaximum() * 1.05);
        c1 -> Print(Form("%s/track_DCA_distance.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_DCA_distance.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_phi_DCA_2D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_phi_DCA_2D.size(); i++)
    {
        track_phi_DCA_2D[i] -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        
        coord_line -> DrawLine(track_phi_DCA_2D[i] -> GetXaxis() -> GetXmin(), DCA_cut.first, track_phi_DCA_2D[i] -> GetXaxis() -> GetXmax(), DCA_cut.first);
        coord_line -> DrawLine(track_phi_DCA_2D[i] -> GetXaxis() -> GetXmin(), DCA_cut.second, track_phi_DCA_2D[i] -> GetXaxis() -> GetXmax(), DCA_cut.second);
        coord_line -> DrawLine(-1 * phi_diff_cut, track_phi_DCA_2D[i] -> GetYaxis() -> GetXmin(), -1 * phi_diff_cut, track_phi_DCA_2D[i] -> GetYaxis() -> GetXmax());
        coord_line -> DrawLine(phi_diff_cut, track_phi_DCA_2D[i] -> GetYaxis() -> GetXmin(), phi_diff_cut, track_phi_DCA_2D[i] -> GetYaxis() -> GetXmax());
        c1 -> Print(Form("%s/track_phi_DCA_2D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_phi_DCA_2D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_ratio_1D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_ratio_1D.size(); i++)
    {
        track_ratio_1D[i] -> SetMinimum(0);
        track_ratio_1D[i] -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/track_ratio_1D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_ratio_1D.pdf)", out_folder_directory.c_str()) );
    
    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/final_track_delta_phi_1D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < final_track_delta_phi_1D.size(); i++)
    {
        for (int i1 = 0; i1 < final_track_delta_phi_1D[i].size(); i1++)
        {   
            double hist_offset = get_dist_offset(final_track_delta_phi_1D[i][i1], 15);
            gaus_pol1_fit->SetParameters( final_track_delta_phi_1D[i][i1] -> GetBinContent(final_track_delta_phi_1D[i][i1] -> GetMaximumBin()) - hist_offset, 0, final_track_delta_phi_1D[i][i1]->GetStdDev()/2., hist_offset, 0);
            final_track_delta_phi_1D[i][i1] -> Fit(gaus_pol1_fit,"NQ");

            // note : par[0] : size
            // note : par[1] : ratio of the two gaussians
            // note : par[2] : mean
            // note : par[3] : width of gaus 1
            // note : par[4] : width of gaus 2
            // note : par[5] : offset
            // note : par[6] : slope
            d_gaus_pol1_fit -> SetParameters(final_track_delta_phi_1D[i][i1] -> GetBinContent(final_track_delta_phi_1D[i][i1] -> GetMaximumBin()) - hist_offset, 0.2, 0, final_track_delta_phi_1D[i][i1]->GetStdDev()/2., final_track_delta_phi_1D[i][i1]->GetStdDev()/2., hist_offset, 0);
            d_gaus_pol1_fit -> SetParLimits(1, 0, 0.5);  // note : the first gaussian is  the main distribution, So it should contain more than 50% of the total distribution
            d_gaus_pol1_fit -> SetParLimits(3, 0, 1000); // note : the width of the gaussian should be positive
            d_gaus_pol1_fit -> SetParLimits(4, 0, 1000); // note : the width of the gaussian should be positive
            final_track_delta_phi_1D[i][i1] -> Fit(d_gaus_pol1_fit,"NQ");

            draw_d_gaus -> SetParameters(d_gaus_pol1_fit->GetParameter(0), d_gaus_pol1_fit->GetParameter(1), d_gaus_pol1_fit->GetParameter(2), d_gaus_pol1_fit->GetParameter(3), d_gaus_pol1_fit->GetParameter(4));
            draw_pol1_line -> SetParameters(d_gaus_pol1_fit -> GetParameter(5), d_gaus_pol1_fit -> GetParameter(6));

            final_track_delta_phi_1D[i][i1] -> SetMinimum(0);
            final_track_delta_phi_1D[i][i1] -> Draw("hist"); 
            d_gaus_pol1_fit -> Draw("lsame");
            draw_d_gaus -> Draw("lsame");
            draw_pol1_line -> Draw("lsame");

            // gaus_pol1_fit -> Draw("lsame");
            // draw_gaus_line -> SetParameters(fabs(gaus_pol1_fit -> GetParameter(0)), gaus_pol1_fit -> GetParameter(1), fabs(gaus_pol1_fit -> GetParameter(2)), 0);
            // draw_gaus_line -> Draw("lsame");
            // draw_pol1_line -> SetParameters(gaus_pol1_fit -> GetParameter(3), gaus_pol1_fit -> GetParameter(4));
            // draw_pol1_line -> Draw("lsame");

            // final_eta_entry[i].push_back((draw_gaus_line -> GetParameter(1) - 3 * draw_gaus_line -> GetParameter(2)), draw_gaus_line -> GetParameter(1) + 3 * draw_gaus_line -> GetParameter(2));
            // cout<<i<<" "<<i1<<" gaus fit par  : "<<fabs(gaus_pol1_fit -> GetParameter(0))<<" "<<(gaus_pol1_fit -> GetParameter(1))<<" "<<fabs(gaus_pol1_fit -> GetParameter(2))<<endl;
            double gaus_integral = fabs(draw_gaus_line -> Integral( (draw_gaus_line -> GetParameter(1) - 3 * fabs(draw_gaus_line -> GetParameter(2))), draw_gaus_line -> GetParameter(1) + 3 * fabs(draw_gaus_line -> GetParameter(2)) )) / final_track_delta_phi_1D[i][i1] -> GetBinWidth(1);
            cout<<i<<" "<<i1<<" gaus integral : "<< gaus_integral <<endl;

            double d_gaus_integral = fabs(draw_d_gaus -> Integral( -0.6, 0.6 )) / final_track_delta_phi_1D[i][i1] -> GetBinWidth(1);
            cout<<i<<" "<<i1<<" D-gaus integral : "<< d_gaus_integral <<endl;
            
            final_dNdeta_1D[i]->SetBinContent(i1 + 1, d_gaus_integral );
            
            ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
            draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s, #eta: %.2f - %.2f",centrality_region[i].c_str(), eta_region_hist -> GetBinCenter(i1+1) - eta_region_hist -> GetBinWidth(i1+1)/2., eta_region_hist -> GetBinCenter(i1+1) + eta_region_hist -> GetBinWidth(i1+1)/2.));
            draw_text -> DrawLatex(0.21, 0.85, Form("Guassian integral : %.2f", gaus_integral));
            draw_text -> DrawLatex(0.21, 0.80, Form("D-Guassian integral : %.2f", d_gaus_integral));
            c1 -> Print(Form("%s/final_track_delta_phi_1D.pdf", out_folder_directory.c_str()));
            c1 -> Clear();
        }
    }
    c1 -> Print( Form("%s/final_track_delta_phi_1D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/final_dNdeta_1D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < final_dNdeta_1D.size(); i++)
    {   
        double N_correction_evt = (i == final_dNdeta_1D_MC.size() - 1) ? N_GoodEvent : N_GoodEvent_vec[i];
        final_dNdeta_1D_MC[i] -> Scale(1./double(final_dNdeta_1D_MC[i] -> GetBinWidth(1) ));
        final_dNdeta_1D_MC[i] -> Scale(1./double(N_correction_evt));
        
        final_dNdeta_1D[i] -> Scale(1./double(final_dNdeta_1D[i] -> GetBinWidth(1) ));
        final_dNdeta_1D[i] -> Scale(1./double(N_correction_evt));
        final_dNdeta_1D[i] -> GetYaxis() -> SetRangeUser(0,800);
        final_dNdeta_1D[i] -> Draw("ep");
        final_dNdeta_1D_MC[i] -> Draw("hist same");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/final_dNdeta_1D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/final_dNdeta_1D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    track_cluster_ratio_multiplicity_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print( Form("%s/track_cluster_ratio_multiplicity_2D.pdf", out_folder_directory.c_str()) );
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    track_correlation_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print( Form("%s/track_correlation_2D.pdf", out_folder_directory.c_str()) );
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    reco_eta_correlation_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print( Form("%s/reco_eta_correlation_2D.pdf", out_folder_directory.c_str()) );
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    reco_eta_diff_reco3P_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print( Form("%s/reco_eta_diff_reco3P_2D.pdf", out_folder_directory.c_str()) );
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    reco_eta_diff_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print( Form("%s/reco_eta_diff_1D.pdf", out_folder_directory.c_str()) );
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_DeltaPhi_eta_2D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_DeltaPhi_eta_2D.size(); i++)
    {
        track_DeltaPhi_eta_2D[i] -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/track_DeltaPhi_eta_2D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_DeltaPhi_eta_2D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_DeltaPhi_DeltaEta_2D.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_DeltaPhi_DeltaEta_2D.size(); i++)
    {
        track_DeltaPhi_DeltaEta_2D[i] -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/track_DeltaPhi_DeltaEta_2D.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_DeltaPhi_DeltaEta_2D.pdf)", out_folder_directory.c_str()) );

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/dNdeta_1D_MC.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < dNdeta_1D_MC.size(); i++)
    {
        // dNdeta_1D_MC[i] -> Scale(1./double(dNdeta_1D_MC[i] -> GetBinWidth(1) ));
        // double N_correction_evt = (i == dNdeta_1D_MC.size() - 1) ? N_GoodEvent : N_GoodEvent_vec[i];
        // dNdeta_1D_MC[i] -> Scale(1./double(N_GoodEvent));
        dNdeta_1D_MC[i] -> GetYaxis() -> SetRangeUser(0,800);
        dNdeta_1D_MC[i] -> Draw("hist");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/dNdeta_1D_MC.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/dNdeta_1D_MC.pdf)", out_folder_directory.c_str()) );


    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    c1 -> Print( Form("%s/track_eta_z_2D_MC.pdf(", out_folder_directory.c_str()) );
    for (int i = 0; i < track_eta_z_2D_MC.size(); i++)
    {
        track_eta_z_2D_MC[i] -> Draw("colz0");
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.90, Form("Centrality : %s",centrality_region[i].c_str()));
        c1 -> Print(Form("%s/track_eta_z_2D_MC.pdf", out_folder_directory.c_str()));
        c1 -> Clear();    
    }
    c1 -> Print( Form("%s/track_eta_z_2D_MC.pdf)", out_folder_directory.c_str()) );
}

void INTTEta::EndRun()
{
    for (int i = 0; i < N_GoodEvent_vec.size(); i++)
    {
        cout<<"N good evt break down : "<< i <<" "<<N_GoodEvent_vec[i]<<endl;
    }

    cout<<"N good evt inclusive : "<<N_GoodEvent<<endl;

    return;
}

double INTTEta::grEY_stddev(TGraphErrors * input_grr)
{
    vector<double> input_vector; input_vector.clear();
    for (int i = 0; i < input_grr -> GetN(); i++) {input_vector.push_back( input_grr -> GetPointY(i) );}

    double average = accumulate( input_vector.begin(), input_vector.end(), 0.0 ) / double(input_vector.size());

    double sum_subt = 0;
	for (int ele : input_vector) {sum_subt+=pow((ele-average),2);}
	
	return sqrt(sum_subt/(input_vector.size()-1));
}	

pair<double, double> INTTEta::mirrorPolynomial(double a, double b) {
    // Interchange 'x' and 'y'
    double mirroredA = 1.0 / a;
    double mirroredB = -b / a;

    return {mirroredA, mirroredB};
}

// note : pair<reduced chi2, eta of the track>
// note : vector : {r,z}
// note : p0 : vertex point {r,z,zerror}
// note : p1 : inner layer
// note : p2 : outer layer
pair<double,double> INTTEta::Get_eta(vector<double>p0, vector<double>p1, vector<double>p2)
{
    // if (track_gr -> GetN() != 0) {cout<<"In INTTEta class, track_gr is not empty, track_gr size : "<<track_gr -> GetN()<<endl; exit(0);}
    
    if (track_gr -> GetN() != 0) {track_gr -> Set(0);}
    
    vector<double> r_vec  = {p0[0], p1[0], p2[0]}; 
    vector<double> z_vec  = {p0[1], p1[1], p2[1]}; 
    vector<double> re_vec = {0, 0, 0}; 
    vector<double> ze_vec = {p0[2], ( fabs( p1[1] ) < 130 ) ? 8. : 10., ( fabs( p2[1] ) < 130 ) ? 8. : 10.}; 

    // note : swap the r and z, to have easier fitting 
    // note : in principle, Z should be in the x axis, R should be in the Y axis
    for (int i = 0; i < 3; i++)
    {
        track_gr -> SetPoint(i,r_vec[i],z_vec[i]);
        track_gr -> SetPointError(i,re_vec[i],ze_vec[i]);

        // cout<<"In INTTEta class, point : "<<i<<" r : "<<r_vec[i]<<" z : "<<z_vec[i]<<" re : "<<re_vec[i]<<" ze : "<<ze_vec[i]<<endl;
    }    
    
    double vertical_line = ( fabs( grEY_stddev(track_gr) ) < 0.00001 ) ? 1 : 0;
    
    if (vertical_line) {return {0,0};}
    else 
    {
        fit_rz -> SetParameters(0,0);

        track_gr -> Fit(fit_rz,"NQ");

        pair<double,double> ax_b = mirrorPolynomial(fit_rz -> GetParameter(1),fit_rz -> GetParameter(0));

        return  {(fit_rz -> GetChisquare() / double(fit_rz -> GetNDF())), -1 * TMath::Log( fabs( tan( atan2(ax_b.first, (ax_b.first > 0) ? 1. : -1. ) / 2 ) ) )};
    }

}

double INTTEta::convertTo360(double radian) {

    double angle = radian * (180./M_PI);
    
    if (fabs(radian) == M_PI) {return 180.;}
    else if (radian == 0)
    {
        return 0;
    }
    else if ( radian > 0 )
    {
        return angle;
    }
    else
    {
        return angle + 360;
    }
}

double INTTEta::get_clu_eta(vector<double> vertex, vector<double> clu_pos)
{
    double correct_x = clu_pos[0] - vertex[0];
    double correct_y = clu_pos[1] - vertex[1];
    double correct_z = clu_pos[2] - vertex[2];
    double clu_r = sqrt(pow(correct_x,2) + pow(correct_y,2));
    // cout<<"correct info : "<<correct_x<<" "<<correct_y<<" "<<correct_z<<" "<<clu_r<<endl;

    return -0.5 * TMath::Log((sqrt(pow(correct_z,2)+pow(clu_r,2))-(correct_z)) / (sqrt(pow(correct_z,2)+pow(clu_r,2))+(correct_z)));
}

#endif