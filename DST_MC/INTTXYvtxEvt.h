#ifndef INTTXYvtxEvt_h
#define INTTXYvtxEvt_h

#include "INTTXYvtx.h"

class INTTXYvtxEvt : public INTTXYvtx{

    public:
        INTTXYvtxEvt(string run_type, string out_folder_directory, pair<double,double> beam_origin, int geo_mode_id, double phi_diff_cut = 0.11, pair<double, double> DCA_cut = {-1,1}, int N_clu_cutl = 20, int N_clu_cut = 10000, bool draw_event_display = true, double peek = 3.32405, double angle_diff_new_l = 0.0, double angle_diff_new_r = 3, bool print_message_opt = true, double window_size = 2.5, int quadrant_trial = 9):
        INTTXYvtx(run_type, out_folder_directory, beam_origin, geo_mode_id, phi_diff_cut, DCA_cut, N_clu_cutl, N_clu_cut, draw_event_display, peek, angle_diff_new_l, angle_diff_new_r, print_message_opt), window_size(window_size), quadrant_trial(quadrant_trial)
        {
            good_tracklet_xy_vec.clear();
            good_tracklet_rz_vec.clear();
            cluster_pair_vec.clear();
            evt_vtx_info.clear();
            InitHist_Evt();
            InitGraphEvt();
            if (draw_event_display == true) {c2 -> Print(Form("%s/temp_event_display.pdf(",out_folder_directory.c_str()));}
        };

        INTTXYvtxEvt(string run_type, string out_folder_directory):INTTXYvtx(run_type, out_folder_directory)
        {
            return;
        };

        // note : the function that new created for this class
        virtual void ProcessEvt(int event_i, vector<clu_info> temp_sPH_inner_nocolumn_vec, vector<clu_info> temp_sPH_outer_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_rz_vec, int NvtxMC, vector<double> TrigvtxMC, bool PhiCheckTag, Long64_t bco_full, pair<double,double> evt_z);
        pair<double,double> GetVtxXYEvt();
        void PrintPlots_Evt();
        void InitHist_Evt();
        vector<pair<double,double>> GetEvtVtxInfo();
        int GetReturnTag();

        // note : the function that inherited from the INTTXYvtx class
        void ClearEvt() override;
        void EndRun() override;

    protected:
        // note : these two functions are from the INTTZvtx class, and I copied them here, not so smart sigh...
        pair<double,double> Get_possible_zvtx(double rvtx, vector<double> p0, vector<double> p1);
        double Get_extrapolation(double given_y, double p0x, double p0y, double p1x, double p1y);

        // note : functions inherited from INTTXYvtx class and would like to modify
        vector<pair<double,double>> MacroVTXSquare(double length, int N_trial, bool draw_plot_opt = true, vector<double> TrigvtxMC = {});
        void subMacroPlotWorking(bool phi_correction, double cos_fit_rangel, double cos_fit_ranger, double guas_fit_range) override;
        
        // note : for the event display
        void temp_bkg(TPad * c1, double peek, pair<double,double> beam_origin, InttConversion * ch_pos_DB);
        double get_z_vertex(clu_info inner_clu, clu_info outer_clu, double target_x, double target_y);
        void InitGraphEvt();
        bool isPointInsideSquare(const std::pair<double, double> point, const std::pair<double, double> square_center, double length);
        void TH2LineFill(TH2F * hist_in, std::pair<double,double> point_1, std::pair<double,double> point_2);
        TCanvas * c2;
        TPad * pad_xy_hist; 
        TPad * pad_xy_hist_bkgrm;
        TPad * pad_xy;
        TPad * pad_rz;
        TPad * pad_dca_inner;
        TPad * pad_dca_outer;
        TPad * pad_delta_phi_inner;
        TPad * pad_delta_phi_outer;
        TPad * pad_inner_outer_phi;
        TPad * pad_track_phi_1D;
        TPad * pad_delta_phi_1D;

        TH2F * evt_display_xy_hist;
        TH2F * evt_display_xy_hist_bkgrm;
        TH2F * evt_dca_inner_2D;
        TH2F * evt_dca_outer_2D;
        TH2F * evt_delta_phi_inner_2D;
        TH2F * evt_delta_phi_outer_2D;
        TH2F * evt_inner_outer_phi_2D;
        TH1F * evt_track_phi_1D;
        TH1F * evt_delta_phi_1D;
        TGraph * evt_display_xy_gr;
        TGraph * evt_display_rz_gr;
        TGraph * true_vertex_gr;
        TGraph * reco_vertex_gr;

        InttConversion * ch_pos_DB;
        TGraph * bkg;
        int return_tag;
        double reco_vtx_x;
        double reco_vtx_y;
        double reco_vtx_xE;
        double reco_vtx_yE;
        double reco_vtx_xStdDev;
        double reco_vtx_yStdDev;


        // note : new created functions only for this class
        // void Init_Hist_Evt();

        // note : new variables declared in this class
        double window_size; 
        int quadrant_trial;
        vector<pair<double,double>> evt_vtx_info;
        vector< pair< pair<double,double>, pair<double,double> > > good_tracklet_xy_vec;
        vector< pair< pair<double,double>, pair<double,double> > > good_tracklet_rz_vec;

        // note : new generated plots only for this class
        TH2F * VTXx_correlation;
        TH2F * VTXy_correlation;

        TH1F * VTXx_diff_1D;
        TH1F * VTXy_diff_1D;

        TH2F * VTXx_diff_Nclus;
        TH2F * VTXy_diff_Nclus;
};



void INTTXYvtxEvt::InitGraphEvt()
{
    c2 = new TCanvas("","",4000,2400);    
    c2 -> cd();

    pad_xy = new TPad(Form("pad_xy"), "", 0.0, 0.66, 0.2, 1.0);
    Characterize_Pad(pad_xy, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_xy -> Draw();

    pad_rz = new TPad(Form("pad_rz"), "", 0.2, 0.66, 0.40, 1.0);
    Characterize_Pad(pad_rz, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_rz -> Draw();

    pad_inner_outer_phi = new TPad(Form("pad_inner_outer_phi"), "", 0.40, 0.66, 0.6, 1.0);
    Characterize_Pad(pad_inner_outer_phi, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_inner_outer_phi -> Draw();
    
    pad_track_phi_1D = new TPad(Form("pad_track_phi_1D"), "", 0.6, 0.66, 0.8, 1.0);
    Characterize_Pad(pad_track_phi_1D, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_track_phi_1D -> Draw();

    pad_delta_phi_1D = new TPad(Form("pad_delta_phi_1D"), "", 0.8, 0.66, 1, 1.0);
    Characterize_Pad(pad_delta_phi_1D, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_delta_phi_1D -> Draw();

    pad_dca_inner = new TPad(Form("pad_dca_inner"), "", 0.0, 0.33, 0.2, 0.66);
    Characterize_Pad(pad_dca_inner, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_dca_inner -> Draw();

    pad_dca_outer = new TPad(Form("pad_dca_outer"), "", 0.2, 0.33, 0.40, 0.66);
    Characterize_Pad(pad_dca_outer, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_dca_outer -> Draw();

    pad_delta_phi_inner = new TPad(Form("pad_delta_phi_inner"), "", 0.4, 0.33, 0.60, 0.66);
    Characterize_Pad(pad_delta_phi_inner, 0.15, 0.1, 0.1, 0.2, 0, 0);
    // pad_delta_phi_inner -> SetLogz(1);
    pad_delta_phi_inner -> Draw();

    pad_delta_phi_outer = new TPad(Form("pad_delta_phi_outer"), "", 0.6, 0.33, 0.80, 0.66);
    Characterize_Pad(pad_delta_phi_outer, 0.15, 0.1, 0.1, 0.2, 0, 0);
    // pad_delta_phi_outer -> SetLogz(1);
    pad_delta_phi_outer -> Draw();

    pad_xy_hist = new TPad(Form("pad_xy_hist"), "", 0.8, 0.33, 1.0, 0.66);
    Characterize_Pad(pad_xy_hist, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_xy_hist -> Draw();

    pad_xy_hist_bkgrm = new TPad(Form("pad_xy_hist_bkgrm"), "", 0.0, 0.0, 0.2, 0.33);
    Characterize_Pad(pad_xy_hist_bkgrm, 0.15, 0.1, 0.1, 0.2, 0, 0);
    pad_xy_hist_bkgrm -> Draw();

    evt_display_xy_gr = new TGraph();
    evt_display_rz_gr = new TGraph();

    true_vertex_gr    = new TGraph();
    true_vertex_gr -> SetMarkerStyle(50);
    true_vertex_gr -> SetMarkerColor(2);
    true_vertex_gr -> SetMarkerSize(0.5);

    reco_vertex_gr    = new TGraph();
    reco_vertex_gr -> SetMarkerStyle(50);
    reco_vertex_gr -> SetMarkerColor(4);
    reco_vertex_gr -> SetMarkerSize(0.5);


    bkg = new TGraph();
    bkg -> SetPoint(0,0,0);

    ch_pos_DB = new InttConversion("ideal", peek); // todo : I make it simplified here, for the setting of the geometry modes.
}

void INTTXYvtxEvt::InitHist_Evt()
{
    VTXx_correlation = new TH2F("","VTXx_correlation",100,-1,0.5,100,-1,0.5);
    VTXx_correlation -> SetStats(0);
    VTXx_correlation -> GetXaxis() -> SetTitle("True VTX_{x} [mm]");
    VTXx_correlation -> GetYaxis() -> SetTitle("Reco VTX_{x} [mm]");
    VTXx_correlation -> GetXaxis() -> SetNdivisions(505);

    VTXy_correlation = new TH2F("","VTXy_correlation",100,1.7,3,100,1.7,3);
    VTXy_correlation -> SetStats(0);
    VTXy_correlation -> GetXaxis() -> SetTitle("True VTX_{y} [mm]");
    VTXy_correlation -> GetYaxis() -> SetTitle("Reco VTX_{y} [mm]");
    VTXy_correlation -> GetXaxis() -> SetNdivisions(505);

    VTXx_diff_1D = new TH1F("","VTXx_diff_1D",100,-0.5,0.5);
    VTXx_diff_1D -> SetStats(0);
    VTXx_diff_1D -> GetXaxis() -> SetTitle("#DeltaX (Reco - True) [mm]");
    VTXx_diff_1D -> GetYaxis() -> SetTitle("Entry");
    VTXx_diff_1D -> GetXaxis() -> SetNdivisions(505);

    VTXy_diff_1D = new TH1F("","VTXy_diff_1D",100,-0.5,0.5);
    VTXy_diff_1D -> SetStats(0);
    VTXy_diff_1D -> GetXaxis() -> SetTitle("#DeltaY (Reco - True) [mm]");
    VTXy_diff_1D -> GetYaxis() -> SetTitle("Entry");
    VTXy_diff_1D -> GetXaxis() -> SetNdivisions(505);

    VTXx_diff_Nclus = new TH2F("","VTXx_diff_Nclus",100,0,8000,100,-1,1);
    VTXx_diff_Nclus -> SetStats(0);
    VTXx_diff_Nclus -> GetXaxis() -> SetTitle("# of clusters");
    VTXx_diff_Nclus -> GetYaxis() -> SetTitle("#DeltaX (Reco - True) [mm]");
    VTXx_diff_Nclus -> GetXaxis() -> SetNdivisions(505);

    VTXy_diff_Nclus = new TH2F("","VTXy_diff_Nclus",100,0,8000,100,-1,1);
    VTXy_diff_Nclus -> SetStats(0);
    VTXy_diff_Nclus -> GetXaxis() -> SetTitle("# of clusters");
    VTXy_diff_Nclus -> GetYaxis() -> SetTitle("#DeltaY (Reco - True) [mm]");
    VTXy_diff_Nclus -> GetXaxis() -> SetNdivisions(505);


    // note : the followings are for the event display
    // note : the event display window is now moved to the nominoal beam position
    evt_display_xy_hist = new TH2F("","evt_display_xy_hist", 100, -1.5 + beam_origin.first, 1.5 + beam_origin.first, 100, -1.5 + beam_origin.second, 1.5 + beam_origin.second);
    evt_display_xy_hist -> SetStats(0);
    evt_display_xy_hist -> GetXaxis() -> SetTitle("X axis [mm]");
    evt_display_xy_hist -> GetYaxis() -> SetTitle("Y axis [mm]");
    evt_display_xy_hist -> GetXaxis() -> SetNdivisions(505);

    evt_display_xy_hist_bkgrm = new TH2F("","evt_display_xy_hist_bkgrm", 100, -1.5 + beam_origin.first, 1.5 + beam_origin.first, 100, -1.5 + beam_origin.second, 1.5 + beam_origin.second);
    evt_display_xy_hist_bkgrm -> SetStats(0);
    evt_display_xy_hist_bkgrm -> GetXaxis() -> SetTitle("X axis [mm]");
    evt_display_xy_hist_bkgrm -> GetYaxis() -> SetTitle("Y axis [mm]");
    evt_display_xy_hist_bkgrm -> GetXaxis() -> SetNdivisions(505);

    evt_dca_inner_2D = new TH2F("","evt_dca_inner_2D",100,0,360,100,-10,10);
    evt_dca_inner_2D -> SetStats(0);
    evt_dca_inner_2D -> GetXaxis() -> SetTitle("Inner cluster #phi [degree]");
    evt_dca_inner_2D -> GetYaxis() -> SetTitle("DCA [mm]");
    evt_dca_inner_2D -> GetXaxis() -> SetNdivisions(505);

    evt_dca_outer_2D = new TH2F("","evt_dca_outer_2D",100,0,360,100,-10,10);
    evt_dca_outer_2D -> SetStats(0);
    evt_dca_outer_2D -> GetXaxis() -> SetTitle("Outer cluster #phi [degree]");
    evt_dca_outer_2D -> GetYaxis() -> SetTitle("DCA [mm]");
    evt_dca_outer_2D -> GetXaxis() -> SetNdivisions(505);

    evt_delta_phi_inner_2D = new TH2F("","evt_delta_phi_inner_2D",361,0,361,100,-1.5,1.5);
    evt_delta_phi_inner_2D -> SetStats(0);
    evt_delta_phi_inner_2D -> GetXaxis() -> SetTitle("Inner cluster #phi [degree]");
    evt_delta_phi_inner_2D -> GetYaxis() -> SetTitle("#Delta#phi (inner - outer) [degree]");
    evt_delta_phi_inner_2D -> GetXaxis() -> SetNdivisions(505);

    evt_delta_phi_outer_2D = new TH2F("","evt_delta_phi_outer_2D",180,0,360,100,-1.5,1.5);
    evt_delta_phi_outer_2D -> SetStats(0);
    evt_delta_phi_outer_2D -> GetXaxis() -> SetTitle("Outer cluster #phi [degree]");
    evt_delta_phi_outer_2D -> GetYaxis() -> SetTitle("#Delta#phi (inner - outer) [degree]");
    evt_delta_phi_outer_2D -> GetXaxis() -> SetNdivisions(505);

    evt_inner_outer_phi_2D = new TH2F("","evt_inner_outer_phi_2D",180,0,360,361,0,361);
    evt_inner_outer_phi_2D -> SetStats(0);
    evt_inner_outer_phi_2D -> GetXaxis() -> SetTitle("Inner cluster #phi [degree]");
    evt_inner_outer_phi_2D -> GetYaxis() -> SetTitle("Outer cluster #phi [degree]");
    evt_inner_outer_phi_2D -> GetXaxis() -> SetNdivisions(505);

    evt_track_phi_1D = new TH1F("","evt_track_phi_1D",180,0,360);
    evt_track_phi_1D -> SetStats(0);
    evt_track_phi_1D -> GetXaxis() -> SetTitle("Tracklet #phi [degree]");
    evt_track_phi_1D -> GetYaxis() -> SetTitle("Entry");
    evt_track_phi_1D -> GetXaxis() -> SetNdivisions(505);

    evt_delta_phi_1D = new TH1F("","evt_delta_phi_1D",100,-10,10);
    evt_delta_phi_1D -> SetStats(0);
    evt_delta_phi_1D -> GetXaxis() -> SetTitle("#Delta#phi (inner - outer) [degree]");
    evt_delta_phi_1D -> GetYaxis() -> SetTitle("Entry");
    evt_delta_phi_1D -> GetXaxis() -> SetNdivisions(505);
}

void INTTXYvtxEvt::PrintPlots_Evt()
{
    c1 -> cd();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    VTXx_correlation -> Draw("colz0");
    VTXx_correlation -> Fit(correlation_fit_MC, "NQ");
    correlation_fit_MC -> Draw("l same");
    draw_text -> DrawLatex(0.21, 0.71, Form("y = %.2f * X +  %.2f",correlation_fit_MC -> GetParameter(1), correlation_fit_MC -> GetParameter(0)));
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print(Form("%s/VTXx_correlation.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    VTXy_correlation -> Draw("colz0");
    VTXy_correlation -> Fit(correlation_fit_MC, "NQ");
    correlation_fit_MC -> Draw("l same");
    draw_text -> DrawLatex(0.21, 0.71, Form("y = %.2f * X +  %.2f",correlation_fit_MC -> GetParameter(1), correlation_fit_MC -> GetParameter(0)));
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print(Form("%s/VTXy_correlation.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    VTXx_diff_1D -> Draw("hist");
    gaus_fit_MC -> SetParameters(VTXx_diff_1D -> GetBinContent( VTXx_diff_1D -> GetMaximumBin() ), VTXx_diff_1D -> GetBinCenter( VTXx_diff_1D -> GetMaximumBin() ), 0.05, 0);
    gaus_fit_MC -> SetParLimits(0,0,100000);  // note : size 
    gaus_fit_MC -> SetParLimits(2,0,10000);   // note : Width
    gaus_fit_MC -> SetParLimits(3,0,10000);   // note : offset
    VTXx_diff_1D -> Fit(gaus_fit_MC, "NQ", "", VTXx_diff_1D -> GetBinCenter( VTXx_diff_1D -> GetMaximumBin() ) - (2 * VTXx_diff_1D -> GetStdDev() ), VTXx_diff_1D -> GetBinCenter( VTXx_diff_1D -> GetMaximumBin() ) + (2 * VTXx_diff_1D -> GetStdDev() ) );
    gaus_fit_MC -> SetRange( gaus_fit_MC->GetParameter(1) - gaus_fit_MC->GetParameter(2) * 2.5, gaus_fit_MC->GetParameter(1) + gaus_fit_MC->GetParameter(2) * 2.5 ); 
    gaus_fit_MC -> Draw("lsame");
    draw_text -> DrawLatex(0.21, 0.71, Form("Gaus mean  : %.2f mm",gaus_fit_MC -> GetParameter(1)));
    draw_text -> DrawLatex(0.21, 0.67, Form("Gaus width : %.2f mm",fabs(gaus_fit_MC -> GetParameter(2))));
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print(Form("%s/VTXx_diff_1D.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    VTXy_diff_1D -> Draw("hist");
    gaus_fit_MC -> SetParameters(VTXy_diff_1D -> GetBinContent( VTXy_diff_1D -> GetMaximumBin() ), VTXy_diff_1D -> GetBinCenter( VTXy_diff_1D -> GetMaximumBin() ), 0.05, 0);
    gaus_fit_MC -> SetParLimits(0,0,100000);  // note : size 
    gaus_fit_MC -> SetParLimits(2,0,10000);   // note : Width
    gaus_fit_MC -> SetParLimits(3,0,10000);   // note : offset
    VTXy_diff_1D -> Fit(gaus_fit_MC, "NQ", "", VTXy_diff_1D -> GetBinCenter( VTXy_diff_1D -> GetMaximumBin() ) - (2 * VTXy_diff_1D -> GetStdDev() ), VTXy_diff_1D -> GetBinCenter( VTXy_diff_1D -> GetMaximumBin() ) + (2 * VTXy_diff_1D -> GetStdDev() ) );
    gaus_fit_MC -> SetRange( gaus_fit_MC->GetParameter(1) - gaus_fit_MC->GetParameter(2) * 2.5, gaus_fit_MC->GetParameter(1) + gaus_fit_MC->GetParameter(2) * 2.5 ); 
    gaus_fit_MC -> Draw("lsame");
    draw_text -> DrawLatex(0.21, 0.71, Form("Gaus mean  : %.2f mm",gaus_fit_MC -> GetParameter(1)));
    draw_text -> DrawLatex(0.21, 0.67, Form("Gaus width : %.2f mm",fabs(gaus_fit_MC -> GetParameter(2))));
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print(Form("%s/VTXy_diff_1D.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    VTXx_diff_Nclus -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print(Form("%s/VTXx_diff_Nclus.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    VTXy_diff_Nclus -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} %s", plot_text.c_str()));
    c1 -> Print(Form("%s/VTXy_diff_Nclus.pdf",out_folder_directory.c_str()));
    c1 -> Clear();
}

// note : this function is for the event by event vertex calculation
void INTTXYvtxEvt::ProcessEvt(int event_i, vector<clu_info> temp_sPH_inner_nocolumn_vec, vector<clu_info> temp_sPH_outer_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_vec, vector<vector<double>> temp_sPH_nocolumn_rz_vec, int NvtxMC, vector<double> TrigvtxMC, bool PhiCheckTag, Long64_t bco_full, pair<double,double> evt_z) // note : evt_z : {z, width} 
{   
    return_tag = 0;

    if (event_i%1 == 0) {cout<<"In INTTXYvtxEvt class, running event : "<<event_i<<endl;}

    total_NClus = temp_sPH_inner_nocolumn_vec.size() + temp_sPH_outer_nocolumn_vec.size();


    if (total_NClus < zvtx_cal_require) {return; cout<<"return confirmation"<<endl;}   
    
    if (run_type == "MC" && NvtxMC != 1) { return; cout<<"In INTTXYvtxEvt class, event : "<<event_i<<" Nvtx : "<<NvtxMC<<" Nvtx more than one "<<endl;}
    if (PhiCheckTag == false)            { return; cout<<"In INTTXYvtxEvt class, event : "<<event_i<<" Nvtx : "<<NvtxMC<<" Not full phi has hits "<<endl;}
    
    if (temp_sPH_inner_nocolumn_vec.size() < 10 || temp_sPH_outer_nocolumn_vec.size() < 10 || total_NClus > N_clu_cut || total_NClus < N_clu_cutl)
    {
        return;
        printf("In INTTXYvtxEvt class, event : %i, low clu continue, NClus : %lu \n", event_i, total_NClus); 
    }
    
    // cout<<"test_1"<<endl;

    for ( int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++ ){ 
        evt_display_xy_gr -> SetPoint(evt_display_xy_gr -> GetN(), temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y);
        evt_display_rz_gr -> SetPoint(evt_display_rz_gr -> GetN(), temp_sPH_inner_nocolumn_vec[inner_i].z, (temp_sPH_inner_nocolumn_vec[inner_i].phi > 180) ? get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x,temp_sPH_inner_nocolumn_vec[inner_i].y) * -1 : get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x,temp_sPH_inner_nocolumn_vec[inner_i].y));
    }   
    // cout<<"test_2"<<endl;
    for ( int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++ ){
        evt_display_xy_gr -> SetPoint(evt_display_xy_gr -> GetN(), temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y);
        evt_display_rz_gr -> SetPoint(evt_display_rz_gr -> GetN(), temp_sPH_outer_nocolumn_vec[outer_i].z, (temp_sPH_outer_nocolumn_vec[outer_i].phi > 180) ? get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x,temp_sPH_outer_nocolumn_vec[outer_i].y) * -1 : get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x,temp_sPH_outer_nocolumn_vec[outer_i].y));
    }
    // cout<<"test_3"<<endl;
    for ( int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++ )
    {    
        // cout<<"test_4"<<endl;
        for ( int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++ )
        {
            // cout<<event_i<<" test_5 "<<inner_i<<" "<<outer_i<<endl;
            // note : calculate the cluster phi after the vertex correction which can enhence the purity of the tracklet selection
            Clus_InnerPhi_Offset = (temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second < 0) ? atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second, temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first) * (180./TMath::Pi());
            Clus_OuterPhi_Offset = (temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second < 0) ? atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi()) + 360 : atan2(temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second, temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first) * (180./TMath::Pi());
            
            // note : before calculating the possible z vertex range of one tracklet, the vertex correction was applied
            pair<double,double> z_range_info = Get_possible_zvtx( 
                0., // get_radius(beam_origin.first,beam_origin.second), 
                {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z}, // note : unsign radius
                {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}  // note : unsign radius
            );
            // cout<<event_i<<" test_6 "<<inner_i<<" "<<outer_i<<endl;
            // note : this is a cut to constraint on the z vertex, only if the tracklets with the range that covers the z vertex can pass this cut
            if (z_range_info.first - z_range_info.second > evt_z.first + evt_z.second || z_range_info.first + z_range_info.second < evt_z.first - evt_z.second) {continue;}
            // cout<<"range test: "<<z_range_info.first - z_range_info.second<<" "<<z_range_info.first + z_range_info.second<<" vertex: "<<evt_z.first<<" "<<evt_z.second<<endl;

            // note : the reason to use the DCA calculation with sign is because that the distribution of 1D DCA will be single peak only
            double DCA_sign = calculateAngleBetweenVectors(
                temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y,
                temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y,
                beam_origin.first, beam_origin.second
            );

            vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
                temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y,
                temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y,
                beam_origin.first, beam_origin.second
            );

            if (DCA_info_vec[0] != fabs(DCA_sign) && fabs( DCA_info_vec[0] - fabs(DCA_sign) ) > 0.1 ){
                cout<<"different DCA : "<<DCA_info_vec[0]<<" "<<DCA_sign<<" diff : "<<DCA_info_vec[0] - fabs(DCA_sign)<<endl;
            }
            // cout<<event_i<<" test_7 "<<inner_i<<" "<<outer_i<<endl;
            evt_dca_inner_2D -> Fill(temp_sPH_inner_nocolumn_vec[inner_i].phi, DCA_sign);
            evt_dca_outer_2D -> Fill(temp_sPH_outer_nocolumn_vec[outer_i].phi, DCA_sign);
            evt_delta_phi_inner_2D -> Fill(temp_sPH_inner_nocolumn_vec[inner_i].phi, temp_sPH_inner_nocolumn_vec[inner_i].phi - temp_sPH_outer_nocolumn_vec[outer_i].phi);
            evt_delta_phi_outer_2D -> Fill(temp_sPH_outer_nocolumn_vec[outer_i].phi, temp_sPH_inner_nocolumn_vec[inner_i].phi - temp_sPH_outer_nocolumn_vec[outer_i].phi);
            evt_inner_outer_phi_2D -> Fill(temp_sPH_inner_nocolumn_vec[inner_i].phi, temp_sPH_outer_nocolumn_vec[outer_i].phi);
            evt_delta_phi_1D -> Fill( temp_sPH_inner_nocolumn_vec[inner_i].phi - temp_sPH_outer_nocolumn_vec[outer_i].phi );

            if (fabs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset) < phi_diff_cut)
            {
                if (DCA_cut.first < DCA_sign && DCA_sign < DCA_cut.second){
                    
                    // // note : before calculating the possible z vertex range of one tracklet, the vertex correction was applied
                    // pair<double,double> z_range_info = Get_possible_zvtx( 
                    //     0., // get_radius(beam_origin.first,beam_origin.second), 
                    //     {get_radius(temp_sPH_inner_nocolumn_vec[inner_i].x - beam_origin.first, temp_sPH_inner_nocolumn_vec[inner_i].y - beam_origin.second), temp_sPH_inner_nocolumn_vec[inner_i].z}, // note : unsign radius
                    //     {get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x - beam_origin.first, temp_sPH_outer_nocolumn_vec[outer_i].y - beam_origin.second), temp_sPH_outer_nocolumn_vec[outer_i].z}  // note : unsign radius
                    // );

                    // // note : this is a cut to constraint on the z vertex, only if the tracklets with the range that covers the z vertex can pass this cut
                    // if (z_range_info.first - z_range_info.second > evt_z.first + evt_z.second || z_range_info.first + z_range_info.second < evt_z.first - evt_z.second) {continue;}
                    // // cout<<"range test: "<<z_range_info.first - z_range_info.second<<" "<<z_range_info.first + z_range_info.second<<" vertex: "<<evt_z.first<<" "<<evt_z.second<<endl;
                    
                    // cout<<event_i<<" test_8 "<<inner_i<<" "<<outer_i<<endl;

                    // note : for the tracklets that pass all the cut listed above can then be used for the vertex XY calculation
                    cluster_pair_vec.push_back({{temp_sPH_inner_nocolumn_vec[inner_i].x,temp_sPH_inner_nocolumn_vec[inner_i].y},{temp_sPH_outer_nocolumn_vec[outer_i].x,temp_sPH_outer_nocolumn_vec[outer_i].y}});
                    // cout<<"we have something here "<<temp_sPH_inner_nocolumn_vec[inner_i].x<<" "<<temp_sPH_inner_nocolumn_vec[inner_i].y<<endl;
                    
                    good_tracklet_xy_vec.push_back({
                        {temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y},
                        {DCA_info_vec[1], DCA_info_vec[2]}
                    });

                    good_tracklet_rz_vec.push_back({
                        {temp_sPH_outer_nocolumn_vec[outer_i].z, (temp_sPH_outer_nocolumn_vec[outer_i].phi > 180) ? get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x,temp_sPH_outer_nocolumn_vec[outer_i].y) * -1 : get_radius(temp_sPH_outer_nocolumn_vec[outer_i].x,temp_sPH_outer_nocolumn_vec[outer_i].y)},
                        {get_z_vertex(temp_sPH_inner_nocolumn_vec[inner_i],temp_sPH_outer_nocolumn_vec[outer_i],DCA_info_vec[1],DCA_info_vec[2]), (temp_sPH_outer_nocolumn_vec[outer_i].phi > 180) ? get_radius(DCA_info_vec[1],DCA_info_vec[2]) * -1 : get_radius(DCA_info_vec[1],DCA_info_vec[2])}
                    });

                    // cout<<event_i<<" test_9 "<<inner_i<<" "<<outer_i<<endl;

                    TH2FSampleLineFill(evt_display_xy_hist, 0.001, {temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y}, {DCA_info_vec[1], DCA_info_vec[2]});
                    // TH2LineFill(evt_display_xy_hist, {temp_sPH_inner_nocolumn_vec[inner_i].x, temp_sPH_inner_nocolumn_vec[inner_i].y}, {temp_sPH_outer_nocolumn_vec[outer_i].x, temp_sPH_outer_nocolumn_vec[outer_i].y});
                    evt_track_phi_1D -> Fill((temp_sPH_inner_nocolumn_vec[inner_i].phi + temp_sPH_outer_nocolumn_vec[outer_i].phi)/2.);

                    // cout<<event_i<<" test_10 "<<inner_i<<" "<<outer_i<<endl;
                }

            }
            
        }
    }

    // note : try to remove the background
    TH2F_FakeClone(evt_display_xy_hist,evt_display_xy_hist_bkgrm);
    TH2F_threshold_advanced_2(evt_display_xy_hist_bkgrm, 0.7);
    
    // cout<<"--- In process event: "<<event_i<<" "<<total_NClus<<" "<<cluster_pair_vec.size()<<endl;

    // evt_vtx_info = MacroVTXSquare(window_size, quadrant_trial, false, TrigvtxMC);
    // evt_vtx_info = {{1,2},{1,2}};
    // double reco_vtx_x = evt_vtx_info[0].first;
    // double reco_vtx_y = evt_vtx_info[0].second;

    // delete DCA_distance_inner_phi_peak_final;
    // delete angle_diff_inner_phi_peak_final;

    reco_vtx_x = evt_display_xy_hist_bkgrm->GetMean(1) + (evt_display_xy_hist_bkgrm -> GetXaxis() -> GetBinWidth(1) / 2.);
    reco_vtx_y = evt_display_xy_hist_bkgrm->GetMean(2) + (evt_display_xy_hist_bkgrm -> GetYaxis() -> GetBinWidth(1) / 2.);
    reco_vtx_xE = evt_display_xy_hist_bkgrm->GetMeanError(1);
    reco_vtx_yE = evt_display_xy_hist_bkgrm->GetMeanError(2);
    reco_vtx_xStdDev = evt_display_xy_hist_bkgrm->GetStdDev(1);
    reco_vtx_yStdDev = evt_display_xy_hist_bkgrm->GetStdDev(2);


    if (total_NClus > 2500)
    {
        VTXx_correlation -> Fill(TrigvtxMC[0]*10., reco_vtx_x);
        VTXy_correlation -> Fill(TrigvtxMC[1]*10., reco_vtx_y);
        VTXx_diff_1D -> Fill(reco_vtx_x - TrigvtxMC[0]*10.);
        VTXy_diff_1D -> Fill(reco_vtx_y - TrigvtxMC[1]*10.);
    }
    
    VTXx_diff_Nclus -> Fill(total_NClus, reco_vtx_x - TrigvtxMC[0]*10.);
    VTXy_diff_Nclus -> Fill(total_NClus, reco_vtx_y - TrigvtxMC[1]*10.);
    
    cout<<"note, event: "<<event_i<<"----------------------- ----------------------- ----------------------- ----------------------- -----------------------"<<endl;
    cout<<"test: "<<reco_vtx_x<<" "<<reco_vtx_y<<endl;
    cout<<"test: "<<TrigvtxMC[0]*10.<<" "<<TrigvtxMC[1]*10.<<endl;
    cout<<"test: deviation "<<reco_vtx_x - TrigvtxMC[0]*10.<<" "<<reco_vtx_y - TrigvtxMC[1]*10.<<endl;
    cout<<"test: distance "<<sqrt(pow(reco_vtx_x - TrigvtxMC[0]*10.,2)+pow(reco_vtx_y - TrigvtxMC[1]*10.,2))<<endl;
    
    return_tag = 1;

    if (draw_event_display == true)
    {
        c2 -> cd();
        // evt_display_xy_gr = new TGraph(temp_sPH_nocolumn_vec[0].size(),&temp_sPH_nocolumn_vec[0][0],&temp_sPH_nocolumn_vec[1][0]);
        evt_display_xy_gr -> SetTitle("INTT event display X-Y plane");
        evt_display_xy_gr -> GetXaxis() -> SetLimits(-150,150);
        evt_display_xy_gr -> GetYaxis() -> SetRangeUser(-150,150);
        evt_display_xy_gr -> GetXaxis() -> SetTitle("X [mm]");
        evt_display_xy_gr -> GetYaxis() -> SetTitle("Y [mm]");
        evt_display_xy_gr -> SetMarkerStyle(20);
        evt_display_xy_gr -> SetMarkerColor(2);
        evt_display_xy_gr -> SetMarkerSize(1);

        // evt_display_rz_gr = new TGraph(temp_sPH_nocolumn_rz_vec[0].size(),&temp_sPH_nocolumn_rz_vec[0][0],&temp_sPH_nocolumn_rz_vec[1][0]);
        evt_display_rz_gr -> SetTitle("INTT event display r-Z plane");
        evt_display_rz_gr -> GetXaxis() -> SetLimits(-500,500);
        evt_display_rz_gr -> GetYaxis() -> SetRangeUser(-150,150);
        evt_display_rz_gr -> GetXaxis() -> SetTitle("Z [mm]");
        evt_display_rz_gr -> GetYaxis() -> SetTitle("Radius [mm]");
        evt_display_rz_gr -> SetMarkerStyle(20);
        evt_display_rz_gr -> SetMarkerColor(2);
        evt_display_rz_gr -> SetMarkerSize(1);

        // note : ----------------------- ----------------------- ----------------------- ----------------------- -----------------------
        pad_xy -> cd();
        temp_bkg(pad_xy, peek, beam_origin, ch_pos_DB);
        evt_display_xy_gr -> Draw("p same");
        for (int track_i = 0; track_i < good_tracklet_xy_vec.size(); track_i++){
            track_line -> DrawLine(
                good_tracklet_xy_vec[track_i].first.first, good_tracklet_xy_vec[track_i].first.second, 
                good_tracklet_xy_vec[track_i].second.first, good_tracklet_xy_vec[track_i].second.second
            );
        }
        // track_line -> Draw("l same");
        draw_text -> DrawLatex(0.2, 0.85, Form("eID : %i, inner Ncluster : %zu, outer Ncluster : %zu",event_i,temp_sPH_inner_nocolumn_vec.size(),temp_sPH_outer_nocolumn_vec.size()));

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_rz -> cd();
        evt_display_rz_gr -> Draw("ap");    
        // eff_sig_range_line -> DrawLine(temp_event_zvtx_info[0],-150,temp_event_zvtx_info[0],150);
        coord_line -> DrawLine(0,-150,0,150);
        coord_line -> DrawLine(-500,0,500,0);
        for (int track_i = 0; track_i < good_tracklet_rz_vec.size(); track_i++){
            track_line -> DrawLine(
                good_tracklet_rz_vec[track_i].first.first, good_tracklet_rz_vec[track_i].first.second, 
                good_tracklet_rz_vec[track_i].second.first, good_tracklet_rz_vec[track_i].second.second
            );
        }
        draw_text -> DrawLatex(0.2, 0.85, Form("Negative radius : Clu_{outer} > 180^{0}"));

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_dca_inner -> cd();
        evt_dca_inner_2D -> Draw("colz0");

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_dca_outer -> cd();
        evt_dca_outer_2D -> Draw("colz0");

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_delta_phi_inner -> cd();
        evt_delta_phi_inner_2D -> Draw("colz0");

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_delta_phi_outer -> cd();
        evt_delta_phi_outer_2D -> Draw("colz0");

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_inner_outer_phi -> cd();
        evt_inner_outer_phi_2D -> Draw("colz0");

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_track_phi_1D -> cd();
        evt_track_phi_1D -> Draw("hist");

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_delta_phi_1D -> cd();
        evt_delta_phi_1D -> Draw("hist");

        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_xy_hist -> cd();
        evt_display_xy_hist -> Draw("colz0");
        true_vertex_gr -> SetPoint(true_vertex_gr->GetN(), TrigvtxMC[0] * 10, TrigvtxMC[1] * 10);
        true_vertex_gr -> Draw("psame");
        // cout<<event_i<<" test_1 "<<endl;
        // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pad_xy_hist_bkgrm -> cd();
        evt_display_xy_hist_bkgrm -> Draw("colz0");
        reco_vertex_gr -> SetPoint(reco_vertex_gr->GetN(), reco_vtx_x , reco_vtx_y);
        true_vertex_gr -> Draw("psame");
        reco_vertex_gr -> Draw("psame");
        // cout<<event_i<<" test_2 "<<endl;

        if(draw_event_display /* && (event_i % print_rate) == 0*/ /*&& good_zvtx_tag == true*/){c2 -> Print(Form("%s/temp_event_display.pdf",out_folder_directory.c_str()));}

        pad_xy -> Clear();
        pad_rz -> Clear();
        pad_xy_hist -> Clear(); 
        pad_dca_inner -> Clear();
        pad_dca_outer -> Clear();
        pad_delta_phi_inner -> Clear();
        pad_delta_phi_outer -> Clear();
        pad_inner_outer_phi -> Clear();
        pad_track_phi_1D -> Clear();
        pad_delta_phi_1D -> Clear();
    }
}

pair<double,double> INTTXYvtxEvt::GetVtxXYEvt() { return evt_vtx_info[0]; }
vector<pair<double,double>> INTTXYvtxEvt::GetEvtVtxInfo() { return {{reco_vtx_x, reco_vtx_y}, {reco_vtx_xE, reco_vtx_yE}, {reco_vtx_xStdDev, reco_vtx_yStdDev}}; }
int INTTXYvtxEvt::GetReturnTag() { return return_tag; }

void INTTXYvtxEvt::ClearEvt()
{
    good_tracklet_xy_vec.clear();
    good_tracklet_rz_vec.clear();
    cluster_pair_vec.clear();
    if (evt_display_xy_gr -> GetN() != 0) {evt_display_xy_gr -> Set(0);}
    if (evt_display_rz_gr -> GetN() != 0) {evt_display_rz_gr -> Set(0);}

    evt_display_xy_hist -> Reset("ICESM");
    evt_dca_inner_2D -> Reset("ICESM");
    evt_dca_outer_2D -> Reset("ICESM");
    evt_delta_phi_inner_2D -> Reset("ICESM");
    evt_delta_phi_outer_2D -> Reset("ICESM");
    evt_inner_outer_phi_2D -> Reset("ICESM");
    evt_track_phi_1D -> Reset("ICESM");
    evt_delta_phi_1D -> Reset("ICESM");
    evt_display_xy_hist_bkgrm -> Reset("ICESM");

    if (true_vertex_gr -> GetN() != 0) {true_vertex_gr -> Set(0);}
    if (reco_vertex_gr -> GetN() != 0) {reco_vertex_gr -> Set(0);}

    return;
}

vector<pair<double,double>> INTTXYvtxEvt::MacroVTXSquare(double length, int N_trial, bool draw_plot_opt, vector<double> TrigvtxMC)
{
    double original_length = length;
    pair<double,double> origin = {0,0};
    // pair<double,double> origin = {beam_origin.first,beam_origin.second};

    vector<pair<double,double>> vtx_vec = Get4vtx(origin,length); // vtx_vec.push_back(origin);
    int small_index;
    vector<double> small_info_vec;
    
    vector<pair<double,double>> vtx_vec_axis = Get4vtxAxis(origin,length);
    int small_index_axis;
    vector<double> small_info_vec_axis[vtx_vec_axis.size()];
    
    vector<double> grr_x; grr_x.clear();
    vector<double> grr_E; grr_E.clear();
    vector<double> grr_y; grr_y.clear();

    if (print_message_opt == true) {cout<<N_trial<<" runs, smart. which gives you the resolution down to "<<length/pow(2,N_trial)<<" mm"<<endl;}

    for (int i = 0; i < N_trial; i++)
    {
        if (print_message_opt == true) {cout<<"~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~"<<" step "<<i<<" ~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~"<<endl;}
        for (int i1 = 0; i1 < vtx_vec.size(); i1++)
        {
            // if (print_message_opt == true) {cout<<"testd vertex : "<<vtx_vec[i1].first<<" "<<vtx_vec[i1].second<<endl;}
            current_vtxX = vtx_vec[i1].first;
            current_vtxY = vtx_vec[i1].second;
            vector<double> info_vec = subMacroVTXxyCorrection(i,i1, draw_plot_opt);
            if (print_message_opt == true) {cout<<"DCA fit E: "<<info_vec[3]<<" delta phi fit E: "<<info_vec[5]<<endl;}
            if (i1 == 0)
            {
                small_info_vec = info_vec;
                small_index = i1;
                
                // DCA_distance_inner_phi_peak_final = (TH2F*)DCA_distance_inner_phi_peak -> Clone();  
                // angle_diff_inner_phi_peak_final = (TH2F*)angle_diff_inner_phi_peak -> Clone();
            }
            else
            {
                // if (info_vec[3] < small_info_vec[3] && info_vec[5] < small_info_vec[5]) // note : the chi2 of the pol0 fit
                if (info_vec[3] + info_vec[5] * 10 < small_info_vec[3] + small_info_vec[5] * 10)
                {
                    small_info_vec = info_vec;
                    small_index = i1;

                    // DCA_distance_inner_phi_peak_final -> Delete();
                    // angle_diff_inner_phi_peak_final -> Delete();
                    // delete DCA_distance_inner_phi_peak_final;
                    // delete angle_diff_inner_phi_peak_final;
                    // DCA_distance_inner_phi_peak_final = (TH2F*)DCA_distance_inner_phi_peak -> Clone();  
                    // angle_diff_inner_phi_peak_final = (TH2F*)angle_diff_inner_phi_peak -> Clone();
                }
            }
            if (print_message_opt == true){cout<<" "<<endl;}

            ClearHist(1);
        }

        // if (print_message_opt == true) {cout<<"the Quadrant "<<small_index<<" won the competition with the vertex: "<<vtx_vec[small_index].first<<" "<<vtx_vec[small_index].second<<endl;}
        
        grr_x.push_back(i);
        grr_y.push_back(small_index);
        grr_E.push_back(0);

        // note : the for loop for the case of the verteix on the axis
        for (int i1 = 0; i1 < vtx_vec_axis.size(); i1++)
        {
            current_vtxX = vtx_vec_axis[i1].first;
            current_vtxY = vtx_vec_axis[i1].second;
            small_info_vec_axis[i1] = subMacroVTXxyCorrection(i,i1, draw_plot_opt);

            ClearHist(1);
        }

        int winner_quadrant_axis = axis_quadrant_map[{
            (small_info_vec_axis[0][3] + small_info_vec_axis[0][5] * 10 < small_info_vec_axis[1][3] + small_info_vec_axis[1][5] * 10) ? 1 : -1,
            (small_info_vec_axis[2][3] + small_info_vec_axis[2][5] * 10 < small_info_vec_axis[3][3] + small_info_vec_axis[3][5] * 10) ? 1 : -1
        }];

        if (true == true)
        {
            // cout<<" "<<endl;
            cout<<"------------------------------------------------------------------------------------"<<endl;
            cout<<"for axis method: "<<endl;
            cout<<"inner DCA, r-chi2: "<<small_info_vec_axis[0][2]<<" fit error: "<<small_info_vec_axis[0][3]<<endl;
            cout<<"inner delta Phi, r-chi2: "<<small_info_vec_axis[0][4]<<" fit error: "<<small_info_vec_axis[0][5]<<endl;
            cout<<"outer DCA, r-chi2: "<<small_info_vec_axis[0][6]<<" fit error: "<<small_info_vec_axis[0][7]<<endl;
            cout<<"outer delta Phi, r-chi2: "<<small_info_vec_axis[0][8]<<" fit error: "<<small_info_vec_axis[0][9]<<endl;
            cout<<"~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~"<<endl;
            cout<<"inner DCA, r-chi2: "<<small_info_vec_axis[1][2]<<" fit error: "<<small_info_vec_axis[1][3]<<endl;
            cout<<"inner delta Phi, r-chi2: "<<small_info_vec_axis[1][4]<<" fit error: "<<small_info_vec_axis[1][5]<<endl;
            cout<<"outer DCA, r-chi2: "<<small_info_vec_axis[1][6]<<" fit error: "<<small_info_vec_axis[1][7]<<endl;
            cout<<"outer delta Phi, r-chi2: "<<small_info_vec_axis[1][8]<<" fit error: "<<small_info_vec_axis[1][9]<<endl;
            cout<<"~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~"<<endl;
            cout<<"inner DCA, r-chi2: "<<small_info_vec_axis[2][2]<<" fit error: "<<small_info_vec_axis[2][3]<<endl;
            cout<<"inner delta Phi, r-chi2: "<<small_info_vec_axis[2][4]<<" fit error: "<<small_info_vec_axis[2][5]<<endl;
            cout<<"outer DCA, r-chi2: "<<small_info_vec_axis[2][6]<<" fit error: "<<small_info_vec_axis[2][7]<<endl;
            cout<<"outer delta Phi, r-chi2: "<<small_info_vec_axis[2][8]<<" fit error: "<<small_info_vec_axis[2][9]<<endl;
            cout<<"~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~      ~~~~"<<endl;
            cout<<"inner DCA, r-chi2: "<<small_info_vec_axis[3][2]<<" fit error: "<<small_info_vec_axis[3][3]<<endl;
            cout<<"inner delta Phi, r-chi2: "<<small_info_vec_axis[3][4]<<" fit error: "<<small_info_vec_axis[3][5]<<endl;
            cout<<"outer DCA, r-chi2: "<<small_info_vec_axis[3][6]<<" fit error: "<<small_info_vec_axis[3][7]<<endl;
            cout<<"outer delta Phi, r-chi2: "<<small_info_vec_axis[3][8]<<" fit error: "<<small_info_vec_axis[3][9]<<endl;
            cout<<" "<<endl;
            cout<<"the trial origin: "<<origin.first<<" "<<origin.second<<" length: "<<length<<endl;
            cout<<"the true vertex: "<<TrigvtxMC[0]*10<<" "<<TrigvtxMC[1]*10<<endl;
            cout<<"the Quadrant "<<small_index<<" won the competition with the vertex: "<<vtx_vec[small_index].first<<" "<<vtx_vec[small_index].second<<endl;
            cout<<"the quadrant, true: "<<find_quadrant({origin},{TrigvtxMC[0]*10, TrigvtxMC[1]*10})<<" the axis method: "<<winner_quadrant_axis<<" the quadrant method: "<<small_index<<endl;
            cout<<"------------------------------------------------------------------------------------"<<endl;
        }

        // if ( winner_quadrant_axis != small_index) {
        //     cout<<" In N trail: "<<i<<", the two methods give different results, "<<"the Quadrant, Axis_method:  "<<winner_quadrant_axis<<" cornor_method: "<<small_index<<endl;
        //     cout<<" The true answer should be : "<<find_quadrant({origin},{TrigvtxMC[0]*10, TrigvtxMC[1]*10})<<endl;
        // }

        // note : generating the new 5 vertex for the next comparison
        // note : start to shrink the square
        if (i != N_trial - 1)
        {
            origin = {(vtx_vec[small_index].first + origin.first)/2., (vtx_vec[small_index].second + origin.second)/2.};
            // cout<<"test : "<<origin.first<<" "<<origin.second<<" length: "<<length<<endl;
            length /= 2.;
            vtx_vec = Get4vtx(origin,length); // vtx_vec.push_back(origin);
            vtx_vec_axis = Get4vtxAxis(origin,length);
        }
    }

    if (draw_plot_opt == true) {DrawTGraphErrors(grr_x, grr_y, grr_E, grr_E, out_folder_directory, {Form("Square_scan_history_%.1fmm_%iTrials", original_length, N_trial), "nth scan", "Winner index", "APL"});}

    // note : final vtx winner, the origin in that test, fit error info.,        fit parameter info.
    return {vtx_vec[small_index],origin, {small_info_vec[3],small_info_vec[5]}, {small_info_vec[10],small_info_vec[11]}};
    // return {{1,2},{1,2}};
}

void INTTXYvtxEvt::subMacroPlotWorking(bool phi_correction, double cos_fit_rangel, double cos_fit_ranger, double guas_fit_range)
{
    for (int i = 0; i < cluster_pair_vec.size(); i++)
    {
        vector<double> DCA_info_vec = calculateDistanceAndClosestPoint(
            cluster_pair_vec[i].first.x, cluster_pair_vec[i].first.y,
            cluster_pair_vec[i].second.x, cluster_pair_vec[i].second.y,
            current_vtxX, current_vtxY
        );

        double DCA_sign = calculateAngleBetweenVectors(
            cluster_pair_vec[i].second.x, cluster_pair_vec[i].second.y,
            cluster_pair_vec[i].first.x, cluster_pair_vec[i].first.y,
            current_vtxX, current_vtxY
        );

        if (phi_correction == true)
        {
            // cout<<"option selected "<<endl;
            Clus_InnerPhi_Offset = (cluster_pair_vec[i].first.y - current_vtxY < 0) ? atan2(cluster_pair_vec[i].first.y - current_vtxY, cluster_pair_vec[i].first.x - current_vtxX) * (180./TMath::Pi()) + 360 : atan2(cluster_pair_vec[i].first.y - current_vtxY, cluster_pair_vec[i].first.x - current_vtxX) * (180./TMath::Pi());
            Clus_OuterPhi_Offset = (cluster_pair_vec[i].second.y - current_vtxY < 0) ? atan2(cluster_pair_vec[i].second.y - current_vtxY, cluster_pair_vec[i].second.x - current_vtxX) * (180./TMath::Pi()) + 360 : atan2(cluster_pair_vec[i].second.y - current_vtxY, cluster_pair_vec[i].second.x - current_vtxX) * (180./TMath::Pi());
        }
        else // note : phi_correction == false 
        {
            Clus_InnerPhi_Offset = (cluster_pair_vec[i].first.y < 0) ? atan2(cluster_pair_vec[i].first.y, cluster_pair_vec[i].first.x) * (180./TMath::Pi()) + 360 : atan2(cluster_pair_vec[i].first.y, cluster_pair_vec[i].first.x) * (180./TMath::Pi()); 
            Clus_OuterPhi_Offset = (cluster_pair_vec[i].second.y < 0) ? atan2(cluster_pair_vec[i].second.y, cluster_pair_vec[i].second.x) * (180./TMath::Pi()) + 360 : atan2(cluster_pair_vec[i].second.y, cluster_pair_vec[i].second.x) * (180./TMath::Pi()); 
        }

        // double phi_test = (cluster_pair_vec[i].first.y < 0) ? atan2(cluster_pair_vec[i].first.y, cluster_pair_vec[i].first.x) * (180./TMath::Pi()) + 360 : atan2(cluster_pair_vec[i].first.y, cluster_pair_vec[i].first.x) * (180./TMath::Pi());
        // double phi_test_offset = (cluster_pair_vec[i].first.y - current_vtxY < 0) ? atan2(cluster_pair_vec[i].first.y - current_vtxY, cluster_pair_vec[i].first.x - current_vtxX) * (180./TMath::Pi()) + 360 : atan2(cluster_pair_vec[i].first.y - current_vtxY, cluster_pair_vec[i].first.x - current_vtxX) * (180./TMath::Pi());
        // cout<<"angle offset test : "<<phi_test<<", "<<phi_test_offset<<endl;
    
        angle_correlation -> Fill(Clus_InnerPhi_Offset, Clus_OuterPhi_Offset);
        angle_diff_DCA_dist -> Fill(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset, DCA_sign);
        angle_diff -> Fill(abs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset));
        angle_diff_inner_phi -> Fill(Clus_InnerPhi_Offset, Clus_InnerPhi_Offset - Clus_OuterPhi_Offset);
        angle_diff_outer_phi -> Fill(Clus_OuterPhi_Offset, Clus_InnerPhi_Offset - Clus_OuterPhi_Offset);
        DCA_point -> Fill(DCA_info_vec[1], DCA_info_vec[2]);
        DCA_distance_inner_phi -> Fill(Clus_InnerPhi_Offset, DCA_sign);
        DCA_distance_outer_phi -> Fill(Clus_OuterPhi_Offset, DCA_sign);
        DCA_distance -> Fill(DCA_sign);
        DCA_distance_inner_X -> Fill(cluster_pair_vec[i].first.x, DCA_sign);
        DCA_distance_inner_Y -> Fill(cluster_pair_vec[i].first.y, DCA_sign);
        DCA_distance_outer_X -> Fill(cluster_pair_vec[i].second.x, DCA_sign);
        DCA_distance_outer_Y -> Fill(cluster_pair_vec[i].second.y, DCA_sign);

        angle_diff_new -> Fill(abs(Clus_InnerPhi_Offset - Clus_OuterPhi_Offset));
    } // note : end of the loop for the cluster pair

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    DCA_distance_inner_phi_peak = (TH2F*)DCA_distance_inner_phi -> Clone();
    TH2F_threshold(DCA_distance_inner_phi_peak, 0.5); // todo : the background cut can be modified, the ratio 0.5
    DCA_distance_inner_phi_peak_profile =  DCA_distance_inner_phi_peak->ProfileX("DCA_distance_inner_phi_peak_profile");
    // TGraph * DCA_distance_inner_phi_peak_profile_graph = new TGraph();
    double point_index = 0;
    vector<double> hist_column_content = SumTH2FColumnContent(DCA_distance_inner_phi_peak);
    for (int i = 0; i < DCA_distance_inner_phi_peak_profile->GetNbinsX(); i++){
        // if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background

        DCA_distance_inner_phi_peak_profile_graph -> SetPoint(point_index, DCA_distance_inner_phi_peak_profile->GetBinCenter(i+1), DCA_distance_inner_phi_peak_profile->GetBinContent(i+1));
        // cout<<"("<<DCA_distance_inner_phi_peak_profile->GetBinCenter(i+1)<<", "<< DCA_distance_inner_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
        point_index += 1;
    }

    // cos_fit -> SetParameters(4, 1.49143e-02,  185, 0.3); // todo : we may have to apply more constaints on the fitting
    // cos_fit -> SetParLimits(0,0,1000); // note : the amplitude has to be positive
    // cos_fit -> SetParLimits(2,0,360); // note : the peak location has to be positive

    // note : here is the test with a gaus fitting to find the peak
    // gaus_fit -> SetParameters(-4.5, cos_fit->GetParameter(2), 50, 0);
    // gaus_fit -> SetParLimits(0,-100,0); // note : the gaus distribution points down
    // DCA_distance_inner_phi_peak_profile -> Fit(gaus_fit, "N","",100, 260);
    // cout<<"test, gaus fit range : "<<gaus_fit->GetParameter(1) - 25<<" "<<gaus_fit->GetParameter(1) + 25<<endl;
       
    horizontal_fit_inner -> SetParameter(0,0);

    // todo : the fit range of the gaussian fit can be modified here 
    DCA_distance_inner_phi_peak_profile_graph -> Fit(horizontal_fit_inner,"NQ","",0,360);
    DCA_distance_inner_phi_peak_profile_graph -> Fit(gaus_fit, "NQ","",cos_fit->GetParameter(2) - guas_fit_range, cos_fit->GetParameter(2) + guas_fit_range); // note : what we want and need is the peak position, so we fit the peak again 
    DCA_distance_inner_phi_peak_profile_graph -> Fit(cos_fit, "NQ","",cos_fit_rangel, cos_fit_ranger);

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    angle_diff_new_bkg_remove = (TH1F*)angle_diff_new -> Clone();
    angle_diff_new_bkg_remove -> SetLineColor(2);
    Hist_1D_bkg_remove(angle_diff_new_bkg_remove, 1.5);

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    angle_diff_inner_phi_peak = (TH2F*)angle_diff_inner_phi -> Clone();
    // TH2F_threshold_advanced(angle_diff_inner_phi_peak, 0.5);
    TH2F_threshold_advanced_2(angle_diff_inner_phi_peak, 0.5); // todo : threshold ratio can be modified here
    hist_column_content = SumTH2FColumnContent(angle_diff_inner_phi_peak);
    angle_diff_inner_phi_peak_profile =  angle_diff_inner_phi_peak->ProfileX("angle_diff_inner_phi_peak_profile");
    // TGraph * angle_diff_inner_phi_peak_profile_graph = new TGraph();
    point_index = 0;
    for (int i = 0; i < angle_diff_inner_phi_peak_profile->GetNbinsX(); i++){
        // if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background
        
        angle_diff_inner_phi_peak_profile_graph -> SetPoint(point_index, angle_diff_inner_phi_peak_profile->GetBinCenter(i+1), angle_diff_inner_phi_peak_profile->GetBinContent(i+1));
        // cout<<"("<<angle_diff_inner_phi_peak_profile->GetBinCenter(i+1)<<", "<< angle_diff_inner_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
        point_index += 1;
    }

    horizontal_fit_angle_diff_inner -> SetParameter(0,0);
    angle_diff_inner_phi_peak_profile_graph -> Fit(horizontal_fit_angle_diff_inner,"NQ","",0,360);

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    DCA_distance_outer_phi_peak = (TH2F*)DCA_distance_outer_phi -> Clone();
    TH2F_threshold(DCA_distance_outer_phi_peak, 0.5); // todo : the background cut can be modified, the ratio 0.5
    DCA_distance_outer_phi_peak_profile =  DCA_distance_outer_phi_peak->ProfileX("DCA_distance_outer_phi_peak_profile");
    // TGraph * DCA_distance_outer_phi_peak_profile_graph = new TGraph();
    point_index = 0;
    hist_column_content = SumTH2FColumnContent(DCA_distance_outer_phi_peak);
    for (int i = 0; i < DCA_distance_outer_phi_peak_profile->GetNbinsX(); i++){
        // if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background

        DCA_distance_outer_phi_peak_profile_graph -> SetPoint(point_index, DCA_distance_outer_phi_peak_profile->GetBinCenter(i+1), DCA_distance_outer_phi_peak_profile->GetBinContent(i+1));
        // cout<<"("<<DCA_distance_outer_phi_peak_profile->GetBinCenter(i+1)<<", "<< DCA_distance_outer_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
        point_index += 1;
    }
       
    horizontal_fit_outer -> SetParameter(0,0);
    // todo : the fit range of the gaussian fit can be modified here 
    DCA_distance_outer_phi_peak_profile_graph -> Fit(horizontal_fit_outer,"NQ","",0,360);

    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------
    angle_diff_outer_phi_peak = (TH2F*)angle_diff_outer_phi -> Clone();
    // TH2F_threshold_advanced(angle_diff_outer_phi_peak, 0.5);
    TH2F_threshold_advanced_2(angle_diff_outer_phi_peak, 0.5); // todo : threshold ratio can be modified here
    hist_column_content = SumTH2FColumnContent(angle_diff_outer_phi_peak);
    angle_diff_outer_phi_peak_profile =  angle_diff_outer_phi_peak->ProfileX("angle_diff_outer_phi_peak_profile");
    // TGraph * angle_diff_outer_phi_peak_profile_graph = new TGraph();
    point_index = 0;
    for (int i = 0; i < angle_diff_outer_phi_peak_profile->GetNbinsX(); i++){
        // if (hist_column_content[i] < 5){continue;} // note : in order to remove some remaining background
        
        angle_diff_outer_phi_peak_profile_graph -> SetPoint(point_index, angle_diff_outer_phi_peak_profile->GetBinCenter(i+1), angle_diff_outer_phi_peak_profile->GetBinContent(i+1));
        // cout<<"("<<angle_diff_outer_phi_peak_profile->GetBinCenter(i+1)<<", "<< angle_diff_outer_phi_peak_profile->GetBinContent(i+1)<<")"<<endl;
        point_index += 1;
    }

    horizontal_fit_angle_diff_outer -> SetParameter(0,0);
    angle_diff_outer_phi_peak_profile_graph -> Fit(horizontal_fit_angle_diff_outer,"NQ","",0,360);
    // note : ----------- ----------- ----------- ---------- ----------- ----------- ---------- ----------- ----------- ----------- -----------



    angle_diff_inner_phi_peak_profile_graph -> Set(0);
    angle_diff_outer_phi_peak_profile_graph -> Set(0);
    DCA_distance_inner_phi_peak_profile_graph -> Set(0);
    DCA_distance_outer_phi_peak_profile_graph -> Set(0);

    if (print_message_opt == true) {cout<<"circle radius : "<<abs(gaus_fit->GetParameter(0) + gaus_fit->GetParameter(3))<<" possible correction angle : "<<gaus_fit->GetParameter(1)<<endl;}
    // return {abs(gaus_fit->GetParameter(0) + gaus_fit->GetParameter(3)), gaus_fit->GetParameter(1)}; // note : gaussian based
    
    // return {
    //     abs(cos_fit->GetParameter(0) - cos_fit->GetParameter(3)), cos_fit->GetParameter(2), // note : cosine wave based
    //     angle_diff_new_bkg_remove -> GetStdDev(), angle_diff_new_bkg_remove -> GetStdDevError() // note : the std dev of the angle diff 1D distribution 
    // }; 
}

pair<double,double> INTTXYvtxEvt::Get_possible_zvtx(double rvtx, vector<double> p0, vector<double> p1) // note : inner p0, outer p1, vector {r,z}, -> {y,x}
{
    vector<double> p0_z_edge = { ( fabs( p0[1] ) < 130 ) ? p0[1] - 8. : p0[1] - 10., ( fabs( p0[1] ) < 130 ) ? p0[1] + 8. : p0[1] + 10.}; // note : vector {left edge, right edge}
    vector<double> p1_z_edge = { ( fabs( p1[1] ) < 130 ) ? p1[1] - 8. : p1[1] - 10., ( fabs( p1[1] ) < 130 ) ? p1[1] + 8. : p1[1] + 10.}; // note : vector {left edge, right edge}

    double edge_first  = Get_extrapolation(rvtx,p0_z_edge[0],p0[0],p1_z_edge[1],p1[0]);
    double edge_second = Get_extrapolation(rvtx,p0_z_edge[1],p0[0],p1_z_edge[0],p1[0]);

    double mid_point = (edge_first + edge_second) / 2.;
    double possible_width = fabs(edge_first - edge_second) / 2.;

    return {mid_point, possible_width}; // note : first : mid point, second : width

}

double INTTXYvtxEvt::Get_extrapolation(double given_y, double p0x, double p0y, double p1x, double p1y) // note : x : z, y : r
{
    if ( fabs(p0x - p1x) < 0.00001 ){ // note : the line is vertical (if z is along the x axis)
        return p0x;
    }
    else {
        double slope = (p1y - p0y) / (p1x - p0x);
        double yIntercept = p0y - slope * p0x;
        double xCoordinate = (given_y - yIntercept) / slope;
        return xCoordinate;
    }
}

void INTTXYvtxEvt::EndRun()
{
    inner_pos_xy -> Reset("ICESM");
    outer_pos_xy -> Reset("ICESM");
    inner_outer_pos_xy -> Reset("ICESM");
    N_cluster_correlation -> Reset("ICESM");
    N_cluster_correlation_close -> Reset("ICESM");
    if (draw_event_display == true) {c2 -> Print(Form("%s/temp_event_display.pdf)",out_folder_directory.c_str()));};

    return;
}


// note : this function shouldn't be used for the official z vertex calculation, it should only be used for the event display
double INTTXYvtxEvt::get_z_vertex(clu_info inner_clu, clu_info outer_clu, double target_x, double target_y)
{
    // note : x = z, 
    // note : y = radius
    double inner_clu_r = sqrt(pow(inner_clu.x,2)+ pow(inner_clu.y,2));
    double outer_clu_r = sqrt(pow(outer_clu.x,2)+ pow(outer_clu.y,2));
    double target_r    = sqrt(pow(target_x,2)   + pow(target_y,2));

    // Use the slope equation (y = ax + b) to calculate the x-coordinate for the target y
    if ( fabs(outer_clu.z - inner_clu.z) < 0.00001 ){
        return outer_clu.z;
    }
    else {
        double slope = (outer_clu_r - inner_clu_r) / (outer_clu.z - inner_clu.z);
        double yIntercept = inner_clu_r - slope * inner_clu.z;
        double xCoordinate = (target_r - yIntercept) / slope;
        return xCoordinate;
    }
    
}

void INTTXYvtxEvt::temp_bkg(TPad * c1, double peek, pair<double,double> beam_origin, InttConversion * ch_pos_DB)
{
    c1 -> cd();

    int N_ladder[4] = {12, 12, 16, 16};
    string ladder_index_string[16] = {"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15"};

    vector<double> x_vec; x_vec.clear();
    vector<double> y_vec; y_vec.clear();

    vector<double> x_vec_2; x_vec_2.clear();
    vector<double> y_vec_2; y_vec_2.clear();

    bkg -> SetTitle("INTT event display X-Y plane");
    bkg -> SetMarkerStyle(20);
    bkg -> SetMarkerSize(0.1);
    // bkg -> SetPoint(1,beam_origin.first,beam_origin.second);
    bkg -> GetXaxis() -> SetLimits(-150,150);
    bkg -> GetYaxis() -> SetRangeUser(-150,150);
    bkg -> GetXaxis() -> SetTitle("X [mm]");
    bkg -> GetYaxis() -> SetTitle("Y [mm]");
    
    bkg -> Draw("ap");

    

    for (int server_i = 0; server_i < 4; server_i++)
    {
        for (int FC_i = 0; FC_i < 14; FC_i++)
        {
            ladder_line -> DrawLine(
                ch_pos_DB -> Get_XY_all(Form("intt%i",server_i),FC_i,14,0).x, ch_pos_DB -> Get_XY_all(Form("intt%i",server_i),FC_i,14,0).y,
                ch_pos_DB -> Get_XY_all(Form("intt%i",server_i),FC_i,1,0).x, ch_pos_DB -> Get_XY_all(Form("intt%i",server_i),FC_i,1,0).y
            );
        }
    }
    
    ladder_line -> Draw("l same");
}


bool INTTXYvtxEvt::isPointInsideSquare(const std::pair<double, double> point, const std::pair<double, double> square_center, double length) {
    // Calculate half-length of the square
    double halfLength = length / 2.0;

    // Calculate boundaries of the square
    double x_min = square_center.first - halfLength;
    double x_max = square_center.first + halfLength;
    double y_min = square_center.second - halfLength;
    double y_max = square_center.second + halfLength;

    // Check if the point lies inside the square
    return (point.first > x_min && point.first < x_max && point.second > y_min && point.second < y_max);
}

void INTTXYvtxEvt::TH2LineFill(TH2F * hist_in, std::pair<double,double> point_1, std::pair<double,double> point_2)
{
    double cell_width_x = hist_in -> GetXaxis() -> GetBinWidth(1);
    double cell_width_y = hist_in -> GetYaxis() -> GetBinWidth(1);
    if (cell_width_x != cell_width_y) {cout<<"the size of the cell is not square!"<<endl; exit(1);}

    for (int xi = 1; xi < hist_in -> GetNbinsX()+1; xi++)
    {
        for (int yi = 1; yi < hist_in -> GetNbinsY()+1; yi++)
        {
            vector<double> DCA_info = calculateDistanceAndClosestPoint(
                point_1.first, point_1.second, 
                point_2.first, point_2.second, 
                hist_in -> GetXaxis() -> GetBinCenter(xi), hist_in -> GetYaxis() -> GetBinCenter(yi)
            ); 

            // cout<<" "<<endl;
            // cout<<xi<<" "<<yi<<" bin center : "<<hist_in -> GetXaxis() -> GetBinCenter(xi)<<" "<< hist_in -> GetYaxis() -> GetBinCenter(yi)<<" DCA: "<<DCA_info[0]<<" DCA point : "<<DCA_info[1]<<" "<<DCA_info[2]<<endl;

            if (isPointInsideSquare({DCA_info[1],DCA_info[2]},{hist_in -> GetXaxis() -> GetBinCenter(xi), hist_in -> GetYaxis() -> GetBinCenter(yi)},cell_width_x))
            {
                // cout<<"------- before check bin content: "<<xi<<" "<<yi<<" "<<hist_in -> GetBinContent(xi,xi)<<endl;
                // cout<<"------- So the cell is filled "<<xi<<" "<<yi<<endl;
                hist_in -> SetBinContent(xi,yi,hist_in -> GetBinContent(xi,yi) + 1);
                // cout<<"------- after check bin content: "<<xi<<" "<<yi<<" "<<hist_in -> GetBinContent(xi,xi)<<endl;
            }
        }
    }
}




#endif