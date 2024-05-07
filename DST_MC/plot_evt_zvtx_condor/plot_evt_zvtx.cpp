#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>

#include "sPhenixStyle.h"
#include "../ana_map_folder/ana_map_v1.h"

double gaus_func(double *x, double *par)
{
    // note : par[0] : size
    // note : par[1] : mean
    // note : par[2] : width
    // note : par[3] : offset 
    return par[0] * TMath::Gaus(x[0],par[1],par[2]) + par[3];
}

int main()
{
    string input_directory = "/sphenix/user/ChengWei/sPH_dNdeta/self_gen_singleMu_v2/complete_file/evt_vtxZ/complete_file";
    string input_file_name = "merged_file.root";
    string out_folder_directory = input_directory + "/merged_result";
    double required_zvtx_diff = 5.; // note : unit : mm
    int zvtx_dist_NClus_cut = 1000;
    string run_type = "MC";
    string plot_text = "simulation";
    
    int eID_in;
    int in_nclu_inner;
    int in_nclu_outer;
    bool good_zvtx_tag;
    int Centrality_bin;
    double reco_zvtx;
    double MC_true_zvtx;

    TFile * file_in = TFile::Open(Form("%s/%s", input_directory.c_str(), input_file_name.c_str()));
    TTree * tree = (TTree*)file_in -> Get("tree_Z");
    tree -> SetBranchStatus("*", 0);
    tree -> SetBranchStatus("nclu_inner",1);
    tree -> SetBranchStatus("nclu_outer",1);
    tree -> SetBranchStatus("LB_Gaus_Mean_mean", 1);
    tree -> SetBranchStatus("good_zvtx_tag", 1);
    tree -> SetBranchStatus("Centrality_bin",1);
    tree -> SetBranchStatus("MC_true_zvtx",1);
    tree -> SetBranchStatus("eID",1);

    tree -> SetBranchAddress("nclu_inner", &in_nclu_inner);
    tree -> SetBranchAddress("nclu_outer", &in_nclu_outer);
    tree -> SetBranchAddress("LB_Gaus_Mean_mean", &reco_zvtx);
    tree -> SetBranchAddress("good_zvtx_tag", &good_zvtx_tag);
    tree -> SetBranchAddress("Centrality_bin", &Centrality_bin);
    tree -> SetBranchAddress("MC_true_zvtx", &MC_true_zvtx);
    tree -> SetBranchAddress("eID", &eID_in);

    cout<<"tree entries : "<<tree->GetEntries()<<endl;

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    ltx->SetTextAlign(31);

    TLatex * draw_text = new TLatex();
    draw_text -> SetNDC();
    draw_text -> SetTextSize(0.025);

    TLegend * leg = new TLegend(0.65,0.75,0.8,0.9);
    leg -> SetTextSize(0.035);

    // todo : change the namespace upon necessity
    map<int,int> centrality_map = ANA_MAP_V2::centrality_map;
    vector<string> centrality_region = ANA_MAP_V2::centrality_region;
    vector<double> z_region = ANA_MAP_V2::z_region;

    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("c1","c1",950,800);
    c1 -> cd();

    // // note : required zvtx tag
    TH2F * true_zvtx_Mbin_Truezpos = new TH2F("","true_zvtx_Mbin_Truezpos;MBin;True zvtx [mm]", centrality_region.size(), 0, centrality_region.size(), z_region.size() - 1, &z_region[0]); 
    TH2F * reco_zvtx_Mbin_Truezpos = new TH2F("","reco_zvtx_Mbin_Truezpos;MBin;True zvtx [mm]", centrality_region.size(), 0, centrality_region.size(), z_region.size() - 1, &z_region[0]); 
    TH2F * ratio_Mbin_Truezpos = new TH2F("","ratio_Mbin_Truezpos;MBin;True zvtx [mm]",         centrality_region.size(), 0, centrality_region.size(), z_region.size() - 1, &z_region[0]); 
    
    TH2F * Z_resolution_Nclu = new TH2F("","Z_resolution_Nclu;# of clusters;#DeltaZ (Reco - True) [mm]",200,0,6000,200,-100,100);
    Z_resolution_Nclu -> GetXaxis() -> SetNdivisions(505);

    TH2F * Z_resolution_pos = new TH2F("","Z_resolution_pos;True Zvtx [mm];#DeltaZ (Reco - True) [mm]",200,-350,350,200,-50,50);
    Z_resolution_pos -> GetXaxis() -> SetNdivisions(505);

    TH2F * Z_resolution_pos_cut = new TH2F("","Z_resolution_pos_cut;True Zvtx [mm];#DeltaZ (Reco - True) [mm]",200,-350,350,200,-50,50);
    Z_resolution_pos_cut -> GetXaxis() -> SetNdivisions(505);

    // TH1F * Z_resolution = new TH1F("","Z_resolution;#DeltaZ (Reco - True) [mm];Entry",200,-30,30);
    // Z_resolution -> GetXaxis() -> SetNdivisions(505);

    TH2F * Z_correlation = new TH2F("","Z_correlation;True Zvtx [mm];Reco Zvtx [mm]",200,-500,100,200,-500,100);
    Z_correlation -> GetXaxis() -> SetNdivisions(505);

    TH1F * reco_Zvtx_dist = new TH1F("","reco_Zvtx_dist;Reco Zvtx [mm];Entry",200,-500,100);
    reco_Zvtx_dist -> GetXaxis() -> SetNdivisions(505);

    TH1F * true_Zvtx_dist = new TH1F("","true_Zvtx_dist;True Zvtx [mm];Entry",200,-500,100);
    true_Zvtx_dist -> GetXaxis() -> SetNdivisions(505);
    true_Zvtx_dist -> SetLineColor(2);


    vector<TH1F *> Z_resolution; Z_resolution.clear();
    for (int i = 0; i < centrality_region.size(); i++)
    {
        Z_resolution.push_back(new TH1F("",Form("Z_resolution Mbin %i;#DeltaZ (Reco - True) [mm];Entry", i),200,-25,25)); // note : unit: mm
        Z_resolution[i] -> GetXaxis() -> SetNdivisions(505);
    }

    TF1 * gaus_fit_2 = new TF1("gaus_fit_2", gaus_func, Z_resolution[0]->GetXaxis()->GetXmin(), Z_resolution[0]->GetXaxis()->GetXmax(),4);
    gaus_fit_2 -> SetLineColor(2);
    gaus_fit_2 -> SetLineWidth(2);
    gaus_fit_2 -> SetNpx(1000);

    TF1 * linear_fit = new TF1("linear_fit","pol1", -500, 500);
    linear_fit -> SetLineColor(2);
    linear_fit -> SetLineWidth(2);

    TGraphErrors * Z_resolution_centrality_gr_fit = new TGraphErrors();
    Z_resolution_centrality_gr_fit -> SetMarkerStyle(20);
    Z_resolution_centrality_gr_fit -> SetMarkerSize(1.5);
    Z_resolution_centrality_gr_fit -> GetXaxis()->SetTitle("Centrality bin"); 
    Z_resolution_centrality_gr_fit -> GetYaxis()->SetTitle("#DeltaZ width");

    TGraphErrors * Z_mean_centrality_gr_fit = new TGraphErrors();
    Z_mean_centrality_gr_fit -> SetMarkerStyle(20);
    Z_mean_centrality_gr_fit -> SetMarkerSize(1.5);
    Z_mean_centrality_gr_fit -> GetXaxis()->SetTitle("Centrality bin"); 
    Z_mean_centrality_gr_fit -> GetYaxis()->SetTitle("Dist. Mean");

    TGraphErrors * Z_resolution_centrality_gr_num = new TGraphErrors();
    Z_resolution_centrality_gr_num -> SetMarkerStyle(20);
    Z_resolution_centrality_gr_num -> SetMarkerSize(1.5);
    Z_resolution_centrality_gr_num -> GetXaxis()->SetTitle("Centrality bin"); 
    Z_resolution_centrality_gr_num -> GetYaxis()->SetTitle("#DeltaZ width");

    TGraphErrors * Z_mean_centrality_gr_num = new TGraphErrors();
    Z_mean_centrality_gr_num -> SetMarkerStyle(20);
    Z_mean_centrality_gr_num -> SetMarkerSize(1.5);
    Z_mean_centrality_gr_num -> GetXaxis()->SetTitle("Centrality bin"); 
    Z_mean_centrality_gr_num -> GetYaxis()->SetTitle("Dist. Mean");

    for (int event_i = 0; event_i < tree -> GetEntries(); event_i++)
    {
        tree -> GetEntry(event_i);

        if (event_i % 4000 == 0) {cout<<"processing, event_i : "<<event_i<<endl;}

        if (good_zvtx_tag == 1 && in_nclu_inner > 0 && in_nclu_outer > 0)
        {
            if ((in_nclu_inner + in_nclu_outer) > zvtx_dist_NClus_cut){
                reco_Zvtx_dist -> Fill(reco_zvtx);
                true_Zvtx_dist -> Fill(MC_true_zvtx);   
            }

            
            if (run_type != "MC") {continue;}

            Z_resolution_Nclu -> Fill((in_nclu_inner + in_nclu_outer), reco_zvtx - MC_true_zvtx);
            Z_resolution_pos  -> Fill(MC_true_zvtx, reco_zvtx - MC_true_zvtx);
            Z_correlation     -> Fill(MC_true_zvtx, reco_zvtx);
            
            if ((in_nclu_inner + in_nclu_outer) > zvtx_dist_NClus_cut) {Z_resolution_pos_cut -> Fill(MC_true_zvtx, reco_zvtx - MC_true_zvtx);}

            Z_resolution[ centrality_map[Centrality_bin] ] -> Fill( reco_zvtx - MC_true_zvtx );
            Z_resolution[ Z_resolution.size() - 1 ]        -> Fill( reco_zvtx - MC_true_zvtx ); // note : the inclusive


            // note : for the passing rate
            true_zvtx_Mbin_Truezpos -> Fill(centrality_map[Centrality_bin], MC_true_zvtx);
            true_zvtx_Mbin_Truezpos -> Fill(centrality_map.size() - 1, MC_true_zvtx);
            if (fabs(reco_zvtx - MC_true_zvtx) < required_zvtx_diff)
            {
                reco_zvtx_Mbin_Truezpos -> Fill(centrality_map[Centrality_bin], MC_true_zvtx);
                reco_zvtx_Mbin_Truezpos -> Fill(centrality_map.size() - 1, MC_true_zvtx);
            }
        }
    }

    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    reco_Zvtx_dist -> Draw("hist"); 
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
    draw_text -> DrawLatex(0.21, 0.87, Form("NClus > %i", zvtx_dist_NClus_cut));
    draw_text -> DrawLatex(0.21, 0.83, Form("Reco. Entries: %.0f", reco_Zvtx_dist -> GetEntries()));
    draw_text -> DrawLatex(0.21, 0.79, Form("Reco. Mean: %.2f mm", reco_Zvtx_dist -> GetMean()));
    draw_text -> DrawLatex(0.21, 0.75, Form("Reco. Width: %.2f mm", reco_Zvtx_dist -> GetStdDev()));

    if (run_type == "MC")
    {
        leg -> AddEntry(reco_Zvtx_dist, "Reco. Z", "f");
        leg -> AddEntry(true_Zvtx_dist, "True Z", "f");
        true_Zvtx_dist -> Draw("hist same");

        draw_text -> DrawLatex(0.21, 0.71, Form("True Entries: %.0f", true_Zvtx_dist -> GetEntries()));
        draw_text -> DrawLatex(0.21, 0.67, Form("True Mean: %.2f mm", true_Zvtx_dist -> GetMean()));
        draw_text -> DrawLatex(0.21, 0.63, Form("True Width: %.2f mm", true_Zvtx_dist -> GetStdDev()));

        leg -> Draw("same");
    }

    c1 -> Print(Form("%s/reco_Zvtx_dist.pdf",out_folder_directory.c_str()));
    c1 -> Clear();

    if (run_type != "MC") {file_in -> Close(); return 0;}


    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    if (run_type == "MC")
    {
        vector<double> fit_width_vec; fit_width_vec.clear();
        vector<double> num_width_vec; num_width_vec.clear();
        vector<double> avg_width_vec; avg_width_vec.clear();

        c1 -> Print(Form("%s/INTT_Z_resolution_centrality.pdf(", out_folder_directory.c_str()));
        for (int i = 0; i < Z_resolution.size(); i++)
        {
            Z_resolution[i] -> SetMinimum(0);

            gaus_fit_2 -> SetParameters(Z_resolution[i] -> GetBinContent( Z_resolution[i] -> GetMaximumBin() ), Z_resolution[i] -> GetBinCenter( Z_resolution[i] -> GetMaximumBin() ), Z_resolution[i] -> GetStdDev(), 0);
            gaus_fit_2 -> SetParLimits(0,0,100000);  // note : size 
            gaus_fit_2 -> SetParLimits(2,0,10000);   // note : Width
            gaus_fit_2 -> SetParLimits(3,0,10000);   // note : offset
            Z_resolution[i] -> Fit(gaus_fit_2, "N", "", Z_resolution[i] -> GetBinCenter( Z_resolution[i] -> GetMaximumBin() ) - (2 * Z_resolution[i] -> GetStdDev() ), Z_resolution[i] -> GetBinCenter( Z_resolution[i] -> GetMaximumBin() ) + (2 * Z_resolution[i] -> GetStdDev() ) );
            cout<<"fit reduced chi2 : "<<i<<" "<<gaus_fit_2 -> GetChisquare() / double(gaus_fit_2 -> GetNDF())<<endl;
            Z_resolution[i] -> Draw("hist"); 
            gaus_fit_2 -> SetRange( gaus_fit_2->GetParameter(1) - gaus_fit_2->GetParameter(2) * 2.5, gaus_fit_2->GetParameter(1) + gaus_fit_2->GetParameter(2) * 2.5 ); 
            gaus_fit_2 -> Draw("lsame");
            draw_text -> DrawLatex(0.21, 0.71, Form("Gaus mean  : %.2f mm",gaus_fit_2 -> GetParameter(1)));
            draw_text -> DrawLatex(0.21, 0.67, Form("Gaus width : %.2f mm",fabs(gaus_fit_2 -> GetParameter(2))));
            ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));

            // Z_resolution_centrality_gr_fit -> SetPoint     (i, i, gaus_fit_2->GetParameter(2));
            // Z_resolution_centrality_gr_fit -> SetPointError(i, 0, gaus_fit_2->GetParError(2));

            // Z_mean_centrality_gr_fit       -> SetPoint     (i, i, gaus_fit_2->GetParameter(1));            
            // Z_mean_centrality_gr_fit       -> SetPointError(i, 0, gaus_fit_2->GetParError(1));

            Z_resolution_centrality_gr_fit -> SetPoint     (i, i, gaus_fit_2->GetParameter(2));
            Z_resolution_centrality_gr_fit -> SetPointError(i, 0, 0);

            Z_mean_centrality_gr_fit       -> SetPoint     (i, i, gaus_fit_2->GetParameter(1));            
            Z_mean_centrality_gr_fit       -> SetPointError(i, 0, 0);


            Z_resolution_centrality_gr_num -> SetPoint     (i, i, Z_resolution[i]->GetStdDev());
            Z_resolution_centrality_gr_num -> SetPointError(i, 0, Z_resolution[i]->GetStdDevError());

            Z_mean_centrality_gr_num       -> SetPoint     (i, i, Z_resolution[i]->GetMean());            
            Z_mean_centrality_gr_num       -> SetPointError(i, 0, Z_resolution[i]->GetMeanError());

            fit_width_vec.push_back(gaus_fit_2->GetParameter(2));
            num_width_vec.push_back(Z_resolution[i]->GetStdDev());
            avg_width_vec.push_back((gaus_fit_2->GetParameter(2) + Z_resolution[i]->GetStdDev())/2.);

            c1 -> Print(Form("%s/INTT_Z_resolution_centrality.pdf", out_folder_directory.c_str()));
            c1 -> Clear();
        }
        c1 -> Print(Form("%s/INTT_Z_resolution_centrality.pdf)", out_folder_directory.c_str()));
        c1 -> Clear();

        for (int i = 0; i < fit_width_vec.size(); i++)
        {
            cout<<i<<" fit: "<<fit_width_vec[i]<<", num: "<<num_width_vec[i]<<", avg: "<<avg_width_vec[i]<<endl;
        }
    }

    if (run_type == "MC")
    {
        c1 -> cd();
        Z_resolution_centrality_gr_fit -> Draw("ap");
        Z_resolution_centrality_gr_fit -> GetXaxis() -> SetLimits(-1,Z_resolution_centrality_gr_fit -> GetN() + 1);
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        c1 -> Print(Form("%s/Z_resolution_centrality_gr_fit.pdf", out_folder_directory.c_str()));
        c1 -> Clear();

        for (int i = 0; i < Z_resolution_centrality_gr_fit->GetN(); i++)
        {
            cout<<i<<" "<<Z_resolution_centrality_gr_fit->GetPointY(i)<<endl;
        }

        c1 -> cd();
        Z_mean_centrality_gr_fit -> Draw("ap");
        Z_mean_centrality_gr_fit -> GetXaxis() -> SetLimits(-1,Z_mean_centrality_gr_fit -> GetN() + 1);
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        c1 -> Print(Form("%s/Z_mean_centrality_gr_fit.pdf", out_folder_directory.c_str()));
        c1 -> Clear();


        c1 -> cd();
        Z_resolution_centrality_gr_num -> Draw("ap");
        Z_resolution_centrality_gr_num -> GetXaxis() -> SetLimits(-1,Z_resolution_centrality_gr_num -> GetN() + 1);
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        c1 -> Print(Form("%s/Z_resolution_centrality_gr_num.pdf", out_folder_directory.c_str()));
        c1 -> Clear();

        for (int i = 0; i < Z_resolution_centrality_gr_num->GetN(); i++)
        {
            cout<<i<<" "<<Z_resolution_centrality_gr_num->GetPointY(i)<<endl;
        }

        c1 -> cd();
        Z_mean_centrality_gr_num -> Draw("ap");
        Z_mean_centrality_gr_num -> GetXaxis() -> SetLimits(-1,Z_mean_centrality_gr_num -> GetN() + 1);
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        c1 -> Print(Form("%s/Z_mean_centrality_gr_num.pdf", out_folder_directory.c_str()));
        c1 -> Clear();

    }


    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (run_type == "MC")
    {    
        Z_resolution_Nclu -> Draw("colz0"); 
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        c1 -> Print(Form("%s/Z_resolution_Nclu.pdf",out_folder_directory.c_str()));
        c1 -> Clear();
    }
    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (run_type == "MC")
    {
        Z_resolution_pos -> Draw("colz0"); 
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        c1 -> Print(Form("%s/Z_resolution_pos.pdf",out_folder_directory.c_str()));
        c1 -> Clear();
    }
    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (run_type == "MC")
    {
        Z_resolution_pos_cut -> Draw("colz0"); 
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        draw_text -> DrawLatex(0.21, 0.87, Form("NClus > %i", zvtx_dist_NClus_cut));
        c1 -> Print(Form("%s/Z_resolution_pos_cut.pdf",out_folder_directory.c_str()));
        c1 -> Clear();
    }

    if (run_type == "MC")
    {
        Z_correlation -> Draw("colz0"); 
        ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
        Z_correlation -> Fit(linear_fit, "NQ");
        linear_fit -> Draw("l same");
        draw_text -> DrawLatex(0.21, 0.87, Form("Y = %.2fx + %.2f", linear_fit -> GetParameter(1), linear_fit -> GetParameter(0)));
        c1 -> Print(Form("%s/Z_correlation.pdf",out_folder_directory.c_str()));
        c1 -> Clear();
    }



    // note : ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // TH2F_division(true_zvtx_Mbin_Truezpos,           reco_zvtx_Mbin_Truezpos, ratio_Mbin_Truezpos);
    // TH2F_division(true_zvtx_Mbin_Truezpos_inclusive, reco_zvtx_Mbin_Truezpos, ratio_Mbin_Truezpos_inclusive);

    ratio_Mbin_Truezpos = (TH2F *) reco_zvtx_Mbin_Truezpos -> Clone("ratio_Mbin_Truezpos");
    ratio_Mbin_Truezpos -> Divide(ratio_Mbin_Truezpos, true_zvtx_Mbin_Truezpos, 1, 1, "B");

    gStyle->SetPaintTextFormat("1.0f");

    c1 -> cd();
    true_zvtx_Mbin_Truezpos -> SetMarkerSize(1.0);
    true_zvtx_Mbin_Truezpos -> Draw("colz0");
    true_zvtx_Mbin_Truezpos -> Draw("HIST TEXT45 SAME");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
    // draw_text->DrawLatex(0.16, 1 - gPad->GetTopMargin() + 0.01, Form("Pass if |Reco - True| < %.0f mm", required_zvtx_diff));
    c1 -> Print(Form("%s/true_zvtx_Mbin_Truezpos.pdf", out_folder_directory.c_str()));
    c1 -> Clear();

    c1 -> cd();
    reco_zvtx_Mbin_Truezpos -> SetMarkerSize(1.0);
    reco_zvtx_Mbin_Truezpos -> Draw("colz0");
    reco_zvtx_Mbin_Truezpos -> Draw("HIST TEXT45 SAME");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
    draw_text->DrawLatex(0.16, 1 - gPad->GetTopMargin() + 0.01, Form("Pass if |Reco - True| < %.0f mm", required_zvtx_diff));
    c1 -> Print(Form("%s/reco_zvtx_Mbin_Truezpos.pdf", out_folder_directory.c_str()));
    c1 -> Clear();

    gStyle->SetPaintTextFormat("1.3f");

    c1 -> cd();
    ratio_Mbin_Truezpos -> SetMarkerSize(0.7);
    ratio_Mbin_Truezpos -> Draw("colz0");
    ratio_Mbin_Truezpos -> Draw("HIST TEXT45E SAME");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", plot_text.c_str()));
    draw_text->DrawLatex(0.16, 1 - gPad->GetTopMargin() + 0.01, Form("Pass if |Reco - True| < %.0f mm", required_zvtx_diff));
    c1 -> Print(Form("%s/ratio_Mbin_Truezpos.pdf", out_folder_directory.c_str()));
    c1 -> Clear();


    true_zvtx_Mbin_Truezpos -> Reset("ICESM");
    reco_zvtx_Mbin_Truezpos -> Reset("ICESM");
    ratio_Mbin_Truezpos -> Reset("ICESM");
    file_in -> Close();

    return 0;
}