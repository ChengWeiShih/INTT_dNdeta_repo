#include <iostream>
#include <fstream>
using namespace std;

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>

#include "sPhenixStyle.h"

vector<string> read_list(string folder_direction, string MC_list_name)
{
    vector<string> file_list;
    string list_buffer;
    ifstream data_list;
	data_list.open((folder_direction + "/" + MC_list_name).c_str());

   file_list.clear();

	while (1)
	{
		data_list >> list_buffer;
		if (!data_list.good())
		{
			break;
		}
		file_list.push_back(list_buffer);
	}
	cout<<"size in the" <<MC_list_name<<": "<<file_list.size()<<endl;

    return file_list;
}

int main()
{
    string input_folder = "/sphenix/user/ChengWei/INTT/INTT_commissioning/ZeroField/20869/condor_runXY_folder";
    string file_list_name = "file_list.txt";
    string output_directory = input_folder + "/" + "merged_result";
    
    TChain *chain_in = new TChain("tree");
    vector<string> file_list = read_list(input_folder, file_list_name);
    for (string file_name : file_list) { chain_in -> Add((file_name).c_str()); }
    cout<<"chain_in -> GetEntries() : "<<chain_in -> GetEntries()<<endl;

    long in_start_evt;
    long in_end_evt;
    double in_quadrant_corner_X;
    double in_quadrant_corner_Y;
    double in_quadrant_center_X;
    double in_quadrant_center_Y;
    double in_line_filled_mean_X;
    double in_line_filled_mean_Y;
    double in_line_filled_stddev_X;
    double in_line_filled_stddev_Y;
    double in_line_filled_variance_X;
    double in_line_filled_variance_Y;
    double in_line_filled_covariance;

    chain_in -> SetBranchAddress("start_evt",&in_start_evt);
    chain_in -> SetBranchAddress("end_evt",&in_end_evt);
    chain_in -> SetBranchAddress("quadrant_corner_X",&in_quadrant_corner_X);
    chain_in -> SetBranchAddress("quadrant_corner_Y",&in_quadrant_corner_Y);
    chain_in -> SetBranchAddress("quadrant_center_X",&in_quadrant_center_X);
    chain_in -> SetBranchAddress("quadrant_center_Y",&in_quadrant_center_Y);
    chain_in -> SetBranchAddress("line_filled_mean_X",&in_line_filled_mean_X);
    chain_in -> SetBranchAddress("line_filled_mean_Y",&in_line_filled_mean_Y);
    chain_in -> SetBranchAddress("line_filled_stddev_X",&in_line_filled_stddev_X);
    chain_in -> SetBranchAddress("line_filled_stddev_Y",&in_line_filled_stddev_Y);
    chain_in -> SetBranchAddress("line_filled_variance_X",&in_line_filled_variance_X);
    chain_in -> SetBranchAddress("line_filled_variance_Y",&in_line_filled_variance_Y);
    chain_in -> SetBranchAddress("line_filled_covariance",&in_line_filled_covariance);

    double X_range_l = -1.;
    double X_range_r = 0.5;
    double Y_range_l = 1.7;
    double Y_range_r = 3.2;

    TH1F * quadrant_VtxX_1D = new TH1F("quadrant_VtxX","quadrant_VtxX;X axis [mm];Entry",100,X_range_l,X_range_r);
    TH1F * quadrant_VtxY_1D = new TH1F("quadrant_VtxY","quadrant_VtxY;Y axis [mm];Entry",100,Y_range_l,Y_range_r);
    TH1F * line_filled_mean_X_1D = new TH1F("line_filled_mean_X","line_filled_mean_X;X axis [mm];Entry",100,X_range_l,X_range_r);
    TH1F * line_filled_mean_Y_1D = new TH1F("line_filled_mean_Y","line_filled_mean_Y;Y axis [mm];Entry",100,Y_range_l,Y_range_r);
    
    TH1F * line_filled_stddev_X_1D = new TH1F("line_filled_stddev_X_1D","line_filled_stddev_X_1D;StdDev, X axis [mm];Entry",100,0,1);
    TH1F * line_filled_stddev_Y_1D = new TH1F("line_filled_stddev_Y_1D","line_filled_stddev_Y_1D;StdDev, Y axis [mm];Entry",100,0,1);
    
    TH1F * line_filled_variance_X_1D = new TH1F("line_filled_variance_X_1D","line_filled_variance_X_1D;Variance, X axis [mm];Entry",100,0,1);    
    TH1F * line_filled_variance_Y_1D = new TH1F("line_filled_variance_Y_1D","line_filled_variance_Y_1D;Variance, Y axis [mm];Entry",100,0,1);
    
    TH1F * line_filled_covariance_1D = new TH1F("line_filled_covariance_1D","line_filled_covariance_1D;Covariance [mm];Entry",100,0,0.12);


    vector<double> evt_index_vec; evt_index_vec.clear();
    vector<double> evt_index_range_vec; evt_index_range_vec.clear();
    vector<double> quadrant_VtxX_vec; quadrant_VtxX_vec.clear();
    vector<double> quadrant_VtxY_vec; quadrant_VtxY_vec.clear();
    vector<double> line_filled_mean_X_vec; line_filled_mean_X_vec.clear();
    vector<double> line_filled_mean_Y_vec; line_filled_mean_Y_vec.clear();
    vector<double> zero_vec; zero_vec.clear();

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    ltx->SetTextAlign(31);

    TLatex * draw_text = new TLatex();
    draw_text -> SetNDC();
    draw_text -> SetTextSize(0.03);

    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("c1","c1",950,800);
    c1 -> cd();

    for (int i = 0; i < chain_in -> GetEntries(); i++)
    {
        chain_in -> GetEntry(i);

        double quadrant_VtxX = (in_quadrant_corner_X + in_quadrant_center_X)/2.;
        double quadrant_VtxY = (in_quadrant_corner_Y + in_quadrant_center_Y)/2.;
        double line_filled_mean_X = in_line_filled_mean_X;
        double line_filled_mean_Y = in_line_filled_mean_Y;

        evt_index_vec.push_back((in_start_evt + in_end_evt)/2.);
        evt_index_range_vec.push_back((in_end_evt - in_start_evt)/2.);
        quadrant_VtxX_vec.push_back(quadrant_VtxX);
        quadrant_VtxY_vec.push_back(quadrant_VtxY);
        line_filled_mean_X_vec.push_back(line_filled_mean_X);
        line_filled_mean_Y_vec.push_back(line_filled_mean_Y);
        zero_vec.push_back(0);

        quadrant_VtxX_1D -> Fill(quadrant_VtxX);
        quadrant_VtxY_1D -> Fill(quadrant_VtxY);
        line_filled_mean_X_1D -> Fill(line_filled_mean_X);
        line_filled_mean_Y_1D -> Fill(line_filled_mean_Y);

        line_filled_stddev_X_1D -> Fill(in_line_filled_stddev_X);
        line_filled_stddev_Y_1D -> Fill(in_line_filled_stddev_Y);
        line_filled_variance_X_1D -> Fill(in_line_filled_variance_X);
        line_filled_variance_Y_1D -> Fill(in_line_filled_variance_Y);
        line_filled_covariance_1D -> Fill(in_line_filled_covariance);
    }

    TGraphErrors * g_quadrant_VtxX_index = new TGraphErrors(chain_in -> GetEntries(), &evt_index_vec[0], &quadrant_VtxX_vec[0], &evt_index_range_vec[0], &zero_vec[0]);
    g_quadrant_VtxX_index -> SetMarkerStyle(20);
    g_quadrant_VtxX_index -> SetMarkerSize(0.7);
    g_quadrant_VtxX_index -> SetMarkerColor(1);
    g_quadrant_VtxX_index -> GetXaxis() -> SetTitle("Event ID");
    g_quadrant_VtxX_index -> GetYaxis() -> SetTitle("Beam Position in X axis [mm]");
    g_quadrant_VtxX_index -> GetXaxis() -> SetLimits(-1000, g_quadrant_VtxX_index -> GetPointX(g_quadrant_VtxX_index -> GetN() - 1) + 10000);
    g_quadrant_VtxX_index -> GetYaxis() -> SetRangeUser(quadrant_VtxX_1D -> GetXaxis() -> GetXmin(),quadrant_VtxX_1D -> GetXaxis() -> GetXmax());

    TGraphErrors * g_quadrant_VtxY_index = new TGraphErrors(chain_in -> GetEntries(), &evt_index_vec[0], &quadrant_VtxY_vec[0], &evt_index_range_vec[0], &zero_vec[0]);
    g_quadrant_VtxY_index -> SetMarkerStyle(20);
    g_quadrant_VtxY_index -> SetMarkerSize(0.7);
    g_quadrant_VtxY_index -> SetMarkerColor(1);
    g_quadrant_VtxY_index -> GetXaxis() -> SetTitle("Event ID");
    g_quadrant_VtxY_index -> GetYaxis() -> SetTitle("Beam Position in Y axis [mm]");
    g_quadrant_VtxY_index -> GetXaxis() -> SetLimits(-1000, g_quadrant_VtxY_index -> GetPointX(g_quadrant_VtxY_index -> GetN() - 1) + 10000);
    g_quadrant_VtxY_index -> GetYaxis() -> SetRangeUser(quadrant_VtxY_1D -> GetXaxis() -> GetXmin(),quadrant_VtxY_1D -> GetXaxis() -> GetXmax());

    TGraphErrors * g_line_filled_mean_X_index = new TGraphErrors(chain_in -> GetEntries(), &evt_index_vec[0], &line_filled_mean_X_vec[0], &evt_index_range_vec[0], &zero_vec[0]);
    g_line_filled_mean_X_index -> SetMarkerStyle(20);
    g_line_filled_mean_X_index -> SetMarkerSize(0.7);
    g_line_filled_mean_X_index -> SetMarkerColor(2);
    g_line_filled_mean_X_index -> SetLineColor(2);
    g_line_filled_mean_X_index -> GetXaxis() -> SetTitle("Event ID");
    g_line_filled_mean_X_index -> GetYaxis() -> SetTitle("Beam Position in X axis [mm]");

    TGraphErrors * g_line_filled_mean_Y_index = new TGraphErrors(chain_in -> GetEntries(), &evt_index_vec[0], &line_filled_mean_Y_vec[0], &evt_index_range_vec[0], &zero_vec[0]);
    g_line_filled_mean_Y_index -> SetMarkerStyle(20);
    g_line_filled_mean_Y_index -> SetMarkerSize(0.7);
    g_line_filled_mean_Y_index -> SetMarkerColor(2);
    g_line_filled_mean_Y_index -> SetLineColor(2);
    g_line_filled_mean_Y_index -> GetXaxis() -> SetTitle("Event ID");
    g_line_filled_mean_Y_index -> GetYaxis() -> SetTitle("Beam Position in Y axis [mm]");

    TLegend * leg = new TLegend(0.25,0.75,0.35,0.9);
    leg -> AddEntry(g_quadrant_VtxX_index, "Quadrant method", "p");
    leg -> AddEntry(g_line_filled_mean_X_index, "Line-filled method", "p");


    c1 -> cd();
    g_quadrant_VtxX_index -> Draw("AP");
    g_line_filled_mean_X_index -> Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    leg -> Draw("same");
    c1 -> Print((output_directory + "/VtxX_index.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    g_quadrant_VtxY_index -> Draw("AP");
    g_line_filled_mean_Y_index -> Draw("p same");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    leg -> Draw("same");
    c1 -> Print((output_directory + "/VtxY_index.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    quadrant_VtxX_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", quadrant_VtxX_1D->GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", quadrant_VtxX_1D->GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", quadrant_VtxX_1D->GetStdDev()));
    c1 -> Print((output_directory + "/quadrant_VtxX_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    quadrant_VtxY_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", quadrant_VtxY_1D->GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", quadrant_VtxY_1D->GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", quadrant_VtxY_1D->GetStdDev()));
    c1 -> Print((output_directory + "/quadrant_VtxY_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    line_filled_mean_X_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", line_filled_mean_X_1D->GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", line_filled_mean_X_1D->GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", line_filled_mean_X_1D->GetStdDev()));
    c1 -> Print((output_directory + "/line_filled_mean_X_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    line_filled_mean_Y_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", line_filled_mean_Y_1D->GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", line_filled_mean_Y_1D->GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", line_filled_mean_Y_1D->GetStdDev()));
    c1 -> Print((output_directory + "/line_filled_mean_Y_1D.pdf").c_str());
    c1 -> Clear();
    

    c1 -> cd();
    line_filled_stddev_X_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", line_filled_stddev_X_1D->GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", line_filled_stddev_X_1D->GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", line_filled_stddev_X_1D->GetStdDev()));
    c1 -> Print((output_directory + "/line_filled_stddev_X_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    line_filled_stddev_Y_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", line_filled_stddev_Y_1D->GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", line_filled_stddev_Y_1D->GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", line_filled_stddev_Y_1D->GetStdDev()));
    c1 -> Print((output_directory + "/line_filled_stddev_Y_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    line_filled_variance_X_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.55, 0.84, Form("Entry: %.0f", line_filled_variance_X_1D->GetEntries()));
    draw_text -> DrawLatex(0.55, 0.80, Form("Mean: %.4f", line_filled_variance_X_1D->GetMean()));
    draw_text -> DrawLatex(0.55, 0.76, Form("StdDev: %.4f", line_filled_variance_X_1D->GetStdDev()));
    c1 -> Print((output_directory + "/line_filled_variance_X_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    line_filled_variance_Y_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.55, 0.84, Form("Entry: %.0f", line_filled_variance_Y_1D->GetEntries()));
    draw_text -> DrawLatex(0.55, 0.80, Form("Mean: %.4f", line_filled_variance_Y_1D->GetMean()));
    draw_text -> DrawLatex(0.55, 0.76, Form("StdDev: %.4f", line_filled_variance_Y_1D->GetStdDev()));
    c1 -> Print((output_directory + "/line_filled_variance_Y_1D.pdf").c_str());
    c1 -> Clear();
    
    c1 -> cd();
    line_filled_covariance_1D -> Draw("hist");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Work-in-progress"));
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", line_filled_covariance_1D->GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", line_filled_covariance_1D->GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", line_filled_covariance_1D->GetStdDev()));
    c1 -> Print((output_directory + "/line_filled_covariance_1D.pdf").c_str());
    c1 -> Clear();

    cout<<"line filled covariance info : "<<line_filled_covariance_1D->GetMean()<<" "<<line_filled_covariance_1D->GetStdDev()<<endl;
    cout<<"line filled vairance X info : "<<line_filled_variance_X_1D->GetMean()<<" "<<line_filled_variance_X_1D->GetStdDev()<<endl;
    cout<<"line filled vairance Y info : "<<line_filled_variance_Y_1D->GetMean()<<" "<<line_filled_variance_Y_1D->GetStdDev()<<endl;
    cout<<"line filled stdDev X   info : " <<line_filled_stddev_X_1D->GetMean()<<" "<<line_filled_stddev_X_1D->GetStdDev()<<endl;
    cout<<"line filled stdDev Y   info : " <<line_filled_stddev_Y_1D->GetMean()<<" "<<line_filled_stddev_Y_1D->GetStdDev()<<endl;

    return 0;
}