#include <iostream>
#include <fstream>
using namespace std;

#include <vector>

#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>

#include "sPhenixStyle.h"


vector<string> read_list(string folder_direction, string file_name_in)
{
    vector<string> file_list;
    string list_buffer;
    ifstream data_list;
	data_list.open((folder_direction + "/" + file_name_in).c_str());

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
	cout<<"size in the" <<file_name_in<<": "<<file_list.size()<<endl;

    return file_list;
}

int main()
{
    string input_folder = "/sphenix/user/ChengWei/sPH_dNdeta/HIJING_ana398_xvtx-0p04cm_yvtx0p24cm_zvtx-20cm_dummyAlignParams/Geo_XY_scan_trial_1/complete_file/full_merged_file";
    string file_list_name = "file_list.txt";
    string output_directory = input_folder + "/" + "merged_result";
    pair<double,double> true_vtx = {-0.4, 2.4}; // note : unit mm

    TChain *chain_in = new TChain("tree_geo_scan");
    vector<string> file_list = read_list(input_folder, file_list_name);
    for (string file_name : file_list) { chain_in -> Add((file_name).c_str()); }
    cout<<"chain_in -> GetEntries() : "<<chain_in -> GetEntries()<<endl;   

    vector<double> * offset_x_vec; offset_x_vec = 0;
    vector<double> * offset_y_vec; offset_y_vec = 0;
    double total_offset_x = 0;
    double total_offset_y = 0;

    double DCA_inner_fitYpos;
    double angle_diff_inner_fitYpos;
    double DCA_inner_fitE;
    double angle_diff_inner_fitE;
    double DCA_outer_fitYpos;
    double angle_diff_outer_fitYpos;
    double DCA_outer_fitE;
    double angle_diff_outer_fitE;
    int    random_seed_out;
    double vtxX;
    double vtxY;
    double trial_originX;
    double trial_originY;
    double angle_diff_mean;
    double angle_diff_stddev;
    double DCA_distance_mean;
    double DCA_distance_stddev;

    chain_in -> SetBranchAddress("offset_x_vec", &offset_x_vec);
    chain_in -> SetBranchAddress("offset_y_vec", &offset_y_vec);
    chain_in -> SetBranchAddress("total_offset_x", &total_offset_x);
    chain_in -> SetBranchAddress("total_offset_y", &total_offset_y);
    chain_in -> SetBranchAddress("DCA_inner_fitYpos", &DCA_inner_fitYpos);
    chain_in -> SetBranchAddress("DCA_inner_fitE", &DCA_inner_fitE);
    chain_in -> SetBranchAddress("angle_diff_inner_fitYpos", &angle_diff_inner_fitYpos);
    chain_in -> SetBranchAddress("angle_diff_inner_fitE", &angle_diff_inner_fitE);
    chain_in -> SetBranchAddress("DCA_outer_fitYpos", &DCA_outer_fitYpos);
    chain_in -> SetBranchAddress("DCA_outer_fitE", &DCA_outer_fitE);
    chain_in -> SetBranchAddress("angle_diff_outer_fitYpos", &angle_diff_outer_fitYpos);
    chain_in -> SetBranchAddress("angle_diff_outer_fitE", &angle_diff_outer_fitE);
    chain_in -> SetBranchAddress("random_seed", &random_seed_out);
    chain_in -> SetBranchAddress("vtxX", &vtxX);
    chain_in -> SetBranchAddress("vtxY", &vtxY);
    chain_in -> SetBranchAddress("trial_originX", &trial_originX);
    chain_in -> SetBranchAddress("trial_originY", &trial_originY);
    chain_in -> SetBranchAddress("angle_diff_mean",     &angle_diff_mean);
    chain_in -> SetBranchAddress("angle_diff_stddev",   &angle_diff_stddev);
    chain_in -> SetBranchAddress("DCA_distance_mean",   &DCA_distance_mean);
    chain_in -> SetBranchAddress("DCA_distance_stddev", &DCA_distance_stddev);

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    ltx->SetTextAlign(31);

    TLatex * draw_text = new TLatex();
    draw_text -> SetNDC();
    draw_text -> SetTextSize(0.03);

    double ladder_offset_l = -0.4;
    double ladder_offset_r = 0.4;
    double single_offset_l = -0.5;
    double single_offset_r = 11.5;
    double total_offset_l = -0.5;
    double total_offset_r = 20.;
    double angle_diff_l = -0.1;
    double angle_diff_r = 0.2;
    double vtx_deviation_l = -0.5;
    double vtx_deviation_r = 3.;
    double angle_diff_fitE_l = -0.05;
    double angle_diff_fitE_r = 0.2;


    TH2F * ladder_offset_XY_2D = new TH2F("ladder_offset_XY_2D", "ladder_offset_XY_2D;Ladder X offset [mm];Ladder Y offset [mm]", 100, ladder_offset_l, ladder_offset_r, 100, ladder_offset_l, ladder_offset_r);
    TH2F * total_offset_XY_2D = new TH2F("total_offset_XY_2D", "total_offset_XY_2D;Total X offset [mm];Total Y offset [mm]", 100, single_offset_l, single_offset_r, 100, single_offset_l, single_offset_r);
    TH2F * total_offset_VtxDeviation_2D = new TH2F("total_offset_VtxDeviation_2D","total_offset_VtxDeviation_2D;Total XY offset [mm];VTX deviation [mm]", 100, total_offset_l, total_offset_r, 100, vtx_deviation_l, vtx_deviation_r);
    double window_size = 1.5;
    TH2F * vtxXY_distribution_2D = new TH2F("vtxXY_distribution_2D","vtxXY_distribution_2D;Vtx X [mm];Vtx Y [mm]", 100, true_vtx.first-window_size, true_vtx.first+window_size, 100, true_vtx.second-window_size, true_vtx.second+window_size);
    TH2F * angle_diff_stddev_deviation_2D = new TH2F("angle_diff_stddev_deviation_2D", "angle_diff_stddev_deviation_2D;#Delta#phi StdDev [deg];VTX deviation [mm]", 100, angle_diff_l, angle_diff_r,100, vtx_deviation_l, vtx_deviation_r);
    
    TH2F * total_offset_angle_diff_stddev_2D = new TH2F("total_offset_angle_diff_stddev_2D", "total_offset_angle_diff_stddev_2D;Total XY offset [mm];#Delta#phi StdDev [deg]", 100, total_offset_l, total_offset_r, 100, angle_diff_l, angle_diff_r);
    TH2F * total_offset_angle_diff_inner_fitE_2D = new TH2F("total_offset_angle_diff_inner_fitE_2D", "total_offset_angle_diff_inner_fitE_2D;Total XY offset [mm];#Delta#phi Inner FitE [deg]", 100, total_offset_l, total_offset_r, 100, angle_diff_fitE_l, angle_diff_fitE_r);


    for (int i = 0; i < chain_in -> GetEntries(); i++)
    {
        chain_in -> GetEntry(i);

        for (int ladder_i = 0; ladder_i < offset_x_vec -> size(); ladder_i++)
        {
            ladder_offset_XY_2D -> Fill(offset_x_vec -> at(ladder_i), offset_y_vec -> at(ladder_i));
        }

        double total_vtx_deviation = sqrt(pow(vtxX - true_vtx.first,2) + pow(vtxY - true_vtx.second,2));

        total_offset_XY_2D                    -> Fill(total_offset_x, total_offset_y);
        total_offset_VtxDeviation_2D          -> Fill(total_offset_x+total_offset_y, total_vtx_deviation);
        vtxXY_distribution_2D                 -> Fill(vtxX,vtxY);
        angle_diff_stddev_deviation_2D        -> Fill(angle_diff_stddev, total_vtx_deviation);
        total_offset_angle_diff_stddev_2D     -> Fill(total_offset_x+total_offset_y, angle_diff_stddev);
        total_offset_angle_diff_inner_fitE_2D -> Fill(total_offset_x+total_offset_y, angle_diff_inner_fitE);
         

    }

    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("c1","c1",950,800);
    c1 -> cd();

    c1 -> cd();
    total_offset_XY_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Simulation"));
    c1 -> Print((output_directory + "/total_offset_XY_2D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    total_offset_VtxDeviation_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Simulation"));
    c1 -> Print((output_directory + "/total_offset_VtxDeviation_2D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    vtxXY_distribution_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Simulation"));
    c1 -> Print((output_directory + "/vtxXY_distribution_2D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    angle_diff_stddev_deviation_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Simulation"));
    c1 -> Print((output_directory + "/angle_diff_stddev_deviation_2D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    total_offset_angle_diff_stddev_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Simulation"));
    c1 -> Print((output_directory + "/total_offset_angle_diff_stddev_2D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    total_offset_angle_diff_inner_fitE_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Simulation"));
    c1 -> Print((output_directory + "/total_offset_angle_diff_inner_fitE_2D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    ladder_offset_XY_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} Simulation"));
    c1 -> Print((output_directory + "/ladder_offset_XY_2D.pdf").c_str());
    c1 -> Clear();

    return 0;
}