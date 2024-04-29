#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

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

int main (/*int argc, char* argv[]*/)
{
    // // note : Check if the correct number of arguments is provided
    // if (argc != 2) {
    //     std::cerr << "Usage: " << argv[0] << " N_event" << std::endl;
    //     return 1;
    // }

    // cout<<"running the folder_"<<argv[1]<<endl;

    // // note : Convert the second argument (argv[1]) to an integer
    // int folder_index = std::atoi(argv[1]);

    string input_directory = "/sphenix/user/ChengWei/sPH_dNdeta/HIJING_ana398_xvtx-0p04cm_yvtx0p24cm_zvtx-20cm_dummyAlignParams/Geo_Z_scan_trial_3/complete_file";
    string input_file_name = "file_list_zvtx_dist.txt"; 
    string output_directory = input_directory + "/final_merged_result";
    pair<double,double> true_vtxXY = {-0.4, 2.4}; // note : [mm]

    double in_reco_avgX;
    double in_reco_avgY;
    double in_zvtx_dist_mean;
    double in_zvtx_dist_width;

    SetsPhenixStyle();
    TCanvas * c1 = new TCanvas("c1","c1",950,800);

    TLatex * ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    ltx->SetTextAlign(31);

    TLatex * draw_text = new TLatex();
    draw_text -> SetNDC();
    draw_text -> SetTextSize(0.025);

    int file_counting = 0; 

    TChain * chain_in = new TChain("tree");
    for (int i = 0; i < 450; i++)
    {
        vector<string> folder_file_list = read_list((input_directory + "/folder_"+ to_string(i)).c_str(), input_file_name);
        cout<<"folder_"<<i<<" N ready files : "<<folder_file_list.size()<<endl;
        file_counting += folder_file_list.size();
        for(string file_name : folder_file_list) {chain_in -> Add(file_name.c_str());}
    }
    
    cout<<"Expected total files : "<<file_counting<<endl;
    cout<<"total events : "<<chain_in -> GetEntries()<<endl;

    chain_in -> SetBranchAddress("reco_avgX", &in_reco_avgX);
    chain_in -> SetBranchAddress("reco_avgY", &in_reco_avgY);
    chain_in -> SetBranchAddress("zvtx_dist_mean", &in_zvtx_dist_mean);
    chain_in -> SetBranchAddress("zvtx_dist_width", &in_zvtx_dist_width);

    TH1F * zvtx_peak_variation_1D = new TH1F("zvtx_peak_variation_1D", "zvtx_peak_variation_1D;zvtx_diff dist peak [mm];Entry", 100, -3,3); // note : unit [mm]
    TH1F * zvtx_width_variation_1D = new TH1F("zvtx_width_variation_1D", "zvtx_width_variation_1D:zvtx_diff dist width [mm];Entry", 100, -6, 6); // note : unit [mm]
    TH2F * recoXY_distance_zvtx_peak_2D = new TH2F("recoXY_distance_zvtx_peak_2D", "recoXY_distance_zvtx_peak_2D;vtxXY deviation [mm];zvtx_diff dist peak [mm]",    100, -0.5, 2, 100, -3, 3);
    TH2F * recoXY_distance_zvtx_width_2D = new TH2F("recoXY_distance_zvtx_width_2D", "recoXY_distance_zvtx_width_2D;vtxXY deviation [mm];zvtx_diff dist width [mm]", 100, -0.5, 2, 100, -6, 6);
    
    // map<int,int> centrality_map = {
    //     {5, 0},
    //     {15, 1},
    //     {25, 2},
    //     {35, 3},
    //     {45, 4},
    //     {55, 5},
    //     {65, 6},
    //     {75, 7},
    //     {85, 8},
    //     {95, 9}
    // };

    // TFile * file_out = new TFile((output_directory + "/combine_file.root").c_str(), "RECREATE");
    // TH1F * zvtx_diff[11];

    // for (int i = 0; i < 11; i++)
    // {
    //     zvtx_diff[i] = new TH1F(("zvtx_diff_M_" + to_string(i)).c_str(), ("zvtx_diff_M_" + to_string(i)+";Reco Z - True Z [mm];Entry").c_str(), 100, -30, 30); // note : unit mm
    // }

    for (int event_i = 0; event_i < chain_in -> GetEntries(); event_i++)
    {
        chain_in -> GetEntry(event_i);

        double vtxXY_deviation = sqrt(pow(in_reco_avgX - true_vtxXY.first, 2) + pow(in_reco_avgY - true_vtxXY.second,2));

        zvtx_peak_variation_1D -> Fill(in_zvtx_dist_mean);
        zvtx_width_variation_1D -> Fill(in_zvtx_dist_width);
        recoXY_distance_zvtx_peak_2D -> Fill(vtxXY_deviation, in_zvtx_dist_mean);
        recoXY_distance_zvtx_width_2D -> Fill(vtxXY_deviation, in_zvtx_dist_width);

        // if (nclu_inner == -1 || nclu_outer == -1) {continue;}
        // if (good_zvtx_tag == false) {continue;}
        // zvtx_diff[centrality_map[centrality_bin]] -> Fill(reco_zvtx - true_zvtx);
        // zvtx_diff[10]                             -> Fill(reco_zvtx - true_zvtx);
    }

    c1 -> cd();
    zvtx_peak_variation_1D -> Draw("hist");
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", zvtx_peak_variation_1D -> GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", zvtx_peak_variation_1D -> GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", zvtx_peak_variation_1D -> GetStdDev()));
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} simulation"));
    c1 -> Print((output_directory + "/zvtx_peak_variation_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    zvtx_width_variation_1D -> Draw("hist");
    draw_text -> DrawLatex(0.2, 0.84, Form("Entry: %.0f", zvtx_width_variation_1D -> GetEntries()));
    draw_text -> DrawLatex(0.2, 0.80, Form("Mean: %.4f", zvtx_width_variation_1D -> GetMean()));
    draw_text -> DrawLatex(0.2, 0.76, Form("StdDev: %.4f", zvtx_width_variation_1D -> GetStdDev()));
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} simulation"));
    c1 -> Print((output_directory + "/zvtx_width_variation_1D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    recoXY_distance_zvtx_peak_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} simulation"));
    c1 -> Print((output_directory + "/recoXY_distance_zvtx_peak_2D.pdf").c_str());
    c1 -> Clear();

    c1 -> cd();
    recoXY_distance_zvtx_width_2D -> Draw("colz0");
    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX INTT}} simulation"));
    c1 -> Print((output_directory + "/recoXY_distance_zvtx_width_2D.pdf").c_str());
    c1 -> Clear();

    // file_out -> cd();
    // for (int i = 0; i < 11; i++)
    // {
    //     zvtx_diff[i] -> Write();
    // }

    // file_out->Close();
}