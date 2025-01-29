#include "FinalResult.h"

FinalResult::FinalResult(
    int runnumber_in,
    int Mbin_in,
    std::string StandardData_directory_in,
    std::string StandardData_file_name_in,
    std::string StandardMC_directory_in,
    std::string StandardMC_file_name_in,
    std::string sPH_label_in,
    std::string Output_directory_in
):
    runnumber(runnumber_in),
    Mbin(Mbin_in),
    StandardData_directory(StandardData_directory_in),
    StandardData_file_name(StandardData_file_name_in),
    StandardMC_directory(StandardMC_directory_in),
    StandardMC_file_name(StandardMC_file_name_in),
    sPH_label(sPH_label_in),
    Output_directory(Output_directory_in)
{
    Output_directory = (Output_directory == "Not_given") ? StandardData_directory : Output_directory;

    PrepareOutputFolderName();
    system(Form("if [ ! -d %s/%s/completed ]; then mkdir -p %s/%s/completed; fi;", Output_directory.c_str(), output_folder_name.c_str(), Output_directory.c_str(), output_folder_name.c_str()));
    final_output_directory = Form("%s/%s", Output_directory.c_str(), output_folder_name.c_str());

    file_in_data_standard = TFile::Open(Form("%s/%s", StandardData_directory.c_str(), StandardData_file_name.c_str()));
    h1D_data_standard = (TH1D*) file_in_data_standard -> Get(StandardData_h1D_name.c_str());

    file_in_MC_standard = TFile::Open(Form("%s/%s", StandardMC_directory.c_str(), StandardMC_file_name.c_str()));
    h1D_truth_standard = (TH1D*) file_in_MC_standard -> Get(StandardTruth_h1D_name.c_str());

    if (
        h1D_data_standard == nullptr ||
        h1D_truth_standard == nullptr || 

        h1D_data_standard -> GetNbinsX()             != h1D_truth_standard -> GetNbinsX() ||
        h1D_data_standard -> GetXaxis() -> GetXmin() != h1D_truth_standard -> GetXaxis() -> GetXmin() ||
        h1D_data_standard -> GetXaxis() -> GetXmax() != h1D_truth_standard -> GetXaxis() -> GetXmax()
    ){
        std::cout << "Error : the standard data or MC is not found" << std::endl;
        exit(1);
    }

    h1D_RunSegmentError_vec.clear();
    h1D_ClusAdcError_vec.clear();
    h1D_GeoOffsetError_vec.clear();
    h1D_DeltaPhiError_vec.clear();

    SetsPhenixStyle();
    c1 = new TCanvas("c1", "c1", 950, 800);

    ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextSize(0.045);
    ltx->SetTextAlign(31);

    draw_text = new TLatex();
    draw_text -> SetNDC();
    draw_text -> SetTextSize(0.03);

    leg_errors = new TLegend(0.31,0.67,0.41,0.87);
    leg_errors -> SetBorderSize(0);
    leg_errors -> SetTextSize(0.03);

    leg_final = new TLegend(0.45,0.8,0.8,0.85);
    leg_final -> SetBorderSize(0);
    leg_final -> SetTextSize(0.03);

    file_out = new TFile(Form("%s/%s.root", final_output_directory.c_str(), output_folder_name.c_str()), "RECREATE");
    
}

void FinalResult::PrepareOutputFolderName()
{
    std::string runnumber_str = std::to_string( runnumber );
    if (runnumber != -1){
        int runnumber_str_len = 8;
        runnumber_str.insert(0, runnumber_str_len - runnumber_str.size(), '0');
    }

    output_folder_name = "Final";
    output_folder_name += "_Mbin" + std::to_string(Mbin);
    
    output_folder_name += Form("_%s",runnumber_str.c_str());
}


void FinalResult::PrepareStatisticalError()
{
    h1D_error_statistic = (TH1D*) h1D_data_standard -> Clone("h1D_error_statistic");
    h1D_error_statistic -> Reset("ICESM");

    for (int i = 1; i <= h1D_data_standard -> GetNbinsX(); i++){
        
        double error = h1D_data_standard -> GetBinError(i) / h1D_data_standard -> GetBinContent(i);
        error = (error == error && error != 1) ? error : 0;
        h1D_error_statistic -> SetBinContent(i, error);
    }

    h1D_error_statistic -> SetMarkerStyle(20);
    h1D_error_statistic -> SetMarkerSize(0.8);
    h1D_error_statistic -> SetLineWidth(0);
    h1D_error_statistic -> SetLineColorAlpha(1,0);
    h1D_error_statistic -> SetMarkerColor(TColor::GetColor(color_code[0].c_str()));
    
}


void FinalResult::PrepareRunSegmentError(std::vector<std::string> file_directory_vec_in)
{
    TFile * temp_file_in;

    for (std::string file_dir : file_directory_vec_in)
    {
        temp_file_in = TFile::Open(file_dir.c_str());

        h1D_RunSegmentError_vec.push_back( (TH1D*) temp_file_in -> Get(StandardData_h1D_name.c_str()) );
        h1D_RunSegmentError_vec.back() -> Divide(h1D_data_standard);
        h1D_to_AbsRatio(h1D_RunSegmentError_vec.back());
    }

    h1D_error_Run_segmentation = h1D_FindLargestOnes("h1D_error_Run_segmentation", h1D_RunSegmentError_vec);

    h1D_error_Run_segmentation -> SetMarkerStyle(20);
    h1D_error_Run_segmentation -> SetMarkerSize(0.8);
    h1D_error_Run_segmentation -> SetLineWidth(0);
    h1D_error_Run_segmentation -> SetLineColorAlpha(1,0);
    h1D_error_Run_segmentation -> SetMarkerColor(TColor::GetColor(color_code[1].c_str()));
}

void FinalResult::PrepareClusAdcError(std::vector<std::string> file_directory_vec_in)
{
    TFile * temp_file_in;

    for (std::string file_dir : file_directory_vec_in)
    {
        temp_file_in = TFile::Open(file_dir.c_str());

        h1D_ClusAdcError_vec.push_back( (TH1D*) temp_file_in -> Get(StandardData_h1D_name.c_str()) );
        h1D_ClusAdcError_vec.back() -> Divide(h1D_data_standard);
        h1D_to_AbsRatio(h1D_ClusAdcError_vec.back());
    }

    h1D_error_ClusAdc = h1D_FindLargestOnes("h1D_error_ClusAdc", h1D_ClusAdcError_vec);

    h1D_error_ClusAdc -> SetMarkerStyle(20);
    h1D_error_ClusAdc -> SetMarkerSize(0.8);
    h1D_error_ClusAdc -> SetLineWidth(0);
    h1D_error_ClusAdc -> SetLineColorAlpha(1,0);
    h1D_error_ClusAdc -> SetMarkerColor(TColor::GetColor(color_code[2].c_str()));

}

void FinalResult::PrepareGeoOffsetError(std::string file_directory_in, std::string alpha_corr_directory_in)
{
    TFile * temp_file_in_GeoUnc = TFile::Open(file_directory_in.c_str());
    TFile * temp_file_in_AlphaCorr = TFile::Open(alpha_corr_directory_in.c_str());

    h1D_error_GeoOffset = (TH1D*) h1D_data_standard -> Clone("h1D_error_GeoOffset");

    auto temp_h1D_GeoUnc = (TH1D*) temp_file_in_GeoUnc -> Get("h1D_RotatedBkg_RecoTrackletEtaPerEvt_VariationMax");
    auto temp_h1D_AlphaCorr = (TH1D*) temp_file_in_AlphaCorr -> Get("h1D_RotatedBkg_alpha_correction");

    h1D_error_GeoOffset -> Divide(
        temp_h1D_GeoUnc,
        temp_h1D_AlphaCorr 
    );

    for (int i = 1; i <= h1D_error_GeoOffset -> GetNbinsX(); i++){
        std::cout<<"GeoOffset: bin: "<< i <<", eta: ["<<h1D_error_GeoOffset->GetXaxis()->GetBinLowEdge(i)<<", "<<h1D_error_GeoOffset->GetXaxis()->GetBinUpEdge(i)<<"]"
        << ", GeoUnc: " << temp_h1D_GeoUnc -> GetBinContent(i)
        << ", AlphaCorr: " << temp_h1D_AlphaCorr -> GetBinContent(i)
        <<", final_error: " << h1D_error_GeoOffset -> GetBinContent(i)
        << std::endl;
    }
    

    h1D_error_GeoOffset -> SetMarkerStyle(20);
    h1D_error_GeoOffset -> SetMarkerSize(0.8);
    h1D_error_GeoOffset -> SetLineWidth(0);
    h1D_error_GeoOffset -> SetLineColorAlpha(1,0);
    h1D_error_GeoOffset -> SetMarkerColor(TColor::GetColor(color_code[3].c_str()));
}


void FinalResult::PrepareDeltaPhiError(std::vector<std::string> file_directory_vec_in)
{
    TFile * temp_file_in;

    for (std::string file_dir : file_directory_vec_in)
    {
        temp_file_in = TFile::Open(file_dir.c_str());

        h1D_DeltaPhiError_vec.push_back( (TH1D*) temp_file_in -> Get(StandardData_h1D_name.c_str()) );
        h1D_DeltaPhiError_vec.back() -> Divide(h1D_data_standard);
        h1D_to_AbsRatio(h1D_DeltaPhiError_vec.back());
    }

    h1D_error_DeltaPhi = h1D_FindLargestOnes("h1D_error_DeltaPhi", h1D_DeltaPhiError_vec);

    h1D_error_DeltaPhi -> SetMarkerStyle(20);
    h1D_error_DeltaPhi -> SetMarkerSize(0.8);
    h1D_error_DeltaPhi -> SetLineWidth(0);
    h1D_error_DeltaPhi -> SetLineColorAlpha(1,0);
    h1D_error_DeltaPhi -> SetMarkerColor(TColor::GetColor(color_code[4].c_str()));
}

void FinalResult::PrepareClusPhiSizeError(std::vector<std::string> file_directory_vec_in)
{
    TFile * temp_file_in;

    for (std::string file_dir : file_directory_vec_in)
    {
        temp_file_in = TFile::Open(file_dir.c_str());

        h1D_ClusPhiSizeError_vec.push_back( (TH1D*) temp_file_in -> Get(StandardData_h1D_name.c_str()) );
        h1D_ClusPhiSizeError_vec.back() -> Divide(h1D_data_standard);
        h1D_to_AbsRatio(h1D_ClusPhiSizeError_vec.back());
    }

    h1D_error_ClusPhiSize = h1D_FindLargestOnes("h1D_error_ClusPhiSize", h1D_ClusPhiSizeError_vec);

    h1D_error_ClusPhiSize -> SetMarkerStyle(20);
    h1D_error_ClusPhiSize -> SetMarkerSize(0.8);
    h1D_error_ClusPhiSize -> SetLineWidth(0);
    h1D_error_ClusPhiSize -> SetLineColorAlpha(1,0);
    h1D_error_ClusPhiSize -> SetMarkerColor(TColor::GetColor(color_code[5].c_str()));

}

void FinalResult::PrepareMCMergedError(std::vector<std::string> file_directory_vec_in)
{
    TFile * temp_file_in;

    for (std::string file_dir : file_directory_vec_in)
    {
        temp_file_in = TFile::Open(file_dir.c_str());

        h1D_MCMergedError_vec.push_back( (TH1D*) temp_file_in -> Get(StandardData_h1D_name.c_str()) );
        h1D_MCMergedError_vec.back() -> Divide(h1D_data_standard);
        h1D_to_AbsRatio(h1D_MCMergedError_vec.back());
    }

    h1D_error_MCMerged = h1D_FindLargestOnes("h1D_error_MCMerged", h1D_MCMergedError_vec);

    h1D_error_MCMerged -> SetMarkerStyle(20);
    h1D_error_MCMerged -> SetMarkerSize(0.8);
    h1D_error_MCMerged -> SetLineWidth(0);
    h1D_error_MCMerged -> SetLineColorAlpha(1,0);
    h1D_error_MCMerged -> SetMarkerColor(TColor::GetColor(color_code[6].c_str()));

}


void FinalResult::PrepareFinalError()
{
    if (h1D_error_statistic == nullptr){
        std::cout << "Error : the statistical error is not found" << std::endl;
        exit(1);
    }

    h1D_error_Final = (TH1D*) h1D_error_statistic -> Clone("h1D_error_Final");
    h1D_error_Final -> Reset("ICESM");

    for (int i = 1; i <= h1D_error_Final -> GetNbinsX(); i++){

        double _error_statistic = (h1D_error_statistic != nullptr) ? h1D_error_statistic -> GetBinContent(i) : 0;
        _error_statistic = (_error_statistic == _error_statistic && _error_statistic != 1) ? _error_statistic : 0;

        double _error_Run_segmentation = (h1D_error_Run_segmentation != nullptr) ? h1D_error_Run_segmentation -> GetBinContent(i) : 0;
        _error_Run_segmentation = (_error_Run_segmentation == _error_Run_segmentation && _error_Run_segmentation != 1) ? _error_Run_segmentation : 0;

        double _error_ClusAdc = (h1D_error_ClusAdc != nullptr) ? h1D_error_ClusAdc -> GetBinContent(i) : 0;
        _error_ClusAdc = (_error_ClusAdc == _error_ClusAdc && _error_ClusAdc != 1) ? _error_ClusAdc : 0;

        double _error_GeoOffset = (h1D_error_GeoOffset != nullptr) ? h1D_error_GeoOffset -> GetBinContent(i) : 0;
        _error_GeoOffset = (_error_GeoOffset == _error_GeoOffset && _error_GeoOffset != 1) ? _error_GeoOffset : 0;

        double _error_DeltaPhi = (h1D_error_DeltaPhi != nullptr) ? h1D_error_DeltaPhi -> GetBinContent(i) : 0;
        _error_DeltaPhi = (_error_DeltaPhi == _error_DeltaPhi && _error_DeltaPhi != 1) ? _error_DeltaPhi : 0;

        double _error_ClusPhiSize = (h1D_error_ClusPhiSize != nullptr) ? h1D_error_ClusPhiSize -> GetBinContent(i) : 0;
        _error_ClusPhiSize = (_error_ClusPhiSize == _error_ClusPhiSize && _error_ClusPhiSize != 1) ? _error_ClusPhiSize : 0;

        double _error_MCMerged = (h1D_error_MCMerged != nullptr) ? h1D_error_MCMerged -> GetBinContent(i) : 0;
        _error_MCMerged = (_error_MCMerged == _error_MCMerged && _error_MCMerged != 1) ? _error_MCMerged : 0;


        double _final_error = std::sqrt(
            pow(_error_statistic, 2) +
            pow(_error_Run_segmentation, 2) +
            pow(_error_ClusAdc, 2) +
            pow(_error_GeoOffset, 2) +
            pow(_error_DeltaPhi, 2) + 
            pow(_error_ClusPhiSize, 2) +
            pow(_error_MCMerged, 2)
        );

        std::cout<<"bin: "<< i <<", eta: ["<<h1D_error_Final->GetXaxis()->GetBinLowEdge(i)<<", "<<h1D_error_Final->GetXaxis()->GetBinUpEdge(i)<<"]"
        << ", statistic: " << _error_statistic 
        << ", Run_segmentation: " << _error_Run_segmentation 
        << ", ClusAdc: " << _error_ClusAdc 
        << ", GeoOffset: " << _error_GeoOffset 
        << ", DeltaPhi: " << _error_DeltaPhi 
        << ", ClusPhiSize: " << _error_ClusPhiSize
        << ", MCMerged: " << _error_MCMerged
        << " final_error: " << _final_error 
        << std::endl;

        h1D_error_Final -> SetBinContent(i, _final_error);
    }

    // todo : the bin removal is here
    for (int i = 1; i <= h1D_error_Final -> GetNbinsX(); i++)
    {
        double lowEdge = h1D_error_Final -> GetXaxis() -> GetBinLowEdge(i);
        double upEdge = h1D_error_Final -> GetXaxis() -> GetBinUpEdge(i);
        double binCenter = h1D_error_Final -> GetXaxis() -> GetBinCenter(i);

        if (binCenter < eta_range.first || binCenter > eta_range.second){
            
            std::cout<<"Rejection, bin: "<< i <<", eta: ["<<lowEdge<<", "<<upEdge<<"]"
            << ", Total error: " << h1D_error_Final -> GetBinContent(i)
            << std::endl;

            h1D_error_Final -> SetBinContent(i, -30);
            h1D_error_Final -> SetBinError(i, 0);
            
            if (h1D_error_statistic != nullptr) {
                h1D_error_statistic -> SetBinContent(i, -30);
                h1D_error_statistic -> SetBinError(i, 0);
            }

            if (h1D_error_Run_segmentation != nullptr) {
                h1D_error_Run_segmentation -> SetBinContent(i, -30);
                h1D_error_Run_segmentation -> SetBinError(i, 0);
            }

            if (h1D_error_ClusAdc != nullptr) {
                h1D_error_ClusAdc -> SetBinContent(i, -30);
                h1D_error_ClusAdc -> SetBinError(i, 0);
            }

            if (h1D_error_GeoOffset != nullptr) {
                h1D_error_GeoOffset -> SetBinContent(i, -30);
                h1D_error_GeoOffset -> SetBinError(i, 0);
            }

            if (h1D_error_DeltaPhi != nullptr) {
                h1D_error_DeltaPhi -> SetBinContent(i, -30);
                h1D_error_DeltaPhi -> SetBinError(i, 0);
            }

            if (h1D_error_ClusPhiSize != nullptr) {
                h1D_error_ClusPhiSize -> SetBinContent(i, -30);
                h1D_error_ClusPhiSize -> SetBinError(i, 0);
            }

            if (h1D_error_MCMerged != nullptr) {
                h1D_error_MCMerged -> SetBinContent(i, -30);
                h1D_error_MCMerged -> SetBinError(i, 0);
            }

        }
    }    


    // Division : ------------------------------------------------------------------------------------------------------------------------------------
    // h1D_error_Final -> SetMaximum(h1D_error_Final -> GetBinContent(h1D_error_Final -> GetMaximumBin()) * 2.);
    h1D_error_Final -> SetMinimum(0);
    h1D_error_Final -> SetMaximum(0.11);
    h1D_error_Final -> SetLineWidth(2);
    h1D_error_Final -> SetLineColor(1);
    h1D_error_Final -> GetXaxis() -> SetTitle("#eta");
    h1D_error_Final -> GetYaxis() -> SetTitle("Relative uncertainty [%]");

    leg_errors -> AddEntry(h1D_error_Final, "Total Uncertainty", "l");
    if (h1D_error_statistic != nullptr) {leg_errors -> AddEntry(h1D_error_statistic, "Stat. Unc.", "p");}
    if (h1D_error_Run_segmentation != nullptr) {leg_errors -> AddEntry(h1D_error_Run_segmentation, "Run segmentation variation", "p");}
    if (h1D_error_ClusAdc != nullptr) {leg_errors -> AddEntry(h1D_error_ClusAdc, "Cluster ADC variation", "p");}
    if (h1D_error_GeoOffset != nullptr) {leg_errors -> AddEntry(h1D_error_GeoOffset, "Geo. misalignment variation", "p");}
    if (h1D_error_DeltaPhi != nullptr) {leg_errors -> AddEntry(h1D_error_DeltaPhi, "#Delta#phi cut variation", "p");}
    if (h1D_error_ClusPhiSize != nullptr) {leg_errors -> AddEntry(h1D_error_ClusPhiSize, "Cluster #phi size variation", "p");}
    if (h1D_error_MCMerged != nullptr) {leg_errors -> AddEntry(h1D_error_MCMerged, "MC merge variation", "p");}

    file_out -> cd();
    
    c1 -> cd();
    h1D_error_Final -> Draw("hist");

    if (h1D_error_statistic != nullptr) {h1D_error_statistic -> Draw("p same");}
    if (h1D_error_Run_segmentation != nullptr) {h1D_error_Run_segmentation -> Draw("p same");}
    if (h1D_error_ClusAdc != nullptr) {h1D_error_ClusAdc -> Draw("p same");}
    if (h1D_error_GeoOffset != nullptr) {h1D_error_GeoOffset -> Draw("p same");}
    if (h1D_error_DeltaPhi != nullptr) {h1D_error_DeltaPhi -> Draw("p same");}
    if (h1D_error_ClusPhiSize != nullptr) {h1D_error_ClusPhiSize -> Draw("p same");}
    if (h1D_error_MCMerged != nullptr) {h1D_error_MCMerged -> Draw("p same");}

    if (AnaDescription.second.size() > 0){
        draw_text -> DrawLatex(AnaDescription.first.first, AnaDescription.first.second, AnaDescription.second.c_str());
    }

    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", sPH_label.c_str()));    

    leg_errors -> Draw("same");

    c1 -> Print(Form("%s/SystUnc_Summary.pdf", final_output_directory.c_str()));
    c1 -> Write("SystUnc_Summary");

    c1 -> Clear();
    c1 -> cd();
    if (h1D_error_statistic != nullptr) {
        h1D_error_statistic -> Draw("p");
        c1 -> Print(Form("%s/SystUnc_Stat.pdf", final_output_directory.c_str()));
    }
}

void FinalResult::PrepareFinalResult()
{
    h1D_data_standard -> SetMarkerStyle(20);
    h1D_data_standard -> SetMarkerSize(0.8);
    h1D_data_standard -> SetMarkerColor(1);
    h1D_data_standard -> SetFillColorAlpha(1,0.5);
    h1D_data_standard -> SetLineColorAlpha(1,0);
    h1D_data_standard -> SetLineWidth(0);

    gE_data_final = new TGraphErrors();
    gE_data_final -> SetMarkerStyle(20);
    gE_data_final -> SetMarkerSize(0.8);
    gE_data_final -> SetMarkerColor(1);
    gE_data_final -> SetFillColorAlpha(1,0.5);
    gE_data_final -> SetLineColorAlpha(1,0);
    gE_data_final -> SetLineWidth(0);

    for (int i = 1; i <= h1D_data_standard -> GetNbinsX(); i++){
        double _content = h1D_data_standard -> GetBinContent(i);
        double _error = h1D_error_Final -> GetBinContent(i);
        double _BinCenter = h1D_data_standard -> GetXaxis() -> GetBinCenter(i);
        double _BinWidth = h1D_data_standard -> GetXaxis() -> GetBinWidth(i);
        
        double final_error = (_content != 0) ? _content * _error : 0.;


        h1D_data_standard -> SetBinError(i, _error * _content);

        if (_BinCenter >= eta_range.first && _BinCenter <= eta_range.second){
            gE_data_final -> SetPoint(gE_data_final -> GetN(), _BinCenter, _content);
            gE_data_final -> SetPointError(gE_data_final -> GetN() - 1, 0, final_error);
        }

        // if (_content == 0) {h1D_data_standard -> SetBinContent(i, -10);}
        if (_BinCenter < eta_range.first || _BinCenter > eta_range.second){
            h1D_data_standard -> SetBinContent(i, -10);
            h1D_data_standard -> SetBinError(i, 0);
        }


        std::cout<<"bin: "<< i <<", eta: ["<<h1D_data_standard->GetXaxis()->GetBinLowEdge(i)<<", "<<h1D_data_standard->GetXaxis()->GetBinUpEdge(i)<<"]"
        << ", content: " << _content
        << ", error: " << _error
        << ", final_error: " << final_error
        << std::endl;
    }

    gE_data_final -> GetYaxis() -> SetRangeUser(0, 400);
    gE_data_final -> GetXaxis() -> SetTitle("#eta");
    gE_data_final -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");

    h1D_data_standard -> SetMinimum(0);
    h1D_data_standard -> SetMaximum(400); // todo: the maximum
    h1D_data_standard -> GetXaxis() -> SetTitle("#eta");
    h1D_data_standard -> GetYaxis() -> SetTitle("dN_{ch}/d#eta");


    h1D_truth_standard -> SetMinimum(0);
    h1D_truth_standard -> SetMaximum(400); // todo: the maximum
    h1D_truth_standard -> SetLineColor(TColor::GetColor("#3288bd"));
    h1D_truth_standard -> SetLineWidth(2);
    h1D_truth_standard -> SetFillColorAlpha(1,0);

    leg_final -> AddEntry(h1D_data_standard, Final_Data_MC_text.first.c_str(), "fp");
    leg_final -> AddEntry(h1D_truth_standard, Final_Data_MC_text.second.c_str(), "l");

    // Division : ------------------------------------------------------------------------------------------------------------------------------------
    file_out -> cd();
    c1 -> cd();
    c1 -> Clear();
    
    // h1D_data_standard -> Draw("E3");
    h1D_truth_standard -> Draw("hist");
    gE_data_final -> Draw("E3");
    gE_data_final -> Draw("p same");
    // h1D_data_standard -> Draw("p same");

    if (AnaDescription.second.size() > 0){
        draw_text -> DrawLatex(AnaDescription.first.first, AnaDescription.first.second, AnaDescription.second.c_str());
    }

    if (Collision_str.second.size() > 0){
        draw_text -> DrawLatex(Collision_str.first.first, Collision_str.first.second, Collision_str.second.c_str());
    }

    ltx->DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("#it{#bf{sPHENIX}} %s", sPH_label.c_str()));
    leg_final -> Draw("same");

    c1 -> Print(Form("%s/FinalResult.pdf", final_output_directory.c_str()));
    c1 -> Write("FinalResult");

}

void FinalResult::EndRun()
{
    gE_data_final -> Write("gE_dNdEta_reco");
    h1D_data_standard -> Write("h1D_dNdEta_reco");
    h1D_truth_standard -> Write("h1D_dNdEta_truth");
    file_out -> Close();
}

void FinalResult::h1D_to_AbsRatio(TH1D * h1D_in)
{
    for (int i = 1; i <= h1D_in -> GetNbinsX(); i++){
        h1D_in -> SetBinContent(i, std::abs( h1D_in -> GetBinContent(i) - 1 ));
    }
}

TH1D * FinalResult::h1D_FindLargestOnes(std::string hist_name, std::vector<TH1D*> h1D_vec_in)
{
    TH1D * h1D_out = (TH1D*) h1D_vec_in.front()->Clone(hist_name.c_str());

    for (int i = 1; i < h1D_vec_in.size(); i++)
    {
        for (int j = 1; j <= h1D_out -> GetNbinsX(); j++)
        {
            if (h1D_vec_in[i] -> GetBinContent(j) > h1D_out -> GetBinContent(j)){
                h1D_out -> SetBinContent(j, h1D_vec_in[i] -> GetBinContent(j));
            }
        }
    }

    return h1D_out;
}