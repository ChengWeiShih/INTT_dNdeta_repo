#include "../RestDist.h"

R__LOAD_LIBRARY(../libRestDist.so)

TH1D * GetReweighting_hist(string input_map_directory, string map_name)
{
  TFile * file_in = TFile::Open(Form("%s", input_map_directory.c_str()));
  TH1D * h1D_INTT_vtxZ_reweighting = (TH1D*)file_in->Get(map_name.c_str()); // todo : the map of the vtxZ reweighting
  return h1D_INTT_vtxZ_reweighting;
}

void Run_RestDist_MC(
  int process_id = 0,
  int run_num = 54280,
  int nevents = -1,
  string input_directory = "/sphenix/user/ChengWei/INTT/INTT/general_codes/CWShih/INTTBcoResolution/macro",
  string input_filename = "file_list_54280_intt.txt",
  string output_directory = "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/seflgendata/run_54280/completed/BCO_check",
  
  // todo : modify here
  std::string output_file_name_suffix = "_MCTrueXY",

  bool Apply_cut = true,
  bool ApplyVtxZReWeighting = true,
  std::pair<bool, int> ApplyEvtBcoFullDiffCut = {false, 61},
  std::pair<bool, std::pair<double,double>> RequireVtxZRange = {true, {-10, 10}},
  std::pair<bool, std::pair<double,double>> isClusQA = {true, {35, 200000}}, // note : adc, phi size

  std::string vtxZReWeighting_input_directory = "/sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/NewCode2024/PrepareAllDist/vtxZDist_comp/plots/Test_20241225_withVtxZQA/INTTvtxZReWeight.root",
  std::string map_name = "HIJING_noZWeight_WithVtxZQA_Inclusive70"
)
{

  string final_output_file_name;

  TH1D * h1D_INTT_vtxZ_reweighting = GetReweighting_hist(vtxZReWeighting_input_directory, map_name);

  // Division : ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // note : with vtxZ range cut, and no ClusQA
  RestDist * RDs1 = new RestDist(
    process_id,
    run_num,
    nevents,
    input_directory,
    input_filename,
    output_directory,

    output_file_name_suffix,

    Apply_cut,
    ApplyVtxZReWeighting,
    ApplyEvtBcoFullDiffCut,
    
    RequireVtxZRange,
    {false,{-10,20000}}
  );
  RDs1->SetINTTvtxZReweighting(h1D_INTT_vtxZ_reweighting);

  final_output_file_name = RDs1->GetOutputFileName();
  cout<<"final_output_file_name: "<<final_output_file_name<<endl;

  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  

  RDs1->PrepareEvent();
  RDs1->EndRun();


  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

  // Division : ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // note : with vtxZ range cut, and w/ ClusQA
  RestDist * RDs2 = new RestDist(
    process_id,
    run_num,
    nevents,
    input_directory,
    input_filename,
    output_directory,

    output_file_name_suffix,

    Apply_cut,
    ApplyVtxZReWeighting,
    ApplyEvtBcoFullDiffCut,
    
    RequireVtxZRange,
    isClusQA
  );
  RDs2->SetINTTvtxZReweighting(h1D_INTT_vtxZ_reweighting);

  final_output_file_name = RDs2->GetOutputFileName();
  cout<<"final_output_file_name: "<<final_output_file_name<<endl;

  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  

  RDs2->PrepareEvent();
  RDs2->EndRun();


  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));


  // Division : ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // note : No vtxZ cut, No ClusQA
  RestDist * RDs3 = new RestDist(
    process_id,
    run_num,
    nevents,
    input_directory,
    input_filename,
    output_directory,

    output_file_name_suffix,

    Apply_cut,
    ApplyVtxZReWeighting,
    ApplyEvtBcoFullDiffCut,
    
    {false, {-1000,1000}},
    {false,{-10,20000}}
  );
  RDs3->SetINTTvtxZReweighting(h1D_INTT_vtxZ_reweighting);

  final_output_file_name = RDs3->GetOutputFileName();
  cout<<"final_output_file_name: "<<final_output_file_name<<endl;

  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  

  RDs3->PrepareEvent();
  RDs3->EndRun();


  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

  return;
}
