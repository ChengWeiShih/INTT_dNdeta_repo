#include "../TrackletHistogramNew.h"

R__LOAD_LIBRARY(../libTrackletHistogramNew.so)

void Run_PrepareHist(
  int process_id = 0,
  int run_num = 54280,
  int nevents = -1,
  string input_directory = "/sphenix/user/ChengWei/INTT/INTT/general_codes/CWShih/INTTBcoResolution/macro",
  string input_filename = "file_list_54280_intt.txt",
  string output_directory = "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/seflgendata/run_54280/completed/BCO_check",
  
  // todo : modify here
  std::string output_file_name_suffix = "_SecondRun",
  std::pair<double, double> vertexXYIncm = {-0.0215087, 0.223197},

  std::pair<bool, TH1D*> vtxZReweight = {false, nullptr},
  bool BcoFullDiffCut = true,
  bool INTT_vtxZ_QA = true,
  std::pair<bool, std::pair<double, double>> isClusQA = {true, {35, 500}}, // note : {adc, phi size}

  bool HaveGeoOffsetTag = false,
  std::pair<bool, int> SetRandomHits = {false, 0}
)
{

  TrackletHistogramNew * TLHN = new TrackletHistogramNew(
    process_id,
    run_num,
    nevents,
    input_directory,
    input_filename,
    output_directory,

    output_file_name_suffix,
    vertexXYIncm,

    vtxZReweight,
    BcoFullDiffCut,
    INTT_vtxZ_QA,
    isClusQA, // note : {adc, phi size}
    HaveGeoOffsetTag,
    SetRandomHits
  );

  string final_output_file_name = TLHN->GetOutputFileName();
  cout<<"final_output_file_name: "<<final_output_file_name<<endl;
  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  

  TLHN -> MainProcess();
  TLHN -> EndRun();

  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

  return;
}
