#include "../TrackletHistogram.h"

R__LOAD_LIBRARY(../libTrackletHistogram.so)

void Run_PrepareHist(
  int process_id = 0,
  int run_num = 54280,
  int nevents = -1,
  string input_directory = "/sphenix/user/ChengWei/INTT/INTT/general_codes/CWShih/INTTBcoResolution/macro",
  string input_filename = "file_list_54280_intt.txt",
  string output_directory = "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/seflgendata/run_54280/completed/BCO_check",
  
  // todo : modify here
  std::string output_file_name_suffix = "_test_N0p07to0p07_nbin140",

  bool vtxZReweight = false,
  bool BcoFullDiffCut = false,
  bool INTT_vtxZ_QA = true,
  bool isWithRotate = true
)
{

  TrackletHistogram * TLH = new TrackletHistogram(
    process_id,
    run_num,
    nevents,
    input_directory,
    input_filename,
    output_directory,

    output_file_name_suffix,

    vtxZReweight,
    BcoFullDiffCut,
    INTT_vtxZ_QA,
    isWithRotate
  );

  string final_output_file_name = TLH->GetOutputFileName();
  cout<<"final_output_file_name: "<<final_output_file_name<<endl;
  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  

  TLH -> MainProcess();
  TLH -> EndRun();

  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

  return;
}
