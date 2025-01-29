#include "../ColumnCheck.h"

R__LOAD_LIBRARY(../libColumnCheck.so)

void Run_ColumnCheck(
  int process_id = 0,
  int run_num = 54280,
  int nevents = -1,
  string input_directory = "/sphenix/user/ChengWei/INTT/INTT/general_codes/CWShih/INTTBcoResolution/macro",
  string input_filename = "file_list_54280_intt.txt",
  string output_directory = "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/seflgendata/run_54280/completed/BCO_check",
  
  // todo : modify here
  std::string output_file_name_suffix = "",
  std::pair<double, double> vertexXYIncm = {-0.0232717, 0.223173},
  double SetMbinFloat = 0.5, // note : 0 - 1
            
  std::pair<double, double> VtxZRange = {-30, 30},
  bool IsZClustering = false,
  bool BcoFullDiffCut = true,
  std::pair<bool, std::pair<double, double>> isClusQA = {true, {35,40}}, // note : {adc, phi size}
 
  bool ColMulMask = false
)
{

  ColumnCheck * CC = new ColumnCheck(
    process_id,
    run_num,
    nevents,
    input_directory,
    input_filename,
    output_directory,

    output_file_name_suffix,
    vertexXYIncm,
    SetMbinFloat,

    VtxZRange,
    IsZClustering,
    BcoFullDiffCut,
    isClusQA, // note : {adc, phi size}

    ColMulMask
  );

  string final_output_file_name = CC->GetOutputFileName();
  cout<<"final_output_file_name: "<<final_output_file_name<<endl;
  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  

  CC -> MainProcess();
  CC -> EndRun();

  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

  return;
}
