#include "../FinalResult.h"
R__LOAD_LIBRARY(../libFinalResult.so)

int Run_Final_new()
{
    int Mbin = 14; std::string Centrality_str = "65-70";
    double Hist_Y_max = 1000;

    int runnumber = 54280;
    std::string sPH_label = "Internal";

    std::string mother_folder = "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/seflgendata/run_54280_HR_Jan172025/Run4/EvtVtxZ/FinalResult";
    std::string range_folder = "vtxZ_-10_10cm_MBin" + std::to_string(Mbin);

    std::string StandardData_directory = mother_folder + "/" + range_folder + "/Folder_BaseLine/Run_0/completed";
    std::string StandardData_file_name = Form("Data_PreparedNdEtaEach_AlphaCorr_AllSensor_VtxZ10_Mbin%d_00054280_00000_dNdEta.root",Mbin);
    std::string StandardMC_directory = StandardData_directory;
    std::string StandardMC_file_name = Form("MC_PreparedNdEtaEach_AlphaCorr_AllSensor_VtxZ10_Mbin%d_00002_dNdEta.root",Mbin);

    FinalResult * final = new FinalResult(
        runnumber,
        Mbin,
        StandardData_directory,
        StandardData_file_name,
        StandardMC_directory,
        StandardMC_file_name,
        sPH_label,
        mother_folder + "/" + range_folder
    );

    final -> SetEtaRange({-1.5, 1.5});
    final -> SetCollisionStr({{0.2, 0.96}, "Au+Au #sqrt{s_{NN}} = 200 GeV"});
    final -> SetAnaDescription({{0.21, 0.9}, Form("Centrality [%s]%%, VtxZ [-10, 10] cm", Centrality_str.c_str())});

    final -> SetFinal_Data_MC_text(
        {
            "Data (PHOBOS approach)",
            "HIJING (generator)"
        }
    );


    std::string output_folder_name = final -> GetOutputFileName();
    std::cout << "output_folder_name = " << output_folder_name << std::endl;

    final -> PrepareStatisticalError();




    // Division : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    std::string RunSegMother_directory = mother_folder + "/" + range_folder + "/Folder_RunSegment"; 
    final -> PrepareRunSegmentError(
        {
            RunSegMother_directory + Form("/Run_0/completed/%s",StandardData_file_name.c_str()),
            RunSegMother_directory + Form("/Run_1/completed/%s",StandardData_file_name.c_str())
        }
    ); // note : run1, run2

    // Division : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    std::string ClusAdcMother_directory = mother_folder + "/" + range_folder + "/Folder_ClusAdcCut";
    final -> PrepareClusAdcError(
        {
            ClusAdcMother_directory + Form("/Run_0/completed/%s",StandardData_file_name.c_str())
        }
    ); 
    
    // Division : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    final -> PrepareGeoOffsetError(
        "/sphenix/user/ChengWei/sPH_dNdeta/Run24AuAuMC/Sim_Ntuple_HIJING_ana443_20241102/GeoOffset_v1/completed/merged_result/FromdNdEta.root",
        StandardData_directory + "/" + Form("MC_PreparedNdEtaEach_AllSensor_VtxZ10_Mbin%d_00001_dNdEta.root",Mbin) // note : for the alpha correction
    );

    // Division : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    std::string DeltaPhiMother_directory = mother_folder + "/" + range_folder + "/Folder_DeltaPhiCut";
    final -> PrepareDeltaPhiError(
        {
            DeltaPhiMother_directory + Form("/Run_0/completed/%s",StandardData_file_name.c_str()),
            DeltaPhiMother_directory + Form("/Run_1/completed/%s",StandardData_file_name.c_str())
        }

    ); // note : 0.018 and 0.024

    // Division : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    std::string ClusPhiSizeMother_directory = mother_folder + "/" + range_folder + "/Folder_ClusPhiCut";
    final -> PrepareClusPhiSizeError(
        {
            ClusPhiSizeMother_directory + Form("/Run_0/completed/%s",StandardData_file_name.c_str())
        }
    );

    // Division : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // std::string MCMergedMother_directory = CW_folder + data_folder_mother + "/EvtVtxZ/TrackHist/completed";
    // final -> PrepareMCMergedError(
    //     {
    //         MCMergedMother_directory + "/dNdEta_MCMergeError1/dNdEta_AllSensor_GeoAccCorr_VtxZ10_Mbin70/completed/Data_PreparedNdEtaEach_AlphaCorr_GeoAccCorr_AllSensor_VtxZ10_Mbin70_00054280_00000_dNdEta.root",
    //         MCMergedMother_directory + "/dNdEta_MCMergeError2/dNdEta_AllSensor_GeoAccCorr_VtxZ10_Mbin70/completed/Data_PreparedNdEtaEach_AlphaCorr_GeoAccCorr_AllSensor_VtxZ10_Mbin70_00054280_00000_dNdEta.root",
    //         MCMergedMother_directory + "/dNdEta_MCMergeError3/dNdEta_AllSensor_GeoAccCorr_VtxZ10_Mbin70/completed/Data_PreparedNdEtaEach_AlphaCorr_GeoAccCorr_AllSensor_VtxZ10_Mbin70_00054280_00000_dNdEta.root",
    //         MCMergedMother_directory + "/dNdEta_MCMergeError4/dNdEta_AllSensor_GeoAccCorr_VtxZ10_Mbin70/completed/Data_PreparedNdEtaEach_AlphaCorr_GeoAccCorr_AllSensor_VtxZ10_Mbin70_00054280_00000_dNdEta.root"
    //     }
    // );

    // Division : -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    final -> PrepareFinalError();
    final -> PrepareFinalResult(Hist_Y_max);
    final -> EndRun();

    return 888;
}