#include "../vtxZDist_comp.h"

R__LOAD_LIBRARY(../libvtxZDist_comp.so)

int Run_vtxZDist_comp()
{
    std::string sPH_labeling = "Internal";
    std::string MC_labeling = "Simulation";

    std::vector<std::pair<std::string, std::pair<std::string,std::string>>> data_input_directory_pair_vec = {
        {
            "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/seflgendata/run_54280_HR_Dec042024/completed/Run24NewCode_EvtVtxZTracklet/completed/vtxZDist/completed/Data_vtxZDist_ApplyCut_EvtBcoFullDiffCut61_DataOldVtxXY_00054280_merged.root",
            {"Data, |MbdVtxZ| < 60 cm, |BcoFullDiff| > 61", "data_OldVtxXY_WithVtxZQA"}
        }
    };

    std::vector<std::pair<std::string, std::pair<std::string,std::string>>> MC_input_directory_pair_vec = {
        {
            "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/sPH_dNdeta/Run24AuAuMC/Sim_Ntuple_HIJING_ana443_20241102/Run24NewCode_EvtVtxZTracklet/completed/vtxZDist/completed/MC_vtxZDist_ApplyCut_DataOldVtxXY_merged.root",
            {"HIJING, |MbdVtxZ| < 60 cm", "HIJING_noZWeight_WithVtxZQA"}
        }
    };

    std::map<std::string, std::vector<std::tuple<double,double,std::string>>> labelling_vec_map = {
        {
            "Data",
            {
                std::make_tuple(0.22, 0.8, "Old vertex XY"),
                std::make_tuple(0.22, 0.76, "With INTT vtxZ QA")
            }
        },

        {
            "MC",
            {
                std::make_tuple(0.22, 0.8, "Truth avg VtxXY"),
                std::make_tuple(0.22, 0.76, "With INTT vtxZ QA")
            }
        },


        {
            "Comp",
            {
                std::make_tuple(0.22, 0.8, "Data, Old vertex XY"),
                std::make_tuple(0.22, 0.76, "MC, Truth avg vtxXY"),
                std::make_tuple(0.22, 0.72, "With INTT vtxZ QA")
            }
        }
    };
    

    std::string output_directory = "/sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/NewCode2024/PrepareAllDist/vtxZDist_comp/plots/Test_20241225_withVtxZQA";

    bool WithVtxZReWeighting = false;

    vtxZDist_comp * VZDC = new vtxZDist_comp(
        sPH_labeling,
        MC_labeling,
        data_input_directory_pair_vec,
        MC_input_directory_pair_vec,
        labelling_vec_map,
        output_directory,
        WithVtxZReWeighting
    );

    VZDC -> MakeDataPlot("hist");
    VZDC -> MakeMCPlot("hist");
    VZDC -> MakeComparisonPlot();
    VZDC -> MakeVtxZCheckPlot();
    VZDC -> PrepareINTTvtxZReWeight();

    return 0;
}