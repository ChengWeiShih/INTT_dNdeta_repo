#include "Constants.h"

namespace Constants{
    // std::vector<double> centrality_edges = {0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<double> centrality_edges = {0, 3, 6, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
    std::vector<std::string> centrality_text = {
        "0-3%",
        "3-6%",
        "6-10%",
        "10-15%",
        "15-20%",
        "20-25%",
        "25-30%",
        "30-35%",
        "35-40%",
        "40-45%",
        "45-50%",
        "50-55%",
        "55-60%",
        "60-65%",
        "65-70%",
        "70-75%",
        "75-80%",
        "80-85%",
        "85-90%",
        "90-95%",
        "95-100%"
    };

    // note : vtxZQA
    std::pair<double, double> cut_vtxZDiff = {-3, 4.}; // note : MBDz - INTTz
    std::pair<double, double> cut_TrapezoidalFitWidth = {1.5, 5.5};
    std::pair<double, double> cut_TrapezoidalFWHM = {2,8};
    std::pair<double, double> cut_INTTvtxZError = {-10000, 100000};

    // note : for analysis
    std::pair<double, double> cut_GlobalMBDvtxZ = {-60, 60};
    std::pair<double, double> cut_AnaVtxZ = {-10, 10}; // note : used by vtxZDist.cpp

    int cut_InttBcoFullDIff_next = 61;

    int Semi_inclusive_bin = 14;
    int Semi_inclusive_interval = 70; // note : 70, for example

    int HighNClus = 500;

    double VtxZEdge_min = -45; // note : cm // note : used by vtxZDist.cpp
    double VtxZEdge_max = 45; // note : cm  // note : used by vtxZDist.cpp
    int nVtxZBin = 18;                      // note : used by vtxZDist.cpp

    double cut_GoodRecoVtxZ = 1.; // note : cm // note : used by vtxZDist.cpp

    // note : for the column multiplicity correction
    // note : this is for the ZID
    int nZbin = 100;
    double Zmin = -25;
    double Zmax = 25;

    // note : almost no change
    double EtaEdge_min = -2.7; // note : used by RestDist.h
    double EtaEdge_max = 2.7;  // note : used by RestDist.h
    int nEtaBin = 27;          // note : used by RestDist.h

    // note : below should never be changed
    double typeA_sensor_half_length_incm = 0.8; // note : [cm]
    double typeB_sensor_half_length_incm = 1.0; // note : [cm] 

    int B0L0_index = 3;
    int B0L1_index = 4;
    int B1L0_index = 5;
    int B1L1_index = 6;
    int nLadder_inner = 12;
    int nLadder_outer = 16;
}