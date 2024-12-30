#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>

namespace Constants{
    extern std::vector<double> centrality_edges;

    // note : vtxZQA
    extern std::pair<double, double> cut_vtxZDiff; // note : MBDz - INTTz
    extern std::pair<double, double> cut_TrapezoidalFitWidth;
    extern std::pair<double, double> cut_TrapezoidalFWHM;
    extern std::pair<double, double> cut_INTTvtxZError;

    // note : for analysis
    extern std::pair<double, double> cut_GlobalMBDvtxZ;
    extern std::pair<double, double> cut_AnaVtxZ;

    extern int cut_InttBcoFullDIff_next;

    extern int Semi_inclusive_bin;

    extern int HighNClus;

    extern double typeA_sensor_half_length_incm; // note : [cm]
    extern double typeB_sensor_half_length_incm; // note : [cm] 

    extern int B0L0_index;
    extern int B0L1_index;
    extern int B1L0_index;
    extern int B1L1_index;
    extern int nLadder_inner;
    extern int nLadder_outer;
}

#endif