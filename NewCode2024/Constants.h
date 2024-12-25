#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants{
    std::vector<double> centrality_edges = {0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    // note : vtxZQA
    std::pair<double, double> cut_vtxZDiff = {-3, 4.}; // note : MBDz - INTTz
    std::pair<double, double> cut_TrapezoidalFitWidth = {1.5, 5.5};
    std::pair<double, double> cut_TrapezoidalFWHM = {2,8};
    std::pair<double, double> cut_INTTvtxZError = {-10000, 100000};

    // note : for analysis
    std::pair<double, double> cut_GlobalMBDvtxZ = {-60, 60};
    std::pair<double, double> cut_AnaVtxZ = {-10, 10};

    int cut_InttBcoFullDIff_next = 61;

    int Semi_inclusive_bin = 7;

    int HighNClus = 500;
}

#endif