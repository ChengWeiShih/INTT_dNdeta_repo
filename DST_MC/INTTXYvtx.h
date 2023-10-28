#ifndef INTTDSTchain_h
#define INTTDSTchain_h

#include "/sphenix/user/ChengWei/INTT/INTT_commissioning/INTT_CW/INTT_commissioning/DAC_Scan/InttConversion_new.h"
#include "/sphenix/user/ChengWei/INTT/INTT_commissioning/INTT_CW/INTT_commissioning/DAC_Scan/InttClustering.h"
#include "../sigmaEff.h"
#include "../sPhenixStyle.C"
#include "../INTTDSTchain.C"

double cos_func(double *x, double *par)
{
    return -1 * par[0] * cos(par[1] * (x[0] + par[2])) + par[3];
}

class INTTXYvtx {
    public :
        INTTXYvtx ();
    
    private :
        double get_radius(double x, double y);        
}

double get_radius(double x, double y)
{
    return sqrt(pow(x,2)+pow(y,2));
}

#endif
