#ifndef GausFunc_h
#define GausFunc_h

double gaus_func(double *x, double *par)
{
    // note : par[0] : size
    // note : par[1] : mean
    // note : par[2] : width
    // note : par[3] : offset 
    return par[0] * TMath::Gaus(x[0],par[1],par[2]) + par[3];
}

double gaus_pol1_func(double *x, double *par)
{
    // note : par[0] : size
    // note : par[1] : mean
    // note : par[2] : width
    // note : par[3] : offset 
    return par[0] * TMath::Gaus(x[0],par[1],par[2]) + par[3] + par[4]*x[0];
}

double d_gaus_pol1_func(double *x, double *par)
{

    // note : par[0] : size
    // note : par[1] : ratio of the two gaussians
    // note : par[2] : mean
    // note : par[3] : width of gaus 1
    // note : par[4] : width of gaus 2
    // note : par[5] : offset
    // note : par[6] : slope
    return par[0] * ( (1. - par[1]) * TMath::Gaus(x[0],par[2],par[3]) + par[1] * TMath::Gaus(x[0],par[2],par[4]) ) + par[5] + par[6] * x[0];
}

double d_gaus_func(double *x, double *par)
{

    // note : par[0] : size
    // note : par[1] : ratio of the two gaussians
    // note : par[2] : mean
    // note : par[3] : width of gaus 1
    // note : par[4] : width of gaus 2
    return par[0] * ( (1. - par[1]) * TMath::Gaus(x[0],par[2],par[3]) + par[1] * TMath::Gaus(x[0],par[2],par[4]) );
}



#endif