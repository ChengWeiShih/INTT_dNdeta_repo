#ifndef INTTReadTree_h
#define INTTReadTree_h

#include "../INTTDSTchain.C"
#include "/sphenix/user/ChengWei/INTT/INTT_commissioning/INTT_CW/INTT_commissioning/DAC_Scan/InttConversion_new.h"
#include "/sphenix/user/ChengWei/INTT/INTT_commissioning/INTT_CW/INTT_commissioning/DAC_Scan/InttClustering.h"

class INTTReadTree
{
    public : 
        vector<clu_info> temp_sPH_inner_nocolumn_vec; 
        vector<clu_info> temp_sPH_outer_nocolumn_vec; 
        vector<vector<double>> temp_sPH_nocolumn_vec;
        vector<vector<double>> temp_sPH_nocolumn_rz_vec;
        
        INTTReadTree(int data_type, string input_directory, string MC_list_name, string tree_name, int clu_size_cut, int clu_sum_adc_cut);
        void EvtInit(long long event_i);
        void EvtSetCluGroup();
        long long GetNEvt();
        unsigned long GetEvtNClus();
        unsigned long GetEvtNClusPost();
        double GetTrigZvtxMC();
        bool GetPhiCheckTag();
        int GetNvtxMC();
        void EvtClear();

    private : 
        string data_type_list[3] = {"MC","data_DST","data_private"};
        long long N_event;

        string input_directory; 
        string MC_list_name;
        string tree_name;
        double clu_sum_adc_cut;
        int clu_size_cut;
        int data_type;   
        int NvtxMC; 
        double TrigZvtxMC;

        TChain * chain_in;
        INTTDSTchain * inttDSTMC;

        unsigned long evt_length;

        void TChainInit_MC();
        double get_radius(double x, double y);
        

};

INTTReadTree::INTTReadTree(int data_type, string input_directory, string MC_list_name, string tree_name, int clu_size_cut, int clu_sum_adc_cut)
:data_type(data_type), input_directory(input_directory), MC_list_name(MC_list_name), tree_name(tree_name), clu_size_cut(clu_size_cut), clu_sum_adc_cut(clu_sum_adc_cut)
{
    temp_sPH_inner_nocolumn_vec.clear();
    temp_sPH_outer_nocolumn_vec.clear();
    
    temp_sPH_nocolumn_vec.clear();
    temp_sPH_nocolumn_vec = vector<vector<double>> (2);

    temp_sPH_nocolumn_rz_vec.clear();
    temp_sPH_nocolumn_rz_vec = vector<vector<double>> (2);
    
    if (data_type_list[data_type] == "MC") 
    {
        std::cout<<"--- INTTReadTree -> input data tupe : MC ---"<<std::endl;
        TChainInit_MC();
        std::cout<<"--- INTTReadTree -> Initialization done ---"<<std::endl;
    }
}

void INTTReadTree::TChainInit_MC()
{
    chain_in = new TChain(tree_name.c_str());
    inttDSTMC = new INTTDSTchain(chain_in, input_directory, MC_list_name);
    N_event = chain_in->GetEntries();
    std::printf("inttDSTMC N event : %lli \n", N_event);
}

void INTTReadTree::EvtInit(long long event_i)
{
    inttDSTMC->LoadTree(event_i);
    inttDSTMC->GetEntry(event_i);

    evt_length = inttDSTMC->NClus;
    NvtxMC     = inttDSTMC->NTruthVtx;
    TrigZvtxMC = inttDSTMC->TruthPV_trig_z;
}

long long INTTReadTree::GetNEvt() { return N_event; }
unsigned long INTTReadTree::GetEvtNClus() { return evt_length; }
double INTTReadTree::GetTrigZvtxMC() {return TrigZvtxMC;}
int INTTReadTree::GetNvtxMC() {return NvtxMC;}

unsigned long INTTReadTree::GetEvtNClusPost() 
{ 
    return temp_sPH_inner_nocolumn_vec.size() + temp_sPH_outer_nocolumn_vec.size(); 
}

void INTTReadTree::EvtSetCluGroup()
{
    if (data_type_list[data_type] == "MC"){
        for (int clu_i = 0; clu_i < evt_length; clu_i++)
        {
            if (int(inttDSTMC -> ClusPhiSize -> at(clu_i)) > clu_size_cut) continue; 
            if (int(inttDSTMC -> ClusAdc -> at(clu_i)) < clu_sum_adc_cut) continue;

            double clu_x = inttDSTMC -> ClusX -> at(clu_i) * 10; // note : change the unit from cm to mm
            double clu_y = inttDSTMC -> ClusY -> at(clu_i) * 10;
            double clu_z = inttDSTMC -> ClusZ -> at(clu_i) * 10;
            double clu_phi = (clu_y < 0) ? atan2(clu_y,clu_x) * (180./TMath::Pi()) + 360 : atan2(clu_y,clu_x) * (180./TMath::Pi());
            int    clu_layer = (inttDSTMC -> ClusLayer -> at(clu_i) == 3 || inttDSTMC -> ClusLayer -> at(clu_i) == 4) ? 0 : 1;
            double clu_radius = get_radius(clu_x, clu_y);

            temp_sPH_nocolumn_vec[0].push_back( clu_x );
            temp_sPH_nocolumn_vec[1].push_back( clu_y );
            
            temp_sPH_nocolumn_rz_vec[0].push_back( clu_z );
            temp_sPH_nocolumn_rz_vec[1].push_back( ( clu_phi > 180 ) ? clu_radius * -1 : clu_radius );
            

            if (clu_layer == 0) {// note : inner
                temp_sPH_inner_nocolumn_vec.push_back({
                    -1, 
                    -1, 
                    int(inttDSTMC -> ClusAdc -> at(clu_i)), 
                    int(inttDSTMC -> ClusAdc -> at(clu_i)), 
                    int(inttDSTMC -> ClusPhiSize -> at(clu_i)), 
                    clu_x, 
                    clu_y, 
                    clu_z, 
                    clu_layer, 
                    clu_phi
                });
            }
            
            if (clu_layer == 1) {// note : outer
                temp_sPH_outer_nocolumn_vec.push_back({
                    -1, 
                    -1, 
                    int(inttDSTMC -> ClusAdc -> at(clu_i)), 
                    int(inttDSTMC -> ClusAdc -> at(clu_i)), 
                    int(inttDSTMC -> ClusPhiSize -> at(clu_i)), 
                    clu_x, 
                    clu_y, 
                    clu_z, 
                    clu_layer, 
                    clu_phi
                });            
            }        
        }
    }
}

double INTTReadTree::get_radius(double x, double y)
{
    return sqrt(pow(x,2)+pow(y,2));
}

bool INTTReadTree::GetPhiCheckTag()
{
    int inner_1_check = 0;
    int inner_2_check = 0;
    int inner_3_check = 0;
    int inner_4_check = 0;
    for (int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++) {
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >= 0 && temp_sPH_inner_nocolumn_vec[inner_i].phi < 90)    inner_1_check = 1;
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >= 90 && temp_sPH_inner_nocolumn_vec[inner_i].phi < 180)  inner_2_check = 1;
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >= 180 && temp_sPH_inner_nocolumn_vec[inner_i].phi < 270) inner_3_check = 1;
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >= 270 && temp_sPH_inner_nocolumn_vec[inner_i].phi < 360) inner_4_check = 1;

        if ( (inner_1_check + inner_2_check + inner_3_check + inner_4_check) == 4 ) break;
    }

    int outer_1_check = 0;
    int outer_2_check = 0;
    int outer_3_check = 0;
    int outer_4_check = 0;
    for (int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++) {
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >= 0 && temp_sPH_outer_nocolumn_vec[outer_i].phi < 90)    outer_1_check = 1;
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >= 90 && temp_sPH_outer_nocolumn_vec[outer_i].phi < 180)  outer_2_check = 1;
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >= 180 && temp_sPH_outer_nocolumn_vec[outer_i].phi < 270) outer_3_check = 1;
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >= 270 && temp_sPH_outer_nocolumn_vec[outer_i].phi < 360) outer_4_check = 1;

        if ( (outer_1_check + outer_2_check + outer_3_check + outer_4_check) == 4 ) break;
    }

    if ( (inner_1_check + inner_2_check + inner_3_check + inner_4_check + outer_1_check + outer_2_check + outer_3_check + outer_4_check) != 8 ) {return false;}
    else { return true; }
}

void INTTReadTree::EvtClear()
{
    temp_sPH_inner_nocolumn_vec.clear();
    temp_sPH_outer_nocolumn_vec.clear();

    temp_sPH_nocolumn_vec.clear();
    temp_sPH_nocolumn_vec = vector<vector<double>> (2);
    
    temp_sPH_nocolumn_rz_vec.clear();
    temp_sPH_nocolumn_rz_vec = vector<vector<double>> (2);   

    evt_length = 0;
    NvtxMC     = 0;
}

#endif