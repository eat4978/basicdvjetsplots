#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include <TLorentzVector.h>

// Helper functions
// Get sample file
std::vector<std::string> get_files(std::string dsid, bool &stop){
    std::vector<std::string> files;
    std::ifstream fin("allntuples.txt");
    std::string l;
    
    // Check if open
    if(!fin.is_open()){
        std::cout << "Error. Unable to open the file allntuples.txt" << std::endl;
        stop = true;
        fin.close();
        return files;
    }
    
    while(std::getline(fin, l)){
        size_t pos = l.find("group.phys-susy." + dsid);
        if(pos != std::string::npos){
            l.erase(0,pos); //delete everything prior to location found, for local file
            files.push_back(l);
        }
    }
    
    fin.close();
    return files;
}

// Get sample information
std::string get_sample_info(std::string dsid, bool &stop, double &lifetime){
    std::ifstream fin("signal_key.txt");
    char* token;
    char buf[1000];
    const char* const dot = ".";
    
    // Check if open
    if(!fin.is_open()){
        std::cout << "Error. Unable to open the file signal_key.txt" << std::endl;
        stop = true;
        return "Error opening file signal_key.txt";
    }
    
    while(fin.good()){
        std::string l;
        getline(fin, l);
        strcpy(buf, l.c_str());
        
        if(buf[0] == 0){
          continue;
        }
        
        if(l.find("mc16_13TeV." + dsid) != std::string::npos){
            // Extract info we want from string
            // parse
            const std::string tmp1(token = strtok(buf, dot));
            const std::string tmp2((token = strtok(0, dot)) ? token : "");
            const std::string tmp3((token = strtok(0, dot)) ? token : "");
            const std::string tmp4((token = strtok(0, "\0")) ? token : "");
            
            //get lifetime
            unsigned long first = tmp3.find_last_of('_');
            unsigned long last = tmp3.find("ns");
            std::string tao = tmp3.substr(first, last - first);
            tao.erase(0, 1);
            
            for(int i = 0; i < tao.length(); i++){
                if(tao[i] == 'p'){
                    tao[i] = '.';
                    break;
                }
            }
            
            fin.close();
            
            lifetime = std::stod(tao);
            
            // return the string without the 3rd token
            return tmp1 + tmp3 + tmp4;
        }
    }
    
    fin.close();
    std::cout << "Sample with DSID = " << dsid << " not found in signal_key.txt" << std::endl;
    stop = true;
    return "Error locating Sample";
}

// Get cross section in pb, using mass of SUSY particles
double get_cross_section(std::string sampleinfo){
    double cross_x = 1;
    
    if (sampleinfo.find("C1N2N1") != std::string::npos){
        if (sampleinfo.find("_100_") != std::string::npos){
            cross_x = 17.125;
        }
        else if (sampleinfo.find("_500_") != std::string::npos){
            cross_x = 0.3003;
        }
        else if (sampleinfo.find("_700_") != std::string::npos){
            cross_x = 0.3003;
        }
        else if (sampleinfo.find("_900_") != std::string::npos){
            cross_x = 0.00738;
        }
        else if (sampleinfo.find("_1100_") != std::string::npos){
            cross_x = 0.0019200;
        }
        else if (sampleinfo.find("_1300_") != std::string::npos){
            cross_x = 0.000570676;
        }
        else if (sampleinfo.find("_1500_") != std::string::npos){
            cross_x = 6.204689e-05;
        }
        else if (sampleinfo.find("_1600_") != std::string::npos){
            cross_x = 3.615468e-05;
        }
        else if (sampleinfo.find("_1700_") != std::string::npos){
            cross_x = 2.14739e-05;
        }
    }
    
    if(cross_x == 1.0){
        std::cout << "Cross section not implemented for this sample (" << sampleinfo << ". Please update get_cross_section() function!" << std::endl;
    }
    
    return cross_x;
}

//std::vector< std::vector<std::string> > read_csv(std::string fileName){
//    std::vector<std::string> row;
//    std::vector< std::vector<std::string> > data;
//    std::string line, cell;
//    std::ifstream fin(fileName.c_str());
//
//    if(!fin.is_open()){
//        std::cout << "Error. Unable to open the file " << fileName << std::endl;
//        fin.close();
//        return data;
//    }
//
//    while(std::getline(fin, line)){
//        row.clear();
//        if(line[0] == '#'){
//            continue;
//        }
//        std::stringstream str(line);
//        while(getline(str, cell, ',')){
//            row.push_back(cell);
//        }
//        data.push_back(row);
//    }
//
//    // delete title line
//    data.erase(data.begin());
//    data.pop_back();
//
////    // testing
////    for(int i = 0; i < data.size(); i++){
////        for(int j = 0; j < data[i].size(); j++){
////            std::cout << data[i][j] << " ";
////        }
////        std::cout<< "\n";
////    }
//
//    fin.close();
//    return data;
//}

//std::vector<double> get_column(std::vector< std::vector<std::string> > data, int col){
//    std::vector<double> column;
//
//    for(int i = 0; i < data.size(); i++){
//        column.push_back(std::stof(data[i][col]));
//    }
//
//    return column;
//}

//selection sort for pair vectors
void sort(std::vector<double> &lifetime, std::vector<std::string> &dsid, std::vector<double> &simulated, std::vector<double> &expected){
    int startScan, minIndex;
    double minValue, minSimulated, minExpected;
    std::string minDsid;

    for(startScan = 0; startScan < (lifetime.size() - 1); startScan++){
        minIndex = startScan;
        minValue = lifetime[startScan];
        minSimulated = simulated[startScan];
        minExpected = expected[startScan];
        minDsid = dsid[startScan];
        
        for(int index = startScan + 1; index < lifetime.size(); index++){
            if(lifetime[index] < minValue){
                minValue = lifetime[index];
                minSimulated = simulated[index];
                minExpected = expected[index];
                minDsid = dsid[index];
                minIndex = index;
            }
        }
        
        lifetime[minIndex] = lifetime[startScan];
        lifetime[startScan] = minValue;
        
        dsid[minIndex] = dsid[startScan];
        dsid[startScan] = minDsid;
        
        simulated[minIndex] = simulated[startScan];
        simulated[startScan] = minSimulated;
        
        expected[minIndex] = expected[startScan];
        expected[startScan] = minExpected;
    }
}

// Get total expected events
void get_n_expected_events(double &events, TTree* t, long nentries, double weight_nomc){
    for(int j = 0; j < nentries; j++){
        t->GetEntry(j);
        float mcEventWeight;
        t->SetBranchAddress("mcEventWeight", &mcEventWeight);
        std::vector<double> * truthJet_Pt;
        t->SetBranchAddress("truthJet_Pt", &truthJet_Pt);
        
        int truthJet_n = (int)truthJet_Pt->size();
        double weight = weight_nomc * mcEventWeight;
        int nDV_passSel = 0;
        
        //flags
        UChar_t DRAW_pass_triggerFlags;
        UChar_t DRAW_pass_DVJETS;
//        UChar_t DRAW_pass_highPtJetFlag;
        UChar_t BaselineSel_HighPtSR;
        UChar_t BaselineSel_TracklessSR;
        UChar_t DRAW_pass_doubleTracklessJetFlag;
        UChar_t DRAW_pass_singleTracklessJetFlag;
        t->SetBranchAddress("DRAW_pass_triggerFlags", &DRAW_pass_triggerFlags);
        t->SetBranchAddress("DRAW_pass_DVJETS", &DRAW_pass_DVJETS);
        t->SetBranchAddress("BaselineSel_HighPtSR", &BaselineSel_HighPtSR);
//        t->SetBranchAddress("DRAW_pass_highPtJetFlag", &DRAW_pass_highPtJetFlag);
        t->SetBranchAddress("BaselineSel_TracklessSR", &BaselineSel_TracklessSR);
        t->SetBranchAddress("DRAW_pass_doubleTracklessJetFlag", &DRAW_pass_doubleTracklessJetFlag);
        t->SetBranchAddress("DRAW_pass_singleTracklessJetFlag", &DRAW_pass_singleTracklessJetFlag);
        std::vector<bool> * DV_passFiducialCut;
        std::vector<bool> * DV_passDistCut;
        std::vector<bool> * DV_passChiSqCut;
        std::vector<bool> * DV_passMaterialVeto_strict;
        std::vector<bool> * DV_passNTrkCut;
        std::vector<bool> * DV_passMassCut;
        std::vector<int> * DV_nTracksSel;
        t->SetBranchAddress("DV_passFiducialCut", &DV_passFiducialCut);
        t->SetBranchAddress("DV_passDistCut", &DV_passDistCut);
        t->SetBranchAddress("DV_passChiSqCut", &DV_passChiSqCut);
        t->SetBranchAddress("DV_passMaterialVeto_strict", &DV_passMaterialVeto_strict);
        t->SetBranchAddress("DV_passNTrkCut", &DV_passNTrkCut);
        t->SetBranchAddress("DV_passMassCut", &DV_passMassCut);
        t->SetBranchAddress("DV_nTracksSel", &DV_nTracksSel);
        
        // trackless
        if (DRAW_pass_triggerFlags != 1){continue;}
        if (!(DRAW_pass_singleTracklessJetFlag == 1 || DRAW_pass_doubleTracklessJetFlag == 1)){continue;}
        if (BaselineSel_HighPtSR != 0){continue;}
        if (BaselineSel_TracklessSR != 1){continue;}
        
        // Loop through displaced vertices (DV) in event
        for (int k = 0; k < truthJet_n; k++){
            if (!(DV_passFiducialCut->at(k) && DV_passDistCut->at(k) && DV_passChiSqCut->at(k) && DV_passMaterialVeto_strict->at(k) && DV_passNTrkCut->at(k) && DV_passMassCut->at(k) && DV_nTracksSel->at(k) >= 2)){
                continue;
            }
            nDV_passSel++;
        }
        
        // At least 1 DV passing all selections
        if (nDV_passSel >= 1){
            events += weight;
        }
    }
}

int main(int argc, const char * argv[]) {
    // Variables
    std::vector<std::string> dsid_list;
    std::vector<double> lifetime_list;
    std::string user_input;
    std::vector<double> simulated;
    std::vector<double> expected;
    
    // get DSID
    while(true){
        std::cout << "Enter DSID (enter q to quit): ";
        std::cin >> user_input;
        if(user_input == "q" || user_input == "Q"){
            break;
        }
        dsid_list.push_back(user_input);
    }
    
    for(int file = 0; file < dsid_list.size(); file++){
        //variables
        std::string tag;
        bool stop = false;
        double lifetime;
        // counters
        int allevents = 0;
        int n_pass_JetSel = 0;
        int pass_R = 0;
        int pass_Dist = 0;
        int pass_d0 = 0;
        int pass_decay = 0;
        int pass_inv_m = 0;
        double n_expected_events = 0;
        
        double pass_JetSel_e = 0;
        double pass_R_e = 0;
        double pass_Dist_e = 0;
        double pass_d0_e = 0;
        double pass_decay_e = 0;
        double pass_inv_m_e = 0;
        
    //    std::vector< std::vector<std::string> > R_1150 = read_csv("HEPData-ins2628398-v1-event_efficiency_HighPt_R_1150_mm.csv");
    //    std::vector<double> sumpt_1150 = get_column(R_1150, 0);
    //    std::vector<double> effiency_1150 = get_column(R_1150, 3);
    //
    //    std::vector< std::vector<std::string> > R_1150_3870 = read_csv("HEPData-ins2628398-v1-event_efficiency_HighPt_R_1150_3870_mm.csv");
    //    std::vector<double> sumpt_1150_3870 = get_column(R_1150_3870, 0);
    //    std::vector<double> effiency_1150_3870 = get_column(R_1150_3870, 3);
    //
    //    std::vector< std::vector<std::string> > R_3870 = read_csv("HEPData-ins2628398-v1-event_efficiency_HighPt_R_3870_mm.csv");
    //    std::vector<double> sumpt_3870 = get_column(R_3870, 0);
    //    std::vector<double> effiency_3870 = get_column(R_3870, 3);
        
        TFile* file_3870 = new TFile("HEPData-ins2628398-v1-event_efficiency_HighPt_R_3870_mm.root");
        TFile* file_1150 = new TFile("HEPData-ins2628398-v1-event_efficiency_HighPt_R_1150_mm.root");
        TFile* file_1150_3870 = new TFile("HEPData-ins2628398-v1-event_efficiency_HighPt_R_1150_3870_mm.root");
        TH1D* h_3870 = (TH1D*)file_3870->Get("event_efficiency_HighPt_R_3870_mm/Hist1D_y1");
        TH1D* h_1150 = (TH1D*)file_1150->Get("event_efficiency_HighPt_R_1150_mm/Hist1D_y1");
        TH1D* h_1150_3870 = (TH1D*)file_1150_3870->Get("event_efficiency_HighPt_R_1150_3870_mm/Hist1D_y1");

                          
        //set dsid
    //    if(argc < 2){
    //        std::cout << "Please specify dsid!" << std::endl;
    //        return 0;
    //    }
    //
    //    std::string dsid = argv[1];
    //    if (argc > 2) {
    //        tag = argv[2];
    //    }
    //    else{
    //        tag = "v1";
    //    }
    //    std::string dsid = "501190";
        
//        std::string dsid;
//        std::cout << "Enter dsid: ";
//        std::cin >> dsid;
        
        std::string dsid = dsid_list.at(file);
        
        // Gather inputs & cross-section
        std::vector<std::string> inputs = get_files(dsid, stop);
        if(stop){ return 1;}
        std::string info = get_sample_info(dsid, stop, lifetime);
        if(stop){ return 1;}
        lifetime_list.push_back(lifetime);
        std::cout << "Analyzing DSID " << dsid << " with sample info = " << info << std::endl;
        double xsec = get_cross_section(info);
        
        // Load tree and loop through events in ntuple
        for(int i  = 0; i < inputs.size(); i++){
            TFile* inFile = TFile::Open(("/Users/catherine/Documents/URAP/" + dsid + "/" + inputs.at(i)).c_str(), "READ");

            // Load Tree
            TTree* t = (TTree*)inFile->Get("trees_SRDV_");
            long nentries = t->GetEntries();
            
            // Get weight information
            TH1D* metadata = (TH1D*)inFile->Get("MetaData_EventCount");
            double total_sum_weights = metadata->GetBinContent(3);
            
            // Determine luminosity of this sample, pb^-1
            UInt_t runNumber;
            t->SetBranchAddress("runNumber", &runNumber);
            t->GetEntry(0);
            double lumi = 1.0;
            if(runNumber == 284500){
                lumi = 36200.; //mc16a (2015/16)
            }
            else if(runNumber == 300000){
                lumi = 44300.; //mc16d (2017)
            }
            else if(runNumber == 310000){
                lumi = 58500.; //mc16e (2018)
            }
            else{
                std::cout << "Run number, " << runNumber << " not recognized. Setting lumi = 1" << std::endl;
            }

            // Calculate event weight
            double weight_nomc = lumi * xsec / total_sum_weights;
            
    //        get_n_expected_events(n_expected_events, t, nentries, weight_nomc);
            
            // Get total expected events
            for(int j = 0; j < nentries; j++){
                t->GetEntry(j);
                float mcEventWeight;
                t->SetBranchAddress("mcEventWeight", &mcEventWeight);
                std::vector<double> * truthJet_Pt;
                t->SetBranchAddress("truthJet_Pt", &truthJet_Pt);

                int DV_n;
                t->SetBranchAddress("DV_n", &DV_n);
    //            int truthJet_n = (int)truthJet_Pt->size();
                double weight = weight_nomc * mcEventWeight;
                int nDV_passSel = 0;

                //flags
                UChar_t DRAW_pass_triggerFlags;
                UChar_t DRAW_pass_DVJETS;
        //        UChar_t DRAW_pass_highPtJetFlag;
                UChar_t BaselineSel_HighPtSR;
                UChar_t BaselineSel_TracklessSR;
                UChar_t DRAW_pass_doubleTracklessJetFlag;
                UChar_t DRAW_pass_singleTracklessJetFlag;
                t->SetBranchAddress("DRAW_pass_triggerFlags", &DRAW_pass_triggerFlags);
                t->SetBranchAddress("DRAW_pass_DVJETS", &DRAW_pass_DVJETS);
                t->SetBranchAddress("BaselineSel_HighPtSR", &BaselineSel_HighPtSR);
        //        t->SetBranchAddress("DRAW_pass_highPtJetFlag", &DRAW_pass_highPtJetFlag);
                t->SetBranchAddress("BaselineSel_TracklessSR", &BaselineSel_TracklessSR);
                t->SetBranchAddress("DRAW_pass_doubleTracklessJetFlag", &DRAW_pass_doubleTracklessJetFlag);
                t->SetBranchAddress("DRAW_pass_singleTracklessJetFlag", &DRAW_pass_singleTracklessJetFlag);
                std::vector<bool> * DV_passFiducialCut;
                std::vector<bool> * DV_passDistCut;
                std::vector<bool> * DV_passChiSqCut;
                std::vector<bool> * DV_passMaterialVeto_strict;
                std::vector<bool> * DV_passNTrkCut;
                std::vector<bool> * DV_passMassCut;
                std::vector<int> * DV_nTracksSel;
                t->SetBranchAddress("DV_passFiducialCut", &DV_passFiducialCut);
                t->SetBranchAddress("DV_passDistCut", &DV_passDistCut);
                t->SetBranchAddress("DV_passChiSqCut", &DV_passChiSqCut);
                t->SetBranchAddress("DV_passMaterialVeto_strict", &DV_passMaterialVeto_strict);
                t->SetBranchAddress("DV_passNTrkCut", &DV_passNTrkCut);
                t->SetBranchAddress("DV_passMassCut", &DV_passMassCut);
                t->SetBranchAddress("DV_nTracksSel", &DV_nTracksSel);

                // trackless
                if (DRAW_pass_triggerFlags != 1){continue;}
                if (!(DRAW_pass_singleTracklessJetFlag == 1 || DRAW_pass_doubleTracklessJetFlag == 1)){continue;}
                if (BaselineSel_HighPtSR != 0){continue;}
                if (BaselineSel_TracklessSR != 1){continue;}

                // Loop through displaced vertices (DV) in event
                for (int k = 0; k < DV_n; k++){
                    if (!(DV_passFiducialCut->at(k) && DV_passDistCut->at(k) && DV_passChiSqCut->at(k) && DV_passMaterialVeto_strict->at(k) && DV_passNTrkCut->at(k) && DV_passMassCut->at(k) && DV_nTracksSel->at(k) >= 2)){
                        continue;
                    }
                    nDV_passSel++;
                }

                // At least 1 DV passing all selections
                if (nDV_passSel >= 1){
                    n_expected_events += weight;
                }
            }
            
            // Loop through entries in TTree
            for(int j = 0; j < nentries; j++){
                t->GetEntry(j);
                //counters
                int Pt137 = 0;
                int Pt101 = 0;
                int Pt83 = 0;
                int Pt55 = 0;
                int LLP_Pt70 = 0;
                int LLP_Pt50 = 0;
                int n_pass_R = 0;
                int n_pass_Dist = 0;
                int n_pass_d0 = 0;
                int n_pass_decay = 0;
                int n_pass_inv_m = 0;
                
                int sum_pt = 0;
                float Rxy_max = 0;
                double efficiency;
                
                TLorentzVector p_tracks;
                
                //other variables
                float mcEventWeight;
                t->SetBranchAddress("mcEventWeight", &mcEventWeight);
                std::vector<double> * truthJet_Pt;
                t->SetBranchAddress("truthJet_Pt", &truthJet_Pt);
                std::vector<float> * x;
                t->SetBranchAddress("truthSparticle_VtxX", &x);
                std::vector<float> * y;
                t->SetBranchAddress("truthSparticle_VtxY", &y);
                std::vector<float> * z;
                t->SetBranchAddress("truthSparticle_VtxZ", &z);
                std::vector<float> * x_LLP;
                t->SetBranchAddress("truthJet_matchedLLP_x", &x_LLP);
                std::vector<float> * y_LLP;
                t->SetBranchAddress("truthJet_matchedLLP_y", &y_LLP);
                std::vector<int> * truthSparticle_Barcode;
                t->SetBranchAddress("truthSparticle_Barcode", &truthSparticle_Barcode);
                std::vector<int> * truth_LLPChild_Parent_barcode;
                t->SetBranchAddress("truth_LLPChild_Parent_barcode", &truth_LLPChild_Parent_barcode);
                std::vector<double> * truth_LLPChild_d0;
                t->SetBranchAddress("truth_LLPChild_d0", &truth_LLPChild_d0);
                std::vector<double> * truth_LLPChild_Charge;
                t->SetBranchAddress("truth_LLPChild_Charge", &truth_LLPChild_Charge);
                std::vector<double> * truth_LLPChild_Pt;
                t->SetBranchAddress("truth_LLPChild_Pt", &truth_LLPChild_Pt);
                std::vector<double> * truth_LLPChild_Eta;
                t->SetBranchAddress("truth_LLPChild_Eta", &truth_LLPChild_Eta);
                std::vector<double> * truth_LLPChild_Phi;
                t->SetBranchAddress("truth_LLPChild_Phi", &truth_LLPChild_Phi);
                
                int truthJet_n = (int)truthJet_Pt->size();
                
                // Set weight
                double weight = weight_nomc * mcEventWeight;
                allevents += weight;
                
                // Jet selection
                for(int k = 0; k < truthJet_n; k++){
                    sum_pt += truthJet_Pt->at(k);
                    
                    // p_T >= 137, n >= 4
                    if(truthJet_Pt->at(k) >= 137){
                        Pt137++;
                    }
                    
                    // p_T >= 101, n >= 5
                    if (truthJet_Pt->at(k) >= 101){
                        Pt101++;
                    }
                    
                    // p_T >= 83, n >= 6
                    if (truthJet_Pt->at(k) >= 83){
                        Pt83++;
                    }
                    
                    // p_T >= 55, n >= 7
                    if (truthJet_Pt->at(k) >= 55){
                        Pt55++;
                    }
                    
                    //LLP decay
                    Rxy_max = sqrt(x_LLP->at(0) * x_LLP->at(0) + y_LLP->at(0) * y_LLP->at(0));
                    if(x_LLP->at(k) != -9999){
                        float Rxy = sqrt(x_LLP->at(k) * x_LLP->at(k) + y_LLP->at(k) * y_LLP->at(k));
                        if(Rxy < 3870){
                            if(truthJet_Pt->at(k) >= 70){
                                LLP_Pt70++;
                            }
                            if(truthJet_Pt->at(k) >= 50){
                                LLP_Pt50++;
                            }
                        }
                        
                        if(Rxy > Rxy_max){
                            Rxy_max = Rxy;
                        }
                    }
                }
                
                //get efficiency
                if(Rxy_max < 1150){
                    efficiency = h_1150->GetBinContent(h_1150->FindBin(sum_pt * 0.001));
                }
                else if (Rxy_max < 3870){
                    efficiency = h_1150_3870->GetBinContent(h_1150_3870->FindBin(sum_pt * 0.001));
                }
                else{
                    efficiency = h_3870->GetBinContent(h_3870->FindBin(sum_pt * 0.001));
                }
                
                if(Pt137 < 4 && Pt101 < 5 && Pt83 < 6 && Pt55 < 7){continue;}
                if(LLP_Pt70 < 1 && LLP_Pt50 < 2){continue;}
                n_pass_JetSel += weight;
                pass_JetSel_e = pass_JetSel_e + (weight * efficiency);
                
                //vertex-level acceptance
                for(int k = 0; k < x->size(); k++){
                    // Rxy < 300 mm and |z| < 300 mm
                    float Rxy = sqrt(x->at(k) * x->at(k) + y->at(k) * y->at(k));
                    if(Rxy >= 300 || abs(z->at(k)) >= 300){continue;}
                    n_pass_R++;
                    
                    // Transverse distance from the primary vertex > 4 mm
                    if(Rxy <= 4){continue;}
                    n_pass_Dist++;
                    
                    // At least 1 charged particle with |d0| > 2 mm
                    for(int l = 0; l < truth_LLPChild_Parent_barcode->size(); l++){
                        if(truth_LLPChild_Parent_barcode->at(l) == truthSparticle_Barcode->at(k)){
                            //At least 1 charged particle with |d0| > 2 mm
                            if(truth_LLPChild_Charge->at(l) != 0 && truth_LLPChild_d0->at(l) > 2){
                                n_pass_d0++;
                            }
                            
                            // n_selected decay products ≥ 5
                            if(truth_LLPChild_Charge->at(l) != 0 && truth_LLPChild_Pt->at(l) / abs(truth_LLPChild_Charge->at(l)) > 1){
                                n_pass_decay++;
                            }
                            
                            TLorentzVector p4_track;
                            p4_track.SetPtEtaPhiM(truth_LLPChild_Pt->at(l), truth_LLPChild_Eta->at(l), truth_LLPChild_Phi->at(l), 0.13957);
                            p_tracks += p4_track;
                        }
                    }
                    
                    // Invariant mass > 10 GeV
                    double invariant_mass = p_tracks.M();
                    if(invariant_mass > 10){
                        n_pass_inv_m++;
                    }
                }
                
                if(n_pass_R < 1){continue;}
                pass_R += weight;
                pass_R_e = pass_R_e + (weight * efficiency);
                
                if(n_pass_Dist < 1){continue;}
                pass_Dist += weight;
                pass_Dist_e = pass_Dist_e + (weight * efficiency);
                
                if(n_pass_d0 < 1){continue;}
                pass_d0 += weight;
                pass_d0_e = pass_d0_e + (weight * efficiency);
                
                if(n_pass_decay < 5){continue;}
                pass_decay += weight;
                pass_decay_e = pass_decay_e + (weight * efficiency);
                
                if(n_pass_inv_m < 1){continue;}
                pass_inv_m += weight;
                pass_inv_m_e = pass_inv_m_e + (weight * efficiency);
            }
        }
            
        //print result
        std::cout << std::endl;
        std::cout << "---------- Acceptance ----------" << std::endl;
        std::cout << "Total events: " << allevents << std::endl;
        std::cout << "Events passing jet selection: " << n_pass_JetSel << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << std::left << std::setw(60) << "Jet selection: " << (double)n_pass_JetSel / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Rxy < 300 mm and |z| < 300 mm: " << (double)pass_R / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Transverse distance from the primary vertex > 4 mm: " << (double)pass_Dist / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Have at least 1 charged particle with |d0| > 2 mm: " << (double)pass_d0 / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(62) << "     n_selected decay products ≥ 5: " << (double)pass_decay / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Invariant mass > 10 GeV: " << (double)pass_inv_m / allevents * 100  << "%" << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        
        // with effiency
        std::cout << std::endl;
        std::cout << "---------- Full Simulation ----------" << std::endl;
        std::cout << "Total events: " << allevents << std::endl;
        std::cout << "Events passing jet selection: " << n_pass_JetSel << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << std::left << std::setw(60) << "Jet selection: " << pass_JetSel_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Rxy < 300 mm and |z| < 300 mm: " << pass_R_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Transverse distance from the primary vertex > 4 mm: " << pass_Dist_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Have at least 1 charged particle with |d0| > 2 mm: " << pass_d0_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(62) << "     n_selected decay products ≥ 5: " << pass_decay_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "     Invariant mass > 10 GeV: " << pass_inv_m_e / allevents * 100  << "%" << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
        
        //expected
        // 0.01 ns: 802.57; 0.032 ns: 1376.99; 0.1 ns: 1735.11; 0.32 ns: 757.08; 1 ns: 244.75; 10 ns: 0
//        std::cout << std::endl;
//        std::cout << std::left << std::setw(15) << "   lifetime" << std::setw(5) << "|" << "Expected effiency * acceptance" << std::endl;
//        std::cout << "------------------------------------------------------" << std::endl;
//        std::cout << std::left << std::setw(15) << "   0.01 ns" << std::setw(10) << "|" << 802.57 / (xsec * 1000 * 139) * 100 << "%" << std::endl;
//        std::cout << std::left << std::setw(15) << "   0.032 ns" << std::setw(10) << "|" << 1376.99 / (xsec * 1000 * 139) * 100 << "%" << std::endl;
//        std::cout << std::left << std::setw(15) << "   0.1 ns" << std::setw(10) << "|" << 1735 / (xsec * 1000 * 139) * 100 << "%" << std::endl;
//        std::cout << std::left << std::setw(15) << "   0.32 ns" << std::setw(10) << "|" << 757.08 / (xsec * 1000 * 139) * 100 << "%" << std::endl;
//        std::cout << std::left << std::setw(15) << "   1 ns" << std::setw(10) << "|" << 244 / (xsec * 1000 * 139) * 100 << "%" << std::endl;
//        std::cout << std::left << std::setw(15) << "   10 ns" << std::setw(10) << "|" << 0 / (xsec * 1000 * 139) * 100 << "%" << std::endl;
//        std::cout << std::endl;
        
        std::cout << "Expected effiency * acceptance = " << n_expected_events / (xsec * 1000 * 139) * 100 << "%" << std::endl;
        std::cout << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl << std::endl;
        
        
        simulated.push_back(pass_inv_m_e / allevents * 100);
        expected.push_back(n_expected_events / (xsec * 1000 * 139) * 100);
        
        // close file
        file_3870->Close();
        file_1150->Close();
        file_1150_3870->Close();
    }
    
    sort(lifetime_list, dsid_list, simulated, expected);
    
    std::cout << std::left << std::setw(10) << "lifetime" << std::setw(10) << "dsid" << std::setw(15) << "simulated" << std::setw(15) << "expected" << std::setw(15) << "percent error" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    for(int i = 0; i < lifetime_list.size(); i++){
        double error = abs(simulated.at(i) - expected.at(i)) / expected.at(i) * 100;
        std::cout << std::left << std::setw(10) << lifetime_list.at(i) << std::setw(10) << dsid_list.at(i) << std::setw(15) << simulated.at(i) << std::setw(15) << expected.at(i) << std::setw(15) << error << std::endl;
    }
    
    std::cout << std::endl;
    return 0;
}
