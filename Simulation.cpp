#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdlib.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include <TLorentzVector.h>

struct File_Info{
    std::string dsid;
    double lifetime;
    double mass;
    double sim_a;
    double sim_full;
    double exp_full;
};

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
std::string get_sample_info(std::string dsid, bool &stop, double &lifetime, int &mass){
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
            
            //get mass
            unsigned long secondLast = tmp3.rfind('_', first - 1);
            std::string m = tmp3.substr(secondLast + 1, first - secondLast - 1);
            mass = std::atoi(m.c_str());

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
        else if (sampleinfo.find("_300_") != std::string::npos){
            cross_x = 0.3003;
        }
        else if (sampleinfo.find("_500_") != std::string::npos){
            cross_x = 0.0361;
        }
        else if (sampleinfo.find("_700_") != std::string::npos){
            cross_x = 0.00738;
        }
        else if (sampleinfo.find("_900_") != std::string::npos){
            cross_x = 0.00192;
        }
        else if (sampleinfo.find("_1100_") != std::string::npos){
            cross_x = 0.000571;
        }
        else if (sampleinfo.find("_1300_") != std::string::npos){
            cross_x = 0.000184;
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

// Get total expected events
void get_n_expected_events(double &events, TTree* t, long nentries, double weight_nomc){
    float mcEventWeight;
    t->SetBranchAddress("mcEventWeight", &mcEventWeight);
    std::vector<double> * truthJet_Pt;
    t->SetBranchAddress("truthJet_Pt", &truthJet_Pt);
    
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
    
    for(int j = 0; j < nentries; j++){
        t->GetEntry(j);
        
        int truthJet_n = (int)truthJet_Pt->size();
        double weight = weight_nomc * mcEventWeight;
        int nDV_passSel = 0;
        
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
    std::string user_input;
    std::vector<File_Info> file_list;
    
    // get DSID
    while(true){
        std::cout << "Enter DSID (enter q to quit): ";
        std::cin >> user_input;
        if(user_input == "q" || user_input == "Q"){
            break;
        }
        File_Info newFile;
        newFile.dsid = user_input;
        file_list.push_back(newFile);
    }
    
    for(int file = 0; file < file_list.size(); file++){
        //variables
        std::string tag;
        bool stop = false;
        double lifetime;
        int mass;
        
        // counters
        double allevents = 0;
        double n_pass_JetSel = 0;
        double pass_R = 0;
        double pass_Dist = 0;
        double pass_d0 = 0;
        double pass_decay = 0;
        double pass_inv_m = 0;
        double n_expected_events = 0;
        
        double pass_JetSel_e = 0;
        double pass_R_e = 0;
        double pass_Dist_e = 0;
        double pass_d0_e = 0;
        double pass_decay_e = 0;
        double pass_inv_m_e = 0;
               
        //event effiency files
        TFile* file_3870 = new TFile("HEPData-ins2628398-v1-event_efficiency_Trackless_R_3870_mm.root");
        TFile* file_1150 = new TFile("HEPData-ins2628398-v1-event_efficiency_Trackless_R_1150_mm.root");
        TFile* file_1150_3870 = new TFile("HEPData-ins2628398-v1-event_efficiency_Trackless_R_1150_3870_mm.root");
        
        TH1D* h_3870 = (TH1D*)file_3870->Get("event_efficiency_Trackless_R_3870_mm/Hist1D_y1");
        TH1D* h_1150 = (TH1D*)file_1150->Get("event_efficiency_Trackless_R_1150_mm/Hist1D_y1");
        TH1D* h_1150_3870 = (TH1D*)file_1150_3870->Get("event_efficiency_Trackless_R_1150_3870_mm/Hist1D_y1");
        
        //vertex effiency files
        TFile* file_eff_vt_22 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_22_mm.root");
        TFile* file_eff_vt_22_25 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_22_25_mm.root");
        TFile* file_eff_vt_25_29 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_25_29_mm.root");
        TFile* file_eff_vt_29_38 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_29_38_mm.root");
        TFile* file_eff_vt_38_46 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_38_46_mm.root");
        TFile* file_eff_vt_46_73 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_46_73_mm.root");
        TFile* file_eff_vt_73_84 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_73_84_mm.root");
        TFile* file_eff_vt_84_111 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_84_111_mm.root");
        TFile* file_eff_vt_111_120 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_111_120_mm.root");
        TFile* file_eff_vt_120_145 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_120_145_mm.root");
        TFile* file_eff_vt_145_180 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_145_180_mm.root");
        TFile* file_eff_vt_180_300 = new TFile("HEPData-ins2628398-v1-vertex_efficiency_R_180_300_mm.root");
        
        TH2D* eff_vt_22 = (TH2D*)file_eff_vt_22->Get("vertex_efficiency_R_22_mm/Hist2D_y1");
        TH2D* eff_vt_22_25 = (TH2D*)file_eff_vt_22_25->Get("vertex_efficiency_R_22_25_mm/Hist2D_y1");
        TH2D* eff_vt_25_29 = (TH2D*)file_eff_vt_25_29->Get("vertex_efficiency_R_25_29_mm/Hist2D_y1");
        TH2D* eff_vt_29_38 = (TH2D*)file_eff_vt_29_38->Get("vertex_efficiency_R_29_38_mm/Hist2D_y1");
        TH2D* eff_vt_38_46 = (TH2D*)file_eff_vt_38_46->Get("vertex_efficiency_R_38_46_mm/Hist2D_y1");
        TH2D* eff_vt_46_73 = (TH2D*)file_eff_vt_46_73->Get("vertex_efficiency_R_46_73_mm/Hist2D_y1");
        TH2D* eff_vt_73_84 = (TH2D*)file_eff_vt_73_84->Get("vertex_efficiency_R_73_84_mm/Hist2D_y1");
        TH2D* eff_vt_84_111 = (TH2D*)file_eff_vt_84_111->Get("vertex_efficiency_R_84_111_mm/Hist2D_y1");
        TH2D* eff_vt_111_120 = (TH2D*)file_eff_vt_111_120->Get("vertex_efficiency_R_111_120_mm/Hist2D_y1");
        TH2D* eff_vt_120_145 = (TH2D*)file_eff_vt_120_145->Get("vertex_efficiency_R_120_145_mm/Hist2D_y1");
        TH2D* eff_vt_145_180 = (TH2D*)file_eff_vt_145_180->Get("vertex_efficiency_R_145_180_mm/Hist2D_y1");
        TH2D* eff_vt_180_300 = (TH2D*)file_eff_vt_180_300->Get("vertex_efficiency_R_180_300_mm/Hist2D_y1");

                          
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
        
        
        std::string dsid = file_list.at(file).dsid;
        
        // Gather inputs & cross-section
        std::vector<std::string> inputs = get_files(dsid, stop);
        if(stop){ return 1;}
        std::string info = get_sample_info(dsid, stop, lifetime, mass);
        if(stop){ return 1;}
        file_list.at(file).mass = mass;
        file_list.at(file).lifetime = lifetime;
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
            
            float mcEventWeight;
            t->SetBranchAddress("mcEventWeight", &mcEventWeight);
            std::vector<double> * truthJet_Pt;
            t->SetBranchAddress("truthJet_Pt", &truthJet_Pt);

            int DV_n;
            t->SetBranchAddress("DV_n", &DV_n);

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
            
            
            /*****************************************************************************************************
             *****************************************************************************************************
             *                                                  Expected                                                                                             *
             *****************************************************************************************************
             *****************************************************************************************************
             */
            
            // Get total expected events
            for(int j = 0; j < nentries; j++){
                t->GetEntry(j);
                
//            int truthJet_n = (int)truthJet_Pt->size();
                double weight = weight_nomc * mcEventWeight;
                int nDV_passSel = 0;
                
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
            
            
            /*****************************************************************************************************
             *****************************************************************************************************
             *                                                  Simulated                                                                                             *
             *****************************************************************************************************
             *****************************************************************************************************
             */
            
            //variables
//                float mcEventWeight;
//                t->SetBranchAddress("mcEventWeight", &mcEventWeight);
//                std::vector<double> * truthJet_Pt;
//                t->SetBranchAddress("truthJet_Pt", &truthJet_Pt);
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
            std::vector<double> * truth_LLPChild_Status;
            t->SetBranchAddress("truth_LLPChild_Status", &truth_LLPChild_Status);
            
            // Loop through entries in TTree
            for(int j = 0; j < nentries; j++){
                t->GetEntry(j);     // ERROR
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
                
                int truthJet_n = (int)truthJet_Pt->size();
                
                // Set weight
                double weight = weight_nomc * mcEventWeight;
                allevents += weight;
                
                // Rxy_max = sqrt(x_LLP->at(0) * x_LLP->at(0) + y_LLP->at(0) * y_LLP->at(0));
                
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
                    // Rxy_max = sqrt(x_LLP->at(0) * x_LLP->at(0) + y_LLP->at(0) * y_LLP->at(0));
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
                float Aeff_vertex = 1.0; // because we compute 1 - mult( 1 - A*eff )
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
                            
                            if(truth_LLPChild_Charge->at(l) != 0 && truth_LLPChild_Pt->at(l) / abs(truth_LLPChild_Charge->at(l)) > 1 && truth_LLPChild_Status->at(l) > 0){
                                n_pass_decay++;
                            }
                            
                            TLorentzVector p4_track;
                            p4_track.SetPtEtaPhiM(truth_LLPChild_Pt->at(l), truth_LLPChild_Eta->at(l), truth_LLPChild_Phi->at(l), 0.13957);
                            p_tracks += p4_track;
                        }
                    }
                    
                    // Invariant mass > 10 GeV
                    double invariant_mass = p_tracks.M();
                    if(invariant_mass < 10) continue;
                    
                    n_pass_inv_m++;
                    
                    // get vertex effiency
                    float eff_vertex;
                    
                    if(Rxy < 22){
                        eff_vertex = eff_vt_22->GetBinContent(eff_vt_22->GetXaxis()->FindBin(invariant_mass), eff_vt_22->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 25){
                        eff_vertex = eff_vt_22_25->GetBinContent(eff_vt_22_25->GetXaxis()->FindBin(invariant_mass), eff_vt_22_25->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 29){
                        eff_vertex = eff_vt_25_29->GetBinContent(eff_vt_25_29->GetXaxis()->FindBin(invariant_mass), eff_vt_25_29->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 38){
                        eff_vertex = eff_vt_29_38->GetBinContent(eff_vt_29_38->GetXaxis()->FindBin(invariant_mass), eff_vt_29_38->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 46){
                        eff_vertex = eff_vt_38_46->GetBinContent(eff_vt_38_46->GetXaxis()->FindBin(invariant_mass), eff_vt_38_46->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 73){
                        eff_vertex = eff_vt_46_73->GetBinContent(eff_vt_46_73->GetXaxis()->FindBin(invariant_mass), eff_vt_46_73->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 84){
                        eff_vertex = eff_vt_73_84->GetBinContent(eff_vt_73_84->GetXaxis()->FindBin(invariant_mass), eff_vt_73_84->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 111){
                        eff_vertex = eff_vt_84_111->GetBinContent(eff_vt_84_111->GetXaxis()->FindBin(invariant_mass), eff_vt_84_111->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 120){
                        eff_vertex = eff_vt_111_120->GetBinContent(eff_vt_111_120->GetXaxis()->FindBin(invariant_mass), eff_vt_111_120->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 145){
                        eff_vertex = eff_vt_120_145->GetBinContent(eff_vt_120_145->GetXaxis()->FindBin(invariant_mass), eff_vt_120_145->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else if (Rxy < 180){
                        eff_vertex = eff_vt_145_180->GetBinContent(eff_vt_145_180->GetXaxis()->FindBin(invariant_mass), eff_vt_145_180->GetYaxis()->FindBin(n_pass_decay));
                    }
                    else{
                        eff_vertex = eff_vt_180_300->GetBinContent(eff_vt_180_300->GetXaxis()->FindBin(invariant_mass), eff_vt_180_300->GetYaxis()->FindBin(n_pass_decay));
                    }
                    
                    Aeff_vertex = Aeff_vertex * ( 1 - eff_vertex);
                }
                
                Aeff_vertex = 1 - Aeff_vertex;
                
                if(n_pass_R < 1){continue;}
                pass_R += weight;
                pass_R_e = pass_R_e + (weight * efficiency * Aeff_vertex);
                
                if(n_pass_Dist < 1){continue;}
                pass_Dist += weight;
                pass_Dist_e = pass_Dist_e + (weight * efficiency * Aeff_vertex);
                
                if(n_pass_d0 < 1){continue;}
                pass_d0 += weight;
                pass_d0_e = pass_d0_e + (weight * efficiency * Aeff_vertex);
                
                if(n_pass_decay < 5){continue;}
                pass_decay += weight;
                pass_decay_e = pass_decay_e + (weight * efficiency * Aeff_vertex);
                
                if(n_pass_inv_m < 1){continue;}
                pass_inv_m += weight;
                pass_inv_m_e = pass_inv_m_e + (weight * efficiency * Aeff_vertex);
            }
        }
            
        //print result
        std::cout << std::endl;
        std::cout << "---------- Acceptance ----------" << std::endl;
        std::cout << "Total events: " << allevents << std::endl;
        std::cout << "Events passing jet selection: " << n_pass_JetSel << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << std::left << std::setw(60) << "Jet selection: " << n_pass_JetSel / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Rxy < 300 mm and |z| < 300 mm: " << pass_R / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Transverse distance from the primary vertex > 4 mm: " << pass_Dist / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Have at least 1 charged particle with |d0| > 2 mm: " << pass_d0 / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(62) << "n_selected decay products ≥ 5: " << pass_decay / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Invariant mass > 10 GeV: " << pass_inv_m / allevents * 100  << "%" << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        
        // with effiency
        std::cout << std::endl;
        std::cout << "---------- Full Simulation ----------" << std::endl;
        std::cout << "Total events: " << allevents << std::endl;
        std::cout << "Events passing jet selection: " << n_pass_JetSel << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << std::left << std::setw(60) << "Jet selection: " << pass_JetSel_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Rxy < 300 mm and |z| < 300 mm: " << pass_R_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Transverse distance from the primary vertex > 4 mm: " << pass_Dist_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Have at least 1 charged particle with |d0| > 2 mm: " << pass_d0_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(62) << "n_selected decay products ≥ 5: " << pass_decay_e / allevents * 100  << "%" << std::endl;
        std::cout << std::left << std::setw(60) << "Invariant mass > 10 GeV: " << pass_inv_m_e / allevents * 100  << "%" << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;

        
        std::cout << "Expected effiency * acceptance = " << n_expected_events / (xsec * 1000 * 139) * 100 << "%" << std::endl;
        std::cout << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl << std::endl;
        
        
        file_list.at(file).sim_a = pass_inv_m / allevents * 100;
        file_list.at(file).sim_full = pass_inv_m_e / allevents * 100;
        file_list.at(file).exp_full = n_expected_events / (xsec * 1000 * 139) * 100;

        
        // close file
        file_3870->Close();
        file_1150->Close();
        file_1150_3870->Close();
        
        file_eff_vt_22->Close();
        file_eff_vt_22_25->Close();
        file_eff_vt_25_29->Close();
        file_eff_vt_29_38->Close();
        file_eff_vt_38_46->Close();
        file_eff_vt_46_73->Close();
        file_eff_vt_73_84->Close();
        file_eff_vt_84_111->Close();
        file_eff_vt_111_120->Close();
        file_eff_vt_120_145->Close();
        file_eff_vt_145_180->Close();
        file_eff_vt_180_300->Close();
        
        delete file_3870;
        delete file_1150;
        delete file_1150_3870;
                
        delete file_eff_vt_22;
        delete file_eff_vt_22_25;
        delete file_eff_vt_25_29;
        delete file_eff_vt_29_38;
        delete file_eff_vt_38_46;
        delete file_eff_vt_46_73;
        delete file_eff_vt_73_84;
        delete file_eff_vt_84_111;
        delete file_eff_vt_111_120;
        delete file_eff_vt_120_145;
        delete file_eff_vt_145_180;
        delete file_eff_vt_180_300;
    }
    
    //sort
    auto compare = [](const File_Info& a, const File_Info& b) {
        if (a.mass == b.mass) {
            // If mass is the same, compare by lifetime
            return a.lifetime < b.lifetime;
        }
        // Otherwise, compare by mass
        return a.mass < b.mass;
    };

    std::sort(file_list.begin(), file_list.end(), compare);
    
    std::cout << std::left << std::setw(9) << "mass" << std::setw(12) << "lifetime" << std::setw(10) << "dsid" << std::setw(15) << "simulated" << std::setw(15) << "expected" << std::setw(15) << "percent error" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    for(int i = 0; i < file_list.size(); i++){
        double error = abs(file_list.at(i).sim_full - file_list.at(i).exp_full) / file_list.at(i).exp_full * 100;
        std::cout << std::left << std::setw(9) << file_list.at(i).mass << std::setw(12) << file_list.at(i).lifetime << std::setw(10) << file_list.at(i).dsid << std::setw(15) << file_list.at(i).sim_full << std::setw(15) << file_list.at(i).exp_full << std::setw(15) << error << std::endl;
    }
    
    std::cout << std::endl;
    return 0;
}
