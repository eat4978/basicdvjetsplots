#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1.h"
#include "atlasstyle-00-04-02/AtlasStyle.h"

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
        if(l.find("group.phys-susy." + dsid) != std::string::npos){
            files.push_back(l);
        }
    }
    
    fin.close();
    return files;
}

// Get sample information
std::string get_sample_info(std::string dsid, bool &stop){
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
            
            fin.close();
            
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



int main(int argc, const char * argv[]) {
    //variables
    std::string tag;
    bool stop = false;
    // counters
    int allevents = 0;
    int passevents = 0;
                      
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
    std::string dsid = "501166";
    
//     //Load Atlas Style
//    gROOT->LoadMacro("atlasstyle-00-04-02/AtlasStyle.C");
//    gROOT->LoadMacro("atlasstyle-00-04-02/AtlasUtils.C");
//    gROOT->LoadMacro("atlasstyle-00-04-02/AtlasLabels.C");
//    SetAtlasStyle();
    
    // Gather inputs & cross-section
    std::vector<std::string> inputs = get_files(dsid, stop);
    if(stop){ return 1;}
    std::string info = get_sample_info(dsid, stop);
    if(stop){ return 1;}
    std::cout << "Analyzing DSID " << dsid << " with sample info = " << info << std::endl;
    double xsec = get_cross_section(info);
    
    // Initialize histograms
    TH1F* h_mDV_all = new TH1F("mDV_all", "mDV_all", 100, 0., 500.);
    TH1F* h_mDV_passPresel = new TH1F("mDV_passPresel", "mDV_passPresel", 100, 0., 500.);
    TH1F* h_mDV_passSel = new TH1F("mDV_passSel", "mDV_passSel", 100, 0., 500.);
    
    // Load tree and loop through events in ntuple
    for(int i  = 0; i < inputs.size(); i++){
//        TFile* inFile = TFile::Open(inputs.at(i).c_str(), "READ");
        TFile* inFile = TFile::Open("/Users/catherine/Documents/URAP/group.phys-susy.501166.e8342_e7400_s3338_r11891_r11748_p4302.28466712._001295.trees.root", "READ");

        // Load Tree
        TTree* t = (TTree*)inFile->Get("trees_SRDV_");
        long nentries = t->GetEntries();
        
        std::cout << "Total number of entries in file: " << nentries << std::endl;
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
            std::cout << "Run number, " << runNumber << " not recognized. Setting lumi=1" << std::endl;
        }
        std::cout << t->GetEntry(0) <<std::endl;

        // Calculate event weight
        double weight_nomc = lumi * xsec / total_sum_weights;
        
        // Loop through entries in TTree
        for(int j = 0; j < nentries; j++){
            float mcEventWeight;
            UChar_t DRAW_pass_triggerFlags;
            UChar_t DRAW_pass_DVJETS;
            UChar_t BaselineSel_HighPtSR;
//            bool BaselineSel_TracklessSR_;
            int DV_n;
            t->SetBranchAddress("mcEventWeight", &mcEventWeight);
            t->SetBranchAddress("DRAW_pass_triggerFlags", &DRAW_pass_triggerFlags);
            t->SetBranchAddress("DRAW_pass_DVJETS", &DRAW_pass_DVJETS);
            t->SetBranchAddress("BaselineSel_HighPtSR", &BaselineSel_HighPtSR);
//            t->SetBranchAddress("BaselineSel_TracklessSR", &BaselineSel_TracklessSR_fromBranch);
            t->SetBranchAddress("DV_n", &DV_n);
            std::vector<int> * DV_m;
            std::vector<bool> * DV_passFiducialCut;
            std::vector<bool> * DV_passDistCut;
            std::vector<bool> * DV_passChiSqCut;
            std::vector<bool> * DV_passMaterialVeto_strict;
            std::vector<bool> * DV_passNTrkCut;
            std::vector<bool> * DV_passMassCut;
            std::vector<int> * DV_nTracksSel;
            t->SetBranchAddress("DV_m", &DV_m);
            t->SetBranchAddress("DV_passFiducialCut", &DV_passFiducialCut);
            t->SetBranchAddress("DV_passDistCut", &DV_passDistCut);
            t->SetBranchAddress("DV_passChiSqCut", &DV_passChiSqCut);
            t->SetBranchAddress("DV_passMaterialVeto_strict", &DV_passMaterialVeto_strict);
            t->SetBranchAddress("DV_passNTrkCut", &DV_passNTrkCut);
            t->SetBranchAddress("DV_passMassCut", &DV_passMassCut);
            t->SetBranchAddress("DV_nTracksSel", &DV_nTracksSel);
            t->GetEntry(j);
            
            // Set weight
            double weight = weight_nomc * mcEventWeight;
            allevents += weight;

            // Apply some selection
            if (DRAW_pass_triggerFlags == 0){
                continue; // Trigger
            }
            if (DRAW_pass_DVJETS == 0){
                continue; // DRAW filter
            }
            if (BaselineSel_HighPtSR == 0){
                continue; // Jet selection in High-pT jet SR
            }
//            if (BaselineSel_TracklessSR == 0){
//                continue; // Jet selection in Trackless jet SR
//            }
            
            // Loop through displaced vertices (DV) in event
            int nDV_passSel = 0;
            for (int k = 0; k < DV_n; k++){

                // Fill histogram
                h_mDV_all->Fill(DV_m->at(k), weight);

                // Preselection
                if (!(DV_passFiducialCut->at(k) && DV_passDistCut->at(k) && DV_passChiSqCut->at(k) && DV_passMaterialVeto_strict->at(k))){
                    continue;
                }
                h_mDV_passPresel->Fill(DV_m->at(k), weight);

                // Full selection
                if (!(DV_passNTrkCut->at(k) && DV_passMassCut->at(k) && DV_nTracksSel->at(k) >= 2)){
                    continue;
                }
                h_mDV_passSel->Fill(DV_m->at(k), weight);
                nDV_passSel++;
            }

            // Remove events that don't have at least 1 DV passing selection
            if (nDV_passSel == 0){
                continue;
            }

            // Count events passing selection
            passevents += weight;
        }
        
    }
    
    std::cout << "Total events: " << allevents << std::endl;
    std::cout << "Events passing selection: " << allevents << std::endl;
    std::cout << "Efficiency = " << passevents/allevents * 100. << " %" << std::endl;
    
    // Save histograms
    TFile* fout = new TFile(("out_" + tag + ".root").c_str(), "RECREATE");
    fout->cd();
    h_mDV_all->Write();
    h_mDV_passPresel->Write();
    h_mDV_passSel->Write();
    fout->Close();
    std::cout << "Writing out_" << tag << ".root" << std::endl;
    
    // Plot histograms
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    h_mDV_all->SetLineWidth(2);
    h_mDV_passPresel->SetLineWidth(2);
    h_mDV_passSel->SetLineWidth(2);
    h_mDV_all->SetLineColor(kBlack);
    h_mDV_passPresel->SetLineColor(kBlue);
    h_mDV_passSel->SetLineColor(kRed);
    h_mDV_all->GetYaxis()->SetTitle("Number of DVs");
    h_mDV_all->GetXaxis()->SetTitle("DV mass [GeV]");
    h_mDV_all->Draw("HIST");
    h_mDV_passPresel->Draw("HIST SAME");
    h_mDV_passSel->Draw("HIST SAME");
    c1->SetLogy(1);

    TLegend* legend = new TLegend(0.4, 0.7, 0.9, 0.9);
    legend->AddEntry(h_mDV_all, "All DVs", "l");
    legend->AddEntry(h_mDV_passPresel, "All DVs passing preselection", "l");
    legend->AddEntry(h_mDV_passSel, "All DVs passing full selection", "l");
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->Draw();

    // Save the plot as a PNG file
    c1->SaveAs(("dvmass_" + tag + ".png").c_str());
    std::cout << "Saved dvmass_" << tag << ".png" << std::endl;
    
    delete legend;
    delete c1;
    delete fout;
    delete h_mDV_passSel;
    delete h_mDV_passPresel;
    delete h_mDV_all;
    
    return 0;
}
