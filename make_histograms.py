import ROOT
import sys
from helper_functions import *

if len(sys.argv)<2:
  print("Please specify dsid!")
  exit(0)

dsid = sys.argv[1]
if (len(sys.argv)>2):
  tag = sys.argv[2]
else: tag = "v1"

# Load atlas style so plots look nice
ROOT.gROOT.LoadMacro("atlasstyle-00-04-02/AtlasStyle.C")
ROOT.gROOT.LoadMacro("atlasstyle-00-04-02/AtlasUtils.C")
ROOT.gROOT.LoadMacro("atlasstyle-00-04-02/AtlasLabels.C")
ROOT.SetAtlasStyle()

# Gather inputs and cross-section
inputs = get_files(dsid)
info = get_sample_info(dsid)
print("Analyzing DSID ", dsid, " with sample info = ", info)
xsec = get_cross_section(info)

# Set counters
allevents = 0
passevents = 0

# Initialize histograms       
h_mDV_all = ROOT.TH1F("mDV_all", "mDV_all", 100, 0., 500.)
h_mDV_passPresel = ROOT.TH1F("mDV_passPresel", "mDV_passPresel", 100, 0., 500.)
h_mDV_passSel = ROOT.TH1F("mDV_passSel", "mDV_passSel", 100, 0., 500.)

# Load tree and loop through events in ntuple
for f in inputs:
   inFile = ROOT.TFile(f, 'READ')

   # Load TTree
   t = inFile.Get("trees_SRDV_")
   nentries = t.GetEntries()

   print("Total number of entries in file: ", nentries)
   # Get weight information
   metadata = inFile.Get("MetaData_EventCount")
   total_sum_weights = metadata.GetBinContent(3)

   # Determine luminosity of this sample, pb^-1
   t.GetEntry(0)
   lumi = 1.0
   if (t.runNumber==284500): lumi=36200. #mc16a (2015/16)
   elif (t.runNumber==300000): lumi=44300. #mc16d (2017)
   elif (t.runNumber==310000): lumi=58500. #mc16e (2018)
   else: print("Run number, ", t.runNumber, " not recognized. Setting lumi=1")

   # Calculate event weight
   weight_nomc = lumi * xsec / total_sum_weights

   # Loop through entries in TTree
   for i in range(0, nentries):
      t.GetEntry(i)

      # Set weight
      weight = weight_nomc * t.mcEventWeight 
      allevents += weight

      # Apply some selection 
      if (t.DRAW_pass_triggerFlags==0): continue # Trigger
      if (t.DRAW_pass_DVJETS==0): continue # DRAW filter
      if (t.BaselineSel_HighPtSR==0): continue # Jet selection in High-pT jet SR
      #if (t.BaselineSel_TracklessSR==0): continue # Jet selection in Trackless jet SR
     
      # Loop through displaced vertices (DV) in event
      nDV_passSel = 0
      for j in range(0, t.DV_n):

         # Fill histogram
         h_mDV_all.Fill(t.DV_m[j], weight)

         # Preselection
         if not (t.DV_passFiducialCut[j] and t.DV_passDistCut[j] and t.DV_passChiSqCut[j] and t.DV_passMaterialVeto_strict[j]): continue
         h_mDV_passPresel.Fill(t.DV_m[j], weight)

         # Full selection
         if not (t.DV_passNTrkCut[j] and t.DV_passMassCut[j] and t.DV_nTracksSel[j]>=2): continue
         h_mDV_passSel.Fill(t.DV_m[j], weight)
        
         nDV_passSel += 1

      # Remove events that don't have at least 1 DV passing selection
      if (nDV_passSel==0): continue
 
      # Count events passing selection
      passevents += weight

print("Total events: ", allevents)
print("Events passing selection: ", passevents)
print("Efficiency = ", passevents/allevents*100., " % ") 

# Save histograms
fout = ROOT.TFile("outputs/out_"+tag+".root", "RECREATE")
h_mDV_all.Write()
h_mDV_passPresel.Write()
h_mDV_passSel.Write()
fout.Close()
print("Writing outputs/out_"+tag+".root")

# Plot histograms
c1 = ROOT.TCanvas("c1", "c1", 800, 600)
h_mDV_all.SetLineWidth(2); h_mDV_passPresel.SetLineWidth(2); h_mDV_passSel.SetLineWidth(2);
h_mDV_all.SetLineColor(ROOT.kBlack); h_mDV_passPresel.SetLineColor(ROOT.kBlue); h_mDV_passSel.SetLineColor(ROOT.kRed);
h_mDV_all.GetYaxis().SetTitle("Number of DVs")
h_mDV_all.GetXaxis().SetTitle("DV mass [GeV]")
h_mDV_all.Draw("HIST")
h_mDV_passPresel.Draw("HIST SAME")
h_mDV_passSel.Draw("HIST SAME")
c1.SetLogy(1)

legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9);
legend.AddEntry(h_mDV_all, "All DVs", "l")
legend.AddEntry(h_mDV_passPresel, "All DVs passing preselection", "l")
legend.AddEntry(h_mDV_passSel, "All DVs passing full selection", "l")
legend.SetTextSize(.04)
legend.SetBorderSize(0)
legend.Draw()

c1.SaveAs("plots/dvmass_"+tag+".png")
print("Saved plots/dvmass_"+tag+".png")
