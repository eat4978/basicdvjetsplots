import re

# Get sample files
def get_files(dsid):
   files = []
   fin = open('allntuples.txt', 'r')
   for l in fin:
      if "group.phys-susy."+str(dsid) in l:
         files += [l.strip()]
   return files

# Get sample information
def get_sample_info(dsid):
   fin = open('signal_key.txt', 'r')
   for l in fin:
      if "mc16_13TeV."+str(dsid) in l:
         # Extract info we want from string
         tmp = l.split('.')[2]
         tmp = re.sub(r'^.*?_', '', tmp)
         tmp = re.sub(r'^.*?_', '', tmp)
         return tmp
   print("Sample with DSID=", dsid, " not found in signal_key.txt")
   return 0

# Get cross section in pb, using mass of SUSY particles
def get_cross_section(sampleinfo):
   cross_x = 1
   if ("C1N2N1" in sampleinfo):
      if ("_100_" in sampleinfo): cross_x = 17.125
      elif ("_300_" in sampleinfo): cross_x = 0.3003
      elif ("_500_" in sampleinfo): cross_x = 0.03605
      elif ("_700_" in sampleinfo): cross_x = 0.00738
      elif ("_900_" in sampleinfo): cross_x = 0.0019200
      elif ("_1100_" in sampleinfo): cross_x = 0.000570676
      elif ("_1300_" in sampleinfo): cross_x = 0.00018362
      elif ("_1500_" in sampleinfo): cross_x = 6.204689e-05
      elif ("_1600_" in sampleinfo): cross_x = 3.615468e-05
      elif ("_1700_" in sampleinfo): cross_x = 2.14739e-05
   
   if (cross_x == 1.0):
      print("Cross section not implemented for this sample (", sampleinfo, ". Please update get_cross_section() function!")
   return cross_x
