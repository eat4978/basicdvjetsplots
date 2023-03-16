# basicdvjetsplots

This repository has instructions for making basic plots of the DV+jets ntuples on the NERSC.

## Environment setup

1. SSH into NERSC computer. See instructions here: https://docs.nersc.gov/connect/ and https://docs.nersc.gov/connect/mfa/#sshproxy

2. Setup ATLAS environment. See instructions here: https://atlaswiki.lbl.gov/software/atlascori

3. Setup root. Type 'lsetup root' to see latest version. 

## Ntuple locations
Information about simulated pp collisions with SUSY particles are stored in ntuples, which are root files. The path to all ntuples we used in the DV+jets analysis on the NERSC are listed in `allntuples.txt`. An examples is:
```
/global/cfs/cdirs/atlas/projecta/atlas/atlaslocalgroupdisk/rucio/group/phys-susy/45/47/group.phys-susy.501166.e8342_e7400_s3338_r11891_r11748_p4302.28466712._001295.trees.root
```
The 6-digit number after 'group.phys-susy.' is the dataset identification number (DSID). The single letter followed by numbers (i.e. e8342) are called tags, and these give us information on how the simulated sample was produced. We can use the DSID and the `signal_key.txt` file to get information on which SUSY particles are simulated in the pp-collisions. In the `signal_key.txt` file, search for a line that contains the DSID of interest. For example, a line with DSID 501166 reads: 
```
mc16_13TeV:mc16_13TeV.501166.MGPy8EG_A14NNPDF23LO_C1N2N1_rpvLF_1500_p1ns.deriv.DAOD_SUSY15.e8342_e7400_s3338_r11891_r11748_p4302
``` 
and the snippet `C1N2N1_rpvLF_1500_p1ns` tells us what in this ntuple:
- `C1N2N1` means that SUSY particles in the simplified EWK RPV model are being simulated
- `rpvLF` means that the SUSY particles decay to light-flavored quarks
- `1500` is the mass of the long-lived SUSY particle in GeV
- `p1ns` means that the lifetime of the long-lived SUSY particle is 0.1 ns

## Making histograms
Run `make_histograms.py`, which takes as arguments a DSID (required) and a tag (optional). Example:
```
python make_histograms.py 501167 v2
```
This script will make a histogram of the reconstructed mass of displaced vertices before and after applying some selection.
