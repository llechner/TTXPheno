''' simple draw script
'''

# Standrard imports
import ROOT
import os

# TTXPheno
from TTXPheno.Tools.user import plot_directory

# Sample
from TTXPheno.samples.gen_fwlite_benchmarks import dim6top_ttZ_ll_LO_currentplane_highStat_scan 

# just 1 file
dim6top_ttZ_ll_LO_currentplane_highStat_scan.files = dim6top_ttZ_ll_LO_currentplane_highStat_scan.files[:1]

# get TChain
c = dim6top_ttZ_ll_LO_currentplane_highStat_scan.chain 

c.Draw("Z_pt>>h_Zpt(50,0,400)") # "weight*(%s)" % weightString(cpt=0.2)

c1 = ROOT.TCanvas()
ROOT.h_Zpt.Draw('hist')

c1.Print(os.path.join(plot_directory, 'myplot.png'))

