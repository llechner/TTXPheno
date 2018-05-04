''' simple draw script
'''

# Standrard imports
import ROOT
import os

# TTXPheno
from TTXPheno.Tools.user import plot_directory
from TTXPheno.Tools.WeightInfo import WeightInfo

# Sample
from TTXPheno.samples.benchmarks import * # dim6top_ttZ_ll_LO_currentplane_highStat_scan 

#lumi = 77 #36+41
#sigmaxBR = 
#Nsim = 5000

# sample = test
sample = dim6top_ttZ_ll_LO_highStat_scan
#sample = dim6top_ttZ_ll_LO_currentplane_highStat_scan

# just 1 file
sample.files = sample.files[:1]

# get TChain
c = sample.chain 

w = WeightInfo(sample.reweight_pkl)
w.set_order( 2 )

c.Draw("Z_pt>>h_Zpt(50,0,500)") # "weight*(%s)" % weightString(cpt=0.2)
c.Draw("Z_pt>>h_Zpt1(50,0,500)", '(' + w.arg_weight_string(cpt=0, cpQM=0, ctZI=-1, ctZ=0) + ')/p_C[0]')

# Remove StatsBox
ROOT.h_Zpt.SetStats(False)
ROOT.h_Zpt1.SetStats(False)

# Line Color
ROOT.h_Zpt.SetLineColor(ROOT.kBlue)
ROOT.h_Zpt1.SetLineColor(ROOT.kRed)

legend = ROOT.TLegend(0.65, 0.75, .9, .9)
#legend->SetNColumns(2)
legend.AddEntry(ROOT.h_Zpt, "p_T(Z) SM", "l")
legend.AddEntry(ROOT.h_Zpt1, "p_T(Z) (cpt, cpQM, ctZI, ctZ)", "l")

h_Zpt_max = ROOT.h_Zpt.GetMaximum()
h_Zpt1_max = ROOT.h_Zpt1.GetMaximum()

if h_Zpt_max < h_Zpt1_max:
    ROOT.h_Zpt1, ROOT.h_Zpt = ROOT.h_Zpt, ROOT.h_Zpt1

# Base Line
line = ROOT.TLine(300,0,300,ROOT.h_Zpt.GetMaximum())
line.SetLineColor(ROOT.kBlack)
line.SetLineWidth(1)

# Plotting
c1 = ROOT.TCanvas()
ROOT.h_Zpt.Draw('HIST')
ROOT.h_Zpt1.Draw('HIST SAME')

#ROOT.gStyle.SetTitleY(1.4)
ROOT.h_Zpt.SetXTitle('p_T(Z)')
ROOT.h_Zpt.SetYTitle('Weighted Entries')
ROOT.h_Zpt.SetTitle('')


#c1.SetLogy()
line.Draw('SAME')
legend.Draw()

c1.Print(os.path.join(plot_directory, 'Zpt.png'))
