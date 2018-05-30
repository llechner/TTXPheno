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

# sample = test
sample = dim6top_ttZ_ll_LO_highStat_scan
#sample = dim6top_ttZ_ll_LO_currentplane_highStat_scan

# just 1 file
sample.files = sample.files[:1]

# get TChain
c = sample.chain 

h = ROOT.TH1F('h_weights','weights',15,0,15)
c.Draw('Iteration$>>h_weights','p_C')
c.Draw('Iteration$>>h_weights_Zpt','p_C*(Z_pt>300)')


c1 = ROOT.TCanvas()
ROOT.h_weights.SetXTitle('Entry')
ROOT.h_weights.SetYTitle('Weight')
ROOT.h_weights.SetTitle('')
ROOT.h_weights.SetStats(False)

# Line Color
ROOT.h_weights.SetLineColor(ROOT.kBlue)
ROOT.h_weights_Zpt.SetLineColor(ROOT.kRed)

# Base Line
line = ROOT.TLine(0,0,15,0)
line.SetLineColor(ROOT.kBlack)
line.SetLineWidth(1)

# Legend
legend = ROOT.TLegend(0.65, 0.75, .9, .9)
legend.AddEntry(ROOT.h_weights, "no cut", "l")
legend.AddEntry(ROOT.h_weights_Zpt, "p_T(Z) > 300", "l")

ROOT.h_weights.Draw('HIST')
ROOT.h_weights_Zpt.Draw('HIST SAME')
line.Draw('SAME')
legend.Draw()

c1.Print(os.path.join(plot_directory, 'weight.png'))
