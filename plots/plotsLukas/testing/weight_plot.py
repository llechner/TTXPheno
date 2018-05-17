''' simple draw script
'''

# Standrard imports
import ROOT
import os

# TTXPheno
from TTXPheno.Tools.user import plot_directory
from TTXPheno.Tools.WeightInfo import WeightInfo
from TTXPheno.Tools.WeightInfo import BinContentToList

# Sample
from TTXPheno.samples.benchmarks import * # dim6top_ttZ_ll_LO_currentplane_highStat_scan 

textsize = 0.04

# sample = test
sample = dim6top_ttZ_ll_LO_highStat_scan
#sample = dim6top_ttZ_ll_LO_currentplane_highStat_scan

# just 1 file
sample.files = sample.files

# get TChain
c = sample.chain 

h = ROOT.TH1F('h_weights','weights',15,0,15)
c.Draw('Iteration$>>h_weights','p_C/p_C[0]')
c.Draw('Iteration$>>h_weights_Zpt','p_C*(Z_pt>300)/p_C[0]')

scale = ROOT.h_weights.GetEntries()/float(ROOT.h_weights.GetNbinsX())
ROOT.h_weights.Scale(1./scale)
#print(scale)
#scale = ROOT.h_weights_Zpt.GetEntries()/float(ROOT.h_weights_Zpt.GetNbinsX())
ROOT.h_weights_Zpt.Scale(1./scale)
#print(scale)

c1 = ROOT.TCanvas()
ROOT.h_weights.SetXTitle('Weight Number')
ROOT.h_weights.SetYTitle('Yield per Weight [c0]')
ROOT.h_weights.SetTitle('')
#ROOT.h_weights.SetTitle('ttZ')
ROOT.h_weights.SetStats(False)

# Line Color
ROOT.h_weights.SetLineColor(ROOT.kBlue)
ROOT.h_weights_Zpt.SetLineColor(ROOT.kRed)

ROOT.h_weights.SetLineWidth(2)
ROOT.h_weights_Zpt.SetLineWidth(2)

ROOT.h_weights.GetXaxis().SetTitleSize(textsize)
ROOT.h_weights.GetYaxis().SetTitleSize(textsize)
ROOT.h_weights.GetXaxis().SetLabelSize(textsize)
ROOT.h_weights.GetYaxis().SetLabelSize(textsize)

ROOT.h_weights_Zpt.GetXaxis().SetLabelSize(textsize)
ROOT.h_weights_Zpt.GetYaxis().SetLabelSize(textsize)
ROOT.h_weights_Zpt.GetXaxis().SetTitleSize(textsize)
ROOT.h_weights_Zpt.GetYaxis().SetTitleSize(textsize)

# Base Line
line = ROOT.TLine(0,0,15,0)
line.SetLineColor(ROOT.kBlack)
line.SetLineWidth(1)

# Legend
legend = ROOT.TLegend(0.62, 0.78, .9, .9)
legend.SetTextSize(textsize)
legend.SetTextFont(1)
legend.AddEntry(ROOT.h_weights, "no cut", "l")
legend.AddEntry(ROOT.h_weights_Zpt, "pT(Z) > 300 GeV", "l")


#labelling
label0 = ROOT.TText()
label0.SetTextAngle(90)
label0.SetNDC()
label0.SetTextFont(1)
label0.SetTextSize(textsize)
label0.SetTextColor(ROOT.kBlack)


ROOT.h_weights.Draw('HIST')
ROOT.h_weights_Zpt.Draw('HIST SAME')


w = WeightInfo(sample.reweight_pkl)
w.set_order( 2 )
for k, el in enumerate(['c0'] + ['*'.join([item.split('rw_')[1] for item in el.split('*')[1:]]) for el in w.weight_string().split('+')][1:]):
    label0.DrawText(0.13+k*.053, 0.35, el)
#label0.Draw('SAME')
line.Draw('SAME')
legend.Draw()

c1.Print(os.path.join(plot_directory, 'weight.png'))
