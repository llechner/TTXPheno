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

pc5thresh = str(0.00005)
vcpt = -7.
vcpQM = 10.
vctZ = 0.
vctZI = 3.
textsize = 0.04
#lumi = 77 #36+41
#sigmaxBR = 
#Nsim = 5000

# sample = test
sample = dim6top_ttZ_ll_LO_highStat_scan
#sample = dim6top_ttZ_ll_LO_currentplane_highStat_scan

# just 1 file
sample.files = sample.files

# get TChain
c = sample.chain 

w = WeightInfo(sample.reweight_pkl)
w.set_order( 2 )

c.Draw("Z_pt>>h_Zpt(50,0,550)") # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("Z_pt>>h_Zpt1(50,0,550)", '(' + w.arg_weight_string(cpt=vcpt, cpQM=vcpQM) + ')/p_C[0]')
c.Draw("Z_pt>>h_Zptarg(50,0,550)", '(' + w.arg_weight_string(cpt=vcpt, cpQM=vcpQM, ctZ=vctZ, ctZI=vctZI) + ')/p_C[0]')
c.Draw("Z_pt>>h_Zptargcut(50,0,550)", '(' + w.arg_weight_string(cpt=vcpt, cpQM=vcpQM, ctZ=vctZ, ctZI=vctZI) + ')/p_C[0]*(p_C[5]<' + pc5thresh + ')')

#num = []
#num_arg = []
#num_argcut = []
#for i in range(1, ROOT.h_Zptarg.GetNbinsX()+1):
#    num.append(ROOT.h_Zpt.GetBinContent(i))
#    num_arg.append(ROOT.h_Zptarg.GetBinContent(i))
#    num_argcut.append(ROOT.h_Zptargcut.GetBinContent(i))

#print(sum(num))
#print(sum(num_arg))
#print(sum(num_argcut))

print(ROOT.h_Zpt.GetEntries())
print(ROOT.h_Zptarg.GetEntries())
print(ROOT.h_Zptargcut.GetEntries())

numarg = ROOT.h_Zptarg.GetEntries()
numargcut = ROOT.h_Zptargcut.GetEntries()
loss_abs = int(numarg - numargcut)
loss_per = loss_abs*100./float(numarg)

h_Zpt_max = ROOT.h_Zpt.GetMaximum()
h_Zptarg_max = ROOT.h_Zptarg.GetMaximum()

# Remove StatsBox
ROOT.h_Zpt.SetStats(False)
ROOT.h_Zptarg.SetStats(False)
ROOT.h_Zptargcut.SetStats(False)

# Line Color
ROOT.h_Zpt.SetLineColor(ROOT.kBlue)
ROOT.h_Zptarg.SetLineColor(ROOT.kRed)
ROOT.h_Zptargcut.SetLineColor(ROOT.kGreen)

ROOT.h_Zpt.SetLineWidth(2)
ROOT.h_Zptarg.SetLineWidth(2)
ROOT.h_Zptargcut.SetLineWidth(2)

#ROOT.h_Zpt1.GetYaxis().SetRangeUser(0,10000)
#ROOT.h_Zpt1.SetAxisRange(0,4, 'Y')

ROOT.h_Zpt.GetXaxis().SetTitleSize(textsize)
ROOT.h_Zpt.GetYaxis().SetTitleSize(textsize)
ROOT.h_Zptarg.GetXaxis().SetTitleSize(textsize)
ROOT.h_Zptarg.GetYaxis().SetTitleSize(textsize)
ROOT.h_Zptargcut.GetXaxis().SetTitleSize(textsize)
ROOT.h_Zptargcut.GetYaxis().SetTitleSize(textsize)

ROOT.h_Zpt.GetXaxis().SetLabelSize(textsize)
ROOT.h_Zpt.GetYaxis().SetLabelSize(textsize)
ROOT.h_Zptarg.GetXaxis().SetLabelSize(textsize)
ROOT.h_Zptarg.GetYaxis().SetLabelSize(textsize)
ROOT.h_Zptargcut.GetXaxis().SetLabelSize(textsize)
ROOT.h_Zptargcut.GetYaxis().SetLabelSize(textsize)

legend = ROOT.TLegend(0.57, 0.72, .9, .9)
legend.SetTextFont(1)
legend.SetTextSize(textsize)
#legend.SetTextSize(0.03)
legend.AddEntry(ROOT.h_Zpt, "pT(Z) SM", "l")
legend.AddEntry(ROOT.h_Zptarg, "pT(Z) dim-6", "l") #}{#splitline{cpt}{#splitline{cPQM}{#splitline{ctZ}{ctZI}}}}", "l")
legend.AddEntry(ROOT.h_Zptargcut, "pT(Z) dim-6, pC5 cut", "l") #}{#splitline{cpt}{#splitline{cPQM}{#splitline{ctZ}{ctZI}}}}", "l")

pt = ROOT.TPaveText(0.1,0.1,.38,.35,"blNDC")
pt.SetTextSize(textsize)
#pt.SetTextSize(0.03)
pt.SetTextFont(1)
pt.SetBorderSize(1)
pt.SetFillColor(ROOT.kWhite)
pt.AddText("cpt = " + str(vcpt).rstrip('0'))
pt.AddText("cpQM = " + str(vcpQM).rstrip('0'))
pt.AddText("ctZ = " + str(vctZ).rstrip('0'))
pt.AddText("ctZI = " + str(vctZI).rstrip('0'))
pt.AddText("pC5 < " + pc5thresh)
pt.AddText("Entries cut loss: %i" %loss_abs)

if h_Zpt_max < h_Zptarg_max:
    ROOT.h_Zptarg, ROOT.h_Zpt = ROOT.h_Zpt, ROOT.h_Zptarg

# Base Line
line1 = ROOT.TLine(340,0,340,ROOT.h_Zpt.GetMaximum())
line2 = ROOT.TLine(360,0,360,ROOT.h_Zpt.GetMaximum())
line = ROOT.TLine(300,0,300,ROOT.h_Zpt.GetMaximum())
line.SetLineColor(ROOT.kBlack)
line.SetLineWidth(2)

# Plotting
c1 = ROOT.TCanvas()
c1.SetCanvasSize(1500, 3000)
c1.SetWindowSize(500, 500)

ROOT.h_Zpt.Draw('HIST')
line.Draw()
#line1.Draw()
#line2.Draw()
#ROOT.h_Zpt.Draw('HIST')
c1.SetLogy()

ROOT.h_Zptarg.Draw('HIST SAME')
ROOT.h_Zptargcut.Draw('HIST SAME')

#ROOT.gStyle.SetTitleY(1.4)
ROOT.h_Zpt.SetXTitle('p_{T}(Z) [GeV]')
ROOT.h_Zpt.SetYTitle('Events')
ROOT.h_Zpt.SetTitle('')
#ROOT.h_Zpt.SetTitle('ttZ')


#c1.SetLogy()
legend.Draw('SAME')
pt.Draw()

c1.Print(os.path.join(plot_directory, 'Zpt_pc5_'+ pc5thresh +'.png'))
