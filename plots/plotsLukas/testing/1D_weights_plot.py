''' simple draw script
'''

# Standrard imports
import ROOT
import os
import numpy as np

# TTXPheno
from TTXPheno.Tools.user import plot_directory
from TTXPheno.Tools.WeightInfo import WeightInfo

# Sample
from TTXPheno.samples.benchmarks import * # dim6top_ttZ_ll_LO_currentplane_highStat_scan 

textsize=0.04
# sample = test
sample = dim6top_ttZ_ll_LO_highStat_scan
#sample = dim6top_ttZ_ll_LO_currentplane_highStat_scan

# just 1 file
sample.files = sample.files[:1]

# get TChain
c = sample.chain 

w = WeightInfo(sample.reweight_pkl)
w.set_order( 2 )

h = ROOT.TH1F('h_weights','weights',15,0,15)
c.Draw('Iteration$>>h_weights','p_C')
#c.Draw('Iteration$>>h_weights','p_C*(Z_pt>300)')

sum_weights = WeightInfo.BinContentToList(ROOT.h_weights)

def plot1Dcoefficient(coeff, range):

    x = np.arange(range[0],range[1],range[2])
    y = np.array([w.Get1DYield(sum_weights, coeff, i) for i in x])/sum_weights[0] #weight to SM
    graph = ROOT.TGraph(len(x), x, y)

    graph.SetLineWidth(2)
    graph.SetLineColor(ROOT.kRed)
    graph.GetXaxis().SetRangeUser(range[0],range[1])
    graph.GetYaxis().SetTitle('Total Yield [SM Yield]')
    graph.GetXaxis().SetTitle(coeff)
#    graph.SetTitle('ttZ')
    graph.SetTitle('')

    graph.GetXaxis().SetTitleSize(textsize)
    graph.GetYaxis().SetTitleSize(textsize)
    graph.GetXaxis().SetLabelSize(textsize)
    graph.GetYaxis().SetLabelSize(textsize)

    # Base Line
    line = ROOT.TLine(range[0],1,range[1],1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kBlue)

    # Legend
    legend = ROOT.TLegend(0.63, 0.78, .9, .9)
    legend.SetTextSize(textsize)
    legend.SetTextFont(1)
    legend.AddEntry(line, "SM Yield", "l")
    legend.AddEntry(graph, "Total Yield", "l")

    # Plotting
    c1 = ROOT.TCanvas()
    graph.Draw()
    line.Draw('SAME')
    legend.Draw()
    c1.Print(os.path.join(plot_directory, '1D%s.png') %(coeff))


def plot1D2coefficient(coeff, range, fixedcoeff, fixedvalue):
    x = np.arange(range[0],range[1],range[2])
    y = []
    for i in x:
        dict = {coeff:i, fixedcoeff:fixedvalue}
        y.append(w.GetNDYield(sum_weights, **dict))
    y = np.array(y)/sum_weights[0] #weight to SM
    graph = ROOT.TGraph(len(x), x, y)

    graph.SetLineWidth(2)
    graph.SetLineColor(ROOT.kRed)
    graph.GetXaxis().SetRangeUser(range[0],range[1])
    graph.GetYaxis().SetTitle('Total Yield [SM Yield]')
    graph.GetXaxis().SetTitle(coeff)
    graph.SetTitle('')
#    graph.SetTitle('ttZ')

    graph.GetXaxis().SetTitleSize(textsize)
    graph.GetYaxis().SetTitleSize(textsize)
    graph.GetXaxis().SetLabelSize(textsize)
    graph.GetYaxis().SetLabelSize(textsize)

    # Base Line
    line = ROOT.TLine(range[0],1,range[1],1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kBlue)

    # Legend
    legend = ROOT.TLegend(0.53, 0.78, .9, .9)
    legend.SetTextSize(textsize)
    legend.SetTextFont(1)
    legend.AddEntry(line, "SM Yield", "l")
    legend.AddEntry(graph, "Total Yield (%s = %s)" %(fixedcoeff, fixedvalue), "l")

    # Plotting
    c1 = ROOT.TCanvas()
    graph.Draw()
    line.Draw('SAME')
    legend.Draw()
    c1.Print(os.path.join(plot_directory, '1D%s_%s%s.png') %(coeff, fixedcoeff,fixedvalue))



for coeff in ['ctZ','ctZI','cpt','cpQM']:
    if coeff in ['ctZ','ctZI']:
        range = [-2,2.1,.1]
        plot1Dcoefficient(coeff, range)
    elif coeff == 'cpt': #fix the cpQM coeff to 10 to get the optimum Graph (Daniels Plot)
        range = [-25,15,.1]
        plot1D2coefficient(coeff, range, 'cpQM', 10)
    elif coeff == 'cpQM': # fix the cpt coeff to +7 to get the optimum Graph (Daniels Plot)
        range = [-10,30,.1]
        plot1D2coefficient(coeff, range, 'cpt', -7)


