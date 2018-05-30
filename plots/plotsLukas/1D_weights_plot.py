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
#c.Draw('Iteration$>>h_weights_Zpt','p_C*(Z_pt>300)')

sum_weights = WeightInfo.BinContentToList(ROOT.h_weights)

def plot1Dcoefficient(coeff, range):

    x = np.arange(range[0],range[1],range[2])
    y = np.array([w.get_weight_yield(sum_weights, **{coeff: i}) for i in x])/sum_weights[0] #weight to SM
    graph = ROOT.TGraph(len(x), x, y)

    # Base Line
    line = ROOT.TLine(range[0],1,range[1],1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(1)

    # Plotting
    c1 = ROOT.TCanvas()
    graph.Draw()
    line.Draw('SAME')
    c1.Print(os.path.join(plot_directory, '1D%s.png') %(coeff))


def plot1D2coefficient(coeff, range, fixedcoeff, fixedvalue):
    x = np.arange(range[0],range[1],range[2])
    y = []
    for i in x:
        dict = {coeff:i, fixedcoeff:fixedvalue}
        y.append(w.get_weight_yield(sum_weights, **dict))
    y = np.array(y)/sum_weights[0] #weight to SM
    graph = ROOT.TGraph(len(x), x, y)

    # Base Line
    line = ROOT.TLine(range[0],1,range[1],1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(1)

    # Plotting
    c1 = ROOT.TCanvas()
    graph.Draw()
    line.Draw('SAME')
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


