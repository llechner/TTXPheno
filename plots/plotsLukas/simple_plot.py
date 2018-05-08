''' simple draw script
'''

# Standrard imports
import ROOT
import os

# TTXPheno
from TTXPheno.Tools.user import plot_directory
from TTXPheno.Tools.WeightInfo import WeightInfo

# Sample
from TTXPheno.samples.gen_fwlite_benchmarks import * # dim6top_ttZ_ll_LO_currentplane_highStat_scan 

sample = dim6top_ttZ_ll_LO_currentplane_highStat_scan

# just 1 file
sample.files = sample.files[:1]

# get TChain
c = sample.chain 

w = WeightInfo(sample.reweight_pkl)
w.set_order( 2 )


#c.Draw("Z_pt>>h_Zpt(50,0,400)") # "weight*(%s)" % weightString(cpt=0.2)
c.Draw("Z_pt>>h_Zpt(50,0,400)",  w.arg_weight_string(cpt=10, cpQM=30) )

c1 = ROOT.TCanvas()
ROOT.h_Zpt.Draw('hist')

#c1.Print(os.path.join(plot_directory, 'myplot.png'))
c1.Print(os.path.join(plot_directory, 'myplot_weight.png'))




#def weightString(cut_list):
#    if type(cut_list)==list: return '&&'.join(['(' + str(item) + ')' for item in cut_list])
#    else: return '(%s)' % cut_list

#print weightString(['cpt=0.2','Zpt<5'])
#exit()
