#!/usr/bin/env python

''' simple draw script
'''

# Standrard imports
import ROOT
import os
import argparse

# TTXPheno
from TTXPheno.Tools.user import plot_directory
from TTXPheno.Tools.WeightInfo import WeightInfo
from TTXPheno.Tools.helpers import create_cut_string

# Sample
from TTXPheno.samples.benchmarks import *

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--sample', action='store', default=fwlite_ttZ_ll_LO_order3_8weights, nargs=1)
argParser.add_argument('--small', action='store_true')
argParser.add_argument('--parameters', action='store', default=['Z_pt'], nargs='*', type = str)
argParser.add_argument('--range', action='store', default=[50,0,500], nargs=3)
argParser.add_argument('--filename', action='store',         default="./reweight_card.dat",  nargs=1,    type = str, help="Output filename")
args = argParser.parse_args()


#######################################################
#######################################################

param = args.parameters #[GenJet_pt, GenJet_eta, GenJet_phi, GenLep_pt, GenLep_eta, GenLep_phi, top_pt, top_eta, top_phi, Z_pt, Z_phi, Z_eta, Z_cosThetaStar]



# Sample Settings
sample = args.sample #fwlite_ttZ_ll_LO_order3_8weights
small = args.small #False
order = 3

# Plot Settings
plot_param = "Z_pt"
log = False
binning = args.range[0] #50
plot_min = args.range[0]0
plot_max = 500
textsize = 0.03

# Naming Settings
shapes = "analysisPlots"
process = "TTZ" #TTgamma #TTW
flavor = "all"
plot = plot_param

# Cuts
lower_cut_dict = {'p_C[5]':0.0005, 'p_C[40]':0.004, 'p_C[20]':0.0004}
upper_cut_dict = {}

# Create Wilson coefficient dictionary
WC_dict = {'cpt':-7., 'cpQM':10., 'ctZ':0., 'ctZI':3.}

#######################################################
#######################################################

def remove_zeros_from_dict(dict):
    return {x:y for x,y in dict.items() if y!=0}

WC_dict = remove_zeros_from_dict(WC_dict)

def get_selection_string(dict, operator):
    if len(dict.keys())==0: return ''

    selection = ''
    addon = ''

    if operator == '<': addon = 'lt'
    elif operator == '>': addon = 'gt'

    for element in dict.keys():
        selection += element.replace('_','')
        selection += '_' + addon
        selection += str(dict[element]).rstrip('0').replace('.','p').replace('-','m').rstrip('p')
        selection += '_'
    return selection.rstrip('_').replace('__','_').replace(']','').replace('[','')

if len(WC_dict.keys())==0: WC = "SM"
else: WC = get_selection_string(WC_dict, '')

if len(lower_cut_dict.keys())==0 and len(upper_cut_dict.keys())==0: selection = "all"
else:
    selection = get_selection_string(upper_cut_dict, '>')
    selection += get_selection_string(lower_cut_dict, '<')

if log: loglin = "log"
else: loglin = "lin"

# Create plot directory
plot_pathdir = '/'.join([shapes, sample.name, process, WC, selection, flavor + '_' + loglin])
if not os.path.exists(os.path.join(plot_directory,plot_pathdir)): os.makedirs(os.path.join(plot_directory,plot_pathdir))
plot_name = plot + '.png'

# Create cut string
cut_string = create_cut_string(lower_cut_dict, '<')
cut_string += create_cut_string(upper_cut_dict, '>')

# Just 1 file
if small: sample.files = sample.files[:10]

# get TChain
c = sample.chain 

w = WeightInfo(sample.reweight_pkl)
w.set_order( order )

c.Draw(plot_param + ">>h_Zpt(" + str(binning) + "," + str(plot_min) + "," + str(plot_max) + ")")
c.Draw(plot_param + ">>h_Zptarg(" + str(binning) + "," + str(plot_min) + "," + str(plot_max) + ")", '(' + w.get_weight_string(**WC_dict) + ')/p_C[0]')


if cut_string != '':
   c.Draw(plot_param + ">>h_Zptargcut(" + str(binning) + "," + str(plot_min) + "," + str(plot_max) + ")", '(' + w.get_weight_string(**WC_dict) + ')/p_C[0]' + cut_string)
   
   numarg = ROOT.h_Zptarg.GetEntries()
   numargcut = ROOT.h_Zptargcut.GetEntries()
   loss_abs = int(numarg - numargcut)


h_Zpt_max = ROOT.h_Zpt.GetMaximum()
h_Zptarg_max = ROOT.h_Zptarg.GetMaximum()

# Remove StatsBox
ROOT.h_Zpt.SetStats(False)
ROOT.h_Zptarg.SetStats(False)

# Line Color
ROOT.h_Zpt.SetLineColor(ROOT.kBlue)
ROOT.h_Zptarg.SetLineColor(ROOT.kRed)

ROOT.h_Zpt.SetLineWidth(2)
ROOT.h_Zptarg.SetLineWidth(2)

ROOT.h_Zpt.GetXaxis().SetTitleSize(textsize)
ROOT.h_Zpt.GetYaxis().SetTitleSize(textsize)
ROOT.h_Zptarg.GetXaxis().SetTitleSize(textsize)
ROOT.h_Zptarg.GetYaxis().SetTitleSize(textsize)

ROOT.h_Zpt.GetXaxis().SetLabelSize(textsize)
ROOT.h_Zpt.GetYaxis().SetLabelSize(textsize)
ROOT.h_Zptarg.GetXaxis().SetLabelSize(textsize)
ROOT.h_Zptarg.GetYaxis().SetLabelSize(textsize)

legend = ROOT.TLegend(0.67, 0.82, .98, .95, '', "NDC")
legend.SetBorderSize(1)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextFont(1)
legend.SetTextSize(textsize)
legend.AddEntry(ROOT.h_Zpt, "pT(Z) SM", "l")
legend.AddEntry(ROOT.h_Zptarg, "pT(Z) dim-6", "l")

elements = (len(WC_dict.keys()) + len(lower_cut_dict.keys()) + len(upper_cut_dict.keys()) + 1)
pt = ROOT.TPaveText(0.16,0.13,.43,.05*elements,"blNDC")
pt.SetTextSize(textsize)
pt.SetTextFont(1)
pt.SetBorderSize(1)
pt.SetFillColor(ROOT.kWhite)

for item in WC_dict.keys():
   pt.AddText(str(item) + " = " + str(WC_dict[item]).rstrip('0'))

if cut_string != '':
   ROOT.h_Zptargcut.SetStats(False)
   ROOT.h_Zptargcut.SetLineColor(ROOT.kGreen)
   ROOT.h_Zptargcut.SetLineWidth(2)
   ROOT.h_Zptargcut.GetXaxis().SetTitleSize(textsize)
   ROOT.h_Zptargcut.GetYaxis().SetTitleSize(textsize)
   ROOT.h_Zptargcut.GetXaxis().SetLabelSize(textsize)
   ROOT.h_Zptargcut.GetYaxis().SetLabelSize(textsize)

   legend.AddEntry(ROOT.h_Zptargcut, "pT(Z) dim-6, cut", "l")

   for item in lower_cut_dict.keys():
      pt.AddText(str(item) + " < " + str(lower_cut_dict[item]).rstrip('0'))
   for item in upper_cut_dict.keys():
      pt.AddText(str(item) + " > " + str(upper_cut_dict[item]).rstrip('0'))

   pt.AddText("Entries cut loss: %i" %loss_abs)


if h_Zpt_max < h_Zptarg_max: ROOT.h_Zptarg, ROOT.h_Zpt = ROOT.h_Zpt, ROOT.h_Zptarg

# Plotting
c1 = ROOT.TCanvas()
c1.SetCanvasSize(800,800)

p1 = ROOT.TPad('', '', 0., 0., 0.95, 0.95)
p1.Draw()
p1.cd()
#c1.SetWindowPosition(-50,50)
#c1.SetWindowSize(800,600)

ROOT.h_Zpt.Draw('HIST')
if loglin == 'log': c1.SetLogy()
ROOT.h_Zptarg.Draw('HIST SAME')
if cut_string != '': ROOT.h_Zptargcut.Draw('HIST SAME')

ROOT.h_Zpt.SetXTitle(str(plot_param) + ' [GeV]')
ROOT.h_Zpt.SetYTitle('Events')
ROOT.h_Zpt.GetXaxis().SetTitleOffset(1.1)
ROOT.h_Zpt.GetYaxis().SetTitleOffset(2.2)


legend.Draw()
pt.Draw()

c1.Update()
c1.Print(os.path.join(plot_directory, plot_pathdir, plot_name))
