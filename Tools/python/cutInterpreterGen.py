''' Class to interpret string based cuts
'''

import logging
logger = logging.getLogger(__name__)

special_cuts = {
    "lepSel1":            "Sum$(genLep_pt>10&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)==1",
#    "lepSel2":            "Sum$(genLep_pt>10&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)==2&&Sum$(genLep_pt>20&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)>=1",
    "lepSel2":            "Sum$(genLep_pt>15&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.4)==2&&Sum$(genLep_pt>25&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.4)>=1",
    "lepSel3":            "Sum$(genLep_pt>10&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)==3&&Sum$(genLep_pt>20&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)>=2&&Sum$(genLep_pt>40&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)>=1",
    "lepSel4":            "Sum$(genLep_pt>10&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)==4&&Sum$(genLep_pt>40&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)>=1",
    "onZ":                "abs(genZ_mass-91.2)<=10",
    "offZ":               "abs(genZ_mass-91.2)>10",
  }

continous_variables = [ ("mll", "genZ_mass"), ("met", "genMet_pt"), ("Zpt","genZ_pt"), ("gammaeta","abs(genPhoton_eta[0])"), ("gammapt","genPhoton_pt[0]"), ("Wpt","genW_pt"), ("relIso","genPhoton_relIso04[0]")]
discrete_variables  = [ ("nlepttgamma2l", "Sum$(genLep_pt>15&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.4)"), ("nlep", "Sum$(genLep_pt>10&&(abs(genLep_pdgId)==11||abs(genLep_pdgId)==13)&&abs(genLep_eta)<2.5)"), ("njet", "Sum$(genJet_pt>30&&abs(genJet_eta)<2.4)"), ("nbjet", "Sum$(genJet_pt>30&&genJet_matchBParton>=1&&abs(genJet_eta)<2.4)") , ]

continous_variables += [ ("genLepZmass", "genLepZ_mass"), ("genLepZpt","genLepZ_pt")]

from TTXPheno.Tools.CutInterpreter import CutInterpreter
cutInterpreter = CutInterpreter( continous_variables, discrete_variables, special_cuts)

#print cutInterpreter.cutString("lepSel2-gammapt20-njet2p-nbjet1p-relIso0to0.4")
