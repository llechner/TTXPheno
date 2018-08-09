''' RECO level cuts 
'''

import logging
logger = logging.getLogger(__name__)

special_cuts = {
    "lepSel1":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==1",
    "lepSel2":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==2&&Sum$(recoLep_pt>20&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=1",
    "lepSel3":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==3&&Sum$(recoLep_pt>20&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=2&&Sum$(recoLep_pt>40&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=1",
    "lepSel4":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==4&&Sum$(recoLep_pt>40&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=1",
    "onZ":                "abs(recoZ_mass-91.2)<=10",
    "offZ":               "abs(recoZ_mass-91.2)>10",
  }

continous_variables = [ ("mll", "recoZ_mass"), ("met", "recoMet_pt"), ("Zpt","recoZ_pt"), ("lepZpt","recoLepZ_pt"), ("gammapt","recoPhoton_pt[0]"), ("Wpt","recoW_pt"), ("relIso","genPhoton_relIso04[0]")]
discrete_variables  = [ ("nlepttgamma2l", "Sum$(recoLep_pt>15&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.4)"), ("nlep", "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)"), ("njet", "Sum$(recoJet_pt>30&&abs(recoJet_eta)<2.4)"), ("nbjet", "Sum$(recoJet_pt>30&&recoJet_bTag>=1&&abs(recoJet_eta)<2.4)") , ]

from TTXPheno.Tools.CutInterpreter import CutInterpreter
cutInterpreter = CutInterpreter( continous_variables, discrete_variables, special_cuts)
