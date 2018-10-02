''' RECO level cuts 
'''

import logging
logger = logging.getLogger(__name__)

special_cuts = {
    "lepSel1":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==1",
    "lepSel2":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==2&&Sum$(recoLep_pt>20&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=1",
#    "lepSel2":            "Sum$(recoLep_pt>15&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.4)==2&&Sum$(recoLep_pt>25&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.4)>=1",
    "lepSel3":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==3&&Sum$(recoLep_pt>20&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=2&&Sum$(recoLep_pt>40&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=1",
    "lepSel4":            "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)==4&&Sum$(recoLep_pt>40&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)>=1",
    "onZ":                "abs(recoZ_mass-91.2)<=10",
    "offZ":               "abs(recoZ_mass-91.2)>10",
    "leptonIso1":         "recoLep_isolationVar[0]<=0.6",
    "leptonIso2":         "recoLep_isolationVar[0]<=0.6&&recoLep_isolationVar[1]<=0.6",
    "leptonIso3":         "recoLep_isolationVar[recoNonZ_l1_index]<=0.06&&recoLep_isolationVar[recoZ_l1_index]<=0.06&&recoLep_isolationVar[recoZ_l2_index]<=0.06",
  }

continous_variables = [ ("mll", "recoZ_mass"), ("met", "recoMet_pt"), ("Zpt","recoZ_pt"), ("lepZpt","recoLepZ_pt"), ("gammapt","recoPhoton_pt[0]"), ("Wpt","recoW_pt"), ("relIso","genPhoton_relIso04[0]")]
discrete_variables  = [ ("nlepttgamma2l", "Sum$(recoLep_pt>15&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<=2.4)"), ("nlep", "Sum$(recoLep_pt>10&&(abs(recoLep_pdgId)==11||abs(recoLep_pdgId)==13)&&abs(recoLep_eta)<2.5)"), ("njet", "Sum$(recoJet_pt>30&&abs(recoJet_eta)<2.4)"), ("nb0jet", "Sum$(recoJet_pt>30&&(recoJet_bTag&(2**0))&&abs(recoJet_eta)<2.4)"), ("nb1jet", "Sum$(recoJet_pt>30&&(recoJet_bTag&(2**1))&&abs(recoJet_eta)<2.4)"),("nb2jet", "Sum$(recoJet_pt>30&&(recoJet_bTag&(2**2))&&abs(recoJet_eta)<2.4)"),("nb3jet", "Sum$(recoJet_pt>30&&(recoJet_bTag&(2**3))&&abs(recoJet_eta)<2.4)"),("nb4jet", "Sum$(recoJet_pt>30&&(recoJet_bTag&(2**4))&&abs(recoJet_eta)<2.4)"),("nb5jet", "Sum$(recoJet_pt>30&&(recoJet_bTag&(2**5))&&abs(recoJet_eta)<2.4)"),("nbjet", "Sum$(recoJet_pt>30&&recoJet_bTag>=1&&abs(recoJet_eta)<2.4)") , ]

from TTXPheno.Tools.CutInterpreter import CutInterpreter
cutInterpreter = CutInterpreter( continous_variables, discrete_variables, special_cuts)
