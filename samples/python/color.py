import ROOT

from TopEFT.samples.helpers import singleton as singleton

@singleton
class color():
  pass

color.data           = ROOT.kBlack
color.nonprompt      = ROOT.kBlue-9
color.singleTop      = 40
color.ttX            = ROOT.kRed-10
color.ttXNoZ         = ROOT.kRed
color.ttH            = ROOT.kRed
color.ttW            = ROOT.kRed
color.ttZ            = ROOT.kOrange
color.ttZtoLLNuNu    = ROOT.kOrange
color.signal         = ROOT.kOrange
color.ttZtoQQ        = ROOT.kBlue
color.ttgamma        = ROOT.kRed
color.tZq            = 9
color.tWZ            = ROOT.kRed-7 #ROOT.kBlue-4
color.WJetsToLNu     = ROOT.kRed-10
color.diBoson        = ROOT.kOrange
color.multiBoson     = ROOT.kOrange
color.ZZ             = ROOT.kGreen+3#kOrange+1
color.WZ             = ROOT.kViolet-1
color.WW             = ROOT.kOrange-7
color.VV             = 30
color.WG             = ROOT.kOrange-5
color.ZG             = ROOT.kGreen#ROOT.kOrange-10
color.triBoson       = ROOT.kYellow
color.rare           = ROOT.kGreen+1
color.rare_noZZ      = ROOT.kGreen+1
color.WZZ            = ROOT.kYellow
color.WWG            = ROOT.kYellow-5
color.QCD            = 46
color.QCD_HT         = 46
color.QCD_Mu5        = 46
color.QCD_EMbcToE    = 46
color.QCD_Mu5EMbcToE = 46

color.other          = 46

color.T2tt_450_0                       = ROOT.kBlack
color.ttbarDMJets_scalar_Mchi1_Mphi200 = ROOT.kBlack
color.ttbarDMJets_scalar_Mchi1_Mphi10  = ROOT.kBlack
color.ttbarDMJets_scalar_Mchi1_Mphi20  = ROOT.kBlack
color.ttbarDMJets_scalar_Mchi1_Mphi100 = ROOT.kBlack
color.ttbarDMJets_pseudoscalar_Mchi1_Mphi100 = ROOT.kRed
color.ttbarDMJets_scalar_Mchi10_Mphi100 = ROOT.kPink
