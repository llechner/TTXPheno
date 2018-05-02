import ROOT

c = ROOT.TChain("Events")
c.Add(" ... .root")
# c.GetListOfBranches().ls()

c.Draw("Z_pt>>h_Zpt(50,0,400)", "weight*(njet>=2)") # "weight*(%s)" % weightString(cpt=0.2)

c1 = ROOT.TCanvas()
ROOT.h_Zpt.Draw('hist')

c1.Print('...plot dir/file.png')

