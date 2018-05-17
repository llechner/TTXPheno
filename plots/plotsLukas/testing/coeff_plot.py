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

index = 5
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
minrange = str(-5)
maxrange = str(5)
c.Draw("p_C[" + str(index) + "]>>h_pc(50,0.,.0003)",'(Z_pt>300)*(p_C[5]<0.002)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[" + str(index) + "]>>h_pcptZ(50,0.158,.16)",'(Z_pt>345)*(Z_pt<355)') # "weight*(%s)" % weightString(cpt=0.2)
c.Draw("p_C[" + str(index) + "]>>h_pcptZ(50,-1,1)",'(Z_pt>345)*(Z_pt<355)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[0]>>h_pc0ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[1]>>h_pc1ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[2]>>h_pc2ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[3]>>h_pc3ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[4]>>h_pc4ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[5]>>h_pc5ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[6]>>h_pc6ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[7]>>h_pc7ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[8]>>h_pc8ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[9]>>h_pc9ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[10]>>h_pc10ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[11]>>h_pc11ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[12]>>h_pc12ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[13]>>h_pc13ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[14]>>h_pc14ptZ(50," + minrange + "," + maxrange + ")",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[" + str(index) + "]>>h_pc1(50,-.00005,.00005)") # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("p_C[" + str(index) + "]>>h_pc1ptZ(50,-.00005,.00005)",'(Z_pt>340)*(Z_pt<360)') # "weight*(%s)" % weightString(cpt=0.2)
#c.Draw("Z_pt>>h_Zpt1(50,0,550)", '(' + w.arg_weight_string(cpt=vcpt, cpQM=vcpQM) + ')/p_C[0]')
#c.Draw("p_C[1]>>h_pc1(50,0,550)", '(' + w.arg_weight_string(cpt=vcpt, cpQM=vcpQM, ctZ=vctZ, ctZI=vctZI) + ')')

# Line Color
ROOT.h_pc.SetLineColor(ROOT.kBlue)
ROOT.h_pcptZ.SetLineColor(ROOT.kRed)
#ROOT.h_pc0ptZ.SetLineColor(0)
#ROOT.h_pc1ptZ.SetLineColor(1)
#ROOT.h_pc2ptZ.SetLineColor(2)
#ROOT.h_pc3ptZ.SetLineColor(3)
#ROOT.h_pc4ptZ.SetLineColor(4)
#ROOT.h_pc5ptZ.SetLineColor(5)
#ROOT.h_pc6ptZ.SetLineColor(6)
#ROOT.h_pc7ptZ.SetLineColor(7)
#ROOT.h_pc8ptZ.SetLineColor(8)
#ROOT.h_pc9ptZ.SetLineColor(9)
#ROOT.h_pc10ptZ.SetLineColor(10)
#ROOT.h_pc11ptZ.SetLineColor(11)
#ROOT.h_pc12ptZ.SetLineColor(12)
#ROOT.h_pc13ptZ.SetLineColor(13)
#ROOT.h_pc14ptZ.SetLineColor(14)

ROOT.h_pc.SetLineWidth(2)
ROOT.h_pcptZ.SetLineWidth(2)

# Plotting
c1 = ROOT.TCanvas()
ROOT.h_pc.Draw('HIST')
#ROOT.h_pcptZ.Draw('HIST SAME')
c1.SetLogy()
#ROOT.h_pc1ptZ.Draw('HIST SAME')
#ROOT.h_pc1ptZ.Draw('HIST SAME')
#ROOT.h_pc2ptZ.Draw('HIST SAME')
#ROOT.h_pc3ptZ.Draw('HIST SAME')
#ROOT.h_pc4ptZ.Draw('HIST SAME')
#ROOT.h_pc5ptZ.Draw('HIST SAME')
#ROOT.h_pc6ptZ.Draw('HIST SAME')
#ROOT.h_pc7ptZ.Draw('HIST SAME')
#ROOT.h_pc8ptZ.Draw('HIST SAME')
#ROOT.h_pc9ptZ.Draw('HIST SAME')
#ROOT.h_pc10ptZ.Draw('HIST SAME')
#ROOT.h_pc11ptZ.Draw('HIST SAME')
#ROOT.h_pc12ptZ.Draw('HIST SAME')
#ROOT.h_pc13ptZ.Draw('HIST SAME')
#ROOT.h_pc14ptZ.Draw('HIST SAME')
c1.Print(os.path.join(plot_directory, 'pc'+str(index)+'.png'))
