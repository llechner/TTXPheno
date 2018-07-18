class GenSearch:

    def __init__(self, genParticles):
        if type( genParticles ) != list:
            self.genParticles = list(genParticles)
        else:
            self.genParticles = genParticles

#    def mother(self, l):
#        ''' Finds the mother by climbing the ladder until there is a mother that has a different pdgId
#        '''
#        if l['nMothers']>1: 
#            return None
#        if l['pdgId']!=self.genParticles[l['motherIndex1']]['pdgId']: return self.genParticles[l['motherIndex1']]
#        else: return self.mother(self.genParticles[l['motherIndex1']])

#    def ancestry(self, l, stop_at_pdgId = [2212]):
#        self.__found = set()
#        return self.__ancestry(l, stop_at_pdgId = stop_at_pdgId )

    @property
    def final_state_particles_no_neutrinos( self ):
        if hasattr(self, "final_state_particles_no_neutrinos_"):
            return self.final_state_particles_no_neutrinos_
        else:
            self.final_state_particles_no_neutrinos_ =  filter( lambda p: p.status()==1 and abs(p.pdgId()) not in [12,14,16], self.genParticles ) 
            return self.final_state_particles_no_neutrinos_


    def daughters(self, p):
        return [p.daughter(i) for i in xrange( p.numberOfDaughters()) ]

    def mothers(self, p):
        return [p.mother(i) for i in xrange( p.numberOfMothers()) ]

    def isLast(self, p ):
        ''' Returns true if the decay products of p do not contain a particle with p.pdgId() '''
        return p.pdgId() not in [d.pdgId() for d in self.daughters( p )]

    def descend(self, p):
        ''' Returns the last particle of the same pdgId in the decay chain started by p
        '''
        cands = filter(lambda q:abs(q.pdgId())==abs(p.pdgId()), self.daughters(p) )
        #if len(cands)>1: print "Warning: Found more than one particle with same pdgId %i in decay chain %r -> %r."%(p.pdgId(), p, cands)
        if len(cands)>0:
            if cands[0]==p:return p
            return self.descend(cands[0])
        return p

    def ascend(self, p):
        ''' Returns the first particle of the same pdgId in the decay chain started by p
        '''
        cands = filter(lambda q:abs(q.pdgId())==abs(p.pdgId()), self.mothers(p) )
        #if len(cands)>1: print "Warning: Found more than one particle with same pdgId %i in decay chain %r -> %r."%(p['pdgId'], p, cands)
        if len(cands)>0:
            if cands[0]==p:return p
            return self.ascend(cands[0])
        return p

    def print_decay( self, p, prefix = ""):
        print prefix+" %s(pt %3.2f, eta %3.2f, phi %3.2f)" % ( pdgToName(p.pdgId()), p.pt(), p.eta(), p.phi())
        for d in self.daughters( self.descend( p ) ):
            self.print_decay( d, prefix = prefix+'--')

    def __ancestry(self, p, stop_at_pdgId):
        ''' Returns indices of all genParticles in the ancestry of l
        '''
        # print "Looking at %s, nMothers %i, already self.__found: %s"% ( pdgToName(l['pdgId']), l['nMothers'], ",".join(str(x) for x in self.__found) )
        for m in self.mothers(p):
            if abs(m.pdgId()) in stop_at_pdgId or m.numberOfMothers()==0:
                # print "Done."
                return self.__found
            else:
                self.__found.add( m )
                self.__found.union( self.__ancestry( m , stop_at_pdgId = stop_at_pdgId ) )
        # print "Done alltogether."
        return self.__found
            

D_mesons = set([ 411, -411, 421, -421, 413, -413, 423, -423, 10411, -10411, 10421, -10421, 10413, -10413, 10423, -10423, 415, -415, 425, -425, 20413, -20413, 20423, -20423, 431, -431, 433, -433, 10431, -10431, 10433, -10433, 435, -435, 20433, -20433])
B_mesons = set([ 511, -511, 521, -521, 513, -513, 523, -523, 10511, -10511, 10521, -10521, 10513, -10513, 10523, -10523, 515, -515, 525, -525, 20513, -20513, 20523, -20523, 531, -531, 533, -533, 10531, -10531, 10533, -10533, 535, -535, 20533, -20533, 541, -541, 543, -543, 10541, -10541, 10543, -10543, 545, -545, 20543, -20543])
massive_vector_bosons = set([23, 24, -24])

D_mesons_abs = set(abs(x) for x in D_mesons)
B_mesons_abs = set(abs(x) for x in B_mesons)
massive_vector_bosons_abs = set([23, 24])
#
pdgToName_ =\
{1:"d",
-1:"dbar",
2:"u",
-2:"ubar",
3:"s",
-3:"sbar",
4:"c",
-4:"cbar",
5:"b",
-5:"bbar",
6:"t",
-6:"tbar",
7:"l",
-7:"anti_l",
8:"h",
-8:"anti_h",
21:"g",
11:"e_minus",
-11:"e_plus",
12:"nu_e",
-12:"anti_nu_e",
13:"mu_minus",
-13:"mu_plus",
14:"nu_mu",
-14:"anti_nu_mu",
15:"tau_minus",
-15:"tau_plus",
16:"nu_tau",
-16:"anti_nu_tau",
17:"L_minus",
-17:"L_plus",
18:"nu_L",
-18:"anti_nu_L",
22:"gamma",
23:"Z0",
24:"W_plus",
-24:"W_minus",
25:"Higgs0",
28:"reggeon",
29:"pomeron",
32:"Z_prime0",
33:"Z_prime_prime0",
34:"W_prime_plus",
-34:"W_prime_minus",
35:"Higgs_prime0",
36:"A0",
37:"Higgs_plus",
-37:"Higgs_minus",
40:"R0",
-40:"anti_R0",
81:"specflav",
82:"rndmflav",
-82:"anti_rndmflav",
83:"phasespa",
84:"c_minushadron",
-84:"anti_c_minushadron",
85:"b_minushadron",
-85:"anti_b_minushadron",
86:"t_minushadron",
-86:"anti_t_minushadron",
89:"Wvirt_plus",
-89:"Wvirt_minus",
90:"diquark",
-90:"anti_diquark",
91:"cluster",
92:"string",
93:"indep",
94:"CMshower",
95:"SPHEaxis",
96:"THRUaxis",
97:"CLUSjet",
98:"CELLjet",
99:"table",
111:"pi0",
211:"pi_plus",
-211:"pi_minus",
210:"pi_diffr_plus",
-210:"pi_diffr_minus",
20111:"pi_2S0",
20211:"pi_2S_plus",
-20211:"pi_2S_minus",
221:"eta",
20221:"eta_2S",
331:"eta_prime",
113:"rho0",
213:"rho_plus",
-213:"rho_minus",
30113:"rho_2S0",
30213:"rho_2S_plus",
-30213:"rho_2S_minus",
40113:"rho_3S0",
40213:"rho_3S_plus",
-40213:"rho_3S_minus",
223:"omega",
30223:"omega_2S",
333:"phi",
10111:"a_00",
10211:"a_0_plus",
-10211:"a_0_minus",
10221:"f_0",
10331:"f_prime_0",
10113:"b_10",
10213:"b_1_plus",
-10213:"b_1_minus",
10223:"h_1",
10333:"h_prime_1",
20113:"a_10",
20213:"a_1_plus",
-20213:"a_1_minus",
20223:"f_1",
20333:"f_prime_1",
115:"a_20",
215:"a_2_plus",
-215:"a_2_minus",
225:"f_2",
335:"f_prime_2",
311:"K0",
-311:"anti_K0",
310:"K_S0",
130:"K_L0",
321:"K_plus",
-321:"K_minus",
313:"K_star0",
-313:"anti_K_star0",
323:"K_star_plus",
-323:"K_star_minus",
10311:"K_0_star0",
-10311:"anti_K_0_star0",
10321:"K_0_star_plus",
-10321:"K_0_star_minus",
10313:"K_10",
-10313:"anti_K_10",
10323:"K_1_plus",
-10323:"K_1_minus",
315:"K_2_star0",
-315:"anti_K_2_star0",
325:"K_2_star_plus",
-325:"K_2_star_minus",
20313:"K_prime_10",
-20313:"anti_K_prime_10",
20323:"K_prime_1_plus",
-20323:"K_prime_1_minus",
411:"D_plus",
-411:"D_minus",
421:"D0",
-421:"anti_D0",
413:"D_star_plus",
-413:"D_star_minus",
423:"D_star0",
-423:"anti_D_star0",
10411:"D_0_star_plus",
-10411:"D_0_star_minus",
10421:"D_0_star0",
-10421:"anti_D_0_star0",
10413:"D_1_plus",
-10413:"D_1_minus",
10423:"D_10",
-10423:"anti_D_10",
415:"D_2_star_plus",
-415:"D_2_star_minus",
425:"D_2_star0",
-425:"anti_D_2_star0",
20413:"D_prime_1_plus",
-20413:"D_prime_1_minus",
20423:"D_prime_10",
-20423:"anti_D_prime_10",
431:"D_s_plus",
-431:"D_s_minus",
433:"D_s_star_plus",
-433:"D_s_star_minus",
10431:"D_s0_star_plus",
-10431:"D_s0_star_minus",
10433:"D_s1_plus",
-10433:"D_s1_minus",
435:"D_s2_star_plus",
-435:"D_s2_star_minus",
20433:"D_prime_s1_plus",
-20433:"D_prime_s1_minus",
511:"B0",
-511:"anti_B0",
521:"B_plus",
-521:"B_minus",
513:"B_star0",
-513:"anti_B_star0",
523:"B_star_plus",
-523:"B_star_minus",
10511:"B_0_star0",
-10511:"anti_B_0_star0",
10521:"B_0_star_plus",
-10521:"B_0_star_minus",
10513:"B_10",
-10513:"anti_B_10",
10523:"B_1_plus",
-10523:"B_1_minus",
515:"B_2_star0",
-515:"anti_B_2_star0",
525:"B_2_star_plus",
-525:"B_2_star_minus",
20513:"B_prime_10",
-20513:"anti_B_prime_10",
20523:"B_prime_1_plus",
-20523:"B_prime_1_minus",
531:"B_s0",
-531:"anti_B_s0",
533:"B_s_star0",
-533:"anti_B_s_star0",
10531:"B_s0_star0",
-10531:"anti_B_s0_star0",
10533:"B_s10",
-10533:"anti_B_s10",
535:"B_s2_star0",
-535:"anti_B_s2_star0",
20533:"B_prime_s10",
-20533:"anti_B_prime_s10",
541:"B_c_plus",
-541:"B_c_minus",
543:"B_c_star_plus",
-543:"B_c_star_minus",
10541:"B_c0_star_plus",
-10541:"B_c0_star_minus",
10543:"B_c1_plus",
-10543:"B_c1_minus",
545:"B_c2_star_plus",
-545:"B_c2_star_minus",
20543:"B_prime_c1_plus",
-20543:"B_prime_c1_minus",
441:"eta_c",
20441:"eta_c_2S",
443:"J_psi",
20443:"psi_2S",
10441:"chi_c0",
10443:"chi_c1",
445:"chi_c2",
20551:"eta_b_2S",
40551:"eta_b_3S",
553:"Upsilon",
20553:"Upsilon_2S",
60553:"Upsilon_3S",
70553:"Upsilon_4S",
80553:"Upsilon_5S",
10553:"h_b",
40553:"h_b_2P",
100553:"h_b_3P",
551:"chi_b0",
20553:"chi_b1",
555:"chi_b2",
30551:"chi_b0_2P",
50553:"chi_b1_2P",
10555:"chi_b2_2P",
50551:"chi_b0_3P",
110553:"chi_b1_3P",
20555:"chi_b2_3P",
40555:"eta_b2_1D",
60555:"eta_b2_2D",
120553:"Upsilon_1_1D",
30555:"Upsilon_2_1D",
557:"Upsilon_3_1D",
130553:"Upsilon_1_2D",
50555:"Upsilon_2_2D",
10557:"Upsilon_3_2D",
1114:"Delta_minus",
-1114:"anti_Delta_plus",
2110:"n_diffr",
-2110:"anti_n_diffr",
2112:"n0",
-2112:"anti_n0",
2114:"Delta0",
-2114:"anti_Delta0",
2210:"p_diffr_plus",
-2210:"anti_p_diffr_minus",
2212:"p_plus",
-2212:"anti_p_minus",
2214:"Delta_plus",
-2214:"anti_Delta_minus",
2224:"Delta_plus_plus",
-2224:"anti_Delta_minus_minus",
3112:"Sigma_minus",
-3112:"anti_Sigma_plus",
3114:"Sigma_star_minus",
-3114:"anti_Sigma_star_plus",
3122:"Lambda0",
-3122:"anti_Lambda0",
3212:"Sigma0",
-3212:"anti_Sigma0",
3214:"Sigma_star0",
-3214:"anti_Sigma_star0",
3222:"Sigma_plus",
-3222:"anti_Sigma_minus",
3224:"Sigma_star_plus",
-3224:"anti_Sigma_star_minus",
3312:"Xi_minus",
-3312:"anti_Xi_plus",
3314:"Xi_star_minus",
-3314:"anti_Xi_star_plus",
3322:"Xi0",
-3322:"anti_Xi0",
3324:"Xi_star0",
-3324:"anti_Xi_star0",
3334:"Omega_minus",
-3334:"anti_Omega_plus",
4112:"Sigma_c0",
-4112:"anti_Sigma_c0",
4114:"Sigma_c_star0",
-4114:"anti_Sigma_c_star0",
4122:"Lambda_c_plus",
-4122:"anti_Lambda_c_minus",
4132:"Xi_c0",
-4132:"anti_Xi_c0",
4212:"Sigma_c_plus",
-4212:"anti_Sigma_c_minus",
4214:"Sigma_c_star_plus",
-4214:"anti_Sigma_c_star_minus",
4222:"Sigma_c_plus_plus",
-4222:"anti_Sigma_c_minus_minus",
4224:"Sigma_c_star_plus_plus",
-4224:"anti_Sigma_c_star_minus_minus",
4322:"Xi_c_plus",
-4322:"anti_Xi_c_minus",
4312:"Xi_prime_c0",
-4312:"Xi_primeanti__c0",
4314:"Xi_c_star0",
-4314:"anti_Xi_c_star0",
4232:"Xi_prime_c_plus",
-4232:"Xi_primeanti__c_minus",
4324:"Xi_c_star_plus",
-4324:"anti_Xi_c_star_minus",
4332:"Omega_c0",
-4332:"anti_Omega_c0",
4334:"Omega_c_star0",
-4334:"anti_Omega_c_star0",
5112:"Sigma_b_minus",
-5112:"anti_Sigma_b_plus",
5114:"Sigma_b_star_minus",
-5114:"anti_Sigma_b_star_plus",
5122:"Lambda_b0",
-5122:"anti_Lambda_b0",
5132:"Xi_b_minus",
-5132:"anti_Xi_b_plus",
5212:"Sigma_b0",
-5212:"anti_Sigma_b0",
5214:"Sigma_b_star0",
-5214:"anti_Sigma_b_star0",
5222:"Sigma_b_plus",
-5222:"anti_Sigma_b_minus",
5224:"Sigma_star_",
-5224:"anti_Sigma_b_star_minus",
5232:"Xi_b0",
-5232:"anti_Xi_b0",
5312:"Xi_prime_b_minus",
-5312:"anti_Xi_prime_b_plus",
5314:"Xi_b_star_minus",
-5314:"anti_Xi_b_star_plus",
5322:"Xi_prime_b0",
-5322:"anti_Xi_prime_b0",
5324:"Xi_b_star0",
-5324:"anti_Xi_b_star0",
5332:"Omega_b_minus",
-5332:"anti_Omega_b_plus",
5334:"Omega_b_star_minus",
-5334:"anti_Omega_b_star_plus",
1101:"dd_0",
-1101:"anti_dd_0",
2101:"ud_0",
-2101:"anti_ud_0",
2201:"uu_0",
-2201:"anti_uu_0",
3101:"sd_0",
-3101:"anti_sd_0",
3201:"su_0",
-3201:"anti_su_0",
3301:"ss_0",
-3301:"anti_ss_0",
4101:"cd_0",
-4101:"anti_cd_0",
4201:"cu_0",
-4201:"anti_cu_0",
4301:"cs_0",
-4301:"anti_cs_0",
4401:"cc_0",
-4401:"anti_cc_0",
5101:"bd_0",
-5101:"anti_bd_0",
5201:"bu_0",
-5201:"anti_bu_0",
5301:"bs_0",
-5301:"anti_bs_0",
5401:"bc_0",
-5401:"anti_bc_0",
5501:"bb_0",
-5501:"anti_bb_0",
1103:"dd_1",
-1103:"anti_dd_1",
2103:"ud_1",
-2103:"anti_ud_1",
2203:"uu_1",
-2203:"anti_uu_1",
3103:"sd_1",
-3103:"anti_sd_1",
3203:"su_1",
-3203:"anti_su_1",
3303:"ss_1",
-3303:"anti_ss_1",
4103:"cd_1",
-4103:"anti_cd_1",
4203:"cu_1",
-4203:"anti_cu_1",
4303:"cs_1",
-4303:"anti_cs_1",
4403:"cc_1",
-4403:"anti_cc_1",
5103:"bd_1",
-5103:"anti_bd_1",
5203:"bu_1",
-5203:"anti_bu_1",
5303:"bs_1",
-5303:"anti_bs_1",
5403:"bc_1",
-5403:"anti_bc_1",
5503:"bb_1",
-5503:"anti_bb_1",
1000011:"s_e_minus_L",
-1000011:"s_e_plus_L",
1000012:"s_nu_e_L",
-1000012:"s_anti_nu_e_L",

1000013:"s_mu_minus_L",
-1000013:"s_mu_plus_L",

1000014:"s_nu_mu_L",
-1000014:"s_anti_nu_mu_L",

1000015:"s_tau_minus_1",
-1000015:"s_tau_plus_1",

1000016:"s_nu_tau_L",
-1000016:"s_anti_nu_tau_L",

2000011:"s_e_minus_R",
-2000011:"s_e_plus_R",

2000013:"s_mu_minus_R",
-2000013:"s_mu_plus_R",

2000015:"s_tau_minus_2",
-2000015:"s_tau_plus_2",

1000021:"s_g",
1000022:"s_chi_0_1",
1000023:"s_chi_0_2",
1000024:"s_chi_plus_1",
1000025:"s_chi_0_3",
1000035:"s_chi_0_4",
1000037:"s_chi_plus_2",
-1000037:"s_chi_minus_2",

1000039:"s_G",
}

def pdgToName( pdg ):
    try:
        return pdgToName_[pdg]
    except KeyError:
        return str(pdg)
