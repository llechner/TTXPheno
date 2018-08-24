''' Profiled log-likelihood fit from card file
'''

# Standard imports
import ROOT
import os
from math import *

import random

# resources
# http://people.na.infn.it/~lista/Statistics/slides/10%20-%20roostats.pdf
# https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsTutorialsJune2013#Exercise_2_Profile_Likelihood_Ca
# https://root.cern.ch/root/html/tutorials/roofit/rf605_profilell.C.html

# Logger
import logging
logger = logging.getLogger(__name__)

# TTXPheno
from TTXPheno.Tools.user import plot_directory

def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def float_cast( string ):
    if string=='-':
        return 1.0
    else:
        return float( string )

def minimize( nll ):
    '''Minimize NLL according to 
       https://root.cern.ch/doc/v606/classRooStats_1_1ProfileLikelihoodCalculator.html#a5fb0beddc24d1bd4bf4346532a517497
    '''
    minimType = ROOT.Math.MinimizerOptions.DefaultMinimizerType()
    minimAlgo = ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()
    strategy  = ROOT.Math.MinimizerOptions.DefaultStrategy()
    level     = ROOT.Math.MinimizerOptions.DefaultPrintLevel() -1
    tolerance = ROOT.Math.MinimizerOptions.DefaultTolerance()

    minimizer = ROOT.RooMinimizer(nll)

    minimizer.setStrategy(strategy)
    minimizer.setEps(tolerance)
    minimizer.setPrintLevel(999)
    minimizer.optimizeConst(2) #to optimize likelihood calculations

    status = minimizer.minimize(minimType,minimAlgo)
    result = minimizer.save()

    status   = -1
    maxtries = 4
    for ntry in range(1, maxtries+1):
       status = minimizer.minimize(minimType,minimAlgo)
       if (status%1000 == 0): # ignore erros from Improve 
          break
       elif ntry < maxtries:
          logger.warning( "minimizer re-try %i: Doing a re-scan first", ntry)
          minimizer.minimize(minimType,"Scan")
          if (ntry == 2):
             if (strategy == 0 ):
                logger.warning( "minimizer re-try %i, trying with strategy = 1", ntry)
                minimizer.setStrategy(1)
             else: 
                ntry+=1 # skip this trial if strategy is already 1 
          if (ntry == 3):
             logger.warning("minimizer re-try %i, trying with migradimproved", ntry)
             minimType = "Minuit"
             minimAlgo = "migradimproved"

    if (status%1000 == 0): return minimizer.save()
        
class ProfiledLoglikelihoodFit:

    def __init__( self, filename):

        self.filename                    = filename
        self.modelname                   = os.path.splitext(os.path.basename(self.filename))[0]
        self.nuisance_names              = []
        self.nuisance_uncertainty_values = {}
        self.seed                        = random.randint(0,1000000)

        # read card filefile
        with open( filename, 'r' ) as file:
            for i_line, line in enumerate(file.readlines()):
                line = line.lstrip().rstrip('\n').split('#')[0]
                #print line, line.split()
                # First line with 'bin' defines the bin names, followd by observations
                if line.startswith('bin') and not hasattr(self, 'bin_names'):
                    self.bin_names = line.split()[1:]
                elif line.startswith('observation'):
                    self.observations = map( int, line.split()[1:] )
                # 2nd line starting with bin lists the process names
                elif line.startswith('bin') and hasattr(self, 'bin_names'):
                    self.bin_name_per_process = line.split()[1:] 
                elif line.startswith( 'process' ):
                    args = line.split()[1:]
                    # one line starting with processes uses the enumeration (0 is signal )
                    if all( map( isInt, args ) ):
                        self.process_number_per_bin = map( int, args )
                        # find the rate positions for the signal for each bin
                        self.signal_position_per_bin = {}
                        for i_number, number in enumerate(self.process_number_per_bin):
                            if number==0:
                                self.signal_position_per_bin[self.bin_name_per_process[i_number]] = i_number
                        if len(self.signal_position_per_bin)!=len(self.bin_names):
                            bin_names_with_signal = [ self.bin_names[i] for i in self.signal_position_per_bin.values() ] 
                            logger.warning( "Signal not present in every bin! No signal in: %s", ",".join( [ name for name in self.bin_names if name not in bin_names_with_signal ]) ) 
                    # the other line starting with process uses the process names
                    else:
                        self.process_name_per_bin   = args
                # the rate is given for each process in each bin
                elif line.startswith( 'rate '):
                    self.process_rate_bin = map(float, line.split()[1:1+len(self.bin_name_per_process)] )
                # Nuisances (need to have >=2 args and 2nd line is 'lnN') 
                elif len(line.split())>2 and line.split()[1]=='lnN':
                    args = line.split()
                    if len(args)>2:
                        self.nuisance_names.append( args[0] )
                        self.nuisance_uncertainty_values[args[0]] = map( float_cast, args[2:2+len(self.bin_name_per_process)] )
                elif len(line.split())==0 or line.startswith('imax') or line.startswith('jmax') or line.startswith('kmax') or line.startswith('---'):
                    pass
                else:
                    raise ValueError( line )

        self.nBins = len(self.bin_names)
        logger.info( "Bins: %i (%s)", self.nBins, ",".join(self.bin_names) )

        self.process_names  = list(set( self.process_name_per_bin ) )
        self.nProcesses     = len( self.process_names )
        logger.info( "Processes: %i (%s)", self.nProcesses, ",".join( self.process_names ) )

        self.nNuisances     = len(self.nuisance_names)  
        logger.info( "Nuisances: %i (%s)", self.nNuisances, ",".join( self.nuisance_names ) )

    def make_workspace( self, rmin=0, rmax=10.):
        self.ws = ROOT.RooWorkspace("workspace_%i"%self.seed)
        ROOT.SetOwnership( self.ws, False )

        # set all observations
        self.observables = ROOT.RooArgSet("observables")
        for i_bin_name, bin_name in enumerate(self.bin_names):
            obs_name = "n_obs_%s"%bin_name
            self.ws.factory( "%s[%i]" %( obs_name, self.observations[i_bin_name] ) )
            self.observables.add( self.ws.var( obs_name ) )
            logger.debug( "Set observation: n_obs_%s to %i", bin_name, self.observations[i_bin_name] )

        # nuisances
        self.nuisances = ROOT.RooArgSet("nuisances")
        self.fixed     = ROOT.RooArgSet("fixed")
        self.glob      = ROOT.RooArgSet("global")

        # signal strength modifier
        self.ws.factory( "r[1,%f,%f]"%(rmin, rmax) )
        self.poi       = ROOT.RooArgSet("poi")
        self.poi.add( self.ws.var("r") )

        for i_nuisance, nuisance in enumerate( self.nuisance_names ):

            # Global value
            self.ws.factory( "global_%s[0,-5,5]" % nuisance )
            self.ws.var("global_%s"%nuisance).setConstant(True)
            self.glob.add( self.ws.var("global_%s"%nuisance) )

            # Nuisance parameter
            self.ws.factory( "beta_%s[0,-5,5]" % nuisance ) # Nuisance parameter
            self.nuisances.add( self.ws.var( "beta_%s"%nuisance ) )

            # Gaussian constraint
            self.ws.factory( "Gaussian::constr_{nuisance}(beta_{nuisance},global_{nuisance},1)".format( nuisance = nuisance ) )

        # nominal rates
        nom_rate_names = {} 
        for i_rate, rate in enumerate(self.process_rate_bin):
            nom_rate_name = "rate_nom_%s_%s"%( self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate] ) 
            nom_rate_names[i_rate] = nom_rate_name
            self.ws.factory( "%s[%f, 0., 10**9]"%( nom_rate_name, rate ) )
            self.fixed.add( self.ws.var(nom_rate_name) )
            self.ws.var(nom_rate_name).setConstant(True)
            logger.debug( "Set %s[%f, 0., 10**9]"%( nom_rate_name, rate ) )

        # quantities for nuisances
        for i_rate, rate in enumerate(self.process_rate_bin):
            factors = []
            for i_nuisance, nuisance in enumerate( self.nuisance_names ):
                # decalare constant uncertainty (e.g. lumi 1.045)
                unc_name  = "unc_%s_%s_%s"%( nuisance, self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate] )
                unc_val   =  self.nuisance_uncertainty_values[nuisance][i_rate]
                if unc_val!=1.:                
                    self.ws.factory( "%s[%f, 0., 10**9]" % (unc_name, self.nuisance_uncertainty_values[nuisance][i_rate]) )
                    self.ws.var(unc_name).setConstant(True)
                    logger.debug( "Nuisance %s. Creating uncertainty %s position %i val %f", nuisance, unc_name, i_rate, self.nuisance_uncertainty_values[nuisance][i_rate]) 
                    # make the alpha factors. We use alpha = unc**beta where unc is the uncertainty and beta is N(global_nuisance[0,-5,5],1)
                    alpha_factor_name = "alpha_%s_%s_%s"%( nuisance, self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate] )
                    self.ws.factory( "expr::{alpha_factor_name}('pow({unc_name},{beta_name})',{unc_name},{beta_name})".format( 
                        alpha_factor_name = alpha_factor_name, 
                        unc_name = unc_name, 
                        beta_name = "beta_%s"%nuisance) )

                    factors.append( alpha_factor_name )

            # Define product of rate with systematic variations, for each yield per process per bin
            self.ws.factory( "prod::{yield_name}({factors})".format(
                yield_name = "yield_%s_%s"%(self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate]), 
                factors = ",".join( (["r"] if self.process_number_per_bin[i_rate]==0 else []) +  [ nom_rate_names[i_rate]] + factors ) ) )

        # Sum up yields for all processes in each bin
        for bin_name in self.bin_names:
            summands = []
            for i_rate, rate in enumerate(self.process_rate_bin):
                if self.bin_name_per_process[i_rate] == bin_name:
                    summands.append( "yield_%s_%s"%(self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate]) )
            self.ws.factory( "sum::{yield_bin}({summands})".format(
                yield_bin = "yield_%s" % bin_name,
                summands  = ",".join(summands)) )
            
            self.ws.factory( "Poisson::poisson_{bin_name}(n_obs_{bin_name},yield_{bin_name})".format(bin_name=bin_name) )

        self.ws.factory( "PROD::model({factors})".format(
            factors = ",".join(
                ["poisson_%s"%bin_name for bin_name in self.bin_names]      # measurements
               +["constr_%s"%nuisance for nuisance in self.nuisance_names]) # nuisances
            ))

        # print the workspace
        self.ws.Print()
      
        # Import data 
        self.data = ROOT.RooDataSet("data_%i"%self.seed, "data_%i"%self.seed, self.observables)
        ROOT.SetOwnership( self.data, False )
        self.data.add( self.observables )
        # import dataset into workspace
        getattr(self.ws, 'import')(self.data)

        # create signal+background Model Config
        self.sbModel = ROOT.RooStats.ModelConfig("sbModel")
        self.sbModel.SetWorkspace( self.ws )
        self.sbModel.SetPdf( self.ws.pdf("model") )
        self.sbModel.SetObservables( self.observables )
        self.sbModel.SetGlobalObservables( self.glob )
        self.sbModel.SetParametersOfInterest( self.poi )
        self.sbModel.SetNuisanceParameters( self.nuisances )
        self.sbModel.SetSnapshot( self.poi )
        getattr(self.ws, 'import')(self.sbModel)

        self.bModel = self.sbModel.Clone()
        self.bModel.SetName("bModel")
        self.poi.first().setVal(0)
        self.bModel.SetSnapshot( self.poi )
        getattr(self.ws, 'import')(self.bModel)

        self.ws.writeToFile('tmp/'+self.modelname+'.root', True);
#        self.ws.writeToFile(self.modelname+'.root', True);

    def calculate_limit( self, calculator = "asymptotic"):

        # asymptotic calculator
        asymptotic_calc = ROOT.RooStats.AsymptoticCalculator(self.data, self.bModel, self.sbModel)
        asymptotic_calc.SetOneSided(True) # for one-side tests (limits)
        ROOT.RooStats.AsymptoticCalculator.SetPrintLevel(-1)
        # frequentist calculator
        frequentist_calc = ROOT.RooStats.FrequentistCalculator(self.data, self.bModel, self.sbModel) #reverse, for exclusion
        frequentist_calc.SetToys(10000,5000)

        if calculator ==    'asymptotic':
            calc = ROOT.RooStats.HypoTestInverter(asymptotic_calc) # for asymptotic 
        elif calculator ==  'frequentist':
            calc = ROOT.RooStats.HypoTestInverter(frequentist_calc) # for asymptotic 
        else:
            raise NotImplementedError( "Calculator %s not known" % calculator )

        # Configure calculator
        calc.SetConfidenceLevel(0.95)
        useCLs = True
        calc.UseCLs(useCLs)
        calc.SetVerbose(False)

        # ToyMCs
        toymcs = calc.GetHypoTestCalculator().GetTestStatSampler()
        profll = ROOT.RooStats.ProfileLikelihoodTestStat(self.sbModel.GetPdf())
        if (useCLs):
            profll.SetOneSided(True)
        toymcs.SetTestStatistic(profll)
        if (not self.sbModel.GetPdf().canBeExtended()):
            toymcs.SetNEventsPerToy(1)

        # Calculate Limit
        npoints = 50  # number of points to scan
        poimin = self.poi.first().getMin()
        poimax = self.poi.first().getMax()
        calc.SetFixedScan(npoints,poimin,poimax)
        r = calc.GetInterval()
        upperLimit = r.UpperLimit()

        # make Brazilian flag plot
        c1 = ROOT.TCanvas() 
        plot = ROOT.RooStats.HypoTestInverterPlot("HTI_Result_Plot_%i"%self.seed,"HypoTest Scan Result_%i"%self.seed,r)
        plot.Draw("goff") 
        c1.SetLogy()
        directory = os.path.join( plot_directory, "limits", self.modelname )
        if not os.path.exists( directory ):
            os.makedirs( directory )
        c1.Print( os.path.join( directory, calculator+"_HTI_Result_Plot.png") )

        # Plots of test statistic 
        n = r.ArraySize()
        if (n> 0 and r.GetResult(0).GetNullDistribution() ):
           if n > 1: 
              ny = ROOT.TMath.CeilNint( sqrt(n) )
              nx = ROOT.TMath.CeilNint(float(n)/ny) 
           for i in range(n): 
              #if (n > 1) c1.cd(i+1)
              pl = plot.MakeTestStatPlot(i)
              pl.SetLogYaxis(True)
              pl.Draw("goff")
              c1.Print(os.path.join( directory, calculator + "teststat_%i.png"%i) )


        plot.IsA().Destructor( plot )
        del plot

        return {i:r.GetExpectedUpperLimit(i) for i in range(-2,3)}

    def cleanup( self ):
        for obj in [ self.data, self.ws]:
            obj.IsA().Destructor( obj )

    def likelihoodTest( self ):
        '''Make likelihood test '''
        # read PDF and data
        pdf =  self.ws.pdf("model")
        data = self.ws.data("data_%i"%self.seed)

        # Here we only test hypothesis
        self.poi.first().setVal(1)
        self.poi.first().setConstant(True)

        # get the fit parameters
        parameters = pdf.getParameters(data)
        ROOT.RooStats.RemoveConstantParameters( parameters )
        
        # create NLL
        self.nll = pdf.createNLL(data, ROOT.RooFit.CloneData(True), ROOT.RooFit.Constrain(parameters))
        ROOT.SetOwnership( self.nll, False )

        # minimize
        fitResult = minimize( self.nll )
        if fitResult is not None:
            logger.info( "Minimum NLL is %f", fitResult.minNll() )
            return fitResult.minNll()
        else:
            logger.warning( "Did not get fitResult!" )

#    def make_profiled_interval( self ):
#
#        # create signal+background Model Config
#        sbModel = ROOT.RooStats.ModelConfig("S_plus_B_model")
#        sbModel.SetWorkspace( self.ws )
#        sbModel.SetPdf( self.ws.pdf("model") )
#        sbModel.SetObservables( self.observables )
#        sbModel.SetGlobalObservables( self.glob )
#        sbModel.SetParametersOfInterest( self.poi )
#        sbModel.SetNuisanceParameters( self.nuisances )
#
#        pl = ROOT.RooStats.ProfileLikelihoodCalculator( self.data, sbModel )
#        ROOT.SetOwnership( pl, False )
#        pl.SetConfidenceLevel(0.683)
#        firstPOI = sbModel.GetParametersOfInterest().first()
#        interval = pl.GetInterval()
#        print firstPOI.GetName(), interval.LowerLimit(firstPOI), interval.UpperLimit(firstPOI)
#        plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
#        plot.Draw("")
#        filename = os.path.splitext(os.path.basename(self.filename))[0]
#        ROOT.c1.Print("/afs/hephy.at/user/r/rschoefbeck/www/etc/%s.png"%filename)

if __name__=="__main__":
    # Logger
    import TTXPheno.Tools.logger as logger
    logger    = logger.get_logger( 'DEBUG', logFile = None)
    
    import sys
    profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( sys.argv[1] )
    profiledLoglikelihoodFit.make_workspace() 

    #expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "frequentist" )
    #expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "asymptotic" )
    #for i in range( -2, 3):
    #    print "expected limit %+i sigma %f"%( i, expected_limit[i] )

    profiledLoglikelihoodFit.likelihoodTest()
