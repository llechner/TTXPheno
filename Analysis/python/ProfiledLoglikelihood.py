''' Profiled log-likelihood fit from card file
'''

# Standard imports
import ROOT
from math import *

# resources
# http://people.na.infn.it/~lista/Statistics/slides/10%20-%20roostats.pdf
# https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsTutorialsJune2013#Exercise_2_Profile_Likelihood_Ca
# https://root.cern.ch/root/html/tutorials/roofit/rf605_profilell.C.html

# Logger
import logging
logger = logging.getLogger(__name__)

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
        
class ProfiledLoglikelihoodFit:

    def __init__( self, filename):

        self.nuisance_names              = []
        self.nuisance_uncertainty_values = {}

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

        self.process_names  = list(set( self.bin_name_per_process ) )
        self.nProcesses     = len( self.process_names )
        logger.info( "Processes: %i (%s)", self.nProcesses, ",".join( self.process_names ) )

        self.nNuisances     = len(self.nuisance_names)  
        logger.info( "Nuisances: %i (%s)", self.nNuisances, ",".join( self.nuisance_names ) )

    def make_workspace( self ):
        self.ws = ROOT.RooWorkspace()
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
        self.ws.factory( "r[1,0,10]" )
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
                    self.ws.factory( "cexpr::{alpha_factor_name}('pow({unc_name},{beta_name})',{unc_name},{beta_name})".format( 
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
        self.data = ROOT.RooDataSet("data", "data", self.observables)
        ROOT.SetOwnership( self.data, False )
        self.data.add( self.observables )
        # import dataset into workspace
        getattr(self.ws, 'import')(self.data)

    def make_profiled_interval( self ):

        # create signal+background Model Config
        sbModel = ROOT.RooStats.ModelConfig("S_plus_B_model")
        sbModel.SetWorkspace( self.ws )
        sbModel.SetPdf( self.ws.pdf("model") )
        sbModel.SetObservables( self.observables )
        sbModel.SetGlobalObservables( self.glob )
        sbModel.SetParametersOfInterest( self.poi )
        sbModel.SetNuisanceParameters( self.nuisances )

        pl = ROOT.RooStats.ProfileLikelihoodCalculator( self.data, sbModel )
        ROOT.SetOwnership( pl, False )
        pl.SetConfidenceLevel(0.683)
        firstPOI = sbModel.GetParametersOfInterest().first()
        interval = pl.GetInterval()
        print firstPOI.GetName(), interval.LowerLimit(firstPOI), interval.UpperLimit(firstPOI)
        plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
        plot.Draw("")
        ROOT.c1.Print("/afs/hephy.at/user/r/rschoefbeck/www/etc/test2.png")

    def make_hypothesis_test( self, r):

        self.poi.first().setVal( r )
        # create signal+background Model Config
        sbModel = ROOT.RooStats.ModelConfig("sbModel")
        sbModel.SetWorkspace( self.ws )
        sbModel.SetPdf( self.ws.pdf("model") )
        sbModel.SetObservables( self.observables )
        sbModel.SetGlobalObservables( self.glob )
        sbModel.SetParametersOfInterest( self.poi )
        sbModel.SetNuisanceParameters( self.nuisances )

        sbModel.SetSnapshot( self.poi )

        bModel = sbModel.Clone()
        bModel.SetName("bHypo")      
        self.poi.first().setVal(0)
        bModel.SetSnapshot( self.poi );

        acl = ROOT.RooStats.AsymptoticCalculator( self.data, sbModel, bModel )
        acl.SetOneSidedDiscovery(True)

        asResult = acl.GetHypoTest()
        asResult.Print()

if __name__=="__main__":
    # Logger
    import TTXPheno.Tools.logger as logger
    logger    = logger.get_logger( 'DEBUG', logFile = None)
    
    import sys
    profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( sys.argv[1] )
    profiledLoglikelihoodFit.make_workspace() 
    #profiledLoglikelihoodFit.make_profiled_interval() 
    profiledLoglikelihoodFit.make_hypothesis_test(r=6) 
