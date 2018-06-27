# Standard imports
import ROOT
from math import *

# resources
# http://people.na.infn.it/~lista/Statistics/slides/10%20-%20roostats.pdf
# https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsTutorialsJune2013#Exercise_2_Profile_Likelihood_Ca
# https://root.cern.ch/root/html/tutorials/roofit/rf605_profilell.C.html

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
                        self.process_number_per_bin = args
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
        ws = ROOT.RooWorkspace()
        ROOT.SetOwnership( ws, False )

        # set all observations
        self.observables = ROOT.RooArgSet("observables")
        for i_bin_name, bin_name in enumerate(self.bin_names):
            obs_name = "n_obs_%s"%bin_name
            ws.factory( "%s[%i]" %( obs_name, self.observations[i_bin_name] ) )
            self.observables.add( ws.var( obs_name ) )
            logger.debug( "Set observation: n_obs_%s to %i", bin_name, self.observations[i_bin_name] )

        # nuisances
        nuisances = ROOT.RooArgSet("nuisances")
        fixed     = ROOT.RooArgSet("fixed")
        poi       = ROOT.RooArgSet("poi")
        for i_nuisance, nuisance in enumerate( self.nuisance_names ):

            # Global value
            ws.factory( "global_%s[0,-5,5]" % nuisance )
            ws.var("global_%s"%nuisance).setConstant(True)
            fixed.add( ws.var("global_%s"%nuisance) )

            # Nuisance parameter
            ws.factory( "beta_%s[0,-5,5]" % nuisance ) # Nuisance parameter
            nuisances.add( ws.var( "beta_%s"%nuisance ) )

            # Gaussian constraint
            ws.factory( "Gaussian::constr_{nuisance}(beta_{nuisance},global_{nuisance},1)".format( nuisance = nuisance ) )

        # nominal rates
        nom_rate_names = {} 
        for i_rate, rate in enumerate(self.process_rate_bin):
            nom_rate_name = "rate_nom_%s_%s"%( self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate] ) 
            nom_rate_names[i_rate] = nom_rate_name
            ws.factory( "%s[%f, 0., 10**9]"%( nom_rate_name, rate ) )
            fixed.add( ws.var(nom_rate_name) )
            ws.var(nom_rate_name).setConstant(True)
            logger.debug( "Set %s[%f, 0., 10**9]"%( nom_rate_name, rate ) )

        # quantities for nuisances
        for i_rate, rate in enumerate(self.process_rate_bin):
            factors = []
            for i_nuisance, nuisance in enumerate( self.nuisance_names ):
                # decalre constant uncertainty (e.g. lumi 1.045)
                unc_name  = "unc_%s_%s_%s"%( nuisance, self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate] )
                unc_val   =  self.nuisance_uncertainty_values[nuisance][i_rate]
                if unc_val!=1.:                
                    ws.factory( "%s[%f, 0., 10**9]" % (unc_name, self.nuisance_uncertainty_values[nuisance][i_rate]) )
                    ws.var(unc_name).setConstant(True)
                    logger.debug( "Nuisance %s. Creating uncertainty %s position %i val %f", nuisance, unc_name, i_rate, self.nuisance_uncertainty_values[nuisance][i_rate]) 
                    
                    alpha_factor_name = "alpha_%s_%s_%s"%( nuisance, self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate] )
                    ws.factory( "cexpr::{alpha_factor_name}('pow({unc_name},{beta_name})',{unc_name},{beta_name})".format( 
                        alpha_factor_name = alpha_factor_name, 
                        unc_name = unc_name, 
                        beta_name = "beta_%s"%nuisance) )

                    factors.append( alpha_factor_name )

            # Define product of rate with systematic variations, for each yield per process per bin
            ws.factory( "prod::{yield_name}({factors})".format(
                yield_name = "yield_%s_%s"%(self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate]), 
                factors = ",".join( [ nom_rate_names[i_rate]] + factors ) ) )

        # Sum up yields for all processes in each bin
        for bin_name in self.bin_names:
            summands = []
            for i_rate, rate in enumerate(self.process_rate_bin):
                if self.bin_name_per_process[i_rate] == bin_name:
                    summands.append( "yield_%s_%s"%(self.bin_name_per_process[i_rate], self.process_name_per_bin[i_rate]) )
            ws.factory( "sum::{yield_bin}({summands})".format(
                yield_bin = "yield_%s" % bin_name,
                summands  = ",".join(summands)) )
            
            ws.factory( "Poisson::poisson_{bin_name}(n_obs_{bin_name},yield_{bin_name})".format(bin_name=bin_name) )

        ws.factory( "PROD::model_core({factors})".format(factors = ",".join(["poisson_%s"%bin_name for bin_name in self.bin_names]+["constr_%s"%nuisance for nuisance in self.nuisance_names])) )
        ws.Print()
      
        # Import data 
        data = ROOT.RooDataSet("data", "data", self.observables)
        ROOT.SetOwnership( data, False )
        data.add( self.observables )
        # import dataset into workspace
        getattr(ws, 'import')(data)


if __name__=="__main__":
    # Logger
    import TTXPheno.Tools.logger as logger
    logger    = logger.get_logger( 'DEBUG', logFile = None)
    
    import sys
    profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( sys.argv[1] )

    profiledLoglikelihoodFit.make_workspace() 

## full event yield
#pWs.factory( "sum::yield(nsig,nbkg)" )
## NOTE: lower-case "sum" create a function. Upper-case "SUM" would create a PDF
## Core model: Poisson probability with mean signal+bkg
## NOTE: "model_core" is a name of the PDF object
#pWs.factory( "Poisson::model_core(n,yield)" )
#
#
## model with systematics
##pWs.factory( "PROD::model(model_core,constr_lumi)" )
##pWs.factory( "PROD::model(model_core,constr_lumi,constr_efficiency)" )
#pWs.factory( "PROD::model (model_core,constr_lumi,constr_efficiency,constr_nbkg)" );
#
## print out the workspace contents
#pWs.Print()
#
#
## create set of observables (will need it for datasets and ModelConfiglater)
#pObs = pWs.var("n") # get the pointer to the observable
#obs  = ROOT.RooArgSet("observables")
#obs.add(pObs)
## create the dataset
#pObs.setVal(15) # this is your observed data: you counted eleven events
#data = ROOT.RooDataSet("data", "data", obs)
#ROOT.SetOwnership( data, False )
#data.add( obs )
## import dataset into workspace
#getattr(pWs, 'import')(data)
## create set of global observables (need to be defined as constants!)
#pWs.var("glob_lumi").setConstant(True)
#pWs.var("glob_efficiency").setConstant(True)
#pWs.var("glob_nbkg").setConstant(True)
#globalObs = ROOT.RooArgSet("global_obs")
##globalObs.add( pWs.var("glob_lumi") )
#globalObs.add( pWs.var("glob_efficiency") )
#globalObs.add( pWs.var("glob_nbkg") )
## create set of parameters of interest (POI)
#poi = ROOT.RooArgSet("poi")
#poi.add( pWs.var("xsec") )
## create set of nuisance parameters
#nuis = ROOT.RooArgSet("nuis")
#nuis.add( pWs.var("alpha_lumi") )
#nuis.add( pWs.var("beta_efficiency") )
#nuis.add( pWs.var("beta_nbkg") ) 
#
## fix all other variables in model:
## everything except observables, POI, and nuisance parameters
## must be constant
#pWs.var("lumi_nom").setConstant(True)
#pWs.var("efficiency_nom").setConstant(True)
#pWs.var("nbkg_nom").setConstant(True)
#pWs.var("lumi_kappa").setConstant(True)
#pWs.var("efficiency_kappa").setConstant(True)
#pWs.var("nbkg_kappa").setConstant(True)
#fixed = ROOT.RooArgSet("fixed")
#fixed.add( pWs.var("lumi_nom") )
#fixed.add( pWs.var("efficiency_nom") )
#fixed.add( pWs.var("nbkg_nom") )
#fixed.add( pWs.var("lumi_kappa") )
#fixed.add( pWs.var("efficiency_kappa") )
#fixed.add( pWs.var("nbkg_kappa") )
## create signal+background Model Config
#sbHypo = ROOT.RooStats.ModelConfig ("SbHypo")
#sbHypo.SetWorkspace( pWs )
#sbHypo.SetPdf( pWs.pdf("model") )
#sbHypo.SetObservables( obs )
#sbHypo.SetGlobalObservables( globalObs )
#sbHypo.SetParametersOfInterest( poi )
#sbHypo.SetNuisanceParameters( nuis )
#
#pl = ROOT.RooStats.ProfileLikelihoodCalculator( data, sbHypo )
#ROOT.SetOwnership( pl, False )
#pl.SetConfidenceLevel(0.683)
#firstPOI = sbHypo.GetParametersOfInterest().first()
#interval = pl.GetInterval()
#print firstPOI.GetName(), interval.LowerLimit(firstPOI), interval.UpperLimit(firstPOI)
#plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
#plot.Draw("")
#ROOT.c1.Print("/afs/hephy.at/user/r/rschoefbeck/www/etc/test.png")
