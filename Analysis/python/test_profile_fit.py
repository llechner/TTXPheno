# Standard imports
import ROOT

# resources
# http://people.na.infn.it/~lista/Statistics/slides/10%20-%20roostats.pdf
# https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsTutorialsJune2013#Exercise_2_Profile_Likelihood_Ca
# https://root.cern.ch/root/html/tutorials/roofit/rf605_profilell.C.html
pWs = ROOT.RooWorkspace()
ROOT.SetOwnership( pWs, False )
# observable: number of events
pWs.factory( "n[10]" )

# signal yield
#pWs.factory( "nsig[0,0,100]" )

# integrated luminosity
#pWs.factory( "lumi[0]" )
# integrated luminosity with systematics
pWs.factory( "lumi_nom[5000.0, 4000.0, 6000.0]" )
pWs.factory( "lumi_kappa[1.045]" )
pWs.factory( "cexpr::alpha_lumi('pow(lumi_kappa,beta_lumi)',lumi_kappa,beta_lumi[0,-5,5])" )
#pWs.factory( "pow::alpha_lumi(lumi_kappa,beta_lumi[0,-5,5])")
pWs.factory( "prod::lumi(lumi_nom,alpha_lumi)" )
pWs.factory( "Gaussian::constr_lumi(beta_lumi,glob_lumi[0,-5,5],1)" ) 

# cross section - parameter of interest
pWs.factory( "xsec[0,0,0.1]" )
# selection efficiency * acceptance
#pWs.factory( "efficiency[0]" )
# selection efficiency * acceptance with systematics
pWs.factory( "efficiency_nom[0.1, 0.05, 0.15]" )
pWs.factory( "efficiency_kappa[1.10]" )
pWs.factory( "cexpr::alpha_efficiency('pow (efficiency_kappa,beta_efficiency)', efficiency_kappa,beta_efficiency[0,-5,5])" )
pWs.factory( "prod::efficiency(efficiency_nom,alpha_efficiency)" )
pWs.factory( "Gaussian::constr_efficiency(beta_efficiency,glob_efficiency[0,-5,5],1)" ) 

# signal yield
pWs.factory( "prod::nsig(lumi,xsec,efficiency)" )

# NOTE: three parameters are "current value", "low bound", "upper bound
# background yield
#pWs.factory( "nbkg[10,0,100]" )
#pWs.factory( "nbkg_nom[10]" )

# background yield with systematics
pWs.factory( "nbkg_nom[10.0, 5.0, 15.0]" )
pWs.factory( "nbkg_kappa[1.10]" )
pWs.factory( "cexpr::alpha_nbkg('pow (nbkg_kappa,beta_nbkg)',nbkg_kappa,beta_nbkg[0,-5,5])" )
pWs.factory( "prod::nbkg(nbkg_nom,alpha_lumi,alpha_nbkg)" )
pWs.factory( "Gaussian::constr_nbkg(beta_nbkg,glob_nbkg[0,-5,5],1)" )

# full event yield
pWs.factory( "sum::yield(nsig,nbkg)" )
# NOTE: lower-case "sum" create a function. Upper-case "SUM" would create a PDF
# Core model: Poisson probability with mean signal+bkg
# NOTE: "model_core" is a name of the PDF object
pWs.factory( "Poisson::model_core(n,yield)" )


# model with systematics
#pWs.factory( "PROD::model(model_core,constr_lumi)" )
#pWs.factory( "PROD::model(model_core,constr_lumi,constr_efficiency)" )
pWs.factory( "PROD::model (model_core,constr_lumi,constr_efficiency,constr_nbkg)" );

# print out the workspace contents
pWs.Print()

# create set of observables (will need it for datasets and ModelConfiglater)
pObs = pWs.var("n") # get the pointer to the observable
obs  = ROOT.RooArgSet("observables")
obs.add(pObs)
# create the dataset
pObs.setVal(15) # this is your observed data: you counted eleven events
data = ROOT.RooDataSet("data", "data", obs)
ROOT.SetOwnership( data, False )
data.add( obs )
# import dataset into workspace
getattr(pWs, 'import')(data)
# create set of global observables (need to be defined as constants!)
pWs.var("glob_lumi").setConstant(True)
pWs.var("glob_efficiency").setConstant(True)
pWs.var("glob_nbkg").setConstant(True)
globalObs = ROOT.RooArgSet("global_obs")
globalObs.add( pWs.var("glob_lumi") )
globalObs.add( pWs.var("glob_efficiency") )
globalObs.add( pWs.var("glob_nbkg") )
# create set of parameters of interest (POI)
poi = ROOT.RooArgSet("poi")
poi.add( pWs.var("xsec") )
# create set of nuisance parameters
nuis = ROOT.RooArgSet("nuis")
nuis.add( pWs.var("beta_lumi") )
nuis.add( pWs.var("beta_efficiency") )
nuis.add( pWs.var("beta_nbkg") ) 

# fix all other variables in model:
# everything except observables, POI, and nuisance parameters
# must be constant
pWs.var("lumi_nom").setConstant(True)
pWs.var("efficiency_nom").setConstant(True)
pWs.var("nbkg_nom").setConstant(True)
pWs.var("lumi_kappa").setConstant(True)
pWs.var("efficiency_kappa").setConstant(True)
pWs.var("nbkg_kappa").setConstant(True)
fixed = ROOT.RooArgSet("fixed")
fixed.add( pWs.var("lumi_nom") )
fixed.add( pWs.var("efficiency_nom") )
fixed.add( pWs.var("nbkg_nom") )
fixed.add( pWs.var("lumi_kappa") )
fixed.add( pWs.var("efficiency_kappa") )
fixed.add( pWs.var("nbkg_kappa") )
# create signal+background Model Config
sbHypo = ROOT.RooStats.ModelConfig ("SbHypo")
sbHypo.SetWorkspace( pWs )
sbHypo.SetPdf( pWs.pdf("model") )
sbHypo.SetObservables( obs )
sbHypo.SetGlobalObservables( globalObs )
sbHypo.SetParametersOfInterest( poi )
sbHypo.SetNuisanceParameters( nuis )

pl = ROOT.RooStats.ProfileLikelihoodCalculator( data, sbHypo )
ROOT.SetOwnership( pl, False )
pl.SetConfidenceLevel(0.683)
firstPOI = sbHypo.GetParametersOfInterest().first()
interval = pl.GetInterval()
print firstPOI.GetName(), interval.LowerLimit(firstPOI), interval.UpperLimit(firstPOI)
plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
plot.Draw("")
ROOT.c1.Print("/afs/hephy.at/user/r/rschoefbeck/www/etc/test.png")

