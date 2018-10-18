#!/usr/bin/env python
import ROOT
import os
import argparse
import shutil

from TTXPheno.Tools.user           import combineReleaseLocation, plot_directory, cardfileLocation


releaseLocation = combineReleaseLocation+"/HiggsAnalysis/CombinedLimit/"

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store', default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],    help="Log level for logging")
argParser.add_argument("--removeDir",      action='store_true',                                                             help="Remove the directory in the combine release after study is done?")
argParser.add_argument("--cores",          action='store', default=6,               nargs='?',                              help="Run on n cores in parallel")
argParser.add_argument("--cardfile",       action='store', default='',              nargs='?',                              help="which cardfile?")
args = argParser.parse_args()


# Logging
import TTXPheno.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )


def wrapper():
    logger.info("Processing impacts")
    #name = "ewkDM_ttZ_ll_DC2A_0p200000_DC2V_0p200000"
    #name = "ewkDM_ttZ_ll_DC2A_0p250000_DC2V_m0p150000"
    name = args.cardfile
    cardFile = name+".txt"
    cardFilePath = cardfileLocation + cardFile
    combineDirname = os.path.join(releaseLocation, "ttZ_YR2018")
    logger.info("Creating %s"%combineDirname)
    if not os.path.isdir(combineDirname): os.makedirs(combineDirname)
    shutil.copyfile(cardFilePath,combineDirname+'/'+cardFile)
    prepWorkspace   = "text2workspace.py %s -m 125"%cardFile
    robustFit       = "combineTool.py -M Impacts -d %s.root -m 125 --doInitialFit "%name
    impactFits      = "combineTool.py -M Impacts -d %s.root -m 125 --doFits --parallel %s "%(name,str(args.cores))
    extractImpact   = "combineTool.py -M Impacts -d %s.root -m 125 -o impacts.json"%name
    plotImpacts     = "plotImpacts.py -i impacts.json -o impacts"
    combineCommand  = "cd %s;eval `scramv1 runtime -sh`;%s;%s;%s;%s;%s"%(combineDirname,prepWorkspace,robustFit,impactFits,extractImpact,plotImpacts)
    logger.info("Will run the following command, might take a few hours:\n%s"%combineCommand)
    
    os.system(combineCommand)

    plotDir = plot_directory + "/impactsYR2018/"
    if not os.path.isdir(plotDir): os.makedirs(plotDir)
    shutil.copyfile(combineDirname+'/impacts.pdf', "%s/%s.pdf"%(plotDir,"ttZ_YR2018"))
    logger.info("Copied result to %s"%plotDir)

    if args.removeDir:
        logger.info("Removing directory in release location")
        shutil.rmtree(combineDirname)


wrapper()
