''' Read Simon Fernbachs HEPMC directories.
Boring but must be done.
'''

hepmc_directory = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/HEPMC"

# RootTools
from RootTools.core.standard import *

# Standard imports
import os
import re

# Logging
import logging
logger = logging.getLogger(__name__)

# Logging
if __name__ == "__main__":
    import TTXPheno.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

## HEPMC directory content:
#PP.hepmc
#1d0_10d0_HH.hepmc
#1d0_10d0_GH.hepmc
#1d0_10d0_HG.hepmc
#1d0_1d0_GH.hepmc
#1d0_1d0_HG.hepmc
#1d0_1d0_HH.hepmc
#1d0_0d1_HG.hepmc
#1d0_0d1_GH.hepmc
#1d0_0d1_HH.hepmc

#1d0_10d0_crosssections.txt
#1d0_1d0_crosssections.txt
#1d0_0d1_crosssections.txt

## xsec file:
#N   PP  GH  HG  HH
#100 0.35(3)e-03 2.3(3)e-03  26(4)e-03   1.1(1)e+00

class HEPMCData(object):

    pdfs      = ["10d0", "0d1", "1d0"]
    processes = ["GH", "HG", "HH"]

    def __init__ ( self, name, directory, isSplit = False, root_directory = None):

        # a .hepmc single file or, if I split it up, a directory of the same name
        self.samples_dict = { 'PP':\
            HEPMCSample.fromFiles( name+"_PP", os.path.join( directory, 'PP.hepmc') ) if not isSplit else  HEPMCSample.fromDirectory( name+"_PP", os.path.join( directory, 'PP') ) }
        if root_directory is not None:
            self.root_samples_dict = {'PP': Sample.fromDirectory( name+"_PP", os.path.join( root_directory, name+"_PP" ) )}
        
        self.name         = name
        for pdf in HEPMCData.pdfs:
            for process in HEPMCData.processes:
                if not isSplit:
                    self.samples_dict[pdf+'_'+process] = HEPMCSample.fromFiles( name+'_'+pdf+'_'+process, os.path.join( directory, '1d0_%s_%s.hepmc'%(pdf, process)) )
                else:
                    self.samples_dict[pdf+'_'+process] = HEPMCSample.fromDirectory( name+'_'+pdf+'_'+process, os.path.join( directory, '1d0_%s_%s'%(pdf, process)) )
                if root_directory is not None:
                    self.root_samples_dict[pdf+'_'+process] = Sample.fromDirectory( name+'_'+pdf+'_'+process, os.path.join( root_directory, '%s_%s_%s'%(name, pdf, process)) )

            nEvents, xsecs = HEPMCData.read_crosssection(  os.path.join( directory, '1d0_%s_crosssections.txt'%pdf) )
            if xsecs is not None:
                if hasattr( self.samples_dict['PP'], "xsec"):
                    if xsecs[0]!=self.samples_dict['PP'].xsec:
                        logger.waring( "Inconsistent PP cross sections! file: %s", os.path.join( directory, '1d0_%s_crosssections.txt'%pdf) )

                self.samples_dict['PP'].xsec = xsecs[0]
                self.samples_dict[pdf+"_GH"].xsec = xsecs[1]
                self.samples_dict[pdf+"_HG"].xsec = xsecs[2]
                self.samples_dict[pdf+"_HH"].xsec = xsecs[3]
                self.samples_dict['PP'].nEvents = nEvents
                self.samples_dict[pdf+"_GH"].nEvents = nEvents
                self.samples_dict[pdf+"_HG"].nEvents = nEvents
                self.samples_dict[pdf+"_HH"].nEvents = nEvents

    def __getitem__( self, key ):
        return self.samples_dict[key]

    @staticmethod
    def read_crosssection( filename ):

        regex = re.compile(r"\([0-9]*\)")

        nEvents = None
        xsecs   = None
        with open(filename) as fp:
            for i, line in enumerate(fp):
                if i == 1: #Simon writes the x-sec in the 2nd line
                    line = line.rstrip().lstrip()
                    line = re.sub( regex, '', line )
                    vals = map( float, line.split() )
                    nEvents, xsecs = vals[0], vals[1:]

        return nEvents, xsecs 

    @property
    def samples( self ):
        return self.samples_dict.values() 

    @property
    def files( self ):
        return sum( [s.files for s in self.samples], [] ) 

    @property
    def root_samples( self ):
        return self.root_samples_dict.values() 

    @property
    def root_files( self ):
        return sum( [s.files for s in self.root_samples], [] ) 
               
hepmc_1event = HEPMCData( "1event", os.path.join( hepmc_directory, "1event" ) )

root_directory = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/RunII_v01/"

#bbar   = HEPMCData( "bbbar", os.path.join( hepmc_directory, "C89_BBBAR" ) )
ttbar  = HEPMCData( "ttbar", os.path.join( hepmc_directory, "C89_TTBAR" ), isSplit = True, root_directory = root_directory)
ttbarZ = HEPMCData( "ttbarZ", os.path.join( hepmc_directory, "C89_TTBARZ" ), isSplit = True, root_directory = root_directory)
