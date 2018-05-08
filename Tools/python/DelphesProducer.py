''' Runs delphes on an edm root file
See https://twiki.cern.ch/twiki/bin/viewauth/CMS/DelphesUPG
'''

# Standard imports
import os
import subprocess

# Logger
import logging
logger = logging.getLogger(__name__)

# cmd line syntax:
#./DelphesCMSFWLite cards/delphes_card_CMS.tcl out.root /afs/hephy.at/data/rschoefbeck02/TTXPheno/small/ttZ0j_rwgt_patch_625_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball_small.root

class DelphesProducer:

    def __init__(self, card = 'cards/delphes_card_CMS.tcl'):
        self.delphes_dir = os.path.expandvars( '$CMSSW_BASE/../delphes' )
        self.card = card

    def produce( self, infiles, outfile ):

        logger.info( "Running Delphes on %i files using %s", len(infiles), self.card )

        # Clean output
        tmp_files = []
        if os.path.exists( outfile ):
            logger.warning( "Found output file %s. Deleting", outfile )
            os.remove( outfile )

        # Create output directory if necessary
        if not os.path.exists( os.path.dirname(outfile) ):
            os.makedirs(os.path.dirname(outfile))

        # process files individually
        for i_infile, infile in enumerate(infiles):
            tmp_file = infile.rstrip('.root')+'_tmp_%i.root'%i_infile
            if os.path.exists( tmp_file ):
                logger.warning( "Found output file %s. Deleting", tmp_file )
                os.remove( tmp_file )
            tmp_files.append( tmp_file )
            subprocess.check_call( [ os.path.join( self.delphes_dir, 'DelphesCMSFWLite'), os.path.join( self.delphes_dir, self.card), tmp_file, infile ] )

        # Hadd
        subprocess.check_call( ['hadd', outfile] + tmp_files ) 

        # Clean up 
        for tmp_file in tmp_files:
            os.remove(tmp_file)
