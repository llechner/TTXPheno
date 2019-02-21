'''Split HEPMC files
'''

# Standard imports
import os

# Logger
import logging
logger = logging.getLogger(__name__)

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def file_open( filename, counter ):
   
    if not filename.endswith('.hepmc'):
        raise RuntimeError( "Need .hepmc file")
    
    dirname = filename.rstrip('.hepmc')
    if not os.path.exists( dirname ):
        os.makedirs( dirname )
    
    fname = os.path.basename( dirname )

    filename = os.path.join( dirname, fname + '_' + str( counter ) + '.hepmc' ) 

    return open( filename, 'w')

def split_file( inputFile, maxEvents):
    event_counter = 0
    with open( inputFile ) as infile:

        file_counter = 0
        outfile = file_open( inputFile, file_counter)

        for iline, line in enumerate(infile):
            if iline==0:
                version = line
            if line.startswith( 'E ' ) and RepresentsInt(line.split()[1]): # new event

                if event_counter>0 and event_counter%maxEvents == 0: # need new file
                    # write last line
                    outfile.write('HepMC::IO_GenEvent-END_EVENT_LISTING\n')
                    outfile.close()
                    logger.info( "After a total of %i events wrote %i events to file %s", event_counter, maxEvents, outfile.name )
                    file_counter += 1
                    outfile = file_open( inputFile, file_counter)
                    logger.info( "Opened new file %s", outfile.name ) 
                    outfile.write(version) 
                    outfile.write('HepMC::IO_GenEvent-START_EVENT_LISTING\n') 
                event_counter += 1            

            outfile.write( line )

        outfile.write('HepMC::IO_GenEvent-END_EVENT_LISTING')
        outfile.close()
        logger.info( "After a total of %i events wrote %i events to file %s", event_counter, maxEvents, outfile )

if __name__ == '__main__':
    #
    # Arguments
    # 
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
    argParser.add_argument('--maxEvents',          action='store',      nargs='?', type=int, default=1000,  help="Maximum number of events")
#    argParser.add_argument('--inputFile',          action='store',      nargs='?', required = True)
    args = argParser.parse_args()

    from TTXPheno.samples.hepmc_samples import *

    #
    # Logger
    #
    import TTXPheno.Tools.logger as _logger
    logger = _logger.get_logger(   args.logLevel, logFile = None)

    for sample in [ttbarZ, bbar]: #ttbar
       for subsample in sample.samples:
            for filename in subsample.files:
                split_file( filename, args.maxEvents )
