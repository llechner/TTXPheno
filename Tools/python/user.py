import os

if os.environ['USER'] in ['schoef', 'rschoefbeck', 'schoefbeck']:
    tmp_directory       = "/afs/hephy.at/data/rschoefbeck02/TTXPheno_tmp/"
    skim_output_directory      = "/afs/hephy.at/data/rschoefbeck02/TTXPheno/skims/"
    skim_directory      = "/afs/hephy.at/data/dspitzbart01/TTXPheno/skims/"
    plot_directory      = "/afs/hephy.at/user/r/rschoefbeck/www/TTXPheno/"

if os.environ['USER'] in ['dspitzbart', 'dspitzba']:
    tmp_directory       = "/afs/hephy.at/data/dspitzbart01/Top_tmp/"
    skim_directory      = "/afs/hephy.at/data/dspitzbart01/TTXPheno/skims/"
    skim_output_directory = "/afs/hephy.at/data/dspitzbart01/TTXPheno/skims/"
    plot_directory      = "/afs/hephy.at/user/d/dspitzbart/www/TTXPheno/"
