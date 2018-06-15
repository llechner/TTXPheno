#!/usr/bin/env python
''' Info about the process specific plot variables. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# Additional plot variables for ttZ
ttZ     = [
#           { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 100000, 0, 1000 ] },
#           { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 20000, 0, 1000 ] },
#           { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 2000, 0, 1000 ] },
#           { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 200, 0, 1000 ] },
           { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 20, 0, 1000 ] },
           { 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 20, -1.2, 1.2 ] },
           { 'plotstring':'m_{ll}',              'var':'Z_mass',         'binning':[ 20, 70, 110 ] },
           { 'plotstring':'#phi(Z)',             'var':'Z_phi',          'binning':[ 20, -np.pi, np.pi ] },
           { 'plotstring':'#eta(Z)',             'var':'Z_eta',          'binning':[ 20, -3, 3 ] },
           { 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',      'binning':[ 20, 0, 1000 ] },
           { 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 20, -np.pi, np.pi ] },
]


# Additional plot variables for ttW
ttW     = [
]

# Additional plot variables for ttgamma
ttgamma = [
]
