#!/usr/bin/env python
''' Info about the process specific plot variables for ROC curves. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# plot variables for ttZ
ttZ     = [
#                  { 'plotstring':'p_{T}(Z) [GeV]',            'var':'Z_pt',           'binning':[ 200, 0, 1000 ],   'plotrange':[20, 0, 500] },

                  { 'plotstring':'m_{ll}',              'var':'Z_mass',         'binning':[ 20, 70, 110 ],    'plotrange':[50, 70, 110] },
                  { 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',  'binning':[ 20, 0, 500 ],    'plotrange':[50, 0, 500] },
#                  { 'plotstring':'#eta(Z)',             'var':'Z_eta',          'binning':[ 20, -3, 3 ],    'plotrange':[50, -3, 3] },
#                  { 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 20, -np.pi, np.pi ],    'plotrange':[50, -np.pi, np.pi] },

#                  { 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 200, -1.2, 1.2 ], 'plotrange':[20, -1.2, 1,2] },
#                  { 'plotstring':'|cos(#theta*)|',        'var':'abs(Z_cosThetaStar)', 'binning':[ 200, 0, 1.2 ], 'plotrange':[20, 0, 1,2] },
]

# plot variables for ttW
ttW     = [
]

# plot variables for ttgamma
ttgamma = [
]
