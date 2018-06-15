#!/usr/bin/env python
''' Info about the process specific plot variables. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# Additional plot variables for ttZ
ttZ     = { '2D':[
#                  { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 100000, 0, 1000 ] },
#                  { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 20000, 0, 1000 ] },
#                  { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 2000, 0, 1000 ] },
                  { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 200, 0, 500 ] },
                  { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 20, 0, 500 ] },
                  { 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 200, -1.2, 1.2 ] },
                  { 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 20, -1.2, 1.2 ] },
                  { 'plotstring':'m_{ll}',              'var':'Z_mass',         'binning':[ 20, 70, 110 ] },
                  { 'plotstring':'#phi(Z)',             'var':'Z_phi',          'binning':[ 20, -np.pi, np.pi ] },
                  { 'plotstring':'#eta(Z)',             'var':'Z_eta',          'binning':[ 20, -3, 3 ] },
                  { 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',      'binning':[ 20, 0, 500 ] },
                  { 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 20, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'plotstring':'p_{T}(Z) : cos(#theta*)',  'var':'Z_cosThetaStar:Z_pt', 'binning':[ 200, 0, 500, 200, -1.2, 1.2 ] },

                  { 'plotstring':'p_{T}(Z) : cos(#theta*)',        'var':'Z_cosThetaStar:Z_pt', 'binning':[ 20, 0, 500, 20, -1.2, 1.2 ] },
                  { 'plotstring':'p_{T}(Z) : #phi(E_{T}^{miss})',  'var':'GenMet_phi:Z_pt',     'binning':[ 20, 0, 500, 20, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(Z) : #phi(Z)',             'var':'Z_phi:Z_pt',          'binning':[ 20, 0, 500, 20, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(Z) : #eta(Z)',             'var':'Z_eta:Z_pt',          'binning':[ 20, 0, 500, 20, -3, 3 ] },

                  { 'plotstring':'cos(#theta*) : #phi(E_{T}^{miss})',  'var':'GenMet_phi:Z_cosThetaStar',     'binning':[ 20, -1.2, 1.2, 20, -np.pi, np.pi ] },
                  { 'plotstring':'cos(#theta*) : #phi(Z)',             'var':'Z_phi:Z_cosThetaStar',          'binning':[ 20, -1.2, 1.2, 20, -np.pi, np.pi ] },
                  { 'plotstring':'cos(#theta*) : #eta(Z)',             'var':'Z_eta:Z_cosThetaStar',          'binning':[ 20, -1.2, 1.2, 20, -3, 3 ] },

                  { 'plotstring':'#phi(E_{T}^{miss}) : #phi(Z)',             'var':'Z_phi:GenMet_phi',          'binning':[ 20, -np.pi, np.pi, 20, -np.pi, np.pi  ] },
                  { 'plotstring':'#phi(E_{T}^{miss}) : #eta(Z)',             'var':'Z_eta:GenMet_phi',          'binning':[ 20, -np.pi, np.pi, 20, -3, 3 ] },

                 ],
}


# Additional plot variables for ttW
ttW     = { '2D':[

                 ],
            '3D':[

                 ],
}

# Additional plot variables for ttgamma
ttgamma = { '2D':[

                 ],
            '3D':[

                 ],
}
