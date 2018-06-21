#!/usr/bin/env python
''' Info about the process specific plot variables. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# Additional plot variables for ttZ
ttZ     = { '2D':[
                  { 'index':1, 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 10, 0, 500 ] },
#                  { 'index':2, 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 5, 0, 500 ] },
                  { 'index':3, 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 10, -1.2, 1.2 ] },
#                  { 'index':4, 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 5, -1.2, 1.2 ] },
#                  { 'index':5, 'plotstring':'m_{ll}',              'var':'Z_mass',         'binning':[ 5, 70, 110 ] },
                  { 'index':6, 'plotstring':'#phi(Z)',             'var':'Z_phi',          'binning':[ 5, -np.pi, np.pi ] },
                  { 'index':7, 'plotstring':'#eta(Z)',             'var':'Z_eta',          'binning':[ 5, -3, 3 ] },
                  { 'index':8, 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',      'binning':[ 5, 0, 500 ] },
                  { 'index':9, 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'index':1, 'plotstring':'p_{T}(Z) : cos(#theta*)',       'var':'Z_cosThetaStar:Z_pt', 'binning':[ 10, 0, 500, 10, -1.2, 1.2 ] },
#                  { 'index':2, 'plotstring':'p_{T}(Z) : cos(#theta*)',       'var':'Z_cosThetaStar:Z_pt', 'binning':[ 5, 0, 500, 5, -1.2, 1.2 ] },
                  { 'index':3, 'plotstring':'p_{T}(Z) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:Z_pt',     'binning':[ 10, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':4, 'plotstring':'p_{T}(Z) : #phi(Z)',            'var':'Z_phi:Z_pt',          'binning':[ 10, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':5, 'plotstring':'p_{T}(Z) : #eta(Z)',            'var':'Z_eta:Z_pt',          'binning':[ 10, 0, 500, 5, -3, 3 ] },

                  { 'index':6, 'plotstring':'cos(#theta*) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:Z_cosThetaStar', 'binning':[ 10, -1.2, 1.2, 5, -np.pi, np.pi ] },
                  { 'index':7, 'plotstring':'cos(#theta*) : #phi(Z)',            'var':'Z_phi:Z_cosThetaStar',      'binning':[ 10, -1.2, 1.2, 5, -np.pi, np.pi ] },
                  { 'index':8, 'plotstring':'cos(#theta*) : #eta(Z)',            'var':'Z_eta:Z_cosThetaStar',      'binning':[ 10, -1.2, 1.2, 5, -3, 3 ] },

#                  { 'index':9, 'plotstring':'#phi(E_{T}^{miss}) : #phi(Z)', 'var':'Z_phi:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
#                  { 'index':10, 'plotstring':'#phi(E_{T}^{miss}) : #eta(Z)', 'var':'Z_eta:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}


# Additional plot variables for ttW
ttW     = { '2D':[
                  { 'index':1, 'plotstring':'p_{T}(W)',            'var':'W_pt',           'binning':[ 10, 0, 500 ] },
                  { 'index':2, 'plotstring':'p_{T}(W)',            'var':'W_pt',           'binning':[ 5, 0, 500 ] },
                  { 'index':3, 'plotstring':'m_{l#nu}',            'var':'W_mass',         'binning':[ 5, 60, 100 ] },
                  { 'index':4, 'plotstring':'#phi(W)',             'var':'W_phi',          'binning':[ 5, -np.pi, np.pi ] },
                  { 'index':5, 'plotstring':'#eta(W)',             'var':'W_eta',          'binning':[ 5, -3, 3 ] },
                  { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',      'binning':[ 5, 0, 500 ] },
                  { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'index':1, 'plotstring':'p_{T}(W) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:W_pt',     'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':2, 'plotstring':'p_{T}(W) : #phi(W)',            'var':'W_phi:W_pt',          'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':3, 'plotstring':'p_{T}(W) : #eta(W)',            'var':'W_eta:W_pt',          'binning':[ 5, 0, 500, 5, -3, 3 ] },

                  { 'index':4, 'plotstring':'#phi(E_{T}^{miss}) : #phi(W)', 'var':'W_phi:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
                  { 'index':5, 'plotstring':'#phi(E_{T}^{miss}) : #eta(W)', 'var':'W_eta:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}

# Additional plot variables for ttgamma
ttgamma = { '2D':[
                  { 'index':1, 'plotstring':'p_{T}(#gamma)',       'var':'gamma_pt',   'binning':[ 10, 0, 500 ] },
                  { 'index':2, 'plotstring':'p_{T}(#gamma)',       'var':'gamma_pt',   'binning':[ 5, 0, 500 ] },
                  { 'index':3, 'plotstring':'m_{#gamma}',          'var':'gamma_mass', 'binning':[ 5, -5, 5 ] },
                  { 'index':4, 'plotstring':'#phi(#gamma)',        'var':'gamma_phi',  'binning':[ 5, -np.pi, np.pi ] },
                  { 'index':5, 'plotstring':'#eta(#gamma)',        'var':'gamma_eta',  'binning':[ 5, -3, 3 ] },
                  { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',  'binning':[ 5, 0, 500 ] },
                  { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi', 'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'index':1, 'plotstring':'p_{T}(#gamma) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:gamma_pt',     'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':2, 'plotstring':'p_{T}(#gamma) : #phi(#gamma)',       'var':'gamma_phi:gamma_pt',          'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':3, 'plotstring':'p_{T}(#gamma) : #eta(#gamma)',       'var':'gamma_eta:gamma_pt',          'binning':[ 5, 0, 500, 5, -3, 3 ] },

                  { 'index':4, 'plotstring':'#phi(E_{T}^{miss}) : #phi(#gamma)', 'var':'gamma_phi:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
                  { 'index':5, 'plotstring':'#phi(E_{T}^{miss}) : #eta(#gamma)', 'var':'gamma_eta:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}
