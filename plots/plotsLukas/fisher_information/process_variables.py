#!/usr/bin/env python
''' Info about the process specific plot variables. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# Additional plot variables for ttZ
ttZ     = { '2D':[
                  { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 10, 0, 500 ] },
                  { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 5, 0, 500 ] },
                  { 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 10, -1.2, 1.2 ] },
                  { 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 5, -1.2, 1.2 ] },
                  { 'plotstring':'m_{ll}',              'var':'Z_mass',         'binning':[ 5, 70, 110 ] },
                  { 'plotstring':'#phi(Z)',             'var':'Z_phi',          'binning':[ 5, -np.pi, np.pi ] },
                  { 'plotstring':'#eta(Z)',             'var':'Z_eta',          'binning':[ 5, -3, 3 ] },
                  { 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',      'binning':[ 5, 0, 500 ] },
                  { 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'plotstring':'p_{T}(Z) : cos(#theta*)',       'var':'Z_cosThetaStar:Z_pt', 'binning':[ 10, 0, 500, 10, -1.2, 1.2 ] },
                  { 'plotstring':'p_{T}(Z) : cos(#theta*)',       'var':'Z_cosThetaStar:Z_pt', 'binning':[ 5, 0, 500, 5, -1.2, 1.2 ] },
                  { 'plotstring':'p_{T}(Z) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:Z_pt',     'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(Z) : #phi(Z)',            'var':'Z_phi:Z_pt',          'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(Z) : #eta(Z)',            'var':'Z_eta:Z_pt',          'binning':[ 5, 0, 500, 5, -3, 3 ] },

                  { 'plotstring':'cos(#theta*) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:Z_cosThetaStar', 'binning':[ 5, -1.2, 1.2, 5, -np.pi, np.pi ] },
                  { 'plotstring':'cos(#theta*) : #phi(Z)',            'var':'Z_phi:Z_cosThetaStar',      'binning':[ 5, -1.2, 1.2, 5, -np.pi, np.pi ] },
                  { 'plotstring':'cos(#theta*) : #eta(Z)',            'var':'Z_eta:Z_cosThetaStar',      'binning':[ 5, -1.2, 1.2, 5, -3, 3 ] },

                  { 'plotstring':'#phi(E_{T}^{miss}) : #phi(Z)', 'var':'Z_phi:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
                  { 'plotstring':'#phi(E_{T}^{miss}) : #eta(Z)', 'var':'Z_eta:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}


# Additional plot variables for ttW
ttW     = { '2D':[
                  { 'plotstring':'p_{T}(W)',            'var':'W_pt',           'binning':[ 10, 0, 500 ] },
                  { 'plotstring':'p_{T}(W)',            'var':'W_pt',           'binning':[ 5, 0, 500 ] },
                  { 'plotstring':'m_{l#nu}',            'var':'W_mass',         'binning':[ 5, 60, 100 ] },
                  { 'plotstring':'#phi(W)',             'var':'W_phi',          'binning':[ 5, -np.pi, np.pi ] },
                  { 'plotstring':'#eta(W)',             'var':'W_eta',          'binning':[ 5, -3, 3 ] },
                  { 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',      'binning':[ 5, 0, 500 ] },
                  { 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'plotstring':'p_{T}(W) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:W_pt',     'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(W) : #phi(W)',            'var':'W_phi:W_pt',          'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(W) : #eta(W)',            'var':'W_eta:W_pt',          'binning':[ 5, 0, 500, 5, -3, 3 ] },

                  { 'plotstring':'#phi(E_{T}^{miss}) : #phi(W)', 'var':'W_phi:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
                  { 'plotstring':'#phi(E_{T}^{miss}) : #eta(W)', 'var':'W_eta:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}

# Additional plot variables for ttgamma
ttgamma = { '2D':[
                  { 'plotstring':'p_{T}(#gamma)',       'var':'gamma_pt',   'binning':[ 10, 0, 500 ] },
                  { 'plotstring':'p_{T}(#gamma)',       'var':'gamma_pt',   'binning':[ 5, 0, 500 ] },
                  { 'plotstring':'m_{#gamma}',          'var':'gamma_mass', 'binning':[ 5, -5, 5 ] },
                  { 'plotstring':'#phi(#gamma)',        'var':'gamma_phi',  'binning':[ 5, -np.pi, np.pi ] },
                  { 'plotstring':'#eta(#gamma)',        'var':'gamma_eta',  'binning':[ 5, -3, 3 ] },
                  { 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',  'binning':[ 5, 0, 500 ] },
                  { 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi', 'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'plotstring':'p_{T}(#gamma) : #phi(E_{T}^{miss})', 'var':'GenMet_phi:gamma_pt',     'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(#gamma) : #phi(#gamma)',       'var':'gamma_phi:gamma_pt',          'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'plotstring':'p_{T}(#gamma) : #eta(#gamma)',       'var':'gamma_eta:gamma_pt',          'binning':[ 5, 0, 500, 5, -3, 3 ] },

                  { 'plotstring':'#phi(E_{T}^{miss}) : #phi(#gamma)', 'var':'gamma_phi:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
                  { 'plotstring':'#phi(E_{T}^{miss}) : #eta(#gamma)', 'var':'gamma_eta:GenMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}
