#!/usr/bin/env python
''' Info about the process specific plot variables. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# Additional plot variables for ttZ
ttZ     = { '2D':[
                  { 'index':1, 'plotstring':'p_{T}(Z)',            'var':'genZ_pt',           'binning':[ 10, 0, 500 ] },
#                  { 'index':2, 'plotstring':'p_{T}(Z)',            'var':'genZ_pt',           'binning':[ 5, 0, 500 ] },
                  { 'index':3, 'plotstring':'cos(#theta*)',        'var':'genZ_cosThetaStar', 'binning':[ 10, -1.2, 1.2 ] },
#                  { 'index':4, 'plotstring':'cos(#theta*)',        'var':'genZ_cosThetaStar', 'binning':[ 5, -1.2, 1.2 ] },
#                  { 'index':5, 'plotstring':'m_{ll}',              'var':'genZ_mass',         'binning':[ 5, 70, 110 ] },
                  { 'index':6, 'plotstring':'#phi(Z)',             'var':'genZ_phi',          'binning':[ 5, -np.pi, np.pi ] },
                  { 'index':7, 'plotstring':'#eta(Z)',             'var':'genZ_eta',          'binning':[ 5, -3, 3 ] },
                  { 'index':8, 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'genMet_pt',      'binning':[ 5, 0, 500 ] },
                  { 'index':9, 'plotstring':'#phi(E_{T}^{miss})',  'var':'genMet_phi',     'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'index':1, 'plotstring':'p_{T}(Z) : cos(#theta*)',       'var':'genZ_cosThetaStar:genZ_pt', 'binning':[ 10, 0, 500, 10, -1.2, 1.2 ] },
#                  { 'index':2, 'plotstring':'p_{T}(Z) : cos(#theta*)',       'var':'genZ_cosThetaStar:genZ_pt', 'binning':[ 5, 0, 500, 5, -1.2, 1.2 ] },
                  { 'index':3, 'plotstring':'p_{T}(Z) : #phi(E_{T}^{miss})', 'var':'genMet_phi:genZ_pt',     'binning':[ 10, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':4, 'plotstring':'p_{T}(Z) : #phi(Z)',            'var':'genZ_phi:genZ_pt',          'binning':[ 10, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':5, 'plotstring':'p_{T}(Z) : #eta(Z)',            'var':'genZ_eta:genZ_pt',          'binning':[ 10, 0, 500, 5, -3, 3 ] },

                  { 'index':6, 'plotstring':'cos(#theta*) : #phi(E_{T}^{miss})', 'var':'genMet_phi:genZ_cosThetaStar', 'binning':[ 10, -1.2, 1.2, 5, -np.pi, np.pi ] },
                  { 'index':7, 'plotstring':'cos(#theta*) : #phi(Z)',            'var':'genZ_phi:genZ_cosThetaStar',      'binning':[ 10, -1.2, 1.2, 5, -np.pi, np.pi ] },
                  { 'index':8, 'plotstring':'cos(#theta*) : #eta(Z)',            'var':'genZ_eta:genZ_cosThetaStar',      'binning':[ 10, -1.2, 1.2, 5, -3, 3 ] },

#                  { 'index':9, 'plotstring':'#phi(E_{T}^{miss}) : #phi(Z)', 'var':'genZ_phi:genMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
#                  { 'index':10, 'plotstring':'#phi(E_{T}^{miss}) : #eta(Z)', 'var':'genZ_eta:genMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}


# Additional plot variables for ttW
ttW     = { '2D':[
                  { 'index':1, 'plotstring':'p_{T}(W)',            'var':'genW_pt',           'binning':[ 10, 0, 500 ] },
                  { 'index':2, 'plotstring':'p_{T}(W)',            'var':'genW_pt',           'binning':[ 5, 0, 500 ] },
                  { 'index':3, 'plotstring':'m_{l#nu}',            'var':'genW_mass',         'binning':[ 5, 60, 100 ] },
                  { 'index':4, 'plotstring':'#phi(W)',             'var':'genW_phi',          'binning':[ 5, -np.pi, np.pi ] },
                  { 'index':5, 'plotstring':'#eta(W)',             'var':'genW_eta',          'binning':[ 5, -3, 3 ] },
                  { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'genMet_pt',      'binning':[ 5, 0, 500 ] },
                  { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',  'var':'genMet_phi',     'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'index':1, 'plotstring':'p_{T}(W) : #phi(E_{T}^{miss})', 'var':'genMet_phi:genW_pt',     'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':2, 'plotstring':'p_{T}(W) : #phi(W)',            'var':'genW_phi:genW_pt',          'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':3, 'plotstring':'p_{T}(W) : #eta(W)',            'var':'genW_eta:genW_pt',          'binning':[ 5, 0, 500, 5, -3, 3 ] },

                  { 'index':4, 'plotstring':'#phi(E_{T}^{miss}) : #phi(W)', 'var':'genW_phi:genMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
                  { 'index':5, 'plotstring':'#phi(E_{T}^{miss}) : #eta(W)', 'var':'genW_eta:genMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}

# Additional plot variables for ttgamma
ttgamma = { '2D':[
                  { 'index':1, 'plotstring':'p_{T}(#gamma)',       'var':'genPhoton_pt',   'binning':[ 10, 0, 500 ] },
                  { 'index':2, 'plotstring':'p_{T}(#gamma)',       'var':'genPhoton_pt',   'binning':[ 5, 0, 500 ] },
                  { 'index':3, 'plotstring':'m_{#gamma}',          'var':'genPhoton_mass', 'binning':[ 5, -5, 5 ] },
                  { 'index':4, 'plotstring':'#phi(#gamma)',        'var':'genPhoton_phi',  'binning':[ 5, -np.pi, np.pi ] },
                  { 'index':5, 'plotstring':'#eta(#gamma)',        'var':'genPhoton_eta',  'binning':[ 5, -3, 3 ] },
                  { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'genMet_pt',  'binning':[ 5, 0, 500 ] },
                  { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',  'var':'genMet_phi', 'binning':[ 5, -np.pi, np.pi ] },
                 ],
            '3D':[
                  { 'index':1, 'plotstring':'p_{T}(#gamma) : #phi(E_{T}^{miss})', 'var':'genMet_phi:genPhoton_pt',     'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':2, 'plotstring':'p_{T}(#gamma) : #phi(#gamma)',       'var':'genPhoton_phi:genPhoton_pt',          'binning':[ 5, 0, 500, 5, -np.pi, np.pi ] },
                  { 'index':3, 'plotstring':'p_{T}(#gamma) : #eta(#gamma)',       'var':'genPhoton_eta:genPhoton_pt',          'binning':[ 5, 0, 500, 5, -3, 3 ] },

                  { 'index':4, 'plotstring':'#phi(E_{T}^{miss}) : #phi(#gamma)', 'var':'genPhoton_phi:genMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -np.pi, np.pi  ] },
                  { 'index':5, 'plotstring':'#phi(E_{T}^{miss}) : #eta(#gamma)', 'var':'genPhoton_eta:genMet_phi', 'binning':[ 5, -np.pi, np.pi, 5, -3, 3 ] },
                 ],
}
