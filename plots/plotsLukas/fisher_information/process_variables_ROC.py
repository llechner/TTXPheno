#!/usr/bin/env python
''' Info about the process specific plot variables for ROC curves. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# plot variables for ttZ
ttZ     = { 'cut':[
                   { 'index':1, 'plotstring':'p_{T}(Z)',             'var':'genZ_pt',                'binning':[ 10, 0, 500 ],        'plotrange':[ 20, 0, 500],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'cos(#theta*)',         'var':'genZ_cosThetaStar',      'binning':[ 5, -1.2, 1.2 ],     'plotrange':[ 20, -1.2, 1,2],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|cos(#theta*)|',       'var':'abs(genZ_cosThetaStar)', 'binning':[ 5, 0, 1.2 ],        'plotrange':[ 20, 0, 1,2],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':4, 'plotstring':'#phi(Z)',              'var':'genZ_phi',               'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#phi(Z)|',            'var':'abs(genZ_phi)',          'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':6, 'plotstring':'#eta(Z)',              'var':'genZ_eta',               'binning':[ 5, -3, 3 ],         'plotrange':[ 20, -3, 3 ],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'|#eta(Z)|',            'var':'abs(genZ_eta)',          'binning':[ 5, 0, 3 ],          'plotrange':[ 20, 0, 3 ],          'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'genMet_pt',           'binning':[ 10, 0, 500 ],        'plotrange':[ 20, 0, 500 ],        'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':9, 'plotstring':'#phi(E_{T}^{miss})',   'var':'genMet_phi',          'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':10, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(genMet_phi)',     'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44]},
                 ],
            'bin':[
                   { 'index':1, 'plotstring':'p_{T}(Z)',             'var':'genZ_pt',                'binning':[ 10, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300,400], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'cos(#theta*)',         'var':'genZ_cosThetaStar',      'binning':[ 5, -1.2, 1.2 ],     'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|cos(#theta*)|',       'var':'abs(genZ_cosThetaStar)', 'binning':[ 5, 0, 1.2 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,0.5,0.8],         'color':[30, 41, 42, 43, 44] },
#                   { 'index':4, 'plotstring':'#phi(Z)',              'var':'genZ_phi',               'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#phi(Z)|',            'var':'abs(genZ_phi)',          'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
#                   { 'index':6, 'plotstring':'#eta(Z)',              'var':'genZ_eta',               'binning':[ 5, -3, 3 ],         'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'|#eta(Z)|',            'var':'abs(genZ_eta)',          'binning':[ 5, 0, 3 ],          'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[1,1.5,2,2.2],       'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'genMet_pt',           'binning':[ 10, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300],     'color':[30, 41, 42, 43, 44] },
#                   { 'index':9, 'plotstring':'#phi(E_{T}^{miss})',   'var':'genMet_phi',          'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
#                   { 'index':10, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(genMet_phi)',     'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44]},
                  ],
}

# plot variables for ttW
ttW     = { 'cut':[
                   { 'index':1, 'plotstring':'p_{T}(W)',             'var':'genW_pt',            'binning':[ 10, 0, 500 ],        'plotrange':[ 20, 0, 500],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':2, 'plotstring':'#phi(W)',              'var':'genW_phi',           'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(W)|',            'var':'abs(genW_phi)',      'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':4, 'plotstring':'#eta(W)',              'var':'genW_eta',           'binning':[ 5, -3, 3 ],         'plotrange':[ 20, -3, 3 ],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(W)|',            'var':'abs(genW_eta)',      'binning':[ 5, 0, 3 ],          'plotrange':[ 20, 0, 3 ],          'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'genMet_pt',       'binning':[ 10, 0, 500 ],        'plotrange':[ 20, 0, 500 ],        'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'genMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(genMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44]},
                  ],
            'bin':[
                   { 'index':1, 'plotstring':'p_{T}(W)',             'var':'genW_pt',            'binning':[ 10, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300,400], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':2, 'plotstring':'#phi(W)',              'var':'genW_phi',           'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(W)|',            'var':'abs(genW_phi)',      'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
#                   { 'index':4, 'plotstring':'#eta(W)',              'var':'genW_eta',           'binning':[ 5, -3, 3 ],         'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(W)|',            'var':'abs(genW_eta)',      'binning':[ 5, 0, 3 ],          'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[1,1.5,2,2.2],       'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'genMet_pt',       'binning':[ 10, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300],     'color':[30, 41, 42, 43, 44] },
#                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'genMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
#                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(genMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                  ],
}

# plot variables for ttgamma
ttgamma = { 'cut':[
                   { 'index':1, 'plotstring':'p_{T}(#gamma)',        'var':'genPhoton_pt',        'binning':[ 10, 0, 500 ],        'plotrange':[ 20, 0, 500],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':2, 'plotstring':'#phi(#gamma)',         'var':'genPhoton_phi',       'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(#gamma)|',       'var':'abs(genPhoton_phi)',  'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':4, 'plotstring':'#eta(#gamma)',         'var':'genPhoton_eta',       'binning':[ 5, -3, 3 ],         'plotrange':[ 20, -3, 3 ],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(#gamma)|',       'var':'abs(genPhoton_eta)',  'binning':[ 5, 0, 3 ],          'plotrange':[ 20, 0, 3 ],          'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'genMet_pt',       'binning':[ 10, 0, 500 ],        'plotrange':[ 20, 0, 500 ],        'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'genMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(genMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                  ],
            'bin':[
                   { 'index':1, 'plotstring':'p_{T}(#gamma)',        'var':'genPhoton_pt',        'binning':[ 10, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300,400], 'color':[30, 41, 42, 43, 44] },
#                   { 'index':2, 'plotstring':'#phi(#gamma)',         'var':'genPhoton_phi',       'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(#gamma)|',       'var':'abs(genPhoton_phi)',  'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
#                   { 'index':4, 'plotstring':'#eta(#gamma)',         'var':'genPhoton_eta',       'binning':[ 5, -3, 3 ],         'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(#gamma)|',       'var':'abs(genPhoton_eta)',  'binning':[ 5, 0, 3 ],          'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[1,1.5,2,2.2],       'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'genMet_pt',       'binning':[ 10, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300],     'color':[30, 41, 42, 43, 44] },
#                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'genMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
#                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(genMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                  ],
}
