#!/usr/bin/env python
''' Info about the process specific plot variables for ROC curves. Define them here!
'''

# Standard imports and batch mode
import numpy as np

# plot variables for ttZ
ttZ     = { 'cut':[
                   { 'index':1, 'plotstring':'p_{T}(Z)',             'var':'Z_pt',                'binning':[ 5, 0, 500 ],        'plotrange':[ 20, 0, 500],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'cos(#theta*)',         'var':'Z_cosThetaStar',      'binning':[ 5, -1.2, 1.2 ],     'plotrange':[ 20, -1.2, 1,2],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|cos(#theta*)|',       'var':'abs(Z_cosThetaStar)', 'binning':[ 5, 0, 1.2 ],        'plotrange':[ 20, 0, 1,2],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':4, 'plotstring':'#phi(Z)',              'var':'Z_phi',               'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#phi(Z)|',            'var':'abs(Z_phi)',          'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'#eta(Z)',              'var':'Z_eta',               'binning':[ 5, -3, 3 ],         'plotrange':[ 20, -3, 3 ],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'|#eta(Z)|',            'var':'abs(Z_eta)',          'binning':[ 5, 0, 3 ],          'plotrange':[ 20, 0, 3 ],          'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'GenMet_pt',           'binning':[ 5, 0, 500 ],        'plotrange':[ 20, 0, 500 ],        'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':9, 'plotstring':'#phi(E_{T}^{miss})',   'var':'GenMet_phi',          'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':10, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(GenMet_phi)',     'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44]},
                 ],
            'bin':[
                   { 'index':1, 'plotstring':'p_{T}(Z)',             'var':'Z_pt',                'binning':[ 5, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300,400], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'cos(#theta*)',         'var':'Z_cosThetaStar',      'binning':[ 5, -1.2, 1.2 ],     'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|cos(#theta*)|',       'var':'abs(Z_cosThetaStar)', 'binning':[ 5, 0, 1.2 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,0.5,0.8],         'color':[30, 41, 42, 43, 44] },
                   { 'index':4, 'plotstring':'#phi(Z)',              'var':'Z_phi',               'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#phi(Z)|',            'var':'abs(Z_phi)',          'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'#eta(Z)',              'var':'Z_eta',               'binning':[ 5, -3, 3 ],         'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'|#eta(Z)|',            'var':'abs(Z_eta)',          'binning':[ 5, 0, 3 ],          'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[1,1.5,2,2.2],       'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'GenMet_pt',           'binning':[ 5, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300],     'color':[30, 41, 42, 43, 44] },
                   { 'index':9, 'plotstring':'#phi(E_{T}^{miss})',   'var':'GenMet_phi',          'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':10, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(GenMet_phi)',     'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44]},
                  ],
}

# plot variables for ttW
ttW     = { 'cut':[
                   { 'index':1, 'plotstring':'p_{T}(W)',             'var':'W_pt',            'binning':[ 5, 0, 500 ],        'plotrange':[ 20, 0, 500],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'#phi(W)',              'var':'W_phi',           'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(W)|',            'var':'abs(W_phi)',      'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':4, 'plotstring':'#eta(W)',              'var':'W_eta',           'binning':[ 5, -3, 3 ],         'plotrange':[ 20, -3, 3 ],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(W)|',            'var':'abs(W_eta)',      'binning':[ 5, 0, 3 ],          'plotrange':[ 20, 0, 3 ],          'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'GenMet_pt',       'binning':[ 5, 0, 500 ],        'plotrange':[ 20, 0, 500 ],        'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'GenMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(GenMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44]},
                  ],
            'bin':[
                   { 'index':1, 'plotstring':'p_{T}(W)',             'var':'W_pt',            'binning':[ 5, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300,400], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'#phi(W)',              'var':'W_phi',           'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(W)|',            'var':'abs(W_phi)',      'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':4, 'plotstring':'#eta(W)',              'var':'W_eta',           'binning':[ 5, -3, 3 ],         'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(W)|',            'var':'abs(W_eta)',      'binning':[ 5, 0, 3 ],          'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[1,1.5,2,2.2],       'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'GenMet_pt',       'binning':[ 5, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300],     'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'GenMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(GenMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                  ],
}

# plot variables for ttgamma
ttgamma = { 'cut':[
                   { 'index':1, 'plotstring':'p_{T}(#gamma)',        'var':'gamma_pt',        'binning':[ 5, 0, 500 ],        'plotrange':[ 20, 0, 500],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'#phi(#gamma)',         'var':'gamma_phi',       'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(#gamma)|',       'var':'abs(gamma_phi)',  'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':4, 'plotstring':'#eta(#gamma)',         'var':'gamma_eta',       'binning':[ 5, -3, 3 ],         'plotrange':[ 20, -3, 3 ],         'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(#gamma)|',       'var':'abs(gamma_eta)',  'binning':[ 5, 0, 3 ],          'plotrange':[ 20, 0, 3 ],          'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'GenMet_pt',       'binning':[ 5, 0, 500 ],        'plotrange':[ 20, 0, 500 ],        'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'GenMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 20, -np.pi, np.pi ], 'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(GenMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 20, 0, np.pi ],      'binningPerGraph':[ 2, 5, 10, 15, 20 ], 'color':[30, 41, 42, 43, 44] },
                  ],
            'bin':[
                   { 'index':1, 'plotstring':'p_{T}(#gamma)',        'var':'gamma_pt',        'binning':[ 5, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300,400], 'color':[30, 41, 42, 43, 44] },
                   { 'index':2, 'plotstring':'#phi(#gamma)',         'var':'gamma_phi',       'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':3, 'plotstring':'|#phi(#gamma)|',       'var':'abs(gamma_phi)',  'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':4, 'plotstring':'#eta(#gamma)',         'var':'gamma_eta',       'binning':[ 5, -3, 3 ],         'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':5, 'plotstring':'|#eta(#gamma)|',       'var':'abs(gamma_eta)',  'binning':[ 5, 0, 3 ],          'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[1,1.5,2,2.2],       'color':[30, 41, 42, 43, 44] },
                   { 'index':6, 'plotstring':'p_{T}(E_{T}^{miss})',  'var':'GenMet_pt',       'binning':[ 5, 0, 500 ],        'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[0,100,200,300],     'color':[30, 41, 42, 43, 44] },
                   { 'index':7, 'plotstring':'#phi(E_{T}^{miss})',   'var':'GenMet_phi',      'binning':[ 5, -np.pi, np.pi ], 'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                   { 'index':8, 'plotstring':'|#phi(E_{T}^{miss})|', 'var':'abs(GenMet_phi)', 'binning':[ 5, 0, np.pi ],      'plotrange':[ 10, 0, 50 ],         'cutPerGraph':[],                  'color':[30, 41, 42, 43, 44] },
                  ],
}
