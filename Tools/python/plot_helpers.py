''' plot helper functions
'''

# RootTools
from RootTools.core.standard             import *

def replace_selectionstrings( selectionList ):
    ''' replace selection string elements with string shown in plots
    '''

    replaced_selectionstrings = []

    for item in selectionList:

        # lepton selection
        if 'lepSel' in item:
            replaced_selectionstrings.append(item[-1] + ' lep' if item[-1]=='1' else item[-1] + ' leps')

        # number of jets/bjets/leps selection
        elif item[0] == 'n':
            num = ''.join(c for c in item if c.isdigit())
            particle = item.split('n')[1].split(num)[0]
            replaced_selectionstrings.append( ' '.join([num, particle] if num=='1' else [num, particle + 's']) )

        # pt selection
        elif 'pt' in item:
            particle = item.split('pt')[0]
            val = item.split('pt')[1]

            if 'to' in val:
                replaced_selectionstrings.append( pre + ' #leq p_{T}(' + particle + ') < ' + post + ' GeV')
            else:
               	replaced_selectionstrings.append( 'p_{T}(' + particle + ') > ' + val + ' GeV')

        # anything
        else:
            replaced_selectionstrings.append( item )

    return replaced_selectionstrings

