import os

def getUncertaintyValue( cardFile, binNumber, process, uncertainty ):

    with open(cardFile, 'r') as f:
        cardFileLines = f.readlines()

    count = 0
    for i, line in enumerate(cardFileLines):
        # read each line of cardfile
        entry = [ item.split('\n')[0] for item in line.split(' ') if item != '' ]
        # rm lnN to get the right column number for the headers and the uncertainties
        if 'lnN' in entry: entry.remove('lnN')
        # 2nd entry of 'bin' is interesting
        if entry[0] == 'bin' and count == 0: count += 1
            # 2nd entry of 'bin' is interesting
        elif entry[0] == 'bin' and count == 1:
            # get the column numbers of the right regions
            columnNumbers = [ j for j, bin in enumerate(entry) if bin == 'Bin' + str(binNumber) or bin == 'Region_' + str(binNumber) ]
            if len(columnNumbers) == 0:
                raise ValueError('BinNumber not found: %i' %binNumber)
            binHeader = entry
        elif entry[0] == 'process' and count == 1:
            # get the column numbers of the right processes for the right regions
            count = None
            columnNumbers = [j for j, item in enumerate(entry) if j in columnNumbers and item == process]
            if len(columnNumbers) == 0:
                raise ValueError('process not found: %s' %process)
            columnNumber = int(columnNumbers[0])
        elif entry[0] == uncertainty:
            # get the uncertainties of the prev. extracted column numbers
            if entry[columnNumber] == '-': return 1.
            else: return float(entry[columnNumber])

    raise ValueError('uncertainty not found: %s' %uncertainty)
    

def getObservationValue( cardFile, binNumber, process ):

    with open(cardFile, 'r') as f:
        cardFileLines = f.readlines()

    count = 0
    for i, line in enumerate(cardFileLines):
        # read each line of cardfile
        entry = [ item for item in line.split(' ') if item != '' ]
        # rm lnN to get the right column number for the headers and the uncertainties
        if 'lnN' in entry: entry.remove('lnN')
        if entry[0] == 'bin' and count == 0: count += 1
            # 2nd entry of 'bin' is interesting
        elif entry[0] == 'bin' and count == 1:
            # get the column numbers of the right regions
            columnNumbers = [ j for j, bin in enumerate(entry) if bin == 'Bin' + str(binNumber) ]
            if len(columnNumbers) == 0:
                raise ValueError('BinNumber not found: %i' %binNumber)
            binHeader = entry
        elif entry[0] == 'process' and count == 1:
            # get the column numbers of the right processes for the right regions
            count = None
            columnNumbers = [j for j, item in enumerate(entry) if j in columnNumbers and item == process]
            if len(columnNumbers) == 0:
                raise ValueError('process not found: %s' %process)
            columnNumber = int(columnNumbers[0])
        elif entry[0] == 'rate':
            # get the uncertainties of the prev. extracted column numbers
            return float(entry[columnNumber])

    raise ValueError('uncertainty not found: %s' %uncertainty)
    

#print getUncertaintyValue( 'TopEFTCardFile.txt', 0, 'TTX', 'ttX')
#print getObservationValue( 'TopEFTCardFile.txt', 0, 'signal', )
