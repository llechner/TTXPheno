from TTXPheno.Analysis.Region import Region

def getRegionsFromThresholds(var, vals, gtLastThreshold = True):
    return [Region(var, (vals[i], vals[i+1])) for i in range(len(vals)-1)]

def getRegions2D(varOne, varOneThresholds, varTwo, varTwoThresholds):
    regions_varOne  = getRegionsFromThresholds(varOne,  varOneThresholds)
    regions_varTwo  = getRegionsFromThresholds(varTwo, varTwoThresholds)

    regions2D = []
    for r1 in regions_varOne:
        for r2 in regions_varTwo:
            regions2D.append(r1+r2)

    return regions2D

#Put all sets of regions that are used in the analysis, closure, tables, etc.

## 3l signal regions
regions = getRegions2D("genZ_pt", [0,100,200,400], "genZ_cosThetaStar", [-1,-0.6, 0.6, 1]) + [Region("genZ_pt", (400, -1))]
genttZRegions = getRegions2D("genZ_pt", [0,100,200,400], "genZ_cosThetaStar", [-1,-0.6, 0.6, 1]) + [Region("genZ_pt", (400, -1))]
genttgammaRegions = getRegionsFromThresholds("genPhoton_pt[0]", [0,100,200,300]) + [Region("genPhoton_pt[0]", (300, -1))]

recottZRegions = getRegions2D("recoZ_pt", [0,100,200,400], "recoZ_cosThetaStar", [-1,-0.6, 0.6, 1]) + [Region("recoZ_pt", (400, -1))]
recottgammaRegions = getRegionsFromThresholds("recoPhoton_pt[0]", [0,100,200,300]) + [Region("recoPhoton_pt[0]", (300, -1))]

genttZRegionsSmall = getRegions2D("genZ_pt", [200,400], "genZ_cosThetaStar", [-1,-0.6])# + [Region("genZ_pt", (400, -1))]
genttgammaRegionsSmall = getRegionsFromThresholds("genPhoton_pt[0]", [200,300])# + [Region("genPhoton_pt[0]", (300, -1))]

recottZRegionsSmall = getRegions2D("recoZ_pt", [200,400], "recoZ_cosThetaStar", [-1,-0.6])# + [Region("recoZ_pt", (400, -1))]
recottgammaRegionsSmall = getRegionsFromThresholds("recoPhoton_pt[0]", [200,300])# + [Region("recoPhoton_pt[0]", (300, -1))]

