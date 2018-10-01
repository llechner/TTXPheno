def readNLLFile(fname):
        import ROOT
        f = ROOT.TFile.Open(fname)
        t = f.Get("limit")
        nll = {}
        t.GetEntry(0) # changed from 1!
        # prefit NLL
        nll["nll0"] = t.nll0
        # delta NLL to prefit (should always be negative since stuff is fitted)
        nll["nll"] = t.nll
        # absolute NLL postfit
        nll["nll_abs"] = t.nll0 + t.nll
        f.Close()
        return nll

import sys
print readNLLFile(sys.argv[1])
