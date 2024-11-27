import ROOT
import sys
from ROOT import gROOT

gROOT.SetBatch(True)

f = ROOT.TFile(sys.argv[1])

c = ROOT.TCanvas()

for key in f.GetListOfKeys():
    c.Clear()
    hist = f.Get(key.GetName())
    hist.Draw()
    if key.GetName() == "hPdgidPos" or key.GetName() == "hPdgidNeg":
        c.SetLogy()
    if len(sys.argv) > 2:
        c.SaveAs(f"{sys.argv[2]}/" + key.GetName() + ".png")
    else:
        c.SaveAs("Plots/" + key.GetName() + ".png")