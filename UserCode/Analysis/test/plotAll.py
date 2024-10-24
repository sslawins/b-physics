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
    c.SaveAs("Plots/" + key.GetName() + ".png")