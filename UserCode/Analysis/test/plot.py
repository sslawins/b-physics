import ROOT
import sys
from ROOT import gROOT

gROOT.SetBatch(True)

filename = sys.argv[1]
file = ROOT.TFile(filename)

c = ROOT.TCanvas()

for key in file.GetListOfKeys():
    hist = file.Get(key.GetName())
    c.Clear()
    hist.Draw()
    c.SaveAs('Plots/' + key.GetName() + '.png')

