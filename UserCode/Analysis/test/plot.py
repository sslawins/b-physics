import ROOT
import sys

filename = sys.argv[1]
file = ROOT.TFile(filename)

for key in file.GetListOfKeys():
    hist = file.Get(key.GetName())

    c = ROOT.TCanvas()
    hist.Draw()
    c.SaveAs('Plots/' + key.GetName() + '.png')

