#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import math
import numpy as np

from ROOT import gStyle, TFile, TCanvas, gROOT, TLatex

gStyle.SetOptStat(1110)
gStyle.SetOptTitle(1)
gStyle.SetStatX(0.85)  # Set X position (right)
gStyle.SetStatY(0.88)  # Set Y position (top)
gStyle.SetStatW(0.2)   # Width of the stats box
gStyle.SetStatH(0.15)  # Height of the stats box


def NormalizeHistogram(histo):
    entries = histo.GetEntries()  # Get the total number of entries in the histogram
    if entries > 0:  # Avoid division by zero
        histo.Scale(1.0 / entries)
    else:
        print(f"Warning: Histogram {histo.GetName()} has zero entries and cannot be normalized.")


def DrawHisto(histo, color, line=1, width=2):
    histo.SetLineColor(line)
    histo.SetLineWidth(width)
    histo.SetFillColor(color)


def SaveCanvas(canvas, histo, name, output_filename="Plots/PhiGamma/{}.png", rebin=1, norm=False, log=False, xmin=None, xmax=None):
    if norm:
        NormalizeHistogram(histo)

    if log:
        canvas.SetLogy(1)

    canvas.SetLeftMargin(0.11)
    canvas.SetRightMargin(0.05)

    canvas.cd()

    # Rebinnowanie histogramu przed rysowaniem
    histo.Rebin(rebin)


    # Rysowanie histogramu
    DrawHisto(histo, 19)
    if xmin is not None and xmax is not None:
        print(f"Setting X-axis range: ({xmin}, {xmax})")  # Debug print
        histo.GetXaxis().SetRangeUser(xmin, xmax)
        
    histo.Draw('hist')

    histo.SetTitleSize(0.02, "t")    
    histo.SetTitleOffset(0.1, "t")  # Move the title up

    cms_label = TLatex()
    cms_label.SetTextFont(42)
    cms_label.SetTextSize(0.04)
    cms_label.SetTextAlign(10)  # Left-align

    canvas.Draw()
    canvas.Print(output_filename.format(name))


print("Hello ROOT")
fileName = "histos_PhiGamma.root"
print('Read data from:', fileName)

f = TFile(fileName)
f.ls()

# Dodano zakresy dla osi X
histograms = {
    'cBPt': ('hBPt', 'BPt', 1, True, False, 0, 50),
    'cBssPt': ('hBssPt', 'BssPt', 1, True, False, 0, 50),
    'cBsToPhiG': ('hBsToPhiG', 'BsToPhiG', 1, False, True, None, None),
    'cBsToB': ('hBsToB', 'BsToB', 1, False, True, None, None),
    'cNinit': ('hNinit', 'Ninit', 1, True, True, None, None),
    'cInitial': ('hInitial', 'Initial', 1, True, True, None, None),
    'cHad': ('hHad', 'Had', 1, True, True, None, None),
    'cBsParents': ('hBsParents', 'BsParents', 1, True, True, None, None),
    'cBsAncestor': ('hBsAncestor', 'BsAncestor', 1, True, True, None, None),
    'cCascade': ('hCascade', 'Cascade', 1, True, True, None, None),
    'cMuMuDeltaR': ('hMuMuDeltaR', 'MuMuDeltaR', 2, True, False, 0, 7),
    'cMuEta_MuEta': ('hMuEta_MuEta', 'MuEta_MuEta', 8, False, False, -10, 10),
    'cMuPt': ('hMuPt', 'MuPt', 1, True, False, 0, 10),
    'cMuPtAll': ('hMuPtAll', 'MuPtAll', 1, False, False, 0, 15),
    'cPhiPt': ('hPhiPt', 'PhiPt', 1, True, False, 0, 30),
    'cPhiEta': ('hPhiEta', 'PhiEta', 2, False, False, -10, 10),
    'cGEta_GPt': ('hGEta_GPt', 'GEta_GPt', 1, False, False, -10, 10),
    'cPhiEta_PhiPt': ('hPhiEta_PhiPt', 'PhiEta_PhiPt', 1, False, False, -10, 10),
    'cGammaPt': ('hGammaPt', 'GammaPt', 1, True, False, 0, 30.),
    'cGammaEta': ('hGammaEta', 'GammaEta', 2, False, False, -10, 10)
}

canvases = {name: TCanvas(name, name, 800, 600) for name in histograms.keys()}
histos = {hname: gROOT.FindObject(hname) for hname, _, _, _, _, _, _ in histograms.values()}

for cname, (hname, title, rebin, norm, log, xmin, xmax) in histograms.items():
    histo = histos.get(hname)
    
    if histo is None:
        print(f"Warning: Histogram '{hname}' not found.")
        continue

    canvas = canvases[cname]
    SaveCanvas(canvas, histo, title, rebin=rebin, norm=norm, log=log, xmin=xmin, xmax=xmax)

input('press enter to exit')
