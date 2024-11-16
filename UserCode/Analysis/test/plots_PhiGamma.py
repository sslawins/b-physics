#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import math
import numpy as np
import cmsstyle as CMS

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

    histo.Rebin(rebin)

    DrawHisto(histo, 19)
    if xmin is not None and xmax is not None:
        print(f"Setting X-axis range: ({xmin}, {xmax})")  # Debug print
        histo.GetXaxis().SetRangeUser(xmin, xmax)
        
    histo.Draw('hist')

    histo.SetTitle( " ")    
    histo.SetTitleOffset(0.1, "t")  # Move the title up


    cms_label = TLatex()
    cms_label.SetTextFont(42)
    cms_label.SetTextSize(0.04)
    cms_label.SetTextAlign(12)  # Left-align
    cms_label.DrawLatexNDC(0.17, 0.92, " #it{Simulation}")

    canvas.Draw()
    canvas.Print(output_filename.format(name))


print("Hello ROOT")
fileName = "histos_PhiGamma.root"
print('Read data from:', fileName)

f = TFile(fileName)
f.ls()

histograms= {
    'cBPt': ('hBPt', 'BPt', 1, True, False, 0, 50),
    'cBssPt': ('hBssPt', 'BssPt', 1, True, False, 0, 50),
    'cBsToPhiG': ('hBsToPhiG', 'BsToPhiG', 1, False, True, None, None),
    'cPhiToMuMu': ('hPhiToMuMu', 'PhiToMuMu', 1, False, True, None, None),
    'cBsToB': ('hBsToB', 'BsToB', 1, False, True, None, None),
    'cNinit': ('hNinit', 'Ninit', 1, True, True, None, None),
    'cInitial': ('hInitial', 'Initial', 1, True, True, None, None),
    'cHad': ('hHad', 'Had', 1, True, True, None, None),
    'cBsParents': ('hBsParents', 'BsParents', 1, True, True, None, None),
    'cBsAncestor': ('hBsAncestor', 'BsAncestor', 1, True, True, None, None),
    'cCascade': ('hCascade', 'Cascade', 1, True, True, None, None),
    'cCascadeGen_Eta_Pt': ('hCascadeGen_Eta_Pt', 'CascadeGen_Eta_Pt', 1, False, False, None, None),
    'cCascadeGen_Eta_Eta': ('hCascadeGen_Eta_Eta', 'CascadeGen_Eta_Eta', 1, False, False, None, None),
    'cCascadeGen_PtAll': ('hCascadeGen_PtAll', 'CascadeGen_PtAll', 1, False, False, 0, 15),
    'cCascadeGen_Pt': ('hCascadeGen_Pt', 'CascadeGen_Pt', 1, True, False, 0, 10),
    'cCascadeMuMuDeltaR': ('hCascadeMuMuDeltaR', 'CascadeMuMuDeltaR', 1, True, False, 0, 7),
    'cCascadeMuMuDeltaEta': ('hCascadeMuMuDeltaEta', 'CascadeMuMuDeltaEta', 2, False, False, None, None),
    'cCascadeMuMuDeltaPhi': ('hCascadeMuMuDeltaPhi', 'CascadeMuMuDeltaPhi', 2, False, False, None, None),
    'cCascadeReco_Pt': ('hCascadeReco_Pt', 'CascadeReco_Pt', 1, True, False, 0, 10),
    'cCascadeReco_Eta_Pt': ('hCascadeReco_Eta_Pt', 'CascadeReco_Eta_Pt', 1, False, False, None, None),
    'cCascade_RecoGenDeltaR': ('hCascade_RecoGenDeltaR', 'Cascade_RecoGenDeltaR', 1, False, False, None, None),
    'cPhiPt': ('hPhiPt', 'PhiPt', 1, True, False, 0, 30),
    'cPhiEta': ('hPhiEta', 'PhiEta', 2, False, False, -10, 10),
    'cPhi_Eta_Pt': ('hPhi_Eta_Pt', 'PhiEta_PhiPt', 1, False, False, -10, 10),
    'cDecayGamma_Eta_Pt': ('hDecayGamma_Eta_Pt', 'DecayGamma_Eta_Pt', 1, False, False, -10, 10),
    'cDecayGamma_Pt': ('hDecayGamma_Pt', 'DecayGamma_Pt', 1, True, False, 0, 30),
    'cDecayGamma_Eta': ('hDecayGamma_Eta', 'DecayGamma_Eta', 2, False, False, -10, 10),
    'cDecayGammaReco_Pt': ('hDecayGammaReco_Pt', 'DecayGammaReco_Pt', 1, True, False, 0, 30),
    'cDecayGammaReco_Eta_Pt': ('hDecayGammaReco_Eta_Pt', 'DecayGammaReco_Eta_Pt', 1, False, False, -10, 10),
    'cDecayGamma_RecoGenDeltaR': ('hDecayGamma_RecoGenDeltaR', 'DecayGamma_RecoGenDeltaR', 1, False, False, None, None),
    'cPhiDecayGen_Pt': ('hPhiDecayGen_Pt', 'PhiDecayGen_Pt', 1, True, False, 0, 10),
    'cPhiDecayGen_Eta_Pt': ('hPhiDecayGen_Eta_Pt', 'PhiDecayGen_Eta_Pt', 1, False, False, -10, 10),
    'cPhiDecayReco_Pt': ('hPhiDecayReco_Pt', 'PhiDecayReco_Pt', 1, True, False, 0, 10),
    'cPhiDecayReco_Eta_Pt': ('hPhiDecayReco_Eta_Pt', 'PhiDecayReco_Eta_Pt', 1, False, False, -10, 10),
    'cPhiDecay_RecoGenDeltaR': ('hPhiDecay_RecoGenDeltaR', 'PhiDecay_RecoGenDeltaR', 1, False, False, None, None)
}

canvases = {name: TCanvas(name, name, 600, 600) for name in histograms.keys()}
histos = {hname: gROOT.FindObject(hname) for hname, _, _, _, _, _, _ in histograms.values()}

for cname, (hname, title, rebin, norm, log, xmin, xmax) in histograms.items():
    histo = histos.get(hname)
    
    if histo is None:
        print(f"Warning: Histogram '{hname}' not found.")
        continue

    canvas = canvases[cname]
    SaveCanvas(canvas, histo, title, rebin=rebin, norm=norm, log=log, xmin=xmin, xmax=xmax)

input('press enter to exit')
