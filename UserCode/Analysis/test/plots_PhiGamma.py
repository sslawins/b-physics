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


def SaveCanvas(canvas, histo, name, output_filename= "Plots/PhiGamma/{}.png", rebin_x=1, rebin_y=1, norm=False, log=False, xmin=None, xmax=None, ymin=None, ymax=None):
    
    if norm:
        NormalizeHistogram(histo)

    if log:
        canvas.SetLogy(1)

    canvas.SetLeftMargin(0.11)
    canvas.SetRightMargin(0.15)

    canvas.cd()

    if rebin_x > 1 or (histo.GetDimension() == 2 and rebin_y > 1):
        if histo.GetDimension() == 2:
            histo.Rebin2D(rebin_x, rebin_y)
        else:
            histo.Rebin(rebin_x)

    DrawHisto(histo, 19)
    if xmin is not None and xmax is not None:
        histo.GetXaxis().SetRangeUser(xmin, xmax)
    if ymin is not None and ymax is not None and histo.GetDimension() == 2:
        histo.GetYaxis().SetRangeUser(ymin, ymax)

        
    if histo.GetDimension() == 2:
        histo.Draw("COLZ")
        colorbar = histo.GetListOfFunctions().FindObject("palette")
        if colorbar:
            colorbar.SetLabelFont(42)
            colorbar.SetLabelSize(0.04)
    else:
        histo.Draw("HIST")

    histo.SetTitle( " ")    
    histo.GetXaxis().SetTitleOffset(1.1)
    histo.GetYaxis().SetTitleOffset(1.4)

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

#canvas, histo, name, output_filename= "Plots/PhiGamma/{}.png", rebin_x=1, rebin_y=1, norm=False, log=False, xmin=None, xmax=None, ymin=None, ymax=None
histograms = {
    #'cBPt': ('hBPt', 'BPt', 1, 1, True, False, 0, 50, None, None),
    #'cBssPt': ('hBssPt', 'BssPt', 1, 1, True, False, 0, 50, None, None),
    #'cBsToPhiG': ('hBsToPhiG', 'BsToPhiG', 1, 1, False, True, None, None, None, None),
    #'cPhiToMuMu': ('hPhiToMuMu', 'PhiToMuMu', 1, 1, False, True, None, None, None, None),
    #'cPhiDecayProducts': ('hPhiDecayProducts', 'PhiDecayProducts', 1, 1, False, True, None, None, None, None),
    'cBsMass': ('hBsMass', 'BsMass', 1, 1, False, False, None, None, None, None),
    'cPhiMass': ('hPhiMass', 'PhiMass', 1, 1, False, False, None, None, None, None),
    #'cBsToB': ('hBsToB', 'BsToB', 1, 1, False, True, None, None, None, None),
    #'cNinit': ('hNinit', 'Ninit', 1, 1, True, True, None, None, None, None),
    #'cInitial': ('hInitial', 'Initial', 1, 1, True, True, None, None, None, None),
    #'cHad': ('hHad', 'Had', 1, 1, True, True, None, None, None, None),
    #'cBsParents': ('hBsParents', 'BsParents', 1, 1, True, True, None, None, None, None),
    #'cBsAncestor': ('hBsAncestor', 'BsAncestor', 1, 1, True, True, None, None, None, None),
    '''
    'cCascade': ('hCascade', 'Cascade', 1, 1, True, True, None, None, None, None),
    'cCascadeGen_Eta_Pt': ('hCascadeGen_Eta_Pt', 'CascadeGen_Eta_Pt', 1, 1, True, True, None, None, None, None),
    'cCascadeGen_Eta_Eta': ('hCascadeGen_Eta_Eta', 'CascadeGen_Eta_Eta', 1, 1, True, True, None, None, None, None),
    'cCascadeGen_PtAll': ('hCascadeGen_PtAll', 'CascadeGen_PtAll', 1, 1, True, True, None, None, None, None),
    'cCascadeGen_Pt': ('hCascadeGen_Pt', 'CascadeGen_Pt', 1, 1, True, True, None, None, None, None),
    'cCascadeMuMuDeltaR': ('hCascadeMuMuDeltaR', 'CascadeMuMuDeltaR', 1, 1, True, True, None, None, None, None),
    'cCascadeMuMuDeltaEta': ('hCascadeMuMuDeltaEta', 'CascadeMuMuDeltaEta', 1, 1, True, True, None, None, None, None),
    'cCascadeMuMuDeltaPhi': ('hCascadeMuMuDeltaPhi', 'CascadeMuMuDeltaPhi', 1, 1, True, True, None, None, None, None),
    'cCascadeReco_Pt': ('hCascadeReco_Pt', 'CascadeReco_Pt', 1, 1, True, True, None, None, None, None),
    'cCascadeReco_Eta_Pt': ('hCascadeReco_Eta_Pt', 'CascadeReco_Eta_Pt', 1, 1, True, True, None, None, None, None),
    'cCascade_RecoGenDeltaR': ('hCascade_RecoGenDeltaR', 'Cascade_RecoGenDeltaR', 1, 1, True, True, None, None, None, None),
    'cCascade_RecoPt_GenPt': ('hCascade_RecoPt_GenPt', 'Cascade_RecoPt_GenPt', 1, 1, True, True, None, None, None, None),
    '''
    'cPhiPt': ('hPhiPt', 'PhiPt', 1, 1, True, False, None, None, None, None),
    'cPhiEta': ('hPhiEta', 'PhiEta', 1, 1, True, False, None, None, None, None),
    'cPhi_Eta_Pt': ('hPhi_Eta_Pt', 'Phi_Eta_Pt', 1, 1, True, False, None, None, None, None),
    'cDecayGamma_Eta_Pt': ('hDecayGamma_Eta_Pt', 'DecayGamma_Eta_Pt', 1, 1, True, False, None, None, None, None),
    'cDecayGamma_Pt': ('hDecayGamma_Pt', 'DecayGamma_Pt', 1, 1, True, False, None, None, None, None),
    'cDecayGamma_Eta': ('hDecayGamma_Eta', 'DecayGamma_Eta', 1, 1, True, False, None, None, None, None),
    'cDecayGammaReco_Pt': ('hDecayGammaReco_Pt', 'DecayGammaReco_Pt', 1, 1, True, False, None, None, None, None),
    'cDecayGammaReco_Eta_Pt': ('hDecayGammaReco_Eta_Pt', 'DecayGammaReco_Eta_Pt', 1, 1, True, False, None, None, None, None),
    'cDecayGamma_RecoGenDeltaR': ('hDecayGamma_RecoGenDeltaR', 'DecayGamma_RecoGenDeltaR', 1, 1, True, False, None, None, None, None),
    'cDecayGamma_RecoPt_GenPt': ('hDecayGamma_RecoPt_GenPt', 'DecayGamma_RecoPt_GenPt', 1, 1, False, False, None, None, None, None),
    'cCutsDecayGammaReco_Pt': ('hCutsDecayGammaReco_Pt', 'CutsDecayGammaReco_Pt', 4, 1, False, False, None, None, None, None),
    'cCutsDecayGammaReco_Eta_Pt': ('hCutsDecayGammaReco_Eta_Pt', 'CutsDecayGammaReco_Eta_Pt', 10, 10, False, False, None, None, None, None),
    'cCutsDecayGamma_RecoGenDeltaR': ('hCutsDecayGamma_RecoGenDeltaR', 'CutsDecayGamma_RecoGenDeltaR', 1, 1, False, False, 0, 0.1, None, None),
    'cCutsDecayGamma_RecoPt_GenPt': ('hCutsDecayGamma_RecoPt_GenPt', 'CutsDecayGamma_RecoPt_GenPt', 5, 5, False, False, None, None, None, None),
    'cCutsDecayGammaGen_Pt': ('hCutsDecayGammaGen_Pt', 'CutsDecayGammaGen_Pt', 4, 1, False, False, 0, 5, None, None),
    'cCutsDecayGammaGen_Eta_Pt': ('hCutsDecayGammaGen_Eta_Pt', 'CutsDecayGammaGen_Eta_Pt', 10, 10, False, False, None, None, None, None),
    'cCutsDecayGammaGenMatched_Pt': ('hCutsDecayGammaGenMatched_Pt', 'CutsDecayGammaGenMatched_Pt', 4, 1, False, False, 0, 5, None, None),
    'cCutsDecayGammaGenMatched_Eta_Pt': ('hCutsDecayGammaGenMatched_Eta_Pt', 'CutsDecayGammaGenMatched_Eta_Pt', 10, 10, True, False, None, None, None, None),
    'cPhiDecayGen_Pt': ('hPhiDecayGen_Pt', 'PhiDecayGen_Pt', 1, 1, True, False, None, None, None, None),
    'cPhiDecayGen_Eta_Pt': ('hPhiDecayGen_Eta_Pt', 'PhiDecayGen_Eta_Pt', 1, 1, False, False, None, None, None, None),
    'cPhiDecayReco_Pt': ('hPhiDecayReco_Pt', 'PhiDecayReco_Pt', 1, 1, False, False, None, None, None, None),
    'cPhiDecayReco_Eta_Pt': ('hPhiDecayReco_Eta_Pt', 'PhiDecayReco_Eta_Pt', 1, 1, True, False, None, None, None, None),
    'cPhiDecay_RecoGenDeltaR': ('hPhiDecay_RecoGenDeltaR', 'PhiDecay_RecoGenDeltaR', 1, 1, False, False, 0., 0.05, None, None),
    'cPhiDecay_RecoPt_GenPt': ('hPhiDecay_RecoPt_GenPt', 'PhiDecay_RecoPt_GenPt', 1, 1, False, False, None, None, None, None),
    'cCutsPhiDecayReco_Pt': ('hCutsPhiDecayReco_Pt', 'CutsPhiDecayReco_Pt', 2, 1, False, False, None, None, None, None),
    'cCutsPhiDecayReco_Eta_Pt': ('hCutsPhiDecayReco_Eta_Pt', 'CutsPhiDecayReco_Eta_Pt', 5, 5, False, False, None, None, None, None),
    'cCutsPhiDecay_RecoGenDeltaR': ('hCutsPhiDecay_RecoGenDeltaR', 'CutsPhiDecay_RecoGenDeltaR', 1, 1, False, False, 0., 0.01, None, None),
    'cCutsPhiDecay_RecoPt_GenPt': ('hCutsPhiDecay_RecoPt_GenPt', 'CutsPhiDecay_RecoPt_GenPt', 1, 1, False, False, None, None, None, None)
}

CMS.SetExtraText("")
CMS.SetCmsText("Simulation")
CMS.SetCmsTextFont(52)
CMS.SetCmsTextSize(0.75*0.76)
CMS.SetLumi("")

canvases = {name: TCanvas(name, name, 600, 600) for name in histograms.keys()}
#canvases = {name: CMS.cmsCanvas('', 0, 1, 0, 1, '', '', square=CMS.kSquare, extraSpace=0.01, iPos=0) for name in histograms.keys()}
histos = {hname: gROOT.FindObject(hname) for hname, _, _, _, _, _, _, _, _, _ in histograms.values()}

#canvas, histo, name, output_filename= "Plots/PhiGamma/{}.png", rebin_x=1, rebin_y=1, norm=False, log=False, xmin=None, xmax=None, ymin=None, ymax=None
for cname, (hname, title, rebin_x, rebin_y, norm, log, xmin, xmax, ymin, ymax) in histograms.items():
    histo = histos.get(hname)
    
    if histo is None:
        print(f"Warning: Histogram '{hname}' not found.")
        continue

    canvas = canvases[cname]
    #
    SaveCanvas(canvas, histo, title, rebin_x=rebin_x, rebin_y = rebin_y, norm=norm, log=log, xmin=xmin, xmax=xmax, ymin =ymin, ymax=ymax)

input('press enter to exit')
