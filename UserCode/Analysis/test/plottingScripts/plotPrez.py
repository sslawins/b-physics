#!/usr/bin/env python3
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

print("Hello ROOT")
fileName = "../rootOutputs/histos_PVdistance4.root"
print('Read data from:', fileName)

f = TFile(fileName)
f.ls()

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
    histo.GetXaxis().SetTitle("Distance [cm]")
    histo.GetYaxis().SetTitle("Events")
    histo.SetTitle(histo.GetName())

hPCAz_reco0SV = gROOT.FindObject("hPCAz_reco0SV")  
hPCAz_genSV = gROOT.FindObject("hPCAz_genSV")   
hPvSvDistance = gROOT.FindObject("hPvSvDistance")
'''
hBsMass = gROOT.FindObject("hBsMass")
hPhiMass = gROOT.FindObject("hPhiMass")
hBssPt = gROOT.FindObject("hBssPt")
hPhiDecayGen_Eta_Pt = gROOT.FindObject("hPhiDecayGen_Eta_Pt")
hDecayGamma_Eta_Pt = gROOT.FindObject("hDecayGamma_Eta_Pt")
hPhiDecayGen_Pt = gROOT.FindObject("hPhiDecayGen_Pt")
hPhiDecay_RecoPt_GenPt = gROOT.FindObject("hPhiDecay_RecoPt_GenPt")
hPhiDecayReco_Eta_Pt = gROOT.FindObject("hPhiDecayReco_Eta_Pt")
hDecayGammaReco_Pt = gROOT.FindObject("hDecayGammaReco_Pt")
hDecayGamma_RecoPt_GenPt = gROOT.FindObject("hDecayGamma_RecoPt_GenPt")
hCutsDecayGamma_RecoPt_GenPt = gROOT.FindObject("hCutsDecayGamma_RecoPt_GenPt")
hCascadeGen_Pt = gROOT.FindObject("hCascadeGen_Pt")
hCascadeReco_Pt = gROOT.FindObject("hCascadeReco_Pt")
hPhiPt = gROOT.FindObject("hPhiPt")


hPhiPt.GetXaxis().SetRangeUser(0, 30)
hBssPt.GetXaxis().SetRangeUser(0, 30)
hBsMass.Rebin(2)

hCascadeReco_Pt.Rebin(4)
hCascadeGen_Pt.Rebin(4)
hCutsDecayGamma_RecoPt_GenPt.Rebin2D(2, 2)

hDecayGamma_RecoPt_GenPt.GetXaxis().SetRangeUser(0, 30)
hDecayGamma_RecoPt_GenPt.GetYaxis().SetRangeUser(0, 30)
'''
hPCAz_reco0SV.GetXaxis().SetRangeUser(0, 0.1)

hPCAz_genSV.GetXaxis().SetRangeUser(0, 0.1)
hPvSvDistance.GetXaxis().SetRangeUser(0, 0.6)

histogram_dict = {
    "cPCAz_reco0SV": (hPCAz_reco0SV, "PCAz_reco0SV.png"),
    "cPCAz_genSV": (hPCAz_genSV, "PCAz_genSV.png"),
    "cPvSvDistance": (hPvSvDistance, "PvSvDistance.png")
}
'''
histogram_dict = {
    "cBsMass": (hBsMass, "Plots/Prez/BsMass.png"),
    "cPhiMass": (hPhiMass, "Plots/Prez/PhiMass.png"),
    "cBssPt": (hBssPt, "Plots/Prez/BssPt.png"),
    "cPhiDecayGen_Eta_Pt": (hPhiDecayGen_Eta_Pt, "Plots/Prez/PhiDecayGen_Eta_Pt.png"),
    "cDecayGamma_Eta_Pt": (hDecayGamma_Eta_Pt, "Plots/Prez/DecayGamma_Eta_Pt.png"),
    "cPhiDecayGen_Pt": (hPhiDecayGen_Pt, "Plots/Prez/PhiDecayGen_Pt.png"),
    "cPhiDecay_RecoPt_GenPt": (hPhiDecay_RecoPt_GenPt, "Plots/Prez/PhiDecay_RecoPt_GenPt.png"),
    "cPhiDecayReco_Eta_Pt": (hPhiDecayReco_Eta_Pt, "Plots/Prez/PhiDecayReco_Eta_Pt.png"),
    "cDecayGammaReco_Pt": (hDecayGammaReco_Pt, "Plots/Prez/DecayGammaReco_Pt.png"),
    "cDecayGamma_RecoPt_GenPt": (hDecayGamma_RecoPt_GenPt, "Plots/Prez/DecayGamma_RecoPt_GenPt.png"),
    "cCutsDecayGamma_RecoPt_GenPt": (hCutsDecayGamma_RecoPt_GenPt, "Plots/Prez/CutsDecayGamma_RecoPt_GenPt.png"),
    "cCascadeGen_Pt": (hCascadeGen_Pt, "Plots/Prez/CascadeGen_Pt.png"),
    "cCascadeReco_Pt": (hCascadeReco_Pt, "Plots/Prez/CascadeReco_Pt.png"),
    "cPhiPt": (hPhiPt, "Plots/Prez/PhiPt.png")
}
'''
for canvas_name, (histogram, output_file) in histogram_dict.items():
    if histogram:
        canvas = TCanvas(canvas_name, canvas_name, 800, 600)

        histogram.SetTitle("")

        canvas.SetLeftMargin(0.14)
        DrawHisto(histogram, 19, 1, 2)
        canvas.cd()
 
        histogram.Draw('hist')
        if histogram.GetDimension() == 2:
            histogram.Draw("COLZ")
            colorbar = histogram.GetListOfFunctions().FindObject("palette")
        cms_label = TLatex()
        cms_label.SetTextFont(42)
        cms_label.SetTextSize(0.04)
        cms_label.SetTextAlign(10)  # Left-align
        cms_label.DrawLatexNDC(0.17, 0.92, " #it{Simulation}")
        canvas.Draw()
        canvas.Print(output_file)
    else:
        print(f"Warning: Histogram for {canvas_name} not found.")

input('press enter to exit')