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
fileName = "histos21.root"
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
    

hGammaWithMuonConditionPt = gROOT.FindObject("hGammaWithMuonConditionPt")
hMatchedGammaWithMuonConditionPt = gROOT.FindObject("hMatchedGammaWithMuonConditionPt")
hRecoGammaWithMuonConditionPt =gROOT.FindObject("hRecoGammaWithMuonConditionPt")
hRecoVsGenGammaWithMuonConditionPt = gROOT.FindObject("hRecoVsGenGammaWithMuonConditionPt")
hGenGammaPtVsEta = gROOT.FindObject("hGenGammaPtVsEta")
hGammaDeltaR = gROOT.FindObject("hGammaDeltaR")

hGenMuPt = gROOT.FindObject("hGenMuPt")
hRecoMuWithMuonConditionPt = gROOT.FindObject("hRecoMuWithMuonConditionPt")
hRecoVsGenMuPtWithMuonCondition = gROOT.FindObject("hRecoVsGenMuPtWithMuonCondition")
hMuDeltaR = gROOT.FindObject("hMuDeltaR")
hMu1VsMu2DeltaR = gROOT.FindObject("hMu1VsMu2DeltaR")

hBs0Pt = gROOT.FindObject("Bs0Pt")
BsInvMass = gROOT.FindObject("BsInvMass")
hGenMuPtVsEta = gROOT.FindObject("hGenMuPtVsEta")
hRecoMuPt = gROOT.FindObject("hRecoMuPt")
hRecoVsGenMuPt = gROOT.FindObject("hRecoVsGenMuPt")
hRecoMuPtVsEta = gROOT.FindObject("hRecoMuPtVsEta")
hRecoGammaPt = gROOT.FindObject("hRecoGammaPt")
hRecoVsGenGammaPt = gROOT.FindObject("hRecoVsGenGammaPt")
hRecoGammaPtVsEta = gROOT.FindObject("hRecoGammaPtVsEta")
hRecoVsGenGammaWithMuonConditionPt = gROOT.FindObject("hRecoVsGenGammaWithMuonConditionPt")

BsInvMass.SetTitle("")
BsInvMass.Rebin(2)
hBs0Pt.GetXaxis().SetRangeUser(0, 50)
hRecoGammaPtVsEta.GetYaxis().SetRangeUser(-3, 3)
hRecoMuPtVsEta.GetYaxis().SetRangeUser(-3, 3)
hGammaWithMuonConditionPt.Rebin(2)

hMatchedGammaWithMuonConditionPt.Rebin(5)

hRecoGammaWithMuonConditionPt.Rebin(4)

hGammaDeltaR.Rebin(2)


hMuDeltaR.GetXaxis().SetRangeUser(0, 0.05)

hMu1VsMu2DeltaR.Rebin(2)
#histo.Rebin2D(rebin_x, rebin_y)



histogram_dict = {
    "cBsInvMass": (BsInvMass, "Plots/Szymon/BsInvMass.png"),
    "cGammaWithMuonConditionPt": (hGammaWithMuonConditionPt, "Plots/Szymon/GammaWithMuonConditionPt.png"),
    "cMatchedGammaWithMuonConditionPt": (hMatchedGammaWithMuonConditionPt, "Plots/Szymon/MatchedGammaWithMuonConditionPt.png"),
    "cRecoGammaWithMuonConditionPt": (hRecoGammaWithMuonConditionPt, "Plots/Szymon/RecoGammaWithMuonConditionPt.png"),
    "cRecoVsGenGammaWithMuonConditionPt": (hRecoVsGenGammaWithMuonConditionPt, "Plots/Szymon/RecoVsGenGammaWithMuonConditionPt.png"),
    "cGenGammaPtVsEta": (hGenGammaPtVsEta, "Plots/Szymon/GenGammaPtVsEta.png"),
    "cGammaDeltaR": (hGammaDeltaR, "Plots/Szymon/GammaDeltaR.png"),
    "cGenMuPt": (hGenMuPt, "Plots/Szymon/GenMuPt.png"),
    "cRecoMuWithMuonConditionPt": (hRecoMuWithMuonConditionPt, "Plots/Szymon/RecoMuWithMuonConditionPt.png"),
    "cRecoVsGenMuPtWithMuonCondition": (hRecoVsGenMuPtWithMuonCondition, "Plots/Szymon/RecoVsGenMuPtWithMuonCondition.png"),
    "cMuDeltaR": (hMuDeltaR, "Plots/Szymon/MuDeltaR.png"),
    "cMu1VsMu2DeltaR": (hMu1VsMu2DeltaR, "Plots/Szymon/Mu1VsMu2DeltaR.png"),
    "cBs0Pt": (hBs0Pt, "Plots/Szymon/Bs0Pt.png"),
    "cGenMuPtVsEta": (hGenMuPtVsEta, "Plots/Szymon/GenMuPtVsEta.png"),
    "cRecoMuPt": (hRecoMuPt, "Plots/Szymon/RecoMuPt.png"),
    "cRecoVsGenMuPt": (hRecoVsGenMuPt, "Plots/Szymon/RecoVsGenMuPt.png"),
    "cRecoMuPtVsEta": (hRecoMuPtVsEta, "Plots/Szymon/RecoMuPtVsEta.png"),
    "cRecoGammaPt": (hRecoGammaPt, "Plots/Szymon/RecoGammaPt.png"),
    "cRecoVsGenGammaPt": (hRecoVsGenGammaPt, "Plots/Szymon/RecoVsGenGammaPt.png"),
    "cRecoGammaPtVsEta": (hRecoGammaPtVsEta, "Plots/Szymon/RecoGammaPtVsEta.png"),
    "cRecoVsGenGammaWithMuonConditionPt": (hRecoVsGenGammaWithMuonConditionPt, "Plots/Szymon/RecoVsGenGammaWithMuonConditionPt.png")
}

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