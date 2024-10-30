#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import math
import numpy as np

from ROOT import gStyle, TPaveText, TFile, TCanvas, gROOT, TLatex, gStyle

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


def DrawHisto ( histo, color, line=1, width=2, stats = 0):
    #histo.SetStats(1)
    
    histo.SetLineColor( line  )
    histo.SetLineWidth( width  )
    histo.SetFillColor( color )

def SaveCanvas(canvas, histo, name, output_filename="Plots/{}.png", rebin = 1, norm = False, log = False):

    if norm :  
        NormalizeHistogram(histo)

    if log :  
        canvas.SetLogy(1)

    canvas.SetLeftMargin(0.11)  
    canvas.SetRightMargin(0.05)  
    #canvas.SetTopMargin(0.08)    
    #canvas.SetBottomMargin(0.9)

    canvas.cd()
    DrawHisto(histo, 19) 
    histo.SetTitleSize(0.02, "t")

    histo.SetTitleOffset(0.1 , "t")  # Move the title up
    histo.Draw('hist')
    histo = histo.Rebin(rebin)

    cms_label = TLatex()
    cms_label.SetTextFont(42)
    cms_label.SetTextSize(0.04)
    cms_label.SetTextAlign(10)  # Left-align
    #cms_label.DrawLatexNDC(0.1, 0.91, "#bf{CMS} #it{Simulation}")
    
    '''
    # Draw the luminosity label
    lumi_label = TLatex()
    lumi_label.SetTextFont(42)
    lumi_label.SetTextSize(0.04)
    lumi_label.SetTextAlign(31)  # Right-align
    lumi_text = "2024 (13.6 TeV)"  # Adjust this text as necessary
    lumi_label.DrawLatexNDC(0.94444, 0.92, lumi_text)
    '''
    canvas.Draw()
    canvas.Print(output_filename.format(name))


def Background( a , b , fun , color ):
    
    fun.SetParameters( a , b)
    fun.SetLineColor( color )  
    fun.SetLineStyle(9)  
    fun.SetLineWidth(2)  


print ("Hello ROOT")
fileName = "histos.root"
print ('Read data from: ', fileName)

f = TFile(fileName);
f.ls();

cBPt = TCanvas( 'cBPt' , 'cBPt' , 800 , 600 )
cBsPt = TCanvas( 'cBsPt' , 'cBsPt' , 800 , 650 )
cBPtFrag = TCanvas( 'cBPtFrag' , 'cBPtFrag' , 800 , 600 )
cBssPtFrag= TCanvas( 'cBssPtFrag' , 'cBssPtFrag' , 800 , 600 )
cBsPtFrag = TCanvas( 'cBsPtFrag' , 'cBsPtFrag' , 800 , 600 )
cAntiBsPtFrag = TCanvas( 'cAntiBsPtFrag' , 'cAntiBsPtFrag' , 800 , 600 )

cBsToMuMuG = TCanvas( 'cBsToMuMuG' , 'cBsToMuMuG' , 800 , 600 )
cBsToBd = TCanvas( 'cBsToBd' , 'cBsToBd' , 800 , 600 )
cBsToB = TCanvas( 'cBsToB' , 'cBsToB' , 800 , 600 )

cHad = TCanvas("cHad", "cHad", 800, 600)
cBsParents = TCanvas("cBsParents", "cBsParents", 800, 600)
cBsAncestor = TCanvas("cBsAncestor", "cBsAncestor", 800, 600)
cBsProduct = TCanvas("cBsProduct", "cBsProduct", 800, 600)

cGBsGDeltaR = TCanvas("cGBsGDeltaR", "cGBsGDeltaR", 800, 600)
cGBsEta_GEta = TCanvas("cGBsEta_GEta", "cGBsEta_GEta", 800, 600)
cGammaBsPt = TCanvas("cGammaBsPt", "cGammaBsPt", 800, 600)

cMuPt = TCanvas("cMuPt", "cMuPt", 800, 600)
cMuEta = TCanvas("cMuEta", "cMuEta", 800, 600)
cMu1Eta_Mu2Eta = TCanvas("cMu1Eta_Mu2Eta", "cMu1Eta_Mu2Eta", 800, 600)
cGammaPt = TCanvas("cGammaPt", "cGammaPt", 800, 600)
cGammaEta = TCanvas("cGammaEta", "cGammaEta", 800, 600)



hBPt = gROOT.FindObject('hBPt')
hBsPt = gROOT.FindObject('hBsPt')
hBPtFrag = gROOT.FindObject('hBPtFrag')
hBssPtFrag = gROOT.FindObject('hBssPtFrag')
hBsPtFrag = gROOT.FindObject('hBsPtFrag')
hAntiBsPtFrag = gROOT.FindObject('hAntiBsPtFrag')
hBsToMuMuG = gROOT.FindObject('hBsToMuMuG')
hBsToBd = gROOT.FindObject('hBsToBd')
hBsToB = gROOT.FindObject('hBsToB')

hHad = gROOT.FindObject('hHad')
hBsParents = gROOT.FindObject('hBsParents')
hBsAncestor = gROOT.FindObject('hBsAncestor')
hBsProduct = gROOT.FindObject('hBsProduct')

hGBsGDeltaR = gROOT.FindObject('hGBsGDeltaR')
hGBsEta_GEta = gROOT.FindObject('hGBsEta_GEta')
hGammaBsPt = gROOT.FindObject('hGammaBsPt')

hMuPt = gROOT.FindObject('hMuPt')
hMuEta = gROOT.FindObject('hMuEta')
hMu1Eta_Mu2Eta = gROOT.FindObject('hMu1Eta_Mu2Eta')
hGammaPt = gROOT.FindObject('hGammaPt')
hGammaEta = gROOT.FindObject('hGammaEta')



SaveCanvas(cBPt, hBPt, 'BPt', rebin = 5,  norm = True)
SaveCanvas(cBsPt, hBsPt, 'BsPt', rebin = 5,  norm = True)
SaveCanvas(cBPtFrag, hBPtFrag, 'BPtFrag', rebin = 5,  norm = True)
SaveCanvas(cBssPtFrag, hBssPtFrag, 'BssPtFrag', rebin = 5,  norm = True)
SaveCanvas(cBsPtFrag, hBsPtFrag, 'BsPtFrag', rebin = 5,  norm = True)
SaveCanvas(cAntiBsPtFrag, hAntiBsPtFrag, 'AntiBsPtFrag', rebin = 5,  norm = True)
SaveCanvas(cBsToMuMuG, hBsToMuMuG, 'BsToMuMuG', log = True)
SaveCanvas(cBsToBd, hBsToBd, 'BsToBd', log = True)
SaveCanvas(cBsToB, hBsToB, 'BsToB', log = True)

SaveCanvas(cHad, hHad, 'Had', rebin=5, norm=True)
SaveCanvas(cBsParents, hBsParents, 'BsParents', norm=True)
SaveCanvas(cBsAncestor, hBsAncestor, 'BsAncestor', norm=True)
SaveCanvas(cBsProduct, hBsProduct, 'BsProduct', norm=True)

SaveCanvas(cGBsGDeltaR, hGBsGDeltaR, 'GBsGDeltaR')
SaveCanvas(cGBsEta_GEta, hGBsEta_GEta, 'GBsEta_GEta')
SaveCanvas(cGammaBsPt, hGammaBsPt, 'GammaBsPt')

SaveCanvas(cMuPt, hMuPt, 'MuPt')
SaveCanvas(cMuEta, hMuEta, 'MuEta')
SaveCanvas(cMu1Eta_Mu2Eta, hMu1Eta_Mu2Eta, 'Mu1Eta_Mu2Eta')
SaveCanvas(cGammaPt, hGammaPt, 'GammaPt')
SaveCanvas(cGammaEta, hGammaEta, 'GammaEta')

input('press enter to exit')