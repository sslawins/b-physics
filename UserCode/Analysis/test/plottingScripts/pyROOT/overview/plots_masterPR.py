#!/usr/bin/env python3
import sys
import cmsstyle as CMS
import ROOT
from ROOT import TCanvas, TH1, TH2, TLatex, TFile

ROOT.gROOT.SetBatch(True)

# Styl globalny
ROOT.gStyle.SetOptStat(1110)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetStatFontSize(0.03)
ROOT.gStyle.SetStatX(0.85)
ROOT.gStyle.SetStatY(0.88)
ROOT.gStyle.SetStatW(0.2)
ROOT.gStyle.SetStatH(0.15)
ROOT.gStyle.SetPalette(ROOT.kViridis)

ROOT.gStyle.SetTitleFontSize(0.05)
ROOT.gStyle.SetLabelSize(0.045, "XYZ")
ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetStatFontSize(0.04)  # box statystyk
ROOT.gStyle.SetLegendFont(42) 

# CMS-style ustawienia
CMS.SetCmsText("Simulation")
CMS.SetExtraText("")
CMS.SetCmsTextFont(52)
CMS.SetCmsTextSize(0.75 * 0.76)
CMS.SetLumi("")

# Plik wejściowy
file = ROOT.TFile.Open("../../../rootOutputs/histos_PhiGamma.root")
f_out = ROOT.TFile("plotted_histograms.root", "RECREATE")
out_plots_path= "plots/"

keys = file.GetListOfKeys()

histos_1D = {}
histos_2D = {}

for key in keys:
    obj = key.ReadObj()
    if obj.InheritsFrom("TH1") and not obj.InheritsFrom("TH2"):
        histos_1D[obj.GetName()] = obj
    elif obj.InheritsFrom("TH2"):
        histos_2D[obj.GetName()] = obj

print(f"Wczytano {len(histos_1D)} histogramów 1D i {len(histos_2D)} histogramów 2D.")

histos_1D.get("hBssPt").GetXaxis().SetRangeUser(0, 35)
histos_1D.get("hBPt").GetXaxis().SetRangeUser(0, 35)

histos_1D.get("hCascade_RecoGenDeltaR").Rebin(2)
histos_1D.get("hCascade_RecoGenDeltaR").GetXaxis().SetRangeUser(0, 1.)

# === Histogramy 1D ===
for name, h in histos_1D.items():
    canv = ROOT.TCanvas(h.GetName() + "_canv", h.GetName(), 800, 600)
    canv.SetLeftMargin(0.15)
    canv.SetBottomMargin(0.15)

    if(h == histos_1D.get("hBsParents") or h == histos_1D.get("hInitial") or h == histos_1D.get("hBsAncestor")):
        canv.SetLogy()
        h.Scale(1.0 / h.Integral())
        h.SetStats(0)  # Wyłączenie statystyk dla histogramów normalizowanych

    if(h ==histos_1D.get("hCascade_RecoGenDeltaR")):
        canv.SetLogy()

    if(h ==histos_1D.get("hCascade")):
        canv.SetLogy()

    if(h == histos_1D.get("hPhiPt") ):
        h.GetXaxis().SetRangeUser(0, 30)
    if(h == histos_1D.get("hPhiDecayReco_Pt") or h == histos_1D.get("hPhiDecayGen_Pt")):
        h.GetXaxis().SetRangeUser(0, 20)

    if(h == histos_1D.get("hPhiDecay_RecoGenDeltaR") ):
        h.GetXaxis().SetRangeUser(0, 0.015)
    h.GetXaxis().SetTitleSize(0.05)  # domyślnie było 0.04
    h.GetYaxis().SetTitleSize(0.05)
    h.GetXaxis().SetLabelSize(0.045)
    h.GetYaxis().SetLabelSize(0.045)
    
    h.SetTitle(" ")
    h.SetLineColor(ROOT.kBlack)
    h.SetFillColorAlpha(19, 0.5)
    h.Draw("HIST")  # pierwszy raz bez statboxa
    ROOT.gPad.Update()

    # teraz narysuj jeszcze raz dla statboxa
    h.Draw("SAME")  # żeby ROOT dodał statbox
    ROOT.gPad.Update()
    stats = h.GetListOfFunctions().FindObject("stats")
    if stats:
        stats.SetTextSize(0.04)  # ✅ Zwiększony rozmiar czcionki

        # ✅ Hack: najpierw przesuń poza ekran, żeby wymusić przeliczenie
        stats.SetX1NDC(1.1)
        stats.SetX2NDC(1.2)
        ROOT.gPad.Modified()
        ROOT.gPad.Update()


        stats.SetX1NDC(0.65)
        stats.SetY1NDC(0.75)
        stats.SetX2NDC(0.90)
        stats.SetY2NDC(0.90)

    # Tekst CMS
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextAlign(11)
    cms_text.SetTextFont(52)
    cms_text.SetTextSize(0.05)
    cms_text.DrawLatex(0.21, 0.905, "Simulation")

    canv.Modified()
    canv.Update()
    canv.Print(f"{out_plots_path + h.GetName().replace(' ', '_')}.png")
    canv.Write()


# === Histogramy 2D ===
for name, h in histos_2D.items():
    if "Eta_Pt" in h.GetName():
        canv = ROOT.TCanvas(h.GetName() + "_canv", h.GetName(), 1200, 600)
    else:
        canv = ROOT.TCanvas(h.GetName() + "_canv", h.GetName(), 800, 600)
    canv.SetLeftMargin(0.15)  # Większy margines z lewej strony
    canv.SetRightMargin(0.15)  # Większy margines z prawej strony
    canv.SetBottomMargin(0.15)
    
    h.Draw("COLZ")  # Narysowanie histogramu
    ROOT.gPad.Update()  # Aktualizacja płótna, aby statbox został utworzony
    canv.Update()
    stats = h.GetListOfFunctions().FindObject("stats")  # Pobranie statboxa
    stats.SetTextSize(0.04)  # ✅ Zwiększony rozmiar czcionki

    # ✅ Hack: najpierw przesuń poza ekran, żeby wymusić przeliczenie
    stats.SetX1NDC(1.1)
    stats.SetX2NDC(1.2)
    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    # ✅ Potem ustaw docelową pozycję
    if "RecoPt_GenPt" in h.GetName() or "Eta_Eta" in h.GetName():
        stats.SetX1NDC(0.60)
        stats.SetY1NDC(0.15)
        stats.SetX2NDC(0.85)
        stats.SetY2NDC(0.40)
    else:
        stats.SetX1NDC(0.60)
        stats.SetY1NDC(0.65)
        stats.SetX2NDC(0.85)
        stats.SetY2NDC(0.90)

    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetZaxis().SetTitleSize(0.05)
    h.GetXaxis().SetLabelSize(0.045)
    h.GetYaxis().SetLabelSize(0.045)
    h.GetZaxis().SetLabelSize(0.045)
    if("hPhiDecay_RecoPt_GenPt" in h.GetName() ):
            h.GetZaxis().SetLabelSize(0.038)
    # Podpisy osi
    h.SetTitle(" ")
    h.GetZaxis().SetTitle("Events")
    # Dodajemy napis CMS ręcznie z wytycznymi
    cms_text = TLatex()
    cms_text.SetNDC()  # Ustawienie na koordynaty jednostkowe
    cms_text.SetTextAlign(11)  # Wyrównanie do lewej, góra

    # Dodatkowy tekst, np. "Simulation"
    cms_text.SetTextFont(52)  # Inny font dla tekstu dodatkowego
    cms_text.SetTextSize(0.05)  # Mniejszy rozmiar
    cms_text.DrawLatex(0.17, 0.905, "Simulation")
    canv.Modified()
    canv.Update()
    canv.Print(f"{out_plots_path + h.GetName().replace(' ', '_')}_2D.png")
    canv.Write()


# === Porównawcze nakładanie histogramów ===

overlay_sets = {
    "gamma_pt": {
        "names": ["hDecayGamma_Pt", "hDecayGammaReco_Pt"],
        "labels": ["Generated", "Reconstructed"],
        "colors": [ROOT.kOrange + 7, ROOT.kAzure + 2],
        "title": "Decay #gamma p_{T}"
    },
    "phi_pt": {
        "names": ["hPhiDecayGen_Pt", "hPhiDecayReco_Pt"],
        "labels": ["Generated", "Reconstructed"],
        "colors": [ROOT.kOrange + 7, ROOT.kAzure + 2],
        "title": "#phi decay p_{T}"
    },
    "cuts_gamma": {
        "names": ["hCutsDecayGammaGen_Pt", "hCutsDecayGammaReco_Pt", "hCutsDecayGammaGenMatched_Pt"],
        "labels": ["Generated", "Reconstructed", "Reconstructed and matched"],
        "colors": [ROOT.kGreen , ROOT.kOrange + 7, ROOT.kAzure + 2],
        "title": "Cuts on decay #gamma p_{T}"
    },
    "cascade_pt": {
        "names": ["hCascadeGen_PtAll","hCascadeGen_Pt" , "hCascadeReco_Pt"],
        "labels": ["Generated (all the events)", "Generated (events with 2 muons)", "Reconstructed (events with 2 muons)"],
        "colors": [ROOT.kGreen , ROOT.kOrange + 7, ROOT.kAzure + 2],
        "title": "Cascade p_{T}"
    }
}

for set_name, info in overlay_sets.items():
    histos = [histos_1D.get(name) for name in info["names"]]
    if not all(histos):
        print(f"Brakuje jednego z histogramów w zestawie: {set_name}")
        continue
    
    canv = ROOT.TCanvas(f"{set_name}_canv", info["title"], 800, 600)
    canv.SetLeftMargin(0.15)
    canv.SetBottomMargin(0.15)
    #canv.SetLeftMargin(0.12)

    positions = {
        "cascade_pt": (0.2, 0.72, 0.88, 0.88),
        "cuts_gamma": (0.4, 0.72, 0.88, 0.88),
        "phi_pt":     (0.62, 0.72, 0.88, 0.88),
        "gamma_pt":   (0.62, 0.72, 0.88, 0.88),
    }
    x1, y1, x2, y2 = positions.get(set_name, (0.62, 0.72, 0.88, 0.88))
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetMargin(0.15)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.05)


    for i, h in enumerate(histos):
        h.GetXaxis().SetTitleSize(0.05)
        h.GetYaxis().SetTitleSize(0.05)
        h.GetZaxis().SetTitleSize(0.05)
        h.GetXaxis().SetLabelSize(0.045)
        h.GetYaxis().SetLabelSize(0.045)
        h.GetZaxis().SetLabelSize(0.045)
        h.SetLineColor(info["colors"][i])
        h.SetFillColorAlpha(info["colors"][i], 0.4)
        h.SetTitle("")
        h.GetXaxis().SetTitle("p_{T} [GeV]")
        h.GetYaxis().SetTitle("Entries")
        h.SetStats(0)

        drawopt = "HIST" if i == 0 else "HIST SAME"
        h.Draw(drawopt)

        # Przykładowy wpis w legendzie (do uzupełnienia)
        leg.AddEntry(h, info["labels"][i], "f")  # <-- Tu możesz zmienić na np. "Generated", "Reconstructed" itd.

    leg.Draw()

    # CMS-style text
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextAlign(11)
    cms_text.SetTextFont(52)
    cms_text.SetTextSize(0.05)
    cms_text.DrawLatex(0.17, 0.905, "Simulation")

    canv.Update()
    canv.Print(f"{out_plots_path}{set_name}_overlay.png")
    canv.Write()

f_out.Close()
input("Press Enter to exit...")
