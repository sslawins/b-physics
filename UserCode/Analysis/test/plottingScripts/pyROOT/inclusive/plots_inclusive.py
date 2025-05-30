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

# CMS-style ustawienia
CMS.SetCmsText("Simulation")
CMS.SetExtraText("")
CMS.SetCmsTextFont(52)
CMS.SetCmsTextSize(0.75 * 0.76)
CMS.SetLumi("")

# Plik wejściowy
file = ROOT.TFile.Open("../../../histos_Phi_inclusive_G_18_03.root")
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



# === Histogramy 1D ===
for name, h in histos_1D.items():
    canv = ROOT.TCanvas(h.GetName() + "_canv", h.GetName(), 800, 600)
    canv.SetLeftMargin(0.12)

    if(h ==histos_1D.get("hCascade_RecoGenDeltaR")):
        canv.SetLogy()

    h.SetTitle(" ")
    h.SetLineColor(ROOT.kBlack)
    h.SetFillColorAlpha(19, 0.5)
    h.Draw()

    # Tekst CMS
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextAlign(11)
    cms_text.SetTextFont(52)
    cms_text.SetTextSize(0.04)
    cms_text.DrawLatex(0.18, 0.905, "Simulation")

    canv.Modified()
    canv.Update()
    canv.Print(f"{out_plots_path + h.GetName().replace(' ', '_')}.png")
    canv.Write()


# === Histogramy 2D ===
for name, h in histos_2D.items():
    canv = ROOT.TCanvas(h.GetName() + "_canv", h.GetName(), 800, 600)
    #canv.SetLeftMargin(0.15)  # Większy margines z lewej strony
    canv.SetRightMargin(0.15)  # Większy margines z prawej strony

    h.Draw("COLZ")  # Narysowanie histogramu
    ROOT.gPad.Update()  # Aktualizacja płótna, aby statbox został utworzony
    canv.Update()
    stats = h.GetListOfFunctions().FindObject("stats")  # Pobranie statboxa
    if "RecoPt_GenPt" in h.GetName():
        stats.SetX1NDC(0.60)  # Lewy dolny róg w układzie NDC (współrzędna X)
        stats.SetY1NDC(0.15)  # Lewy dolny róg w układzie NDC (współrzędna Y)
        stats.SetX2NDC(0.8)   # Prawy górny róg w układzie NDC (współrzędna X)
        stats.SetY2NDC(0.3)  # Dolna krawędź w układzie NDC
        ROOT.gPad.Modified()  # Zaktualizowanie płótna po zmianach
    
    # Podpisy osi
    h.SetTitle(" ")
    h.GetZaxis().SetTitle("Events")
    # Dodajemy napis CMS ręcznie z wytycznymi
    cms_text = TLatex()
    cms_text.SetNDC()  # Ustawienie na koordynaty jednostkowe
    cms_text.SetTextAlign(11)  # Wyrównanie do lewej, góra

    # Dodatkowy tekst, np. "Simulation"
    cms_text.SetTextFont(52)  # Inny font dla tekstu dodatkowego
    cms_text.SetTextSize(0.04)  # Mniejszy rozmiar
    cms_text.DrawLatex(0.16, 0.905, "Simulation")
    canv.Modified()
    canv.Update()
    canv.Print(f"{out_plots_path + h.GetName().replace(' ', '_')}_2D.png")
    canv.Write()

#TEST
# Tworzymy dwa histogramy
h1 = ROOT.TH1F("h1", "Histogram 1", 100, -10, 10)
h2 = ROOT.TH1F("h2", "Histogram 2", 500, 0, 10)

# Wypełniamy losowo
h1.FillRandom("gaus", 10000)
h2.FillRandom("gaus", 800000)

# Ustawienia histogramu 1 (szary, półprzezroczysty)
h1.SetLineColor(ROOT.kBlack)
h1.SetFillColorAlpha(ROOT.kYellow, 0.5)  # szary
h1.SetTitle("")
h1.GetXaxis().SetTitle("X")
h1.GetYaxis().SetTitle("Entries")

# Ustawienia histogramu 2 (niebieski, półprzezroczysty)
h2.SetLineColor(ROOT.kBlue + 1)
h2.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)

# Canvas
c = ROOT.TCanvas("c", "Two Histograms", 800, 600)
c.SetLeftMargin(0.1)

# Rysujemy histogramy
h1.Draw("HIST")
h2.Draw("HIST SAME")

# Legenda
leg = ROOT.TLegend(0.65, 0.7, 0.88, 0.85)
leg.AddEntry(h1, "Histogram 1 (gray)", "f")
leg.AddEntry(h2, "Histogram 2 (blue)", "f")
leg.Draw()

c.Update()
c.Print(f"{out_plots_path }two_histograms.png")
phiMass_info = {
    "names": ["hGenPhotonFourMomentaMass", "hRecoPhotonFourMomentaMass", "hRecoKaons_RecoPhotonDir_GenPhotonMag", "hRecoKaons_GenPhotonDir_RecoPhotonMag"],
    "labels": ["Fully generated #gamma", "Fully reconstructed #gamma", "Reconstructed direction,", "Reconstructed energy,", "generated energy", "generated direction"],
    "colors": [ROOT.kOrange + 7, ROOT.kAzure + 2, ROOT.kGreen, ROOT.kBlack],
    "title": "Phi Mass"
}

histos = [histos_1D.get(name) for name in phiMass_info["names"]]
if not all(histos):
    print(f"Brakuje jednego z histogramów w zestawie: phiMass")
    print(f"Zestaw: phiMass, brakujące histogramy: {[name for name, h in zip(phiMass_info['names'], histos) if h is None]}")
else:
    canv = ROOT.TCanvas("phiMass_canv", phiMass_info["title"], 800, 600)
    canv.SetLeftMargin(0.15)
    canv.SetBottomMargin(0.15)

    leg = ROOT.TLegend(0.49, 0.55, 0.88, 0.88)
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
        h.SetLineColor(phiMass_info["colors"][i])
        h.SetFillColorAlpha(phiMass_info["colors"][i], 0.4)
        h.SetTitle("")
        h.GetXaxis().SetTitle("M_{inv} [GeV]")
        h.GetYaxis().SetTitle("Events")
        h.SetStats(0)

        # Ustaw zakres tylko dla pierwszego histogramu przed pierwszym Draw
        if i == 0:
            h.GetXaxis().SetRangeUser(4.75, 6.5)

        drawopt = "HIST" if i == 0 else "HIST SAME"
        h.Draw(drawopt)

        leg.AddEntry(h, phiMass_info["labels"][i], "f")
        if i >= 2:
            leg.AddEntry(h, phiMass_info["labels"][i+2], "")

    leg.Draw()

    # CMS-style text
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextAlign(11)
    cms_text.SetTextFont(52)
    cms_text.SetTextSize(0.05)
    cms_text.DrawLatex(0.17, 0.905, "Simulation")

    canv.Update()
    canv.Print(f"{out_plots_path}phiMass_overlay.png")
    canv.Write()
f_out.Close()
input("Press Enter to exit...")
