#!/usr/bin/env python3
import sys
import math
import numpy as np
import cmsstyle as CMS

from ROOT import gStyle, TFile, TCanvas, gROOT, TLatex

gStyle.SetOptStat(1110)
gStyle.SetOptTitle(0)  # Title off
gStyle.SetStatX(0.85)  # Set X position (right)
gStyle.SetStatY(0.88)  # Set Y position (top)
gStyle.SetStatW(0.2)   # Width of the stats box
gStyle.SetStatH(0.15)  # Height of the stats box

import uproot
import matplotlib.pyplot as plt
import mplhep as hep

hep.style.use("CMS")

# === 1. Wczytanie wszystkich histogramów ===
file = uproot.open("../../rootOutputs/histos_PhiGamma.root")
print(file.keys())
histograms_1D = {}
histograms_2D = {}

for key in file.keys():
    obj = file[key]
    if obj.classname.startswith("TH1"):
        histograms_1D[key] = obj
    elif obj.classname.startswith("TH2"):
        histograms_2D[key] = obj

print(f"Wczytano {len(histograms_1D)} histogramów 1D i {len(histograms_2D)} histogramów 2D.")

# === 2. Funkcja do parsowania tytułów z ROOT ===
def parse_labels(title):
    parts = title.split(";")
    if len(parts) >= 3:
        return parts[1], parts[2]  # xlabel, ylabel
    elif len(parts) == 2:
        return parts[1], "Events"
    return "", "Events"

# === 3. Rysowanie histogramów 1D ===
for name, hist in histograms_1D.items():
    values, edges = hist.to_numpy()
    title = hist.title or name
    xlabel, ylabel = parse_labels(title)  # Fetch labels
    
    fig, ax = plt.subplots()
    hep.histplot(values, bins=edges, ax=ax, label=title, histtype='fill', alpha=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    # Dodanie boxa z statystykami
    mean = np.mean(values)
    stddev = np.std(values)
    n_entries = len(values)
    
    stats_text = f"Mean: {mean:.2f}\nStd: {stddev:.2f}\nEntries: {n_entries}"
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=12, verticalalignment='top', 
            horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'))
    
    # CMS label
    hep.cms.label("Simulation", loc=2)
    
    fig.savefig(f"{name.replace(' ', '_')}.png")
    plt.close()

# === 4. Rysowanie histogramów 2D ===
for name, hist in histograms_2D.items():
    values, xedges, yedges = hist.to_numpy()
    title = hist.title or name
    xlabel, ylabel = parse_labels(title)  # Fetch labels
    
    fig, ax = plt.subplots()
    pcm = ax.pcolormesh(xedges, yedges, values.T, shading='auto', cmap='viridis')
    fig.colorbar(pcm, ax=ax, label="Counts")
    
    # Dodanie boxa z statystykami dla 2D
    mean = np.mean(values)
    stddev = np.std(values)
    n_entries = np.sum(values)
    
    stats_text = f"Mean: {mean:.2f}\nStd: {stddev:.2f}\nEntries: {n_entries:.0f}"
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=12, verticalalignment='top',
            horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'))
    
    # CMS label
    hep.cms.label("Simulation", loc=0)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    fig.savefig(f"{name.replace(' ', '_')}_2D.png")
    plt.close()

input('press enter to exit')
