#include "TFile.h"
#include "TH1D.h"

void myPlot()
{
    gROOT->SetBatch();
    gStyle->SetOptStat(0);
    TFile* f = new TFile("histosMarcinsIdea.root", "READ");

    TProfile* Profile = (TProfile*)f->Get("hRecoVsGenEnergyProfile");
    Profile->Approximate();

    TCanvas* c = new TCanvas("c", "c", 800, 600);
    TLine* line = new TLine(0, 0, 35, 35);

    TFitResultPtr fr = Profile->Fit("pol1", "S");

    double a = fr->Parameter(1);
    double b = fr->Parameter(0);
    std::cout << "Fit parameters: a = " << a << ", b = " << b << std::endl;
    TGraph* g = new TGraph();

    for(int i = 0; i < Profile->GetNbinsX(); i++)
    {
        double x = Profile->GetBinCenter(i + 1);
        double y = (Profile->GetBinContent(i + 1) - b) / a; // Adjusted to match the fit line
        g->SetPoint(i, x, y);
    }

    Profile->Draw();
    line->Draw("same");
    g->SetLineColor(kGreen);
    g->Draw("same *");


    c->SaveAs("PlotsMarcin/EnergyProfile.png");
}
