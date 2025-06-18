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

    TFitResultPtr fr = Profile->Fit("pol1", "S", "", 10, 35);

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

    Profile->SetTitle("Photon Energy Profile;gen E^{#gamma} [GeV];reco E^{#gamma} [GeV]");
    Profile->Draw();
    line->Draw("same");
    g->SetLineColor(kGreen);
    g->Draw("same *");


    c->SaveAs("PlotsMarcin/EnergyProfile.png");
    //
    //
    c->Clear();
    TProfile* Profile2 = (TProfile*)f->Get("hRecoVsGenPtBsFrameMuonProfile");
    TLine* line2 = new TLine(0, 0, 3, 3);
    Profile2->Draw();
    line2->Draw("same");
    c->SaveAs("PlotsMarcin/PtBsFrameMuonProfile.png");
    //
    c->Clear();
    TProfile* Profile3 = (TProfile*)f->Get("hRecoVsGenPtBsFramePhotonProfile");
    TLine* line3 = new TLine(0, 0, 3, 3);
    Profile3->Draw();
    line3->Draw("same");
    c->SaveAs("PlotsMarcin/PtBsFramePhotonProfile.png");
    //
    c->Clear();
    TProfile* Profile4 = (TProfile*)f->Get("hRecoVsGenPtBsFramePhotonProfileCorrected");
    TLine* line4 = new TLine(0, 0, 3, 3);
    Profile4->Draw();
    line4->Draw("same");
    c->SaveAs("PlotsMarcin/PtBsFramePhotonProfileCorrected.png");
    //
    //
    TPaveText *pt = new TPaveText(.78, .6, .98, .77, "NDC");
    TH1D* hDimuonVertexXResidual = (TH1D*)f->Get("hDimuonVertexXResidual");
    fr = hDimuonVertexXResidual->Fit("gaus", "S");
    hDimuonVertexXResidual->SetTitle("X Residual Dimuon; x [cm]; counts");
    hDimuonVertexXResidual->Draw();
    pt = new TPaveText(.78, .6, .98, .77, "NDC");
    pt->AddText(TString::Format("#mu: %.3f", fr->Parameter(fr->Index("Mean"))));
    pt->AddText(TString::Format("#sigma: %.3f", fr->Parameter(fr->Index("Sigma"))));
    pt->Draw();
    c->SaveAs("PlotsMarcin/DimuonXResidual.png");
    c->Clear();
    //
    TH1D* hDimuonVertexZResidual = (TH1D*)f->Get("hDimuonVertexZResidual");
    fr = hDimuonVertexZResidual->Fit("gaus", "S");
    hDimuonVertexZResidual->SetTitle("Z Residual Dimuon; z [cm]; counts");
    hDimuonVertexZResidual->Draw();
    pt = new TPaveText(.78, .6, .98, .77, "NDC");
    pt->AddText(TString::Format("#mu: %.3f", fr->Parameter(fr->Index("Mean"))));
    pt->AddText(TString::Format("#sigma: %.3f", fr->Parameter(fr->Index("Sigma"))));
    pt->Draw();
    c->SaveAs("PlotsMarcin/DimuonZResidual.png");
    c->Clear();
    
}
