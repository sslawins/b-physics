void fit()
{
    gROOT->SetBatch(kTRUE); // Run in batch mode to avoid GUI pop-ups
    TFile *f = TFile::Open("histosRecoPhoton.root");
    TH1D *h = (TH1D*)f->Get("hBsMass");

    h->Fit("gaus");
    TF1 *f1 = h->GetFunction("gaus");
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->SetRange(5.0, 6.0);
    f1->SetTitle("Gaussian Fit to Bs Mass Distribution");
    h->Draw();

    TCanvas *c1 = new TCanvas("c1", "Fit Result", 800, 600);
    c1->cd();
    h->Draw();
    f1->Draw("same");
    c1->SaveAs("fit_result.png");
    f->Close();
}