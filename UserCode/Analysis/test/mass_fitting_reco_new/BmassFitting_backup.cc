#include <iostream>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TFile.h"
#include "RooWorkspace.h"

#include "RooExponential.h"
#include "RooJohnson.h"

using namespace std;
using namespace RooFit;


void BmassFitting()
{
    gROOT->SetBatch();

/*
      // exp background
    RooRealVar conts("conts", "conts", -0.3, -0.5, -0.1);
    RooExponential expoBg("expoBg", "expoBG", x, conts);
    // pol background
    // Polynomial background model
    RooRealVar m1par("m1par", "m1par", 60, -100, 100);
    RooRealVar m2par("m2par", "m2par", 0, -1, 1);
    RooRealVar m3par("m3par", "m3par", 0, -0.5, 0.5);
    RooGenericPdf polyBg("backgroundPdf", "@0 + @1*@3 + @2*@3*@3", RooArgSet(m1par, m2par, m3par, x));
*/
    //gROOT->SetBatch();

    // signal
   /* RooRealVar mu("mu", "mu", 5.36, 4, 6);
    RooRealVar lambda("lambda", "lambda", 0.2, 0.0, 20);
    RooRealVar gamma("gamma", "gamma", 0.02, 0.0, 1.);
    RooRealVar delta("delta", "delta", 0.733407, 0, 20);
    RooJohnson john("signalPdf", "john", x, mu, lambda, gamma, delta);


    
    RooRealVar slope("slope", "Decay constant", -0.3, -10, 10);
    RooExponential expo("expo", "Exponential PDF", x, slope);

    RooRealVar mean("mean", "Gaussian Mean", 5.36, 4, 6);
    RooRealVar sigma("sigma", "Gaussian Sigma", 1, 0.1, 10);
    RooGaussian gauss("gauss", "Gaussian PDF", x, mean, sigma);
    RooRealVar frac("frac", "Fraction", 0.1, 0.0, 1.0);
    RooAddPdf backgroundPdf("backgroundPdf", "Background Model", RooArgList(expo, gauss), frac);



    //RooProdPdf backgroundPdf("backgroundPdf", "Background Model", RooArgSet(expo, gauss));

  
    // signal + background
    RooRealVar nsig("N_S", "nsig", 1000, 0, 10000);
    RooRealVar nbkg("N_B", "nbkg", 5000, 0, 50000);
    //RooAddPdf model("totalPdf", "model", RooArgList(john, polyBg), RooArgList(nsig, nbkg));
    RooAddPdf model("totalPdf", "model", RooArgList(john, backgroundPdf), RooArgList(nsig, nbkg));

    // model.plotOn(frame);

    model.fitTo(data);

    model.plotOn(frame);
    model.plotOn(frame, Components(backgroundPdf), LineStyle(kDashed));
    model.plotOn(frame, Components(john), LineStyle(kDashed), LineColor(kRed));
    double integral = backgroundPdf.createIntegral(x)->getVal();
    std::cout << "Integral of backgroundPdf over x: " << integral << std::endl;

    RooWorkspace *w = new RooWorkspace("workspace", "workspace");

    w->import(data);
    w->import(model);*/

    TFile* f = new TFile("histosRecoPhoton.root");
    TH1D* h = (TH1D*)f->Get("hBsMass");

    RooRealVar x("Massdistribution", "M_{B^{0}_{s}} (GeV/c^{2})",4.0, 7.0);
    RooDataHist data("dataset", "data", x, Import(*h));

    RooPlot *frame = x.frame(Title("M_{B^{0}_{s}} (GeV/c^{2})"));
    data.plotOn(frame, DataError(RooAbsData::SumW2));


    RooRealVar mu("mu", "mu", 5.2, 5.5);
    RooRealVar lambda("lambda", "lambda", 0.4, 0.2, 1.1);
    RooRealVar gamma("gamma", "gamma", 0.02, 0.0, 0.5);
    RooRealVar delta("delta", "delta", 5, 5, 100);
    RooJohnson john("signalPdf", "john", x, mu, lambda, gamma, delta);
    RooRealVar slope("slope", "Decay constant", -0.6, -20, 20);
    RooExponential expo("expo", "Exponential PDF", x, slope);
    RooRealVar mean("mean", "Gaussian Mean", 5.36, 4, 6);
    RooRealVar sigma("sigma", "Gaussian Sigma", 1, 0.1, 10);
    RooGaussian gauss("gauss", "Gaussian PDF", x, mean, sigma);
    RooRealVar frac("frac", "Fraction", 0.1, 0.0, 1.0);
    RooAddPdf backgroundPdf("backgroundPdf", "Background Model", RooArgList(expo), RooArgList(frac));

    
    RooRealVar nbkg("N_B", "nbkg", 5000, 0, 50000);
    RooExtendPdf extendedBackgroundPdf("extendedBackgroundPdf", "Extended Background Model", expo, nbkg);
    RooRealVar nsig("N_S", "nsig", 1000, 0, 10000);
    RooExtendPdf extendedSignalPdf("extendedSignalPdf", "Extended Signal Model", john, nsig);
    
    RooAddPdf model("totalPdf", "Total Model", RooArgList(extendedSignalPdf, extendedBackgroundPdf));
    
    

    RooFitResult* fitRes = model.fitTo(data,Save(),NumCPU(8));


    RooArgList finalPars = fitRes->floatParsFinal();
    for (int i = 0; i < finalPars.getSize(); ++i)
    {
        RooRealVar* par = (RooRealVar*)finalPars.at(i);
        std::cout << par->GetName() << " = " << par->getVal() 
        << " ± " << par->getError() << std::endl;
    }
    RooRealVar* nsig_fit = (RooRealVar*)fitRes->floatParsFinal().find("N_S");
    double nsigValue = nsig_fit->getVal();
    double nsigError = nsig_fit->getError();
    std::cout << "Fitted Nsig: " << nsigValue << " +- " << nsigError << std::endl;


    double L_eq_1M = 1.323e-5; // from production
    double cms_lumi = 200.0;  // in fb⁻¹
    double fraction = 100000.0 / 1000000.0;  // Scaling factor for 100k events
    double effective_lumi = fraction * L_eq_1M;
    double MC_N_S =  nsigValue / effective_lumi;
    double BF = 5.51e-9;  // Theoretical branching fraction for BsToMuMuGamma decay
    double data_N_S = MC_N_S * cms_lumi *BF ;
    std::cout<< " data N_S = " << data_N_S << "\n";
    std::cout << "Effective Luminosity: " << effective_lumi/BF << " fb^-1" << std::endl;


    fitRes->Print("v");
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    TH2* hcorr = fitRes->correlationHist() ;
    TCanvas* ccor = new TCanvas("Corr Matrix","mass par corr",800,400) ;
    gPad->SetLeftMargin(0.15) ; hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz") ;
    model.plotOn(frame);
    model.plotOn(frame, Components(extendedSignalPdf),  LineColor(kGreen), RooFit::Name("sig"), LineWidth(2), LineStyle(4));
    model.plotOn(frame, Components(extendedBackgroundPdf), LineColor(kRed), RooFit::Name("bkg"), LineWidth(2), LineStyle(2));
    
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(frame->findObject("sig"),"B_{s}#rightarrow#mu#mu#gamma","l");
    leg->AddEntry(frame->findObject("bkg"),"Combinatorial","l");

    double integral = extendedBackgroundPdf.createIntegral(x)->getVal();
    


    Double_t chisquare_mass = frame->chiSquare();
    std::cout << "Chi2/NdF - chisqure :  = " << chisquare_mass <<"\n";
    double checkIntegral = extendedBackgroundPdf.createIntegral(x)->getVal();
    std::cout << "Integral of extendedBackgroundPdf over x: " << checkIntegral << std::endl;


    TCanvas *c = new TCanvas();
    frame->Draw();
    leg->Draw("same");
    c->SaveAs("BmassFitting.png");



    
        /*RooRealVar mu("mu", "Signal Mean", 5.36, 4.0, 6.0);
        RooRealVar sigma("sigma", "Signal Width", 0.02, 0.1, 2.5);
        RooRealVar alpha("alpha", "Alpha", 1.5, -2.5, 10.0);
        RooRealVar n("n", "n", 5, 0.1, 20);
        RooCBShape signalPdf("signalPdf", "Crystal Ball Signal PDF", x, mu, sigma, alpha, n);

        RooRealVar slope("slope", "Decay constant", -0.3, -10, 10);
        RooExponential expo("expo", "Exponential PDF", x, slope);

        RooRealVar mean_bg("mean_bg", "Background Gaussian Mean", 5.0, 4.0, 6.0);
        RooRealVar sigma_bg("sigma_bg", "Background Gaussian Sigma", 1.0, 0.1, 10);
        RooGaussian gauss_bg("gauss_bg", "Background Gaussian PDF", x, mean_bg, sigma_bg);

        RooRealVar frac("frac", "Fraction", 0.1, 0.0, 1.0);
        RooAddPdf backgroundPdf("backgroundPdf", "Background Model", RooArgList(expo, gauss_bg), RooArgList(frac));

        
        RooRealVar nbkg("N_B", "Number of Background Events", 5000, 0, 50000);
        RooRealVar nsig("N_S", "Number of Signal Events", 1000, 0, 10000);

        
        RooExtendPdf extendedSignalPdf("extendedSignalPdf", "Extended Signal Model", signalPdf, nsig);
        RooExtendPdf extendedBackgroundPdf("extendedBackgroundPdf", "Extended Background Model", backgroundPdf, nbkg);

        
        RooAddPdf model("totalPdf", "Total Model", RooArgList(extendedSignalPdf, extendedBackgroundPdf));
        model.fitTo(data);
        
        std::cout << "Fitted Signal Yield (N_S): " << nsig.getVal() << " ± " << nsig.getError() << std::endl;

       
        
        model.plotOn(frame);
        model.plotOn(frame, Components(extendedBackgroundPdf), LineStyle(kDashed));
        model.plotOn(frame, Components(extendedSignalPdf), LineStyle(kDashed), LineColor(kRed));

        
        double integral_bkg = extendedBackgroundPdf.createIntegral(x)->getVal();
        double integral_sig = extendedSignalPdf.createIntegral(x)->getVal();

        std::cout << "Integral of backgroundPdf over x: " << integral_bkg << std::endl;
        std::cout << "Integral of signalPdf over x: " << integral_sig << std::endl;

        
        RooWorkspace w("workspace", "workspace");
        w.import(data);
        w.import(model);
        w.writeToFile("workspace.root");

        
        TCanvas *c = new TCanvas();
        frame->Draw();
        c->SaveAs("BmassFitting.png");*/

        
}
