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

    TFile* f = new TFile("histosConvertedPhotons.root");
    TH1D* h = (TH1D*)f->Get("hBsMass");

    // TFile* f = new TFile("toyData.root");
    // TH1D* h = (TH1D*)f->Get("hBsMass__Massdistribution");

    RooRealVar x("Massdistribution", "M_{B^{0}_{s}} (GeV/c^{2})",4.0, 6.0);
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
    RooAddPdf backgroundPdf("backgroundPdf", "Background Model", RooArgList(expo, gauss), RooArgList(frac));

    
    RooRealVar nbkg("N_B", "nbkg", 5000, 0, 50000);
    RooExtendPdf extendedBackgroundPdf("extendedBackgroundPdf", "Extended Background Model", backgroundPdf, nbkg);
    RooRealVar nsig("N_S", "nsig", 1000, 0, 10000);
    RooExtendPdf extendedSignalPdf("extendedSignalPdf", "Extended Signal Model", john, nsig);
    
    RooAddPdf model("totalPdf", "Total Model", RooArgList(extendedSignalPdf, extendedBackgroundPdf));
    
    

    RooFitResult* fitRes = model.fitTo(data,Save(),NumCPU(8));

    // nsig.setVal(8);

    // // generate a toy dataset from the model
    // RooDataSet* toyData = model.generate(x, 100000);
    // // fill a histogram of the toy dataset
    // TH1* hToy = toyData->createHistogram("hBsMass", x, Binning(100));

    // // plot the toy dataset
    // TCanvas* cToy = new TCanvas("Toy Data","Toy Data",800,400);
    // hToy->Draw();
    // cToy->SaveAs("toyData.png");

    // // save the toy dataset to a file
    // TFile* fToy = new TFile("toyData.root", "RECREATE");
    // hToy->Write();
    // fToy->Close();


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

    RooWorkspace *w = new RooWorkspace("workspace", "workspace");

    RooStats::ModelConfig modelConfigSB("S+B Model");
    modelConfigSB.SetWorkspace(*w);
    modelConfigSB.SetPdf(model);
    modelConfigSB.SetParametersOfInterest(nsig);
    modelConfigSB.SetObservables(x);
    modelConfigSB.SetSnapshot(nsig);

    
    RooStats::ModelConfig modelConfigB("B-only Model");
    modelConfigB.SetWorkspace(*w);
    modelConfigB.SetPdf(extendedBackgroundPdf);
    modelConfigB.SetObservables(x);
    modelConfigB.SetParametersOfInterest(nsig);
    modelConfigB.SetNuisanceParameters(nbkg);  // N_B is a nuisance parameter
    modelConfigB.SetSnapshot(nsig);  // Snapshot for B-only (N_S=0)
    modelConfigB.SetSnapshot(nbkg);


    
    w->import(data);
    w->import(modelConfigSB);
    w->import(modelConfigB);
    w->import(model);
    w->writeToFile("workspace.root");

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
