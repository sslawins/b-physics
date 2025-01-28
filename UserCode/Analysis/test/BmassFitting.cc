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

#include "RooExponential.h"
#include "RooJohnson.h"

using namespace std;
using namespace RooFit;


void BmassFitting()
{
    TFile* f = new TFile("histosConvertedPhotons_copy.root");
    TH1D* h = (TH1D*)f->Get("hBsMass");

    RooRealVar x("x", "x", 3., 7.);
    RooDataHist data("data", "data", x, Import(*h));

    RooPlot *frame = x.frame(Title("Bs mass"));
    data.plotOn(frame, DataError(RooAbsData::SumW2));

    // signal
    RooRealVar mu("mu", "mu", 5.36, 4, 6);
    RooRealVar lambda("lambda", "lambda", 0.2, 0.0, 20);
    RooRealVar gamma("gamma", "gamma", 0.02, 0.0, 0.1);
    RooRealVar delta("delta", "delta", 0.733407, 0, 2);
    RooJohnson john("john", "john", x, mu, lambda, gamma, delta);

    // exp background
    RooRealVar conts("conts", "conts", -0.3, -0.5, -0.1);
    RooExponential expoBg("expoBg", "expoBG", x, conts);
    // pol background
    // Polynomial background model
    RooRealVar m1par("m1par", "m1par", 0.3, 0.6);
    RooRealVar m2par("m2par", "m2par", 0.5, 0.7);
    RooRealVar m3par("m3par", "m3par", 0.1, 0.5);
    RooGenericPdf polyBg("polyBg", "@0 + @1*@3 + @2*@3*@3", RooArgSet(m1par, m2par, m3par, x));

    // signal + background
    RooRealVar nsig("nsig", "nsig", 1000, 0, 10000);
    RooRealVar nbkg("nbkg", "nbkg", 5000, 0, 50000);
    RooAddPdf model("model", "model", RooArgList(john, expoBg), RooArgList(nsig, nbkg));

    // model.plotOn(frame);

    model.fitTo(data);

    model.plotOn(frame);
    frame->Draw();
}