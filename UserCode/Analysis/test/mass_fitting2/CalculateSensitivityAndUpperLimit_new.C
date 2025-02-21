#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "RooRandom.h"
#include "TFile.h"


/*#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSystem.h"
#include "TROOT.h"
 
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
 
#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
 
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"*/ 

//https://root.cern/doc/master/StandardHypoTestDemo_8C.html
//https://root.cern/doc/master/StandardHypoTestInvDemo_8C.html


void CalculateSensitivityAndUpperLimit_new() {
    
    TFile file("workspace.root");
    RooWorkspace* w = (RooWorkspace*)file.Get("workspace");

    if (!w) {
        std::cerr << "Error: Workspace not found in file!" <<"\n";
        return;
    }

    
    RooAbsPdf* signalPdf = w->pdf("extendedSignalPdf");
    RooAbsPdf* extendedBackgroundPdf = w->pdf("extendedBackgroundPdf");
    RooRealVar* m_mumugamma = w->var("Massdistribution");
    RooRealVar* N_S = w->var("N_S");
    RooRealVar* N_B = w->var("N_B");

    RooDataHist* data = (RooDataHist*)w->data("dataset");

    if (!signalPdf || !extendedBackgroundPdf || !m_mumugamma || !N_S || !N_B) {
        std::cerr << "Error: Missing components in the workspace!" << std::endl;
        return;
    }

   
    RooStats::ModelConfig modelConfigSB("S+B Model");
    modelConfigSB.SetWorkspace(*w);
    modelConfigSB.SetPdf(*w->pdf("totalPdf"));
    modelConfigSB.SetParametersOfInterest(*N_S);
    modelConfigSB.SetObservables(*m_mumugamma);
    modelConfigSB.SetSnapshot(*N_S);

    
    RooStats::ModelConfig modelConfigB("B-only Model");
    modelConfigB.SetWorkspace(*w);
    modelConfigB.SetPdf(*extendedBackgroundPdf);
    modelConfigB.SetObservables(*m_mumugamma);
    modelConfigB.SetParametersOfInterest(*N_S);
    modelConfigB.SetNuisanceParameters(*N_B);  // N_B is a nuisance parameter
    N_S->setVal(0);
    modelConfigB.SetSnapshot(*N_S);  // Snapshot for B-only (N_S=0)
    modelConfigB.SetSnapshot(*N_B);
     

   
    //RooArgSet poi(*N_S);
    //modelConfigSB.SetSnapshot(poi);
    //modelConfigB.SetSnapshot(poi);
    //modelConfigB.SetSnapshot(NuisanceParameter);
   

    // N_S->setVal(0);  // Fix N_S = 0 for B hypothesis
    // N_S->setConstant(true);
    // RooAbsData* asimovData = RooStats::AsymptoticCalculator::GenerateAsimovData(*w->pdf("totalPdf"), *m_mumugamma);
    RooAbsData* asimovData = RooStats::AsymptoticCalculator::GenerateAsimovData(*w->pdf("totalPdf"), *m_mumugamma);
    // asimovData->Print("V");
    RooStats::ProfileLikelihoodCalculator plc(*data, modelConfigB);
    RooStats::SimpleLikelihoodRatioTestStat slrts (*w->pdf("extendedBackgroundPdf"), *w->pdf("totalPdf"));
    RooArgSet nullPOI(*N_S);
    double eval = slrts.Evaluate(*asimovData, nullPOI);
    //double eval = slrts.Evaluate (*asimovData, *N_S);
    
    RooStats::HypoTestResult* nullHypothesisTest = plc.GetHypoTest();
    if (!nullHypothesisTest) {
        std::cerr << "Error: Null hypothesis test failed!" <<"\n";
        return;
    }
    double CLb = nullHypothesisTest->CLb();
    std::cout << "Confidence Level for Background-only hypothesis (CLb): " << CLb <<"\n";

     
   
    cout << "-------------------------------------------------" <<  "\n";
    cout << "The p-value for the null is " << nullHypothesisTest->NullPValue() <<  "\n";
    cout << "Corresponding to a significance of " << nullHypothesisTest->Significance() << "\n";
    cout<< " the value of the test stat : likelihood ration: "<< eval << "\n";
    cout << "-------------------------------------------------\n\n" << "\n";
    
   
   RooFitResult* fitResult = w->pdf("totalPdf")->fitTo(*asimovData, RooFit::Save()); //, RooFit::PrintLevel(0));
    if (fitResult->status() != 0) {
        std::cerr << "Warning: Fit did not converge properly!" <<"\n";
    }

    double fitted_N_S = N_S->getVal();
    double fitted_N_B = N_B->getVal();

    std::cout << "Fitted N_S from Asimov dataset: " << fitted_N_S <<"\n";
    std::cout << "Fitted N_B from Asimov dataset: " << fitted_N_B <<"\n";
    //N_S->setMax(10000);


    RooStats::AsymptoticCalculator ac(*data, modelConfigB, modelConfigSB);
    ac.SetOneSided(true);
    // ac.SetPrintLevel(0);
    RooStats::HypoTestInverter inverter(ac, N_S);
    inverter.SetConfidenceLevel(0.95);
    // inverter.UseCLs(true);
    inverter.SetFixedScan(20, 0, 1500);
     RooStats::HypoTestInverterResult* result = inverter.GetInterval();
     if (!result) {
        std::cerr << "Error: HypoTest failed!" << std::endl;
        return;
    }
    std::cout<<" upper limit "<< result->UpperLimit()<< "\n";
    std::cout<<" CLs "<< result->ConfidenceLevel()<< "\n";


    TCanvas* c1 = new TCanvas("c1","c1");
    RooStats::HypoTestInverterPlot* plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",result);
    plot->Draw("CLb 2CL");  // plot also CLb and CLs+b
    c1->Draw();
    c1->SaveAs("HTI_Result_Plot.png");
    
    //delete result;
    // delete asimovData;
    file.Close();

   
    
    /*

    
    
    I checked that this is 0.5. A confidence level of 0.5 means that under the 
    background-only hypothesis, the observed data (or Asimov data) falls right in the middle of the expected background distribution.
    This is not an error—it's actually a good sign that the background model is being correctly evaluated. If CL_b =0 , 
    it would mean the background hypothesis is completely incompatible with the observed data.
    

    double expectedUpperLimit = result->UpperLimit();
    std::cout << "Expected Upper Limit on N_S at 95% CL: " << expectedUpperLimit << "\n";

    // Convert to branching fraction:
    double efficiency = 0.1;
    double N_Bs = 1.2e9; // More realistic estimate of Bs production in CMS
    double luminosity = 200; // in fb⁻¹
    double branchingFractionUpperLimit = expectedUpperLimit / (efficiency * N_Bs * luminosity);

    std::cout << "Expected Upper Limit on Branching Fraction at 95% CL: " << branchingFractionUpperLimit << "\n";

    
     double theoretical_BF = 1e-8;  // Theoretical branching fraction for BsToMuMuGamma decay
     double fraction = 100000.0 / 1000000.0;  // Scaling factor for 100k events
     double expected_N_S = fraction * sigma_MC * effectiveLumi;
     double upperLimit_BF = (upperLimit_N_S / expected_N_S) * theoretical_BF;
    // Effective luminosity calculation:
    double sigma_MC = 7.557e7; // pb
    double L_eq_1M = 1.323e-5; // 1/fb for 1M events
    double effectiveLumi = L_eq_1M * 1000.0; // Convert from 1/fb to 1/pb
    double expected_N_S = sigma_MC * effectiveLumi;
    delete result;
    delete asimovData;
    file.Close();
     
     
     
     RooStats::HypoTestInverter calc(ac);
     calc.SetConfidenceLevel(0.95);
     calc.UseCLs(true);
     calc.SetVerbose(true);
     RooStats::HypoTestInverterResult* result = calc.GetInterval();
     std::cout<<" upper limit "<< result->UpperLimit()<< "\n";
     
     */
}

int main() {
    CalculateSensitivityAndUpperLimit_new();
    return 0;
}
