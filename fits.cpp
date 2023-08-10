#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPDF.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
#include <TMath.h>
#include <TLegend.h>

void fits() {
  TFile *_file0 = TFile::Open("hist.root");
  TTree *t = (TTree *)_file0->Get("data");
  
  const int numBins = 1;
  double binEdges[numBins + 1];
  binEdges[0] = 0;
  binEdges[1] = 60;

  TPDF *pdfOutput = new TPDF("histograms.pdf");
  
  TH1D *histograms[numBins];
  TF1 *fg1mu = new TF1("fg1mu", "gaus", 0, 20);
  TF1 *fg2mu = new TF1("fg2mu", "gaus", 0, 20);
  TF1 *fg3mu = new TF1("fg3mu", "gaus", 0, 20);
  TF1 *fh2s = new TF1("fh2s", "gaus + gaus(3) + gaus(6)", 0, 20);
  
  // Define Erfc functions
  TF1 *erffg1mu = new TF1("erffg1mu", "0.5 * TMath::Erfc((x[0] - [0]) / (TMath::Sqrt(2) * [1]))", 0, 20);
  TF1 *erffg2mu = new TF1("erffg2mu", "0.5 * TMath::Erfc((x[0] - [0]) / (TMath::Sqrt(2) * [1]))", 0, 20);
  
  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600); // Declare canvas before the loop
  
  for (int i = 0; i < numBins; ++i) {
    TString histName = Form("h_%d", i);
    TString histTitle = Form("ecalE Range %d-%d", static_cast<int>(binEdges[i]), static_cast<int>(binEdges[i+1]));
    histograms[i] = new TH1D(histName, histTitle, 150, 0, 20);
    
    TString drawOption = Form("(%.1f<ecalE&&ecalE<%.1f)", binEdges[i], binEdges[i + 1]);
    t->Draw("ehcal2E>>" + histName, drawOption);
    
    int numHistogramBins = histograms[i]->GetNbinsX();
    double totalEnergyRange = 20.0; // Total energy range of the histogram
    double energyWidthPerBin = totalEnergyRange / numHistogramBins;
    
    std::cout << "Energy width per bin: " << energyWidthPerBin << " units" << std::endl;
    
    histograms[i]->Fit(fg2mu, "", "", 2, 10);
    
    std::cout << "Parameters of fg2mu after fit:" << std::endl;
    std::cout << "Mean: " << fg2mu->GetParameter(1) << std::endl;
    std::cout << "Sigma: " << fg2mu->GetParameter(2) << std::endl;
    
    fh2s->SetParLimits(7, 8, 10);
    fh2s->SetParLimits(8, fg2mu->GetParameter(2) * 0.7, fg2mu->GetParameter(2) * 1.3);
    fh2s->SetParameters(fg2mu->GetParameter(0), fg2mu->GetParameter(1), fg2mu->GetParameter(2),
                        fg2mu->GetParameter(0) * 0.3, 3., fg2mu->GetParameter(2),
                        fg2mu->GetParameter(0) * 0.1, 9, fg2mu->GetParameter(2));
    
    histograms[i]->Fit(fh2s, "BS", "", 0., 12.);
    
    fg1mu->SetParameters(fh2s->GetParameter(3), fh2s->GetParameter(4),  std::abs(fh2s->GetParameter(5)));
    fg2mu->SetParameters(fh2s->GetParameter(0), fh2s->GetParameter(1), fh2s->GetParameter(2));
    fg3mu->SetParameters(fh2s->GetParameter(6), fh2s->GetParameter(7), fh2s->GetParameter(8));
    
    std::cout << "Parameters of fg1mu:" << std::endl;
    std::cout << "Mean: " << fg1mu->GetParameter(1) << std::endl;
    std::cout << "Sigma: " << fg1mu->GetParameter(2) << std::endl;
    std::cout << "Sigma from fh2s: " << fh2s->GetParameter(5) << std::endl;
    
    std::cout << "Parameters of fg3mu:" << std::endl;
    std::cout << "Mean: " << fg3mu->GetParameter(1) << std::endl;
    std::cout << "Sigma: " << fg3mu->GetParameter(2) << std::endl;
    
    fg1mu->SetLineColor(kBlue + 2);
    fg2mu->SetLineColor(kGreen + 2);
    fg3mu->SetLineColor(kMagenta);
    
    
    histograms[i]->GetXaxis()->SetTitle("Energy (units)");
    histograms[i]->GetYaxis()->SetTitle("Counts");
    histograms[i]->Draw();
    fg1mu->Draw("same");
    fg2mu->Draw("same");
    fg3mu->Draw("same");
    
    gPad->Update();
    // Create a new canvas for the Erfc fits
    TCanvas *canvasErf = new TCanvas("canvasErf", "Erfc Fits", 800, 600);
    TString erfTitle = Form("Erfc fits for ecalE Range %d-%d", static_cast<int>(binEdges[i]), static_cast<int>(binEdges[i+1]));
    canvasErf->SetTitle(erfTitle);
    
    // Set parameters for Erfc
    erffg1mu->SetParameters(fg1mu->GetParameter(1), fg1mu->GetParameter(2));
    erffg2mu->SetParameters(fg2mu->GetParameter(1), fg2mu->GetParameter(2));
    
    erffg1mu->SetLineColor(kRed);
    erffg2mu->SetLineColor(kOrange);
    
    erffg1mu->Draw();
    erffg2mu->Draw("same");
    
    // legend
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(erffg1mu, "erffg1mu", "l");
    legend->AddEntry(erffg2mu, "erffg2mu", "l");
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0);
    legend->SetTextSize(0.03);
    legend->Draw();
    
    // Calculate area under fg1mu and fg2mu
    double a_fg1mu = -std::numeric_limits<double>::infinity(); // Lower limit for integration
    double b_fg1mu = std::numeric_limits<double>::infinity();  // Upper limit for integration
    double area_fg1mu = fg1mu->Integral(a_fg1mu, b_fg1mu);
    std::cout << "Area under fg1mu: " << area_fg1mu << std::endl;
    
    double a_fg2mu = -std::numeric_limits<double>::infinity(); // Lower limit for integration
    double b_fg2mu = std::numeric_limits<double>::infinity();  // Upper limit for integration
    double area_fg2mu = fg2mu->Integral(a_fg2mu, b_fg2mu);
    std::cout << "Area under fg2mu: " << area_fg2mu << std::endl;
    
    
    // Calculate the energy represented by the areas
    //  double energy_fg1mu = area_fg1mu * energyWidthPerBin;
    // double energy_fg2mu = area_fg2mu * energyWidthPerBin;
    
    // std::cout << "Energy represented by fg1mu: " << energy_fg1mu << " units" << std::endl;
    //std::cout << "Energy represented by fg2mu: " << energy_fg2mu << " units" << std::endl;
    
    // calculate quantity of events represented by the areas
    double totalEvents = histograms[i]->GetEntries();
    double events_fg1mu = area_fg1mu / energyWidthPerBin;
    double events_fg2mu = area_fg2mu / energyWidthPerBin;
    
    std::cout << "Number of events under fg1mu: " << events_fg1mu << std::endl;
    std::cout << "Number of events under fg2mu: " << events_fg2mu << std::endl;
    
    // Calculate the signal-to-noise ratio (SNR)
    double snr = area_fg2mu / area_fg1mu;
    std::cout << "Signal-to-Noise Ratio (SNR): " << snr << std::endl;
    
    
    
    // Calculate the optimal x-axis value using Erfc functions
    double optimalX = -1;
    double maxAreaDifference = -std::numeric_limits<double>::infinity();
    for (double x = 0; x <= 20; x += energyWidthPerBin) {
      double a_erffg1mu_x = erffg1mu->Eval(x);
      double a_erffg2mu_x = erffg2mu->Eval(x);
      double areaDifference_x = a_erffg2mu_x - a_erffg1mu_x;
      
      if (areaDifference_x > maxAreaDifference) {
	maxAreaDifference = areaDifference_x;
	optimalX = x;
      }
    }
    
    std::cout << "Optimal x for maximum area difference: " << optimalX << " units" << std::endl;
    
    
    canvasErf->Print("erf_fits.pdf", "pdf");
    delete legend;
    delete canvasErf;
    
    if (i < numBins - 1) {
      canvas->Print("histograms.pdf(", "pdf");
    } else {
      canvas->Print("histograms.pdf)", "pdf");
    }
  }
  
  delete canvas;
  delete pdfOutput;
  delete fg1mu;
  delete fg2mu;
  delete fg3mu;
  delete fh2s;
  for (int i = 0; i < numBins; ++i)
    delete histograms[i];
  
  _file0->Close();
}

