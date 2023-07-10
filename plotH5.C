#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TObject.h>
#include <TString.h>
#include <iostream>

const TString inputFilename = "/home/rocio/Desktop/dimu2/hist.root"; // read file
const TString oFilename = "outputsorted5.pdf";

// Function to sort detector type using histogram name
TString GetDetectorType(const TString& histName) {
  if (histName.Contains("ecal", TString::kIgnoreCase)) {
    return "ECAL";
  } else if (histName.Contains("hcal", TString::kIgnoreCase)) {
    return "HCAL";
  } else if (histName.Contains("srd", TString::kIgnoreCase)) {
    return "SRD";
  } else if (histName.Contains("wcal", TString::kIgnoreCase)) {
    return "WCAL";
  } else if (histName.Contains("gem", TString::kIgnoreCase)) {
    return "GEM";
  } else if (histName.Contains("momentum", TString::kIgnoreCase)) {
    return "momentum";
  } else if (histName.Contains("mm", TString::kIgnoreCase)) {
    return "micromega";
  } else if (histName.Contains("veto", TString::kIgnoreCase)) {
    return "veto";
  } else if (histName.Contains("cellEnergyHisto", TString::kIgnoreCase)) { //all of these histos have form cellEnergyHisto_x_y in rccode.cc 
    return "cellEnergyHistograms";
  } else {
    return "Other";
  }
}

void plotH5(const char* filename = inputFilename, const char* outputFilename = oFilename) {
  gStyle->SetPalette(kRainBow);
  TApplication* app = new TApplication("app", NULL, NULL);

  // Open root file
  TFile* file = new TFile(filename, "READ");
  if (file->IsZombie()) {
    printf("Error opening ROOT file\n");
    return;
  }

  // Create the output PDF
  TCanvas canvas;
  canvas.Print(Form("%s[", outputFilename)); // Start of PDF

  // Group histograms by detector type
  std::map<TString, std::vector<TH1*>> histogramMap;
  histogramMap["Other"];

  // Iteration over objects
  TIter nextkey(file->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)nextkey())) {
    // Reading object
    TObject* obj = key->ReadObj();

    // Check for TH1 and TH2 classes
    if (obj->InheritsFrom(TH1::Class())) {
      if (TH1* hist1D = dynamic_cast<TH1*>(obj)) {
        // Handle 1D histograms
        // Extract detector type from histogram name
        TString histName = hist1D->GetName();
        TString detectorType = GetDetectorType(histName); // Call GetDetectorType function

        // Convert detectorType to lowercase for consistent sorting
        detectorType.ToLower();

        // Group histograms by detector type
        histogramMap[detectorType].push_back(hist1D);
      } else if (TH2* hist2D = dynamic_cast<TH2*>(obj)) {
        // Handle 2D histograms
        // Extract detector type from histogram name
        TString histName = hist2D->GetName();
        TString detectorType = GetDetectorType(histName);

        //Convert detectorType to lowercase for sorting
	 detectorType.ToLower();

        // Group histograms by detector type
        histogramMap[detectorType].push_back(hist2D);
      }
    }
  }

  // Print histogramMap
  for (const auto& pair : histogramMap) {
    const TString& detectorType = pair.first;
    const std::vector<TH1*>& histograms = pair.second;

    std::cout << "Detector Type: " << detectorType << std::endl;
    for (const auto& hist : histograms) {
      std::cout << "Histogram Name: " << hist->GetName() << std::endl;
    }
    std::cout << "-----------------------" << std::endl;
  }

  // Check if "cellEnergyHistograms" exist in histogramMap
  if (histogramMap.count("cellenergyhistograms") > 0) {
    const std::vector<TH1*>& cellEnergyHistograms = histogramMap["cellenergyhistograms"];
    const int numCellHistograms = cellEnergyHistograms.size();

    // Calculate the number of pages needed to display cellEnergyHistograms
    int numPages = (numCellHistograms + 35) / 36;  // Display 36 histograms per page (6x6 grid)
    int canvasSize = 1800;  // Size of each sub-canvas
    int subCanvasSize = canvasSize / 6;  // Size of each histogram sub-canvas

    for (int page = 0; page < numPages; ++page) {
      TCanvas cellCanvas("cellCanvas", "Cell Energy Histograms", canvasSize, canvasSize); //creation of different canvas for cellEnergyHistograms case
      cellCanvas.Divide(6, 6, 0, 0); 

      int startHistIndex = page * 36;
      int endHistIndex = std::min(startHistIndex + 36, numCellHistograms);

      for (int i = startHistIndex; i < endHistIndex; ++i) {
        cellCanvas.cd(i % 36 + 1);

        // Draw histogram
        if (TH1* hist1D = dynamic_cast<TH1*>(cellEnergyHistograms[i])) {
          hist1D->Draw("colz");
        }
      }

      // Save the cellCanvas to the output PDF
      cellCanvas.Print(outputFilename);
    }
  }

  // Draw and save other histograms of different case  grouped by detector type
  for (const auto& pair : histogramMap) {
    const TString& detectorType = pair.first;
    const std::vector<TH1*>& histograms = pair.second;

    // Skip cellEnergyHistograms since it has already been handled
    if (detectorType == "cellenergyhistograms") {
      continue;
    }

    int numHistograms = histograms.size();

    // Calculate the number of pages needed to display rest of histograms 
    int numPages = (numHistograms + 3) / 4;  // display 4 histograms per page (2x2 grid)
    int canvasSize = 800;  // Size of each histogram canvas

    for (int page = 0; page < numPages; ++page) {
      TCanvas histCanvas("histCanvas", "Histogram", canvasSize, canvasSize);
      histCanvas.Divide(2, 2);

      int startHistIndex = page * 4;
      int endHistIndex = std::min(startHistIndex + 4, numHistograms);

      for (int i = startHistIndex; i < endHistIndex; ++i) {
        histCanvas.cd(i % 4 + 1);

        // Draw histogram
        if (TH1* hist1D = dynamic_cast<TH1*>(histograms[i])) {
          hist1D->Draw("colz");
        } else if (TH2* hist2D = dynamic_cast<TH2*>(histograms[i])) {
          hist2D->Draw("colz");
        }
      }

      // Save the histCanvas to the output PDF
      histCanvas.Print(outputFilename);
    }
  }

  canvas.Print(Form("%s]", outputFilename));

  // Close root file
  delete file;
  delete app;
}
