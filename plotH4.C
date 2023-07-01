#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TObject.h>
#include <TString.h>
#include <iostream>

const TString inputFilename = "/home/rocio/Desktop/dimu2/example.root"; // read file
const TString oFilename = "outputsorted4.pdf";


// Function to sort detector type using histogram name

TString GetDetectorType(const TString& histName) {
  if (histName.Contains("ecal")) {
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
  } else {
    return "Other";
  }
}

void plotH4(const char* filename = inputFilename, const char* outputFilename = oFilename) {
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
        TString detectorType = GetDetectorType(histName); //call GetDetectorType function

        // Convert detectorType to lowercase for consistent sorting
        detectorType.ToLower();

        // Group histograms by detector type
        histogramMap[detectorType].push_back(hist1D);
      }
      else if (TH2* hist2D = dynamic_cast<TH2*>(obj)) {
        // Handle 2D histograms
        // Extract detector type from histogram name
        TString histName = hist2D->GetName();
        TString detectorType = GetDetectorType(histName);

        // Convert detectorType to lowercase for sorting
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

  // maximum number of histograms per page
  const int histogramsPerPage = 4;

  // Draw and save histograms grouped by detector type
  for (const auto& pair : histogramMap) {
    const TString& detectorType = pair.first;
    const std::vector<TH1*>& histograms = pair.second;

    // Calculate pages as a function of histo number
    int numPages = histograms.size() / histogramsPerPage;
    if (histograms.size() % histogramsPerPage != 0)
      numPages++; // add page if required

    for (int page = 0; page < numPages; page++) {
      // create  page in PDF for each detector type
      canvas.Clear();
      canvas.Divide(2, 2);

      int startIdx = page * histogramsPerPage;
      int endIdx = startIdx + histogramsPerPage;

      //
      const std::vector<TH1*>& otherHistograms = histogramMap["Other"];
      
      // Draw histograms for the current page
      for (int i = startIdx; i < endIdx && i < histograms.size(); ++i) {
        canvas.cd(i % histogramsPerPage + 1);

        // Draw histogram
        if (TH1* hist1D = dynamic_cast<TH1*>(histograms[i])) {
          hist1D->Draw();
        } else if (TH2* hist2D = dynamic_cast<TH2*>(histograms[i])) {
          hist2D->Draw("colz"); // NOT WORKING
        }
      }

          // Save page to PDF

      canvas.Print(outputFilename);
    }
     
  }

  canvas.Print(Form("%s]", outputFilename));

  // Close root file
  // file->Close(); /DEALLOATE MEMORY

  delete file;
  delete app;
}
