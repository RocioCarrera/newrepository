// DaqDataDecoding
#include "DaqEventsManager.h"
#include "DaqEvent.h"
using namespace CS;

// ROOT
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include "TStyle.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TCanvas.h"

// c++
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

// P348 reco
#include "p348reco.h"
#include "shower.h"
#include "badburst.h"
#include "tracking.h"
#ifndef TRACKINGTOOLS
#include "MManager.h"
#endif

// timing
#include "timecut.h"

// tracking tools
#ifdef TRACKINGTOOLS
#include "tracking_new.h"
#endif

// for MC with tracking-tools
#include "simu.h"

int procev(RecoEvent &e);
void endrun();


int ievproc=-1;
int ifmc=0;
int ISetup;
int maxstep=0;
string stepname[20];

DaqEventsManager manager;

ofstream* foutd;

TH1D* momentum;

TH2I* ehplot;
TH1I* srdplot;
TH1I* srdplot0;
TH1I* srdplot1;
TH1I* srdplot2;
TH1I* mcecalentries;


//new histos
TH2F* cellEnergyHist; //2d array
TH1D*** cellEnergyHistograms;
 
TH1D* ecalratio;
TH1D* ecalcm;
TH1D* mcecaly;

TH1D* ecalrad;
TH1D* ecalrad1;
TH1D* chi2profile;

TH1D* eecal;
TH1D* eecalnoedge;
TH1D* eecalcentral;
TH1D* eecalcentraltot;
TH1D* eecalweighted;
TH1D* etotecal;
TH1D* epileup;
TH1D* epileup33;
TH1D* eprs;
TH1D* prstoecal;
TH1D* eveto;
TH1D* evhcal;
TH1D* eecalbad;

TH1D* etotwcal;
TH1D* wpileup;
TH1D* wpileup0;
TH1D* wpileup1;
TH1D* wpileup2;
TH1D* ewecal;
TH1D* ewprs;
TH1D* ewecalmain;
TH1D* ecatcher;
TH1D* ewecalecal;
TH2D* eweplot;
TH2D* eweplot1;
TH1D* eweplot2;
TH1D* eweplot3;
//TH2D* ewprsmain;
TH1D* evtwc;
TH1D* evtec;
TH1D* es2;
TH1D* ew2;
TH1D* ew22;
TH1D* ew23;
TH1D* ew2all;
TH1D* ev2;
TH1D* es4;
TH1D* edm;
TH2D* ev2w2;

TH1D* flow;

TH1D* ts0;
TH1D* tsrd;
TH1D* tdiffsrd;
TH1D* tecal;
TH1D* tdiffecal;
TH1D* thcal;
TH1D* tdiffhcal;
TH1D* tdet;
TH1D* tdiffdet;

TH1D* ehcal;
TH1D* ehcal0;
TH1D* ehcal1;
TH1D* ehcal2;
TH1D* ehcal3;
TH1D* ehcalcatcher;
TH2D* ehcal0veto;
TH1D* Rhcal;
TH1D* Rhcal12;
TH1D* Rhcal1;
TH1D* etot;

TH1D* ehcell0;
TH1D* ehcell1;
TH1D* ehcell2;
TH1D* ehcell3;
TH1D* ehcell4;
TH1D* ehcell5;
TH1D* ehcell6;
TH1D* ehcell7;
TH1D* ehcell8;

TH1F* cm0;
TH1F* cm5;
TH1F* cm10;
TH1F* cm15;

TH1D* ebgo0;
TH1D* ebgo3;
TH1D* ebgo4;
TH1D* ebgo7;

TH1D* hdata;
TH1D* hMC;
TFile* myfile_data;
TFile* myfile_MC;



// time histos
#include "timehisto.inc"

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: ./rccode.exe <file name>" << std::endl;
    return 1;
  }
  
  int NevBadBursts = 0;

  int NPassed = 0;

  std::string inpfile(argv[1]);
  int inputtype = 0;
  if( inpfile.find(std::string(".d")) == inpfile.size()-2 ) inputtype = 1;

  // ROOT calls
  TFile* hOutputFile = new TFile("hist.root", "RECREATE");

  if(inputtype == 0) { // Reading data file(s) ---------------------------

    std::cout << "Reading data files, first is " << argv[1] << std::endl;

//    DaqEventsManager manager;
    manager.SetMapsDir("../maps");
    for (int i = 1; i < argc; ++i)
      manager.AddDataSource(argv[i]);
    manager.Print();

    int maxev = 1000000;
    int nev = 0;

    // DST file
    foutd = new std::ofstream("mydsti.d");
    if (!foutd) {
      cerr << "ERROR: failed to open DST file" << endl;
      return 1;
    }

    // Output file
    ofstream fout("output.dat");
    if (!fout) {
      cerr << "ERROR: failed to open events output file" << endl;
      return 1;
    }


    //while (manager.ReadEvent() && nev < maxev) { // event loop, read event by event
    while (manager.ReadEvent()) { // event loop, read event by event
 
      nev++;

      const int nevt = manager.GetEventsCounter();

      // print progress
      if (nevt % 1000 == 1)
        cout << "===> Event #" << manager.GetEventsCounter() << "  N passed = " << NPassed << endl;

      // decode event (prepare digis)
      const bool decoded = manager.DecodeEvent();

      // skip events with decoding problems
      if (!decoded) {
        cout << "WARNING: fail to decode event #" << nevt << endl;
        continue;
      }

      // run reconstruction
      RecoEvent e = RunP348Reco(manager);

      // process only "physics" events (no random or calibration trigger)
      if (!e.isPhysics) continue;

      // skip known "bad" spills
      if (IsBadBurst(e.run, e.spill)) {NevBadBursts++; continue;}

      int ifproc = procev(e);
      NPassed += ifproc;

      // Output file
//      if(ifproc) fout.write((const char*) manager.GetEvent().GetBuffer(), manager.GetEvent().GetLength());

    }
    std::cout << "Number of events in bad bursts = " << NevBadBursts << std::endl;
  }  // End of data file(s) reading -------------------------------------

  if(inputtype == 1) { // Reading MC file -------------------------------

    std::cout << "Reading MC file " << argv[1] << std::endl;

    std::ifstream inFile(argv[1]);

    if (!inFile) {
      std::cerr << "ERROR: can't open file " << argv[2] << std::endl;
      return 1;
    }
  
    const int Nevents = 100000000;
    int NevMCRead = 0;

    // event loop, read event by event
    for (int iev = 0; iev < Nevents; iev++) {

      RecoEvent e;
      //if(inputtype == 1) e = RunP348RecoMC(inFile); // same as in exampletrackingmc in master
      if(inputtype == 1) e = RunP348RecoMC(inFile, 1, 1); // (&inFile, bool doDigi = 1, bool doSmear = 0)

      if (!e.mc) break;

      NevMCRead++;

      // print progress
      if (iev % 100 == 1)
        std::cout << "===> Event #" << iev << "  N passed = " << NPassed << std::endl;

      int ifproc = procev(e);
      NPassed += ifproc;

    }
    std::cout << std::endl;
    std::cout << "Number of MC events read in = " << NevMCRead << std::endl;
  } // End of MC file reading -------------------------------------------

  endrun();

  hOutputFile->Write();

  return 0;
}


int procev(RecoEvent &e0)
{
  ievproc++;
  if(ievproc == 0) { // initialization

    momentum = new TH1D("momentum", "momentum from MM", 400, 0., 200.);

    ehplot = new TH2I("ehplot", "ECAL vs HCAL;ECAL, GeV;HCAL, GeV;#nevents", 120, 0, 120, 120, 0, 120);
    srdplot = new TH1I("srdplot", "SRD;Energy, MeV;#nevents", 250, 0, 250);
    srdplot0 = new TH1I("srdplot0", "SRD;Energy, MeV;#nevents", 250, 0, 250);
    srdplot1 = new TH1I("srdplot1", "SRD;Energy, MeV;#nevents", 250, 0, 250);
    srdplot2 = new TH1I("srdplot2", "SRD;Energy, MeV;#nevents", 250, 0, 250);

    mcecalentries = new TH1I("MC_ecalentries", "MC ECAL entry X coordinates in mm;X, mm;#nevents", 100, -200., 200.);

    ecalratio = new TH1D("ecalratio", "E1/E2-1 vs coord", 39, -19.5, 19.5);
    ecalratio->Sumw2();
    ecalcm = new TH1D("ecalcm", "ECAL CM vs coord", 39, -19.5, 19.5);
    ecalcm->Sumw2();
    mcecaly = new TH1D("MC_ecalY", "MC ECAL entry Y coordinates in mm;X, mm;#nevents", 39, -19.5, 19.5);
    mcecaly->Sumw2();

    //new histos
    cellEnergyHist = new TH2F("cellEnergyHist", "Cell Energy Distribution", 7, -1, 6, 7, -1, 6);

    cellEnergyHistograms = new TH1D**[6];

    for (int x = 0; x < 6; ++x) {
      cellEnergyHistograms[x] = new TH1D*[6];

      for (int y = 0; y < 6; ++y) {
	std::string histName = "CellEnergyHist_" + std::to_string(x) + "_" + std::to_string(y);
	cellEnergyHistograms[x][y] = new TH1D(histName.c_str(), histName.c_str(), 100, 0, 100);
      }
    }




    ecalrad = new TH1D("ecalrad", "Shower radius", 50, 0., 100.);
    ecalrad1 = new TH1D("ecalrad1", "Shower radius 1", 50, 0., 100.);
    chi2profile = new TH1D("chi2profile", "ECAL profile chi2", 100, 0, 10.);

    eecal = new TH1D("eecal", "ECAL energy", 1000, 0., 200.);
    eecalnoedge = new TH1D("eecalnoedge", "ECAL 5*5 energy", 1000, 0., 200.);
    eecalcentral = new TH1D("eecalcentral", "ECAL 3*3 energy", 1000, 0., 200.);
    eecalcentraltot = new TH1D("eecalcentraltot", "ECAL 3*3 total energy", 1000, 0., 200.);
    eecalweighted = new TH1D("eecalweighted", "ECAL energy, weight used", 200, 0., 200.);
	etotecal = new TH1D("etotecal", "ECAL total energy", 1000, 0., 200.);
    epileup = new TH1D("epileup", "ECAL energy pileup", 200, 0., 40.);
	epileup33 = new TH1D("epileupcentral", "ECAL energy pileup 3*3", 200, 0., 40.);
    eprs = new TH1D("eprs", "PRS energy", 800, 0., 20.);
    prstoecal = new TH1D("prstoecal", "PRS/ECAL", 100, 0., 0.5);
    eveto = new TH1D("eveto", "VETO energy", 100, 0., 0.1);
    evhcal = new TH1D("evhcal", "VHCAL energy", 100, 0., 10.);
    eecalbad = new TH1D("eecalbad", "ECAL bad energy", 100, 0., 100.);

    etotwcal = new TH1D("etotwcal", "WCAL total energy", 1000, 0., 200.);
    wpileup = new TH1D("wpileup", "WCAL energy pileup", 200, 0., 40.);
    wpileup0 = new TH1D("wpileup0", "WCAL0 energy pileup", 200, 0., 40.);
    wpileup1 = new TH1D("wpileup1", "WCAL1 energy pileup", 200, 0., 40.);
    wpileup2 = new TH1D("wpileup2", "WCAL2 energy pileup", 200, 0., 40.);

    ewecal = new TH1D("ewecal", "WECAL energy", 1000, 0., 200.);
    ewprs = new TH1D("ewprs", "WECAL-PRS energy", 1000, 0., 20.);
    ewecalmain = new TH1D("ewecalmain", "WECAL main energy", 1000, 0., 200.);
    ecatcher = new TH1D("ecatcher", "W catcher energy", 200, 0., 10.);
    ewecalecal = new TH1D("ewecal+ecal", "WECAL+ECAL energy", 200, 0., 200.);
    eweplot = new TH2D("eweplot", "ECAL vs WCAL", 80, 0., 160., 80, 0., 160.);
    eweplot1 = new TH2D("eweplot1", "ECAL+HCAL vs WCAL", 90, 0., 180., 90, 0., 180.);
    eweplot2 = new TH1D("eweplot2", "WECAL energy wgt", 50, 0., 100.);
    eweplot3 = new TH1D("eweplot3", "WECAL energy", 50, 0., 100.);
    //ewprsmain = new TH2D("ewprsmain", "WPRS vs WMAIN", 100, 0., 200., 40, 0., 40.);
    evtwc = new TH1D("evtwc", "VTWC energy", 100, 0., 0.02);
    evtec = new TH1D("evtec", "VTEC energy", 100, 0., 0.01);
    es2 = new TH1D("es2", "S2 energy", 100, 0., 0.01);
    ew2 = new TH1D("ew2", "wcal2 energy", 100, 0., 0.01);
    ew22 = new TH1D("ew22", "wcal2 energy W<50", 100, 0., 0.01);
    ew23 = new TH1D("ew23", "wcal2 energy WPRS>8", 100, 0., 0.01);
    ew2all = new TH1D("ew2all", "wcal2 all events", 100, 0., 0.01);
    ev2 = new TH1D("ev2", "V2 energy", 100, 0., 0.02);
    es4 = new TH1D("es4", "S4 energy", 100, 0., 0.002);
    edm = new TH1D("edm", "DM energy", 100, 0., 0.01);
    ev2w2 = new TH2D("ev2w2", "WCAL[2] vs V2", 50, 0., 0.02, 50, 0., 0.01);

    flow = new TH1D("flow", "cut flow", 20, -0.5, 19.5);

    ts0 = new TH1D("ts0", "Time S0", 400, -100., 300.);
    tsrd = new TH1D("tsrd", "Time SRD0", 400, -100., 300.);
    tdiffsrd = new TH1D("tsrd-ts0", "Time diff SRD", 400, -200., 200.);
    tecal = new TH1D("tecal", "Time ECAL cell", 400, -200., 200.);
    tdiffecal = new TH1D("tecal-ts0", "Time diff ecal", 400, -200., 200.);
    thcal = new TH1D("thcal", "Time HCAL cell", 400, -200., 200.);
    tdiffhcal = new TH1D("thcal-ts0", "Time diff hcal", 400, -200., 200.);
    tdet = new TH1D("tdet", "Time some det", 400, -100., 300.);
    tdiffdet = new TH1D("tdet-ts0", "Time some det", 400, -200., 200.);

    ehcal = new TH1D("ehcal", "HCAL energy", 400, 0., 200.);
    ehcal0 = new TH1D("ehcal0", "HCAL module 0 energy", 200, 0., 40.);
    ehcal1 = new TH1D("ehcal1", "HCAL module 1 energy", 40, 0., 20.);
    ehcal2 = new TH1D("ehcal2", "HCAL module 2 energy", 40, 0., 20.);
    ehcal3 = new TH1D("ehcal3", "HCAL module 3 energy", 80, 0., 40.);
    ehcalcatcher = new TH1D("ehcalcatcher", "HCAL + catcher energy", 400, 0., 200.);
    ehcal0veto = new TH2D("ehcal0veto", "VETO vs ehcal0", 100, 0., 20., 20, 0., 40.);
    Rhcal  = new TH1D("Rhcal", "HCAL R value", 50, 0., 1.);
    Rhcal12 = new TH1D("Rhcal12", "HCAL 1-2 R value", 50, 0., 1.);
    Rhcal1 = new TH1D("Rhcal1", "HCAL1 R value", 50, 0., 1.);
    etot = new TH1D("etot", "Energy WCAL+ECAL+HCAL", 200, 0., 200.);

    ehcell0 = new TH1D("ehcell0", "HCAL1 cell 0 energy", 400, 0., 200.);
    ehcell1 = new TH1D("ehcell1", "HCAL1 cell 1 energy", 400, 0., 200.);
    ehcell2 = new TH1D("ehcell2", "HCAL1 cell 2 energy", 400, 0., 200.);
    ehcell3 = new TH1D("ehcell3", "HCAL1 cell 3 energy", 400, 0., 200.);
    ehcell4 = new TH1D("ehcell4", "HCAL1 cell 4 energy", 400, 0., 200.);
    ehcell5 = new TH1D("ehcell5", "HCAL1 cell 5 energy", 400, 0., 200.);
    ehcell6 = new TH1D("ehcell6", "HCAL1 cell 6 energy", 400, 0., 200.);
    ehcell7 = new TH1D("ehcell7", "HCAL1 cell 7 energy", 400, 0., 200.);
    ehcell8 = new TH1D("ehcell8", "HCAL1 cell 8 energy", 400, 0., 200.);

    cm0 = new TH1F("cm0", "cm for coord 0 pm 0.5", 10, -10., 10.);
    cm5 = new TH1F("cm5", "cm for coord 5 pm 0.5", 20, 0., 20.);
    cm10 = new TH1F("cm10", "cm for coord 10 pm 0.5", 20, 5., 25.);
    cm15 = new TH1F("cm15", "cm for coord 15 pm 0.5", 20, 10., 30.);

    ebgo0 = new TH1D("ebgo0", "BGO cristal 0 energy", 100, 0., 200.);
    ebgo3 = new TH1D("ebgo3", "BGO cristal 3 energy", 100, 0., 200.);
    ebgo4 = new TH1D("ebgo4", "BGO cristal 4 energy", 100, 0., 200.);
    ebgo7 = new TH1D("ebgo7", "BGO cristal 7 energy", 100, 0., 200.);
	
    if (e0.mc) {
      NominalBeamEnergy = 100.;
      ISetup = -1.;
      if(MCRunInfo.find(string("setup_number")) != MCRunInfo.end()) ISetup = cast<int>(MCRunInfo["setup_number"]);
      if(MCRunInfo.find(string("SetupNumber")) != MCRunInfo.end()) ISetup = cast<int>(MCRunInfo["SetupNumber"]);
      if(ISetup < 0) std::cout << "Setup number not found in the Run Info, assumed -1" << std::endl;
      if(ISetup > 27) // MC for the visible mode 2018 setup
        NominalBeamEnergy = 150.;
    }

  }  // end of initialization

  //crear histogramas aqui 

  int NCellsHCALX = 0;
  if(e0.mc) {
    NCellsHCALX = mcrunoption<int>("HCAL:NCellsHCALX", 3);
  } else {
    NCellsHCALX = HCAL_pos.Nx;
  }

#include "mkvars.inc"

//#include "mkcuts_hcontamination1.inc"
//#include "mkcuts_hcontamination_control.inc"
//#include "mkcuts2021_calib_e.inc"
//#include "mkcuts2021_invis.inc"
#include "mkcuts_calib_e.inc"
//#include "mkcuts_calib_e_sign.inc"
//#include "mkcuts_calib_pi.inc"
//#include "mkcuts_calib_pi_50.inc"
//#include "mkcuts_calib_pi_2021mu.inc"
//#include "mkcuts_tracking_2021mu.inc"
//#include "mkcuts_calib_pi_hadrinvis50.inc"
//#include "mkcuts_hadrinvis50.inc"
//#include "mkcuts_calib_mu.inc"
//#include "mkcuts_calib_mu_sign.inc"
//#include "mkcuts_invis.inc"
//#include "mkcuts_invis_nosrd.inc"
//#include "mkcuts_invis_nuage.inc"
//#include "mkcuts_invis_herm.inc"
//#include "mkcuts_invis_loose.inc"
//#include "mkcuts_invis_loose_srd.inc"
//#include "mkcuts_invis_avis.inc"
//#include "mkcuts_invis_avis_1.inc"
//#include "mkcuts_dimu_100.inc"
#include "mkcuts_dimu_100_2022.inc"
//#include "mkcuts_dimu_100_sign.inc"
//#include "mkcuts_dimu_100_singlemu.inc"
//#include "mkcuts_W_2021.inc"
//#include "mkcuts_W_100.inc"
//#include "mkcuts_W_150.inc"
//#include "mkcuts_W_150_new.inc"
//#include "mkcuts_W_150_soft.inc"
//#include "mkcuts_W_150_noW2.inc"
//#include "mkcuts_W_pitest_150.inc"
//#include "mkcuts_W_pitest_150_soft.inc"
//#include "mkcuts_W_pitest_150_noW2_soft.inc"
//#include "mkcuts_W0_100.inc"
//#include "mkcuts_W0_150.inc"
//#include "mkcuts_W0_pitest_150.inc"
//#include "mkcuts_W0_pitest_150_hard.inc"
//#include "mkcuts_W0_150_noW2.inc"
//#include "mkcuts_W0_150_soft.inc"
//#include "mkcuts_W_calib_mu.inc"
//#include "mkcuts_W_calib_e.inc"
//#include "mkcuts_W_calib_e_150.inc"
//#include "mkcuts_W_calibecal_pi.inc"
//#include "mkcuts_W_dimu.inc"
//#include "mkcuts_W_dimu_150.inc"
//#include "mkcuts_W_K0.inc"
//#include "mkcuts_W_K0_150.inc"
//#include "mkcuts_W_K0_150_noW2.inc"
//#include "mkcuts_W0_K0.inc"
//#include "mkcuts_W0_K0_150.inc"
//#include "mkcuts_W0_K0_150_hard.inc"

//  std::cout << "Event #" << ievproc;

  if(e0.mc) ifmc = 1;

//  timehisto(e0);

  int ipass = 1;

  if(ipass) {
    if(e0.mc) {
      flow->Fill(-1., e0.mc->Weight);
    } else {
      flow->Fill(-1.);
    }
  }

  int itimingok = 1;

  // apply time cuts
  RecoEvent e = e0;
  double ecalbad = 0.;
  double hcalbad = 0.;
  if (!e.mc) { // timecut is relevant only for real data events
    if(!e0.hasMasterTime()) itimingok = 0;
    timecut(e, TimeCut);
    const double ETot = e0.ecalTotalEnergy();
    const double EGood = e.ecalTotalEnergy();
    if(e.run < 2467 || (e.run >= 3855 && e.run < 4224) || (e.run >= 6500 && e.run <= 8395) ) { // Use bad energy cut only for invisible mode
      //if (EGood/ETot < 0.7) itimingok = 0; // old BadEnergy cut
      if (ETot - EGood > 30.) itimingok = 0; // new BadEnergy cut
    }
    ecalbad = ETot - EGood;

    double ETotH = e0.hcalTotalEnergy();
    double EGoodH = e.hcalTotalEnergy();
    if(e.run >= 5317 && e.run <= 5324) { // check bad HCAL energy in hadron runs of M2 run 2021
      if (ETotH - EGoodH > 5.) itimingok = 0;
    }

    ETotH = 0.; EGoodH = 0.; // Check HCAL bad energy (only in modules along the bent beam)
    for (int d = 0; d < 3; ++d) {
      for (int x = 0; x < NCellsHCALX; ++x) {
        for (int y = 0; y < 3; ++y) {
          ETotH += e0.HCAL[d][x][y].energy;
          EGoodH += e.HCAL[d][x][y].energy;
        }
      }
    }
    hcalbad  = ETotH - EGoodH;
    if(hcalbad > 50.) itimingok = 0;

  }

  if(ECALBadEnergyCut > 0.5) {
    if(!itimingok) ipass = 0;
  }
  
  int istep = 0;
  stepname[istep] = "initial";
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  // Straw cuts
  int mpStraw4X = 0; // Straw hits
  int mpStraw4Y = 0;
  int mpStraw3X = 0;
  int mpStraw3Y = 0;
  if (!e.mc) { // data
    for (auto &d_it : manager.GetEventDigits()) {
      const Chip::Digit *digit = d_it.second;
      const ChipNA64TDC::Digit *tdc = dynamic_cast<const ChipNA64TDC::Digit*>(digit);
      //if (tdc) std::cout << " tdc: " << tdc << std::endl; 
      if (tdc) {
        if (tdc->GetDetID().GetName() == "ST04X" && tdc->GetTime() < 920 && tdc->GetTime() > 840 && tdc->GetWire() != 64) {
          mpStraw4X++;
        }
        if (tdc->GetDetID().GetName() == "ST04Y" && tdc->GetTime() < 920 && tdc->GetTime() > 840 && tdc->GetWire() != 64) {
          mpStraw4Y++;
        }
        if (tdc->GetDetID().GetName() == "ST03X" && tdc->GetTime() < 920 && tdc->GetTime() > 840 && tdc->GetWire() != 64) {
          mpStraw3X++;
        }
        if (tdc->GetDetID().GetName() == "ST03Y" && tdc->GetTime() < 920 && tdc->GetTime() > 840 && tdc->GetWire() != 64) {
          mpStraw3Y++;
        }
      }
    }
  } else { // MC
    for (auto strawHit : e.mc->ST3truth) {
      if(strawHit.IProj == 0) mpStraw3X++;
      if(strawHit.IProj == 1) mpStraw3Y++;
    }
    //mpStraw3X = (e.mc->ST3Xtruth).size();
    //mpStraw3Y = (e.mc->ST3Ytruth).size();
  }
  if(StrawCut > 0.5) {if(mpStraw3X < 1 || mpStraw3X > 5) ipass = 0; if(mpStraw3Y < 1 || mpStraw3Y > 5) ipass = 0;}

  istep = 1;
  stepname[istep] = "straw3 ";  // -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  int NCellsECALX = 0;
  int NCellsECALY = 0;
  if(e.mc) {
    NCellsECALX = mcrunoption<int>("ECAL:NCellsECALX", 6);
    NCellsECALY = mcrunoption<int>("ECAL:NCellsECAL", 6);
  } else {
    NCellsECALX = ECAL_pos.Nx;
    NCellsECALY = ECAL_pos.Ny;
  }
  int imaxcellx=-1, imaxcelly=-1;
  double emaxcell=0.;
 
 


  for (int x = 0; x < NCellsECALX; ++x) {
 
   for (int y = 0; y < 6; ++y) {

      if((e.ECAL[0][x][y].energy+e.ECAL[1][x][y].energy) > emaxcell) {
        emaxcell = e.ECAL[0][x][y].energy+e.ECAL[1][x][y].energy;
        imaxcellx = x;
        imaxcelly = y;
      }

      if ((e.ECAL[0][x][y].energy + e.ECAL[1][x][y].energy) > 0.0) { // fill if energy > 0
	cellEnergyHist->Fill(x, y, e.ECAL[0][x][y].energy + e.ECAL[1][x][y].energy); //cell 6x6 matrix
	cellEnergyHistograms[x][y]->Fill(e.ECAL[0][x][y].energy + e.ECAL[1][x][y].energy);  

   }
      
    }
  }

 


  if(ECALMaxCell > -1) { // check max cell
    if(e.mc || ECALMaxCell <= 29) { // set needed cell numbers from cuts
      if(e.mc && ECALMaxCell > 29) {
        std::cout << "Automatic guess of needed ecalmaxcell does not work for MC, exiting" << std::endl;
        exit(1);
      }
      if(NCellsECALX == 6 && NCellsECALY == 6) {
        if(imaxcellx != ECALMaxCell || imaxcelly != ECALMaxCell) ipass = 0;
      } else {
	if(NCellsECALY == 5) {
          if(imaxcellx != ECALMaxCell || imaxcelly != 2) ipass = 0;
        }
	if(NCellsECALX == 5) {
          if(imaxcellx != 2 || imaxcelly != 3) ipass = 0;
        }
      }
    } else { // set needed cell numbers from run number
      if(e.run <= 4309) {
        if(imaxcellx != 3 || imaxcelly != 3) ipass = 0;
      }
      if(e.run >= 4310 && e.run <= 5261) {
        if(imaxcellx != 2 || imaxcelly != 2) ipass = 0;
      }
      if(e.run >= 5262) {
        if(imaxcellx != 2 || imaxcelly != 3) ipass = 0;
      }
    }
  }

  istep = 2;
  stepname[istep] = "maxcell"; 	// -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }


  // SRD ---------------------------------------------------------------------------------------

  // Find out configuration:
  int NSRDModules = 3;
  int ISRDDesignYear = 2016;
  if(e.mc) { // MC
    NSRDModules = mcrunoption<int>("SRD:NModules", 3);
    ISRDDesignYear = mcrunoption<int>("SRD:DesignYear", 2016);
    if(!mcrunoption<int>("SRD:NModules", 0)) { // old samples, there is no "SRD:NModules" in the run info
      if(ISetup > 27) NSRDModules = 2; // visible mode 2018 setup, 150 GeV 2 SRD modules
    }
  }
  if(!e.mc) { // data
    if(e.run >= 4224 && e.run <= 4309) { // visible mode 2018 setup, 150 GeV 2 SRD modules
      NSRDModules = 2;
    }
  }

  // Rescale cuts for new design:
  double ScaleSRDCutModule = 1.;
  double ScaleSRDCut = 1.;
  if(ISRDDesignYear == 2021) {
    ScaleSRDCutModule = 0.886;
    if(NSRDModules == 4) ScaleSRDCut = 0.662; // for 3 modules
  }

  // Sum calculation:
  double srd=0.;
  for(int i=0; i < NSRDModules; i++) { srd += (e.SRD[i].energy) / MeV; }

  // One specific period of datataking
  if(!e.mc && e.run >= 2484 && e.run < 2505) { // bad positioning of SRD, no shield
    srd = (e.SRD[0].energy + e.SRD[1].energy) / MeV;
    SRDLowerCut = 6.6;
    SRDUpperCut = 100.;
  }

  // Cutting:
#if 0
  // Standard cuts
  if(srd < SRDLowerCut * ScaleSRDCut || srd > SRDUpperCut * ScaleSRDCut) ipass = 0;
  if(NSRDModules <= 3) for(int i=0; i < NSRDModules; i++) { if(e.SRD[i].energy/MeV < SRDEachLowerCut * ScaleSRDCutModule) ipass = 0; }
  int NSRDHit = 0;
  for(int i=0; i < NSRDModules; i++) { if(e.SRD[i].energy/MeV >= SRDEachLowerCut * ScaleSRDCutModule) NSRDHit++; }
  if(NSRDModules == 4 && NSRDHit < 3) ipass = 0;
#endif

  // Toropin cuts 2022
  if(srd < SRDLowerCut * ScaleSRDCut || srd > SRDUpperCut * ScaleSRDCut) ipass = 0;
  if(SRDEachLowerCut > -0.5) {
    int isrdok[3] = {0,0,0};
    int isrdok1[3] = {0,0,0};
    double cutcorr[3] = {1.,1.,1.};
    // if(e.mc) cutcorr[2] = 0.7; // srd[2] is overcalibrated in 2022 periods 0 and 1
    for (int i=0; i < 3; i++) {
      if(e.SRD[i].energy/MeV > cutcorr[i]*SRDEachLowerCut && e.SRD[i].energy/MeV < 80.) isrdok[i] = 1;
      if(e.SRD[i].energy/MeV > cutcorr[i]*SRDEachLowerCut) isrdok1[i] = 1;
    }
    int ipass1 = 0;
    if(isrdok[0] && isrdok[1] && isrdok1[2]) ipass1 = 1;
    if(isrdok[1] && isrdok[2] && isrdok1[0]) ipass1 = 1;
    if(isrdok[0] && isrdok[2] && isrdok1[1]) ipass1 = 1;
    if(!ipass1) ipass = 0;
  }

  istep = 3;
  stepname[istep] = "SRD    "; 	// -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }


#ifndef TRACKINGTOOLS 
   //preparing MM manager
    MManager man;
  
    //initialize magnet from the geometry of the run
    man.GenfitInit(conddbpath() + "NA64Geom2018.root");
    SimpleTrack track;
    double mymom = -0.5;
    if (!e.mc) { // data

      //invisible mode 2018 and earlier
      if(e.run >= 3855 && e.run <= 4223) {
        man.AddEvent(e,NOCUT);
        track = man.SimpleTrackingMM();	
      }
      //invisible mode 2021 electron run
      if(e.run >= 5021 && e.run <= 5186) {
        man.AddEvent(e,NOCUT);
        track = man.SimpleTrackingMM();
      }
      //invisible mode 2022 electron run
      if(e.run >= 6035) {
        track = simple_tracking_4MM(e);
      }

      //visible mode
      if(e.run >= 4224 && e.run <= 4309) {
        // simple track reconstruction
        man.AddMicromega(e.MM3,NOCUT);
        man.AddMicromega(e.MM4,NOCUT);
        man.AddMicromega(e.MM5,NOCUT);
        man.AddMicromega(e.MM6,NOCUT);
        track = man.PointDeflectionTracking(); //simple_tracking_4MM(e);
      }

      if (track) {
        //cout << "mom = " << track.momentum << endl;
        mymom = track.momentum;
      }
    } 
    else {                             // MC:
      if(MCFileFormatVersion >= 8) {
        //if(ievproc == 0) {
        //  Reset_Calib();                   // This block no more needed explicitly, it is done in RunP348RecoMC
        //  Init_Calib_2022B_MM(1664928923);
        //  LoadMCGeometryFromMetaData();
        //}
        //DoDigitization(e);
        man.Reset();
        man.AddMicromega(e.MM1);    
        man.AddMicromega(e.MM2);    
        man.AddMicromega(e.MM3);
        man.AddMicromega(e.MM4);
     
        track = simple_tracking_4MM_2016(e);
        if(track) mymom = track.momentum;
      }
    }
#endif
#ifdef TRACKINGTOOLS
    // track reconstruction with Tracking Tools
    double mymom;
    double mychisq;
    if(!e.mc) {
      mymom = -0.5;
      mychisq = -0.5;
      bool ifupstream = true;
      // to be set to false for downstream momentum in the muon runs
      if(MomentumLowerCut > -10.) {
        const std::vector<Track_t> mytracks = GetReconstructedTracks(e, ifupstream);
        for (unsigned int p = 0; p < mytracks.size(); p++) {
          Track_t track = mytracks.at(p);
          mymom = track.momentum;
          mychisq = track.chisquare;        
          if(mymom > 159.9999 && mymom < 160.0001) mychisq = 500.; // temporary
          if(mymom > 99.9999 && mymom < 100.0001) mychisq = 500.; // temporary

          //std::cout << "ECAL extrapolation: " << (track.getECALHit())[0] << " " << (track.getECALHit())[1] << std::endl;

        }
        if(mytracks.size() > 1) mychisq = 500.;
        if(mytracks.size() > 1) std::cout << "More than 1 track, N = " << mytracks.size() << std::endl;
      }

      //detect bad values of chi squared and apply chi squared cut
      if ( isnan((float)mychisq) || isinf((float)mychisq) ) { ipass = 0; }
      else {
        if(MomentumLowerCut > -0.5) {
          if(mychisq < MomentumChiSqLowerCut) ipass = 0;
          if(mychisq > MomentumChiSqUpperCut) ipass = 0;
        }
      }
    }

    if(e.mc) {
      mymom = -0.5;
      mychisq = -0.5;
      const bool debug = false;
      std::vector<std::pair<TVector3, double> > hits;
      FillAbsMMHits(e.MM1, hits);
      FillAbsMMHits(e.MM2, hits);
      FillAbsMMHits(e.MM3, hits);
      FillAbsMMHits(e.MM4, hits);
      if(debug)
        for(auto hit : hits)
          cout << "Reco point in MM: (" << hit.first.x() << ", " << hit.first.y() << ", " << hit.first.z() << ")" << endl;
      if(MomentumLowerCut > -10.) {
        const std::vector<Track_t> tracks = GetReconstructedTracksFromHits(hits);
        if(!tracks.empty()) {
          double _best_chi2(9999);
          double _best_mom(9999);
          // Look for the best track with the lowest chi2 value
          for (auto elem : tracks) {
            if (elem.chisquare < _best_chi2 && elem.chisquare > -1.) {
              _best_chi2 = elem.chisquare;
              _best_mom = elem.momentum;
            }
          }
          mychisq = _best_chi2;
          mymom = _best_mom;
          if (debug) std::cout << "mymom = " << mymom <<std::endl;
          if (debug) std::cout << "mychisq = " << mychisq <<std::endl;
        }
      }
    }
#endif
  
  if(mymom < MomentumLowerCut) ipass = 0;
  if(mymom > MomentumUpperCut) ipass = 0;

  istep = 4;
  stepname[istep] = "trkmom "; 	// -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  double ecal = 0.;
  double ecalcentraltot = 0.;
  double ecaltot = 0.;
  double epileupall = 0.;
  double epileupcentral = 0.;  
  double ecalnoedge = 0.;
  double ecalcentral = 0.;
  double ecaltailcentral = 0.;
  double ecaltailnoedge = 0.;
  double prs = 0.;
  for (int d = 0; d < 2; ++d) {
    for (int x = 0; x < NCellsECALX; ++x) {
      for (int y = 0; y < 6; ++y) {
        double totECal = e.ECAL[d][x][y].energy;
        double tail = 0.;
        ecaltot += totECal;
        if(!e.mc) { // study tails from previous peaks in ECAL
          tail = e.ECAL[d][x][y].pileup();
          if(tail < totECal) e.ECAL[d][x][y].energy = e.ECAL[d][x][y].energy_nopileup();
       	  if(tail >= totECal) e.ECAL[d][x][y].energy = 0.;
          epileupall += tail;
        }
        if(abs(x-3) <= 1 && abs(y-3) <= 1 ) {
          ecalcentraltot += totECal;
          epileupcentral += tail; 
          ecaltailcentral += tail;
          ecalcentral += e.ECAL[d][x][y].energy;
        }
        ecal += e.ECAL[d][x][y].energy;
        if(d == 0) prs += e.ECAL[d][x][y].energy;
        if( x * y != 0 ) {
          ecalnoedge += e.ECAL[d][x][y].energy;
          ecaltailnoedge += tail;
        }
      }
    }
  }

  double wpileupall = 0.; 
  double wpileupall0 = 0.;
  double wpileupall1 = 0.;
  double wpileupall2 = 0.;
  double wcaltot = 0.;
  for (int d = 0; d < 3; ++d) {
    double totWCal = e.WCAL[d].energy;
    double tail = 0.;
    if(!e.mc) { // study tails from previous peaks in WCAL
      tail = e.WCAL[d].pileup();
      wpileupall += tail;
      if (d == 0) wpileupall0 += tail;
      if (d == 1) wpileupall1 += tail;
      if (d == 2) wpileupall2 += tail;
      if(tail < totWCal) e.WCAL[d].energy = e.WCAL[d].energy_nopileup();
      if(tail >= totWCal) e.WCAL[d].energy = 0.;
    }
    wcaltot += totWCal;
  }

  if(e.mc) { // smear ECAL
    for (int d = 0; d < 2; ++d) {
      for (int x = 0; x < NCellsECALX; ++x) {
        for (int y = 0; y < 6; ++y) {
          e.ECAL[d][x][y].energy *= 1. + gRandom->Gaus(0., 0.04);
        }
      }
    }
    ecal = 0.;
    prs = 0.;
    for (int d = 0; d < 2; ++d) {
      for (int x = 0; x < NCellsECALX; ++x) {
        for (int y = 0; y < 6; ++y) {
          ecal += e.ECAL[d][x][y].energy;
          if(d == 0) prs += e.ECAL[d][x][y].energy;
        }
      }
    }
  }
  if(e.mc) { // WCAL[2] inefficiency
    if(e.WCAL[2].energy < 0.0009) {
      double eff = pow(e.WCAL[2].energy/0.0009, 6.);
      if(gRandom->Rndm() > eff) e.WCAL[2].energy = 0.;
    }
  }
  if(e.mc) { // smear S4 and WCAL[2]
    double adds4 = gRandom->Gaus(0., 0.000048);
    if(e.S4.energy < 0.00024) adds4 *= e.S4.energy/0.00024;
    e.S4.energy += adds4;

    double phemip = 10.;
    double oldval = e.WCAL[2].energy;
    double newval = 0.001 * gRandom->Poisson(phemip*e.WCAL[2].energy/0.001) / phemip;
    e.WCAL[2].energy = newval;
    double addw2 = gRandom->Gaus(0., 0.00015);
    addw2 = 0.;
    if(e.WCAL[2].energy > 0.0001) addw2 = 0.0008*(e.WCAL[1].energy/100.) * fabs(gRandom->Gaus(0., 1.));
    if(addw2 > e.WCAL[2].energy) e.WCAL[2].energy += addw2;
    addw2 = gRandom->Gaus(0., 0.0001);
    if(e.WCAL[2].energy > 0.0009) addw2 *= 1.5;
    if(e.WCAL[2].energy > 0.0001) e.WCAL[2].energy += addw2;

    if(e.WCAL[2].energy < 0.) e.WCAL[2].energy = 0.;

    if(e.WCAL[2].energy < 0.00025 && gRandom->Rndm() < 0.5) e.WCAL[2].energy = 0.;
    if(e.WCAL[2].energy < 0.0001 && gRandom->Rndm() < 0.5) e.WCAL[2].energy = 0.;
  }

  if(e.mc) { // smear WCAL[0] and WCAL[1]
//    double addw0 = gRandom->Gaus(0., 0.05*e.WCAL[0].energy);
//    e.WCAL[0].energy += addw0;
//    double addw1 = gRandom->Gaus(0., 0.04*e.WCAL[1].energy);
//    e.WCAL[1].energy += addw1;
  }

  //if(ecaltailcentral > 0.007) ipass = 0; // tail cut for muon calibration
  //if(ecaltailnoedge > 0.1) ipass = 0; // tail cut for hadron calibration

  //if(sumtailcentral > 0.015) ipass = 0; // sum tail cut for subtraction  muon calibration
  //if(sumtailnoedge > 0.2) ipass = 0; // sum tail cut for subtraction pion calibration

  //if(sumtailW > 0.05) ipass = 0;

  int narrowsh=0;
  if(imaxcellx && imaxcelly) {
    if(e.ECAL[1][imaxcellx-1][imaxcelly-1].energy < 0.01*emaxcell &&
       e.ECAL[1][imaxcellx][imaxcelly-1].energy < 0.01*emaxcell &&
       e.ECAL[1][imaxcellx+1][imaxcelly-1].energy < 0.01*emaxcell &&
       e.ECAL[1][imaxcellx-1][imaxcelly].energy < 0.01*emaxcell &&
       e.ECAL[1][imaxcellx+1][imaxcelly].energy < 0.01*emaxcell &&
       e.ECAL[1][imaxcellx-1][imaxcelly+1].energy < 0.01*emaxcell &&
       e.ECAL[1][imaxcellx][imaxcelly+1].energy < 0.01*emaxcell &&
       e.ECAL[1][imaxcellx+1][imaxcelly+1].energy < 0.01*emaxcell) narrowsh=1;
  }

  if(e.WCAL[0].energy < WCAL0LowerCut) ipass = 0;
  if(e.WCAL[0].energy > WCAL0UpperCut) ipass = 0;        
  if(e.WCAL[1].energy < WCAL1LowerCut) ipass = 0;        
  if(e.WCAL[1].energy > WCAL1UpperCut) ipass = 0;        
  if((e.WCAL[0].energy + e.WCAL[1].energy) < WCALSumLowerCut) ipass = 0;        
  if((e.WCAL[0].energy + e.WCAL[1].energy) > WCALSumUpperCut) ipass = 0;

  if(!e.mc) {
    if(e.run >= 3466 && e.run <= 3560) { // use catcher as a veto for the second part of visible mode 2017
      if(e.WCAL[2].energy > WCATVetoUpperCut) ipass = 0;
    }
    if(e.run >= 4224 && e.run <= 4309) { // use catcher as a veto for the visible mode 2018
      if(e.WCAT.energy > WCATVetoUpperCut) ipass = 0;
    }
  }

  if(e.WCAT.energy < WCATLowerCut) ipass = 0;
  if(e.WCAT.energy > WCATUpperCut) ipass = 0;

  istep = 5;
  stepname[istep] = "WCAT   ";  // -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if(prs < PRSLowerCut) ipass = 0;
  if(prs > PRSUpperCut) ipass = 0;

  istep = 6;
  stepname[istep] = "PRS    "; 	// -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if(ECALNarrowShower > 0 && narrowsh) ipass = 0;
  if(ecal < ECALLowerCut) ipass = 0;
  if(ecal > ECALUpperCut) ipass = 0;

  istep = 7;
  stepname[istep] = "EECAL  ";  // -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if(e.WCAL[2].energy < WCAL2LowerCut) ipass = 0;
  if(e.WCAL[2].energy > WCAL2UpperCut) ipass = 0;
  if(e.V2.energy < V2LowerCut) ipass = 0;
  if(e.V2.energy > V2UpperCut) ipass = 0;
  //if(sumtailW2 > 0.0001) ipass = 0;

  istep = 8;
  stepname[istep] = "WCAL+V2";  // -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if(e.S4.energy < S4LowerCut) ipass = 0;
  if(e.S4.energy > S4UpperCut) ipass = 0;

  istep = 9;
  stepname[istep] = "S4     ";  // -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  double veto = 0.;
  double a1, a2, aveto;
  for (int iv = 0; iv < 3; ++iv) {
    a1 = e.VETO[2*iv].energy;
    a2 = e.VETO[2*iv+1].energy;
    aveto = sqrt(a1*a2);
    if(a1 < 0.3*a2) aveto = a2;
    if(a2 < 0.3*a1) aveto = a1;
    veto += aveto;
  }

  double hcalraw=0, hcal0=0, hcal1=0, hcal2=0, hcal3=0, hcal=0, hcalcentral=0, hcal12raw=0, hcal1raw=0, hcal12=0, hcal12central=0, hcal1central=0;
  for (int d = 0; d < 4; ++d)
    for (int x = 0; x < NCellsHCALX; ++x)
      for (int y = 0; y < 3; ++y) {
        if(d < 3) {
          hcalraw += e.HCAL[d][x][y].energy;
          if(x == 1 && y == 1) hcalcentral += e.HCAL[d][x][y].energy;
        }
        if(d == 1) {
          hcal12raw += e.HCAL[d][x][y].energy;
          hcal1raw += e.HCAL[d][x][y].energy;
          if(x == 1 && y == 1) {
            hcal12central += e.HCAL[d][x][y].energy;
            hcal1central += e.HCAL[d][x][y].energy;
          }
        }
        if(d == 2) {
          hcal12raw += e.HCAL[d][x][y].energy;
          if(!e.mc) { // DATA
            if(e.run <= 4309) {
              if((x == 1 && y == 1) || (x == 2 && y == 1)) hcal12central += e.HCAL[d][x][y].energy; // in 2018 HCAL2 was not properly centered
            } else {
              if((x == 1 && y == 1)) hcal12central += e.HCAL[d][x][y].energy;
            }
          } else { // MC: FIXME : the code is now for >= 2021, check hcal2 X coordinate?
            if((x == 1 && y == 1)) hcal12central += e.HCAL[d][x][y].energy;
          }
        }
	if(d == 0) hcal0 += e.HCAL[d][x][y].energy;
        if(d == 1) hcal1 += e.HCAL[d][x][y].energy;
        if(d == 2) hcal2 += e.HCAL[d][x][y].energy;
        if(d == 3) hcal3 += e.HCAL[d][x][y].energy;
      }

  double hcalrawforR = 0.;
  if(NCellsHCALX == 6) {
    hcalcentral = 0.;
    for (int d = 0; d < 4; ++d)
      for (int x = 0; x < NCellsHCALX; ++x)
        for (int y = 0; y < 3; ++y) {
          hcalrawforR += e.HCAL[d][x][y].energy;
          if((x == 2 && y == 1) || (x == 3 && y == 1)) hcalcentral += e.HCAL[d][x][y].energy;
        }
  }

  if(e.mc) {
    double addhcal0=0., addhcal1, addhcal2;
    if(hcal0 < 7.) {
      addhcal0 = gRandom->Gaus(0., 0.4*hcal0/4.5); // smearing for MC
    } else {
      addhcal0 = gRandom->Gaus(0., 0.505); // smearing for MC
    }
    hcal0 += addhcal0;
    if(hcal1 < 7.) {
      addhcal1 = gRandom->Gaus(0., 0.35*hcal1/4.5); // smearing for MC
    } else {
      addhcal1 = gRandom->Gaus(0., 0.505); // smearing for MC
    }
    hcal1 += addhcal1;
    if(hcal2 < 7.) {
      addhcal2 = gRandom->Gaus(0., 0.35*hcal2/4.5); // smearing for MC
    } else {
      addhcal2 = gRandom->Gaus(0., 0.505); // smearing for MC
    }
    hcal2 += addhcal2;
    hcal = hcal0 + hcal1 + hcal2;
  } else {
    hcal = hcalraw;
  }
  double Rhcalval = (hcalraw-hcalcentral)/hcalraw;
  double Rhcal12val = (hcal12raw-hcal12central)/hcal12raw;
  double Rhcal1val = (hcal1raw-hcal1central)/hcal1raw;
  if(NCellsHCALX == 6) {
    Rhcalval = (hcalrawforR-hcalcentral)/hcalrawforR;
  }

#if 0
  if(!e.mc) {
    if(ipass && veto < 0.03) {
      *foutd << (e.WCAL[0].energy+e.WCAL[1].energy) << " " << e.WCAL[2].energy << " " << ecal << " "
             << veto << " " << hcal0 << " " << imaxcellx << " " << imaxcelly << " " << e.run << std::endl;
    }
  }
#endif

  double vhcalsum = 0.;
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < 4; ++x) {
      vhcalsum += e.VHCAL0[x][y].energy;
    }
  }
  if(vhcalsum > VHCALUpperCut) ipass = 0;

  if(veto < VETOLowerCut) ipass = 0;
  if(veto > VETOUpperCut) ipass = 0; // Total energy cut

  for (int iv = 0; iv < 3; ++iv) { // Each cell cut
    a1 = e.VETO[2*iv].energy;
    a2 = e.VETO[2*iv+1].energy;
    aveto = sqrt(a1*a2);
    if(a1 < 0.3*a2) aveto = a2;
    if(a2 < 0.3*a1) aveto = a1;
    if(iv != 1) { // periphery
      if(!e.mc) {if(aveto > VETOEachUpperCut) ipass = 0;}
      if(e.mc) {if(aveto > 0.85*VETOEachUpperCut) ipass = 0;} // temporary because veto in MC is slightly under-noised (or in data overcalibrated)
    } else { // central
      if(!e.mc) {if(aveto > VETOCentralUpperCut) ipass = 0;}
      if(e.mc) {if(aveto > 0.85*VETOCentralUpperCut) ipass = 0;} // temporary because veto in MC is slightly under-noised (or in data overcalibrated)
    }
  }


  if((e.ZDCAL[0].energy + e.ZDCAL[1].energy) > ZDCALUpperCut) ipass = 0;

  istep = 10;
  stepname[istep] = "VETO+ZD";  // -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if(hcal0 < HCAL0LowerCut) ipass = 0;
  if(hcal0 > HCAL0UpperCut) ipass = 0;
  if(hcal1 < HCAL1LowerCut) ipass = 0;
  if(hcal1 > HCAL1UpperCut) ipass = 0;
  if(hcal2 < HCAL2LowerCut) ipass = 0;
  if(hcal2 > HCAL2UpperCut) ipass = 0;
  if(hcal3 < HCAL3LowerCut) ipass = 0;
  if(hcal3 > HCAL3UpperCut) ipass = 0;
  if(hcal < HCALSumLowerCut) ipass = 0;
  if(hcal > HCALSumUpperCut) ipass = 0;
  if(Rhcal12val > Rhcal12UpperCut) ipass = 0;
  if(Rhcalval < RhcalLowerCut) ipass = 0;
  if(Rhcal1to0UpperCut < 999.) {
    if(hcal1/hcal0 > Rhcal1to0UpperCut) ipass = 0;
  }
  int ifmuon = 0;
  if(hcal2 > 2.) ifmuon = 1.;
  if(e.HCAL[0][1][1].energy > 1.5 && e.HCAL[1][1][1].energy > 1.5 && hcal2 > 1.3) ifmuon = 1.;
  if(HCALMuonCut > 0.5 && ifmuon) ipass = 0;
  int ihctopo = 0;
  int nperiph = 0;
  double eperiph = 0.;
  for (int x = 0; x < NCellsHCALX; ++x) {
    for (int y = 0; y < 3; ++y) {
      if((x != 1 || y != 1) && e.HCAL[0][x][y].energy > 0.5) nperiph++;
      if(x != 1 || y != 1) eperiph += e.HCAL[0][x][y].energy;
    }
  }
  if(nperiph && eperiph > e.HCAL[0][1][1].energy) ihctopo = 1;
  if(HCALTopoCut && ihctopo) ipass = 0;

  if(e.WCAL[0].energy+e.WCAL[1].energy + ecal + hcal0 + hcal1 < ETotLowerCut) ipass = 0; // Total energy in calorimeters

  istep = 11;
  stepname[istep] = "HCAL   "; 	// -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  double mychi2 = 0.;
  //mychi2 = calcShowerChi2(e, 2, 2, 3, 3); // very old call
  //methods to get entry coordinate: SP_MM34, SP_ECAL 
  mychi2 = calcShowerChi2(e, imaxcellx, imaxcelly, 2, SP_ECAL, true, kCorrSigma);
  //mychi2 = calcShowerChi2_sqrt(e); // short default call

  if(ECALChi2UpperCut < 101. && mychi2 > ECALChi2UpperCut) ipass = 0;

  istep = 12;
  stepname[istep] = "chi2   "; 	// -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if( (e.WCAL[0].energy+e.WCAL[1].energy+ecal) < WCALECALLowerCut) ipass = 0.;
  if( (e.WCAL[0].energy+e.WCAL[1].energy+ecal) > WCALECALUpperCut) ipass = 0.;

  istep = 13;
  stepname[istep] = "W+ECAL ";  // -------------------------- step ------------------------------------
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if(e.mc && MM56Cut > 0.5) {
    if((e.mc->GEM1truth).size() > 2) ipass = 0;
    if((e.mc->GEM2truth).size() > 2) ipass = 0;
  }

  istep = 14;
  stepname[istep] = "MM56   ";  // -------------------------- step ------------------------------------
  maxstep = istep;
  if(ipass) {
    if(e.mc) {
      flow->Fill((double)istep, e.mc->Weight);
    } else {
      flow->Fill((double)istep);
    }
  }

  if(ecal < ECALLowerCut2) ipass = 0;

  if(!ipass) {
//    std::cout << std::endl;
    return 0;
  }

  // Event passed, analysing and histogramming ------------------------

  timehisto(e0);


  if(IPrint) std::cout << "event passed, id = " << e.id() << std::endl;

//  std::cout << "tail = " << ecaltailcentral << std::endl;

  eecalbad->Fill(ecalbad);

  for (int y = 0; y < NCellsECALY; ++y) {
    for (int x = 0; x < NCellsECALX; ++x) {
//      std::cout << e.ECAL[0][x][y].energy+e.ECAL[1][x][y].energy << " ";
    }
//    std::cout << std::endl;
  }

  double ecalup = 0., ecaldown = 0.;
  double ecalcmvalx = 0., ecalcmvaly = 0.;
  double eecalcut = 0.;
  for (int d = 0; d < 2; ++d) {
    for (int x = 1; x < NCellsECALX; ++x) {
      for (int y = 1; y < 6; ++y) {
        eecalcut += e.ECAL[d][x][y].energy;
        if(y < 3) {ecaldown += e.ECAL[d][x][y].energy;}
        else {ecalup += e.ECAL[d][x][y].energy;}
        ecalcmvalx += e.ECAL[d][x][y].energy * 38.2 * ( (double)x - 3. );
        ecalcmvaly += e.ECAL[d][x][y].energy * 38.2 * ( (double)y - 3. );
      }
    }
  }
  ecalcmvalx = ecalcmvalx / eecalcut;
  ecalcmvaly = ecalcmvaly / eecalcut;
  double ecalr = 0.;
  double ecalr1 = 0.;
  for (int d = 0; d < 2; ++d) {
    for (int x = 1; x < NCellsECALX; ++x) {
      for (int y = 1; y < 6; ++y) {
        double xcoorddiff = 38.2 * ( (double)x - 3. ) - ecalcmvalx;
        double ycoorddiff = 38.2 * ( (double)y - 3. ) - ecalcmvaly;
        double rad = sqrt(xcoorddiff*xcoorddiff + ycoorddiff*ycoorddiff);
        ecalr += e.ECAL[d][x][y].energy * rad;
        if(e.mc) {
          xcoorddiff = 38.2 * ( (double)x - 3. ) - e.mc->ECALEntryX;
          ycoorddiff = 38.2 * ( (double)y - 3. ) - e.mc->ECALEntryY;
          rad = sqrt(xcoorddiff*xcoorddiff + ycoorddiff*ycoorddiff);
          ecalr1 += e.ECAL[d][x][y].energy * rad;
        }
      }
    }
  }
  ecalr /= eecalcut;
  ecalr1 /= eecalcut;


  eveto->Fill(veto);
  evhcal->Fill(vhcalsum);

  if(IPrint) {
    std::cout << " Subdet: "
              << " SRD = " << srd
              << " ECAL = " << ecal
              << " HCAL0-2 = " << hcal0+hcal1+hcal2
              << " HCAL3 = " << hcal3
              << " Veto = " << veto
              << std::endl;
  }
  
  momentum->Fill(mymom);   //track.momentum

  ehplot->Fill(ecal, hcal);
  srdplot->Fill(srd);
  srdplot0->Fill(e.SRD[0].energy/MeV);
  srdplot1->Fill(e.SRD[1].energy/MeV);
  srdplot2->Fill(e.SRD[2].energy/MeV);

  ebgo0->Fill(e.BGO[0].energy/MeV);
  ebgo3->Fill(e.BGO[3].energy/MeV);
  ebgo4->Fill(e.BGO[4].energy/MeV);
  ebgo7->Fill(e.BGO[7].energy/MeV);

//  double mychi2 = calcShowerChi2(e, 2, 2, 3, 3);

  if(ecal > 2.) {
    if(e.mc) {
      ecalrad->Fill(ecalr, e.mc->Weight);
      chi2profile->Fill(mychi2, e.mc->Weight);
    } else {
      ecalrad->Fill(ecalr);
      chi2profile->Fill(mychi2);
    }
  }

  if(IPrint && Mode == 1) {
    std::cout << " WCAL = " << e.WCAL[0].energy+e.WCAL[1].energy << "  WCAL[2] = " << e.WCAL[2].energy
              << " no timing WCAL[2] = " << e0.WCAL[2].energy << std::endl;
    std::cout << " S4 = " << e.S4.energy << " no timing S4 = " << e0.S4.energy << std::endl;
    std::cout << " HCAL0 = " << hcal0 << "  ECAL chi2 = " << mychi2 << std::endl;
    std::cout << " V2 = " << e.V2.energy << " no timing V2 = " << e0.V2.energy << std::endl;
  }

  if(e.mc) {
    mcecalentries->Fill(e.mc->ECALEntryX);
    mcecaly->Fill(e.mc->ECALEntryY);

    ecalratio->Fill( e.mc->ECALEntryY, (ecalup/ecaldown - 1.) );

    ecalcm->Fill( e.mc->ECALEntryY, ecalcmvaly );

    if(ecal > 2.) {
      //if(fabs(e.mc->ECALEntryX) < 50. && fabs(e.mc->ECALEntryY) < 50.)ecalrad->Fill(ecalr);
      if(fabs(e.mc->ECALEntryX) < 50. && fabs(e.mc->ECALEntryY) < 50.)ecalrad1->Fill(ecalr1);
    }
  }

  eecal->Fill(ecal);
  etotecal->Fill(ecaltot);
  epileup->Fill(epileupall);
  epileup33->Fill(epileupcentral); //(ecaltailcentral);
  eecalnoedge->Fill(ecalnoedge);
  eecalcentral->Fill(ecalcentral);
  eecalcentraltot->Fill(ecalcentraltot);
  
  
  etotwcal->Fill(wcaltot);
  wpileup->Fill(wpileupall);
  wpileup0->Fill(wpileupall0);
  wpileup1->Fill(wpileupall1);
  wpileup2->Fill(wpileupall2);

  
  if(e.mc) {
    eecalweighted->Fill(ecal, e.mc->Weight);
  } else {
    eecalweighted->Fill(ecal);
  }
  eprs->Fill(prs);
  prstoecal->Fill(prs/ecal);

  if(e.WCAL[0].energy+e.WCAL[1].energy + ecal+hcal0+hcal1 > 50.) {
    eweplot2->Fill((e.WCAL[0].energy+e.WCAL[1].energy), ecal+hcal0+hcal1);
    eweplot3->Fill((e.WCAL[0].energy+e.WCAL[1].energy), 1.);
  }

  if(e.mc) {
    ewecal->Fill((e.WCAL[0].energy+e.WCAL[1].energy), e.mc->Weight);
    ewprs->Fill(e.WCAL[0].energy, e.mc->Weight);
    ewecalmain->Fill(e.WCAL[1].energy, e.mc->Weight);
    //ecatcher->Fill(e.WCAL[2].energy, e.mc->Weight);
    ecatcher->Fill(e.WCAT.energy, e.mc->Weight);
    ewecalecal->Fill((e.WCAL[0].energy+e.WCAL[1].energy) + ecal, e.mc->Weight);
    eweplot->Fill((e.WCAL[0].energy+e.WCAL[1].energy), ecal, e.mc->Weight);
    eweplot1->Fill((e.WCAL[0].energy+e.WCAL[1].energy), ecal+hcal0+hcal1, e.mc->Weight);
    //ewprsmain->Fill(e.WCAL[1].energy, e.WCAL[0].energy, e.mc->Weight);
    evtwc->Fill(e.VTWC.energy, e.mc->Weight);
    evtec->Fill(e.VTEC.energy, e.mc->Weight);
    es2->Fill(e.S2.energy, e.mc->Weight);

    ew2->Fill(e.WCAL[2].energy, e.mc->Weight);

    if((e.WCAL[0].energy+e.WCAL[1].energy) < 50.) ew22->Fill(e.WCAL[2].energy, e.mc->Weight);
    if(e.WCAL[0].energy > 8.) ew23->Fill(e.WCAL[2].energy, e.mc->Weight);
    ev2->Fill(e.V2.energy, e.mc->Weight);
    es4->Fill(e.S4.energy, e.mc->Weight);
    edm->Fill(e.DM[0].energy, e.mc->Weight);
    ev2w2->Fill(e.V2.energy, e.WCAL[2].energy, e.mc->Weight);
  } else {
    ewecal->Fill((e.WCAL[0].energy+e.WCAL[1].energy));
    eweplot1->Fill((e.WCAL[0].energy+e.WCAL[1].energy), ecal+hcal0+hcal1);
    ewprs->Fill(e.WCAL[0].energy);
    ewecalmain->Fill(+e.WCAL[1].energy);
    if(e.run >= 3466 && e.run <= 3560) { 
      ecatcher->Fill(e.WCAL[2].energy);
    }
    //if(e.run >= 4224 && e.run <= 4309) { 
    //  ecatcher->Fill(e.WCAT.energy);
    //}
    if(e.run >= 4224) {
      ecatcher->Fill(e.WCAT.energy);
    }
    ewecalecal->Fill((e.WCAL[0].energy+e.WCAL[1].energy) + ecal);
    eweplot->Fill((e.WCAL[0].energy+e.WCAL[1].energy), ecal);
    //ewprsmain->Fill(e.WCAL[1].energy, e.WCAL[0].energy);
    evtwc->Fill(e.VTWC.energy);
    evtec->Fill(e.VTEC.energy);
    es2->Fill(e.S2.energy);
    ew2->Fill(e.WCAL[2].energy);
    if((e.WCAL[0].energy+e.WCAL[1].energy) < 50.) ew22->Fill(e.WCAL[2].energy);
    if(e.WCAL[0].energy > 8.) ew23->Fill(e.WCAL[2].energy);
    ev2->Fill(e.V2.energy);
    es4->Fill(e.S4.energy);
    edm->Fill(e.DM[0].energy);
    ev2w2->Fill(e.V2.energy, e.WCAL[2].energy);
  }

  if(e.mc) {
    ehcal->Fill(hcal, e.mc->Weight);
  } else {
    ehcal->Fill(hcal);
  }

  ehcal0->Fill(hcal0);
  ehcal1->Fill(hcal1);
  ehcal2->Fill(hcal2);
  ehcal3->Fill(hcal3);

  ehcalcatcher->Fill(hcal + e.WCAT.energy);

  ehcal0veto->Fill(hcal0, veto*1000.);

  if(e.mc) {
    if(ecal < 80.) Rhcal->Fill(Rhcalval, e.mc->Weight);
    if(ecal < 80.) Rhcal12->Fill(Rhcal12val, e.mc->Weight);
    if(ecal < 80.) Rhcal1->Fill(Rhcal1val, e.mc->Weight);
//    if(ecal < 30.) Rhcal->Fill(Rhcalval, e.mc->Weight);
  } else {
    if(ecal < 80.) Rhcal->Fill(Rhcalval);
    if(ecal < 80.) Rhcal12->Fill(Rhcal12val);
    if(ecal < 80.) Rhcal1->Fill(Rhcal1val);
//    if(ecal < 30.) Rhcal->Fill(Rhcalval);
  }

  if(e.mc) {
    etot->Fill(e.WCAL[0].energy + e.WCAL[1].energy + ecal + hcal0 + hcal1 + hcal2, e.mc->Weight);
  } else {
    etot->Fill(e.WCAL[0].energy + e.WCAL[1].energy + ecal + hcal0 + hcal1 + hcal2);
  }

  int ix0 = 0;
  if(NCellsHCALX == 6) ix0 = 1;
  ehcell0->Fill(e.HCAL[0][ix0+0][0].energy);
  ehcell1->Fill(e.HCAL[0][ix0+1][0].energy);
  ehcell2->Fill(e.HCAL[0][ix0+2][0].energy);
  ehcell3->Fill(e.HCAL[0][ix0+0][1].energy);
  ehcell4->Fill(e.HCAL[0][ix0+1][1].energy);
  ehcell5->Fill(e.HCAL[0][ix0+2][1].energy);
  ehcell6->Fill(e.HCAL[0][ix0+0][2].energy);
  ehcell7->Fill(e.HCAL[0][ix0+1][2].energy);
  ehcell8->Fill(e.HCAL[0][ix0+2][2].energy);

  if(e.mc) {
    if(e.mc->ECALEntryY > -0.5 && e.mc->ECALEntryY < 0.5) cm0->Fill(ecalcmvaly);
    if(e.mc->ECALEntryY > 4.5 && e.mc->ECALEntryY < 5.5) cm5->Fill(ecalcmvaly);
    if(e.mc->ECALEntryY > 9.5 && e.mc->ECALEntryY < 10.5) cm10->Fill(ecalcmvaly);
    if(e.mc->ECALEntryY > 14.5 && e.mc->ECALEntryY < 15.5) cm15->Fill(ecalcmvaly);
  }

  return 1;
}


void endrun()
{
  for(int ibin=1; ibin <= 39; ibin++) {
    mcecaly->SetBinError(ibin, 0.);
  }

  ecalcm->Divide(mcecaly);
  ecalratio->Divide(mcecaly);

  eweplot2->Divide(eweplot3);

  std::cout << "Integral in ECAL          = " << eecal->Integral() << std::endl;
  std::cout << "Integral in ECAL weighted = " << eecalweighted->Integral() << std::endl;
  //std::cout << "Integral in WCAL          = " << ewecal->Integral() << std::endl;
  std::cout << "ECAL, HCAL bad energy: " << 100. * (flow->GetBinContent(0) - flow->GetBinContent(1)) / flow->GetBinContent(0) << "%" << std::endl;

  std::cout << endl;
  for(int istep=1; istep <= 14; istep++) { // cutflow printout
    if(flow->GetBinContent(istep) > 0.) {
      double drop = (flow->GetBinContent(istep) - flow->GetBinContent(istep+1)) / flow->GetBinContent(istep);
      if(istep < 10) {
        std::cout << "step  " << istep << "   " << stepname[istep] << "  " << drop*100. << "%" << std::endl;
      } else {
        std::cout << "step  " << istep << "  " << stepname[istep] << "  " << drop*100. << "%" << std::endl;
      }
    }
  }
  double drop = flow->GetBinContent(maxstep) / flow->GetBinContent(0);
  std::cout << endl;
  std::cout << "efficiency" << "  " << drop*100. << "%" << std::endl;

  std::cout << endl;

}
