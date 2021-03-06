#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
#include "../KMCDetectorFwd.h"
#include "../KMCProbeFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "MeasurementUtils.h"
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;
double zTOF = 1.0;
double vX = 0, vY = 0, vZ = 0; // event vertex

TDatime dt;

void runBkgVT(Int_t nevents = 5, 
	      double Eint = 40.,//"../setups/setup_EHN1-H8_short_10pixel_1.5T_BB.txt"
	      const char *setup = "../setups/setup-EHN1_BetheBloch.txt",
        TString suffix = "_layer5_E40_longB",
	      bool simulateBg=kFALSE,
        int minITSHits = 5)
{

  int refreshBg = 100;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);
  //gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8210/xmldoc")); // check if pythia8 path is set correctly !!!!
  
  TH1D *hEff[3];
  for(int i=0; i<3; i++){
    hEff[i] = new TH1D(Form("hEff_%s",all_particle_name[i].Data()),Form("%s;efficiency;counts",all_particle_name[i].Data()),20,0,1);
  }

  TH2F *hResPy = new TH2F("hResPy", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 25, -0.001, 0.001, 30, 0, 3);
  TH2F *hResPyPos = new TH2F("hResPyPos", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 25, -0.001, 0.001, 30, 0, 3);
  TH2F *hResPyNeg = new TH2F("hResPyNeg", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 25, -0.001, 0.001, 30, 0, 3);

  TH3F *h3DPiBkg = new TH3F("h3DPiBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH3F *h3DKBkg = new TH3F("h3DKBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH3F *h3DPBkg = new TH3F("h3DPBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH1F *hNevents = new TH1F("hNevents", "", 1, 0, 1);
  TH1F* hGenStat = new TH1F("hGenStat","",18,0.5,18.5);
  TH1F *hnSigmaTOF = new TH1F("hnSigmaTOF", ";n#sigma_{#beta};counts", 100, -5., 5.);
  TH2F *hTOF = new TH2F("hTOF", ";#it{p}(GeV/#it{c});TOF #beta;counts", 200, 0, 8, 300, 0.8, 1.15);

  
  hGenStat->GetXaxis()->SetBinLabel(1,"#pi to gen");
  hGenStat->GetXaxis()->SetBinLabel(2,"#pi bad SolveSingleTrack");
  hGenStat->GetXaxis()->SetBinLabel(3,"#pi bad GetWinnerMCTrack");
  hGenStat->GetXaxis()->SetBinLabel(4,"#pi bad chi2");
  hGenStat->GetXaxis()->SetBinLabel(5,"#pi in tree");
  hGenStat->GetXaxis()->SetBinLabel(6,"#pi with fake clusters");
  hGenStat->GetXaxis()->SetBinLabel(7,"K to gen");
  hGenStat->GetXaxis()->SetBinLabel(8,"K bad SolveSingleTrack");
  hGenStat->GetXaxis()->SetBinLabel(9,"K bad GetWinnerMCTrack");
  hGenStat->GetXaxis()->SetBinLabel(10,"K bad chi2");
  hGenStat->GetXaxis()->SetBinLabel(11,"K in tree");
  hGenStat->GetXaxis()->SetBinLabel(12,"K with fake clusters");
  hGenStat->GetXaxis()->SetBinLabel(13,"p to gen");
  hGenStat->GetXaxis()->SetBinLabel(14,"p bad SolveSingleTrack");
  hGenStat->GetXaxis()->SetBinLabel(15,"p bad GetWinnerMCTrack");
  hGenStat->GetXaxis()->SetBinLabel(16,"p bad chi2");
  hGenStat->GetXaxis()->SetBinLabel(17,"p in tree");
  hGenStat->GetXaxis()->SetBinLabel(18,"p with fake clusters");
  
  
  //Magnetic field and detector parameters
  MagField *mag = new MagField(1);
  int BNreg = mag->GetNReg();
  const double *BzMin = mag->GetZMin();
  const double *BzMax = mag->GetZMax();
  const double *BVal;
  printf("*************************************\n");
  printf("number of magnetic field regions = %d\n", BNreg);
  for (int i = 0; i < BNreg; i++){
    BVal = mag->GetBVals(i);
    printf("*** Field region %d ***\n", i);
    if (i == 0){
      printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    }else if (i == 1){
      printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    }
  }

  KMCDetectorFwd *det = new KMCDetectorFwd();
  printf("Setup file = %s\n",setup);
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT
  det->SetMinITSHits(minITSHits); //NA60+
  // we don't need MS part here, even if it is in the setup
  //det->SetMinITSHits(det->GetNumberOfActiveLayersITS()-1); //NA60
  det->SetMinMSHits(0); //NA60+
  //det->SetMinMSHits(det->GetNumberOfActiveLayersMS()-1); //NA60
  det->SetMinTRHits(0);
  //
  // max number of seeds on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000);
  //
  // set chi2 cuts
  det->SetMaxChi2Cl(10.);  // max track to cluster chi2
  det->SetMaxChi2NDF(3.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(-20);  // fiducial cut on chi2 of convergence to vtx
  
  // IMPORTANT FOR NON-UNIFORM FIELDS
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); //NA60+
  //det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  //
  det->BookControlHistos();

  // Get Pi, K, P spectral shapes
  TF1* fdNdYPi=det->GetdNdYPi();
  TF1* fdNdYK=det->GetdNdYK();
  TF1* fdNdYP=det->GetdNdYP();
  TF1* fdNdPtPi=det->GetdNdPtPi();
  TF1* fdNdPtK=det->GetdNdPtK();
  TF1* fdNdPtP=det->GetdNdPtP();





  TFile *f = new TFile(Form("treeBkgEvents%s.root",suffix.Data()), "RECREATE");
  TTree *tree = new TTree("tree", "tree Bkg");
  TClonesArray *arrtr = new TClonesArray("KMCProbeFwd");
  TClonesArray &aarrtr = *arrtr;
  tree->Branch("tracks", &arrtr);
  int old_count=0;
  for (Int_t iev = 0; iev < nevents; iev++){
    aarrtr.Clear();
    printf(" ***************  ev = %d \n", iev);
    double T0 = smearT(0);
    double pxyz[3];
    hNevents->Fill(0.5);
    if (simulateBg && (iev % refreshBg) == 0)
      det->GenBgEvent(vX, vY, vZ);
    
    //gRandom->SetSeed(newseed++);//seed);
    double ntrPi = gRandom->Poisson(det->GetNChPi());
    //printf("fNChPi=%f ntrPi=%f\n", det->GetNChPi(), ntrPi);

    double yrap, pt, phi;
    int charge;
    double mass;
    double L;
    Int_t icount = 0;
    for (int itr = 0; itr < ntrPi; itr++){
      //gRandom->SetSeed(newseed++);//seed);
      yrap = fdNdYPi->GetRandom();
      pt = fdNdPtPi->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.52 ? 1 : -1;
      mass = KMCDetectorFwd::kMassPi;
      h3DPiBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      hGenStat->Fill(1);

      TLorentzVector *ppi = new TLorentzVector(0., 0., 0., 0.);
      ppi->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      if (!det->SolveSingleTrack(ppi->Pt(), ppi->Rapidity(), ppi->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99)){
        hGenStat->Fill(2);
        continue;
      }
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw){
        hGenStat->Fill(3);
        continue;
      }
      if (trw->GetNormChi2(kTRUE) > ChiTot){
        hGenStat->Fill(4);
        continue;
      }
      trw->GetPXYZ(pxyz);
      double tHit= smearT(hasTOF(pxyz, mass, L, zTOF));
      trw->SetTOF(tHit-T0);
      trw->SetL(L);
      double beta = betaTOF(L, tHit-T0);
      hTOF->Fill(TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]),beta);
      hnSigmaTOF->Fill((beta-betaGen(pxyz, mass))/sigmaBeta(zTOF));
      //printf("charge = %d, %f \n", charge, trw->GetCharge());
      new (aarrtr[icount]) KMCProbeFwd(*trw);
      hGenStat->Fill(5);
      if(trw->GetNFakeITSHits()>0) hGenStat->Fill(6);
      icount++;
    }
    hEff[0]->Fill(icount/ntrPi);
    // kaons
    //gRandom->SetSeed(newseed++);
    double ntrK = gRandom->Poisson(det->GetNChK());
    //printf("fNChK=%f ntrK=%f\n", det->GetNChK(), ntrK);
    old_count = icount;
    for (int itr = 0; itr < ntrK; itr++){
      //gRandom->SetSeed(newseed++);
      yrap = fdNdYK->GetRandom();
      pt = fdNdPtK->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.3 ? 1 : -1;
      mass = KMCDetectorFwd::kMassK;
      h3DKBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      hGenStat->Fill(7);
      
      TLorentzVector *pk = new TLorentzVector(0., 0., 0., 0.);
      pk->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      if (!det->SolveSingleTrack(pk->Pt(), pk->Rapidity(), pk->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99)){
        hGenStat->Fill(8);
        continue;
      }
      KMCProbeFwd *trw2 = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw2){
        hGenStat->Fill(9);
        continue;
      }
      if (trw2->GetNormChi2(kTRUE) > ChiTot){
        hGenStat->Fill(10);
        continue;
      }

      double tHit= smearT(hasTOF(pxyz, mass, L, zTOF));
      trw2->SetTOF(tHit-T0);
      trw2->SetL(L);
      trw2->GetPXYZ(pxyz);
      double beta = betaTOF(L, tHit-T0);
      hTOF->Fill(TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]),beta);
      hnSigmaTOF->Fill((beta-betaGen(pxyz, mass))/sigmaBeta(zTOF));
      //hResPy->Fill(pxyz[1]-py,pk->Pt());
      //printf("charge = %d, %f \n", charge, trw2->GetCharge());
      new (aarrtr[icount]) KMCProbeFwd(*trw2);
      hGenStat->Fill(11);
      if(trw2->GetNFakeITSHits()>0) hGenStat->Fill(12);
      icount++;
    }
    
    hEff[1]->Fill((icount-old_count)/ntrK);
    // protons
    //gRandom->SetSeed(newseed++);
    double ntrP = gRandom->Poisson(det->GetNChP());
    old_count = icount;
    //printf("fNChP=%f ntrP=%f\n", det->GetNChP(), ntrP);
    for (int itr = 0; itr < ntrP; itr++){
      //gRandom->SetSeed(newseed++);
      yrap = fdNdYP->GetRandom();
      pt = fdNdPtP->GetRandom();
      //std::cout<<"pt: "<<pt<<" y: "<<yrap<<std::endl;
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = 1;
      mass = KMCDetectorFwd::kMassP;
      h3DPBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      hGenStat->Fill(13);
      
      TLorentzVector *pp = new TLorentzVector(0., 0., 0., 0.);
      pp->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);

      if (!det->SolveSingleTrack(pp->Pt(), pp->Rapidity(), pp->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99)){
        hGenStat->Fill(14);
        continue;
      }
      KMCProbeFwd *trw3 = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw3){
        hGenStat->Fill(15);
        continue;
      }
      if (trw3->GetNormChi2(kTRUE) > ChiTot){
        hGenStat->Fill(16);
        continue;
      }

      double tHit= smearT(hasTOF(pxyz, mass, L, zTOF));
      trw3->SetTOF(tHit-T0);
      trw3->SetL(L);
      trw3->GetPXYZ(pxyz);
      double beta = betaTOF(L, tHit-T0);
      hTOF->Fill(TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]),beta);
      hnSigmaTOF->Fill((beta-betaGen(pxyz, mass))/sigmaBeta(zTOF));
      //hResPy->Fill(pxyz[1]-py,pp->Pt());
      //printf("charge = %d, %f \n", charge, trw3->GetCharge());
      new (aarrtr[icount]) KMCProbeFwd(*trw3);
      hGenStat->Fill(17);
      if(trw3->GetNFakeITSHits()>0) hGenStat->Fill(18);
      icount++;
    }
    hEff[2]->Fill((icount-old_count)/ntrP);
    printf("Pions+Kaons+Protons in array = %d out of %.0f \n",icount,ntrPi+ntrK+ntrP);
    tree->Fill();
  }
  f->cd();
  tree->Write();
  f->Close();
  
  TFile *outfile = new TFile(Form("bkgdistributions%s.root",suffix.Data()), "recreate");
  for(int i=0; i<3; i++){
    hEff[i]->Write();
  }
  outfile->cd();
  hResPy->Write();
  hResPyNeg->Write();
  hResPyPos->Write();
  hNevents->Write();
  hGenStat->Write();
  h3DPiBkg->Write();
  h3DKBkg->Write();
  h3DPBkg->Write();
  hTOF->Write();
  hnSigmaTOF->Write();
  outfile->Close();
}
