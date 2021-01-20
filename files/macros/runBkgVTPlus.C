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


#include "AliAODRecoDecay.h"
#include "AliDecayer.h"
#include "AliDecayerEvtGen.h"
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;
// rapidity range
double yminSG = -2.; 
double ymaxSG = 8.;
// pT range
double ptminSG = 0.;
double ptmaxSG = 3.;
// event vertex
double vX = 0, vY = 0, vZ = 0;
//number of energy of the colliding beam
const int NEnergy = 5;
//number of particles 
const int NParticles = 5;
//array of the energy
double Elab[NEnergy] = {20,30,40,80,158};
//                                                    phi                        K               Lambda                Omega                    Csi
double Tslope[2][NParticles][NEnergy] = {{{196.8,237.4,244.6,239.8,298.7},{0,0,228.9, 223.1, 228.9},{244,249,258,265,301},{0,0,218,0,267},{221,233,222,227,277}},//matter e rapidity distribution [matter/antimatter][particle][beam energy]
                                         {{196.8,237.4,244.6,239.8,298.7},{0,0,226,217,226},{339,284,301,292,303},{0,0,218,0,259},{311,277,255,321,0}}};//antimatter
//sigma parameter of the gaussians of th                    phi                        K                       Lambda                Omega                    Csi
double sigma_rapidity[2][NParticles][NEnergy] = {{{0.425,0.538,0.696,0.658,1.451},{0,0,0.674, 0.743, 0.84},{0.51,0.66,0.91,0.87,1.0},{0,0,0.6,0,1.2},{0.45,0.56,0.76,0.71,1.18}},//matter
                                                 {{0.425,0.538,0.696,0.658,1.451},{0,0,0,0,0},{0,0,0.71,0.85,0.95},{0,0,0.6,0,1.0},{0,0,0,0,1.2}}};//antimatter
//anti omega sigma  e lambda
//mu paramter of the gaussian of the rapidity distribution [matter/antimatter][particle][beam energy]
//                                                       phi                        K                   Lambda                   Omega                    Csi
double y0_rapidity[2][NParticles][NEnergy] = {{{0.425,0.538,0.487,0.682,0.},{0,0,0.619, 0.701, 0.775},{0.49,0.59,0.65,0.94,0},{0,0,0,0,0},{0.45,0.47,0.54,0.68,0}},//matter
                                              {{0.425,0.538,0.487,0.682,0.},{0,0,0.569,0.668,0.727},{0,0,0,0,0},{0,0,0.0,0},{0,0,0,0,0}}};//antimatter
//multiplicity for event [matter/antimatter][particle][beam energy]
//                                                       phi                        K                   Lambda                   Omega                    Csi
double multiplicity[2][NParticles][NEnergy] = {{{1.89,1.84,2.55,4.04,8.46},{0,0,59.1,76.9,103.0},{27.1,36.9,43.1,50.1,44.9},{0,0,0.14,0,0.43},{1.50,2.42,2.96,3.80,4.04}},//matter
                                               {{1.89,1.84,2.55,4.04,8.46},{0,0,19.2,32.4,51.9},{0.16,0.39,0.68,1.82,3.07},{0,0,0.14,0,0.19},{0.12,0.13,0.58,0.66,0}}};//antimatter
//name of the particles
TString particle_name[NParticles] = {"phi","K0s","Lambda","Omega","Csi"};
//pdg code of the particles
double pdg_mother[NParticles] = {333,310,3122,3334,3312};
//pdg code of the daughters
//                                     K+    K-   pi+  pi-    p    pi-  Lambda  K-  Lambda pi-
double pdg_daughter[NParticles][2] = {{321,-321},{211,-211},{2212,-211},{3122,-321},{3122,-211}};


double GetTslope(int pdgParticle, double Eint, bool matter);
double GetSigmaRapidity(int pdgParticle, double Eint, bool matter);
double GetY0Rapidity(int pdgParticle, double Eint, bool matter);
double GetMultiplicity(int pdgParticle, double Eint, bool matter);
double GetY0(double Eint);
void GetPDGDaughters(int pdgParticle,int pdgDaughters[], bool matter);

TDatime dt;

void runBkgVT(Int_t nevents = 5, 
	      double Eint = 158.,//"../setups/setup_EHN1-H8_short_10pixel_1.5T_BB.txt",/
	      const char *setup = "../setups/setup-EHN1_BetheBloch.txt",
        TString suffix = "_layer5",
	      bool simulateBg=kTRUE,
        int pdg_target = 310,
        int minITSHits = 5)
{

  int refreshBg = 100;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");

  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); //read the default decay table DECAY.DEC and particle table
  printf("-- Use existing decay modes in aliroot\n");
  fDecayer->SetForceDecay(kHadronicD); 
  
  /*
  bool privTab=kFALSE;
  if (strlen(privateDecayTable)>0){
    if(gSystem->Exec(Form("ls -l %s",privateDecayTable))==0){
      fDecayer->SetDecayTablePath((char*)privateDecayTable);
      fDecayer->ReadDecayTable();
      printf("-- Use decay table from file %s\n",privateDecayTable);
      privTab=kTRUE;
    }
  }
  if(!privTab){
  }
  */
  
  //gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8210/xmldoc")); // check if pythia8 path is set correctly !!!!
  
  TH2F *hResPy = new TH2F("hResPy", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 25, -0.001, 0.001, 30, 0, 3);
  TH2F *hResPyPos = new TH2F("hResPyPos", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 25, -0.001, 0.001, 30, 0, 3);
  TH2F *hResPyNeg = new TH2F("hResPyNeg", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 25, -0.001, 0.001, 30, 0, 3);

  TH3F *h3DPiBkg = new TH3F("h3DPiBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH3F *h3DKBkg = new TH3F("h3DKBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH3F *h3DPBkg = new TH3F("h3DPBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH1F *hNevents = new TH1F("hNevents", "", 1, 0, 1);
  TH1F* hGenStat = new TH1F("hGenStat","",6,0.5,6.5);
  hGenStat->GetXaxis()->SetBinLabel(1,"gen");
  hGenStat->GetXaxis()->SetBinLabel(2,"bad SolveSingleTrack");
  hGenStat->GetXaxis()->SetBinLabel(3,"bad GetWinnerMCTrack");
  hGenStat->GetXaxis()->SetBinLabel(4,"bad chi2");
  hGenStat->GetXaxis()->SetBinLabel(5,"in tree");
  hGenStat->GetXaxis()->SetBinLabel(6,"with fake clusters");
  
  
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

  //define the pT and rapidity probability function
  TF1 *fpt = new TF1("fpt","x*exp(-TMath::Sqrt(x**2+[0]**2)/[1])",ptminSG,ptmaxSG);
  TF1 *fy = new TF1("fy","exp(-0.5*((x-[0])/[1])**2)+exp(-0.5*((x+[0])/[1])**2) ",-10,10);

  TFile *f = new TFile(Form("NewtreeBkgEvents%s.root",suffix.Data()), "RECREATE");
  TTree *tree = new TTree("tree", "tree Bkg");
  TClonesArray *arrtr = new TClonesArray("KMCProbeFwd");
  TClonesArray &aarrtr = *arrtr;
  tree->Branch("tracks", &arrtr);
  

  Double_t charge = 0;
  Double_t yrap = 0;
  Double_t pt = 0;
  Double_t phi = 0;
  Double_t pxGen = 0;
  Double_t pyGen = 0;
  Double_t mt = 0;
  Double_t pzGen = 0;
  Double_t en = 0;
  Double_t mass = 0;
  Double_t multiplicity = 0;
  Double_t ntrack = 0;
  
  for (Int_t iev = 0; iev < nevents; iev++){
    aarrtr.Clear();
    printf(" ***************  ev = %d \n", iev);
    hNevents->Fill(0.5);

    if (simulateBg && (iev % refreshBg) == 0)
      det->GenBgEvent(vX, vY, vZ);

    Int_t icount = 0;
    printf("ok");
    double ntrPi = gRandom->Poisson(det->GetNChPi());
    printf("fNChPi=%f ntrPi=%f\n", det->GetNChPi(), ntrPi);

    printf(" ***************  pdg = pion \n");

    for (int itr = 0; itr < ntrPi; itr++){	
      yrap = fdNdYPi->GetRandom();
      pt = fdNdPtPi->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.52 ? 1 : -1;
      mass = KMCDetectorFwd::kMassPi;
      h3DPiBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};


      TLorentzVector *ppi = new TLorentzVector(0., 0., 0., 0.);
      ppi->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      double py= pxyz[1];
      if (!det->SolveSingleTrack(ppi->Pt(), ppi->Rapidity(), ppi->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99))
        continue;
        
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw){
        continue;
      }
      if (trw->GetNormChi2(kTRUE) > ChiTot){
        continue;
      }
      new (aarrtr[icount]) KMCProbeFwd(*trw);
      icount++;
    }

    // kaons
    double ntrK = gRandom->Poisson(det->GetNChK());
    //printf("fNChK=%f ntrK=%f\n", det->GetNChK(), ntrK);
    
    printf(" ***************  pdg = kaon \n");
    for (int itr = 0; itr < ntrK; itr++){
      yrap = fdNdYK->GetRandom();
      pt = fdNdPtK->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.3 ? 1 : -1;
      mass = KMCDetectorFwd::kMassK;
      h3DKBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      
      TLorentzVector *pk = new TLorentzVector(0., 0., 0., 0.);
      pk->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      double py= pxyz[1];
      if (!det->SolveSingleTrack(pk->Pt(), pk->Rapidity(), pk->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99))
        continue;
      
      KMCProbeFwd *trw2 = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw2){
        continue;
      }
      if (trw2->GetNormChi2(kTRUE) > ChiTot){
        continue;
      }
      new (aarrtr[icount]) KMCProbeFwd(*trw2);
      icount++;
    }
    
    // protons
    double ntrP = gRandom->Poisson(det->GetNChP());

    printf(" ***************  pdg = proton \n");
    //printf("fNChP=%f ntrP=%f\n", det->GetNChP(), ntrP);
    for (int itr = 0; itr < ntrP; itr++){
      yrap = fdNdYP->GetRandom();
      pt = fdNdPtP->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = 1;
      mass = KMCDetectorFwd::kMassP;
      h3DPBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      
      TLorentzVector *pp = new TLorentzVector(0., 0., 0., 0.);
      pp->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      double py= pxyz[1];

      if (!det->SolveSingleTrack(pp->Pt(), pp->Rapidity(), pp->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99))
        continue;
      
      KMCProbeFwd *trw3 = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw3){
        continue;
      }
      if (trw3->GetNormChi2(kTRUE) > ChiTot){
        continue;
      }
      new (aarrtr[icount]) KMCProbeFwd(*trw3);
      icount++;
    }
    //if(false)
    for (int im = 0; im < NParticles; im++){
      int pdg_mom = pdg_mother[im];
      if(pdg_target==pdg_mom)
        continue;
      printf(" ***************  pdg = %d \n", pdg_mom);


      mass = TDatabasePDG::Instance()->GetParticle(pdg_mom)->Mass();
      multiplicity = GetMultiplicity(pdg_mom,Eint,true)+GetMultiplicity(pdg_mom,Eint,false);
      ntrack = gRandom->Poisson(multiplicity);
      //std::cout<<"ntrack: "<<ntrack<<std::endl;


      for (int itr = 0; itr < ntrack; itr++){
        //printf(" ***************  trk = %d \n", itr);
        charge = gRandom->Rndm() > GetMultiplicity(pdg_mom,Eint,false)/multiplicity ? 1 : -1;
        fpt->SetParameter(0,mass);
        fpt->SetParameter(1,GetTslope(pdg_mom,Eint,charge==1)/1000);
        fy->SetParameter(0,GetY0Rapidity(pdg_mom,Eint,charge==1));
        fy->SetParameter(1,GetSigmaRapidity(pdg_mom,Eint,charge==1));
        yrap = fpt->GetRandom();
        pt = fy->GetRandom();
        phi = gRandom->Rndm() * TMath::Pi() * 2;
        pxGen = pt * TMath::Cos(phi);
        pyGen = pt * TMath::Sin(phi);
        mt = TMath::Sqrt(pt * pt + mass * mass);
        pzGen = mt * TMath::SinH(yrap);
        en = mt * TMath::CosH(yrap);


        TClonesArray *particles = new TClonesArray("TParticle", 1000);
        TLorentzVector *mom = new TLorentzVector();
        mom->SetPxPyPzE(pxGen, pyGen, pzGen, en);
        Int_t np;
        do{
          fDecayer->Decay(pdg_mom, mom);
          np = fDecayer->ImportParticles(particles);
        } while (np < 0);

        h3DPiBkg->Fill(pt, yrap, phi);
        hGenStat->Fill(1);

        // loop on decay products
        //printf(" ***************  nptot = %d \n", np);
        for (int i = 0; i < np; i++) {

          TParticle *iparticle = (TParticle *)particles->At(i);
          Int_t kf = iparticle->GetPdgCode();
          if(kf == 22) continue;
          //printf(" ***************  np = %d pdg = %d\n", i, kf);
          vX = iparticle->Vx();
          vY = iparticle->Vy();
          vZ = iparticle->Vz();
          Int_t crg= (iparticle->GetPdgCode()>0) ? 1 : -1;
          TLorentzVector partVector;
          partVector.SetXYZM(iparticle->Px(), iparticle->Py(), iparticle->Pz(), iparticle->GetMass());
          if (!det->SolveSingleTrack(partVector.Pt(), partVector.Rapidity(), partVector.Phi(), iparticle->GetMass(), crg, vX, vY, vZ, 0, 1, 99)){  
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
          new (aarrtr[icount]) KMCProbeFwd(*trw);
          hGenStat->Fill(5);
          if(trw->GetNFakeITSHits()>0) hGenStat->Fill(6);
          icount++;
        }
      }
    }
    
    
  
    printf("Background in array = %d track p %d pi %d k %d \n",icount,ntrP,ntrPi,ntrK);
    tree->Fill();
  }
  f->cd();
  tree->Write();
  f->Close();
  
  TFile *outfile = new TFile(Form("Newbkgdistributions%s.root",suffix.Data()), "recreate");
  outfile->cd();
  hResPy->Write();
  hResPyNeg->Write();
  hResPyPos->Write();
  hNevents->Write();
  hGenStat->Write();
  h3DPiBkg->Write();
  h3DKBkg->Write();
  h3DPBkg->Write();
  outfile->Close();
}

double GetTslope(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  return Tslope[matter ? 0 : 1][index_pdg][index_E];
}

double GetSigmaRapidity(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  return sigma_rapidity[matter ? 0 : 1][index_pdg][index_E];
}

double GetY0Rapidity(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  return y0_rapidity[matter ? 0 : 1][index_pdg][index_E];
}

double GetMultiplicity(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  return multiplicity[matter ? 0 : 1][index_pdg][index_E];
}

double GetY0(double Eint){
  double y0BG = 2.22;
  if (Eint == 20)
    y0BG = 1.9;   // gaussian y mean - 20 GeV
  else if (Eint == 40)
    y0BG = 2.22; // gaussian y mean - 40 GeV
  else if (Eint == 60)
    y0BG = 2.42;  // gaussian y mean - 60 GeV
  else if (Eint == 80)
    y0BG = 2.57;  // gaussian y mean - 80 GeV
  else if (Eint == 158)
    y0BG = 2.9;   // gaussian y mean - 160 GeV
  else if (Eint == 400)
    y0BG = 3.37;  // gaussian y mean - 400 GeV
  return y0BG;
}


void GetPDGDaughters(int pdgParticle, int pdgDaughters[], bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  int sign = matter ? 1 : -1;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      for(int i = 0; i < 2; i++)
        pdgDaughters[i] = pdg_daughter[index_pdg][i]*sign;
    }
    counter++;
  }
}