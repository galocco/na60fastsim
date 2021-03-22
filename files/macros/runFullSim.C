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
double ptmaxSG = 4.;
// event vertex
double vX = 0, vY = 0, vZ = 0;
//number of energy of the colliding beam
const int NEnergy = 5;
//number of particles 
const int NParticles = 5;
//array of the energy in GeV
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
double multiplicity[2][NParticles][NEnergy] = {{{1.89,1.84,2.55,4.04,8.46},{0,0,(59.1+19.2)/2.,76.9,77.3},{27.1,36.9,43.1,50.1,44.9},{0,0,0.14,0,0.43},{1.50,2.42,2.96,3.80,4.04}},//matter
                                               {{0,0,0,0,0},{0,0,0,0,0},{0.16,0.39,0.68,1.82,3.07},{0,0,0.14,0,0.19},{0.12,0.13,0.58,0.66,0}}};//antimatter
//name of the particles
TString particle_name[NParticles] = {"phi","K0s","Lambda","Omega","Csi"};
TString all_particle_name[NParticles+3] = {"pion","kaon","proton","phi","K0s","Lambda","Omega","Csi"};

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
int GetArrayPosition(int pdg);

TDatime dt;

void runFullSim(Int_t nevents = 5, 
	      double Eint = 158.,//"../setups/setup_EHN1-H8_short_10pixel_1.5T_BB.txt",/
	      const char *setup = "../setups/setup-EHN1_BetheBloch.txt",
        TString suffix = "_layer5_new",
	      bool simulateBg=kFALSE,
        int minITSHits = 5,
        bool only_prompt = false)
{

  double y0 = GetY0(Eint);
  int refreshBg = 100;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");

  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); //read the default decay table DECAY.DEC and particle table
  printf("-- Use existing decay modes in aliroot\n");
  //fDecayer->SetDecayTablePath("../decaytables/USERTABPHI.DEC");
  //fDecayer->ReadDecayTable();
  fDecayer->SetForceDecay(kHadronicD); 
  TH1F *hNevents = new TH1F("hNevents", "", 1, 0, 1);
  TH1D *hPDGSaved = new TH1D("hPDGSaved",";PDG code; entries",20000,-10000.5,10000.5);
  TH1D *hPDGSurv1 = new TH1D("hPDGSurv1",";PDG code; entries",20000,-10000.5,10000.5);
  TH1D *hPDGSurv2 = new TH1D("hPDGSurv2",";PDG code; entries",20000,-10000.5,10000.5);
  TH1D *hPDGGen = new TH1D("hPDGGen",";PDG code; entries",20000,-10000.5,10000.5);
  TH1F *hNPartRec = new TH1F("hNPartRec","",1,0,1);
  TH1F *hNPartGen = new TH1F("hNPartGen","",1,0,1);

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


  TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class()) {
    MagField* mag = (MagField*) fld;
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
  }

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

  TF1 *fm = new TF1("fm","1/((x*x-[0]*[0])*(x*x-[0]*[0])*[2]+([0]*[1])*([0]*[1]))");
  fm->SetParameter(2,TMath::Power(TMath::C(),4));

  TFile *f = new TFile(Form("treeBkgEvents%s.root",suffix.Data()), "RECREATE");
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
  Double_t massbw = 0;
  Double_t lifetime = 0;
  Double_t multiplicity = 0;
  Int_t ntrack = 0;
  TLorentzVector *mom = new TLorentzVector();
  Int_t index_list[NParticles+4];
  Int_t particle_index = 0;
  Double_t eff = 0;
  for (Int_t iev = 0; iev < nevents; iev++){
    // montecarlo ID
    for(int i=0; i<NParticles+3; i++)
      index_list[i] = -1;

    aarrtr.Clear();
    printf(" ***************  ev = %d \n", iev);
    hNevents->Fill(0.5);
    if (simulateBg && (iev % refreshBg) == 0){
      det->GenBgEvent(0, 0, 0);
    }

    Int_t icount = 0;
    //pions
    double ntr = gRandom->Poisson(det->GetNChPi());
    for (int itr = 0; itr < ntr; itr++){	
      yrap = fdNdYPi->GetRandom();
      pt = fdNdPtPi->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.52 ? 1 : -1;
      mass = KMCDetectorFwd::kMassPi;
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};


      TLorentzVector *ppi = new TLorentzVector(0., 0., 0., 0.);
      ppi->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      double py= pxyz[1];
      index_list[0]++;

      if (!det->SolveSingleTrack(ppi->Pt(), ppi->Rapidity(), ppi->Phi(), mass, charge, 0, 0, 0, 0, 1, 99)){
        continue;
      }
        
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw){
        continue;
      }
      if (trw->GetNormChi2(kTRUE) > ChiTot){
        continue;
      }

      trw->SetIndex(index_list[0]);
      trw->SetPdg(211*charge);
      trw->SetPdgMother(0);
      new (aarrtr[icount]) KMCProbeFwd(*trw);
      icount++;
    }

    // kaons
    ntr = gRandom->Poisson(det->GetNChK());
    for (int itr = 0; itr < ntr; itr++){
      yrap = fdNdYK->GetRandom();
      pt = fdNdPtK->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.3 ? 1 : -1;
      mass = KMCDetectorFwd::kMassK;
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      
      TLorentzVector *pk = new TLorentzVector(0., 0., 0., 0.);
      pk->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      double py= pxyz[1];

      index_list[1]++;
      if (!det->SolveSingleTrack(pk->Pt(), pk->Rapidity(), pk->Phi(), mass, charge, 0, 0, 0, 0, 1, 99))
        continue;
      
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw){
        continue;
      }
      if (trw->GetNormChi2(kTRUE) > ChiTot){
        continue;
      }
      trw->SetIndex(index_list[1]);
      trw->SetPdg(321*charge);
      trw->SetPdgMother(0);
      new (aarrtr[icount]) KMCProbeFwd(*trw);
      icount++;
    }
    // protons
    ntr = gRandom->Poisson(det->GetNChP());
    for (int itr = 0; itr < ntr; itr++){
      yrap = fdNdYP->GetRandom();
      pt = fdNdPtP->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = 1;
      mass = KMCDetectorFwd::kMassP;
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      
      TLorentzVector *pp = new TLorentzVector(0., 0., 0., 0.);
      pp->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      double py= pxyz[1];

      index_list[2]++;
      if (!det->SolveSingleTrack(pp->Pt(), pp->Rapidity(), pp->Phi(), mass, charge, 0, 0, 0, 0, 1, 99))
        continue;
      
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw){
        continue;
      }
      if (trw->GetNormChi2(kTRUE) > ChiTot){
        continue;
      }
      trw->SetIndex(index_list[2]);
      trw->SetPdg(2212);
      trw->SetPdgMother(0);
      new (aarrtr[icount]) KMCProbeFwd(*trw);
      icount++;
    }
    if(!only_prompt)
    for (int im = 0; im < NParticles; im++){
      int pdg_mom = pdg_mother[im];
      mass = TDatabasePDG::Instance()->GetParticle(pdg_mom)->Mass();
      lifetime = TDatabasePDG::Instance()->GetParticle(pdg_mom)->Lifetime();
      if(lifetime==0)
        lifetime = TMath::Power(8.59,-11);
      fm->SetRange(mass*0.9,mass*1.1);
      fm->SetParameter(0,mass);
      fm->SetParameter(1,TMath::Hbar()*6.242*TMath::Power(10,9)/(lifetime));

      multiplicity = GetMultiplicity(pdg_mom,Eint,true)+GetMultiplicity(pdg_mom,Eint,false);
      ntrack = gRandom->Poisson(multiplicity);
      //if(pdg_mom==310)
      //  std::cout<<multiplicity<<" "<<ntrack<<std::endl;
      for (int itr = 0; itr < ntrack; itr++){
        charge = gRandom->Rndm() > GetMultiplicity(pdg_mom,Eint,false)/multiplicity ? 1 : -1;
        fpt->SetParameter(0,mass);
        fpt->SetParameter(1,GetTslope(pdg_mom,Eint,charge==1)/1000);
        fy->SetParameter(0,GetY0Rapidity(pdg_mom,Eint,charge==1));
        fy->SetParameter(1,GetSigmaRapidity(pdg_mom,Eint,charge==1));
        yrap = fy->GetRandom()+y0;
        pt = fpt->GetRandom();
        massbw = fm->GetRandom();
        phi = gRandom->Rndm() * TMath::Pi() * 2;
        pxGen = pt * TMath::Cos(phi);
        pyGen = pt * TMath::Sin(phi);
        mt = TMath::Sqrt(pt * pt + massbw * massbw);
        pzGen = mt * TMath::SinH(yrap);
        en = mt * TMath::CosH(yrap);


        TClonesArray *particles = new TClonesArray("TParticle", 1000);
        mom->SetPxPyPzE(pxGen, pyGen, pzGen, en);
        Int_t np;
        do{
          fDecayer->Decay(pdg_mom, mom);
          np = fDecayer->ImportParticles(particles);
        } while (np < 0);
        // loop on decay products
        for (int i = 0; i < np; i++) {

          TParticle *iparticle = (TParticle *)particles->At(i);
          Int_t kf = iparticle->GetPdgCode();
          index_list[GetArrayPosition(kf)]++;
          hPDGGen->Fill(kf);
          if(kf==89 || TDatabasePDG::Instance()->GetParticle(kf)->Charge()==0) continue;
          vX = iparticle->Vx();
          vY = iparticle->Vy();
          vZ = iparticle->Vz();
          Int_t crg= (iparticle->GetPdgCode()>0) ? 1 : -1;
          TLorentzVector *partVector = new TLorentzVector(0., 0., 0., 0.);
          partVector->SetXYZM(iparticle->Px(), iparticle->Py(), iparticle->Pz(), iparticle->GetMass());
          //if(pdg_mom==333)
          //std::cout<<"vX: "<<vX<<" vY: "<<vY<<" vZ: "<<vZ<<std::endl;
          if (!det->SolveSingleTrack(partVector->Pt(), partVector->Rapidity(), partVector->Phi(), iparticle->GetMass(), crg, vX, vY, vZ, 0, 1, 99)){  
            continue;
          }

          KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
          if (!trw){
            continue;
          }

          if (trw->GetNormChi2(kTRUE) > ChiTot){
            continue;
          }
          trw->SetIndex(index_list[GetArrayPosition(kf)]);
          //std::cout<<"setindex: "<<index_list[GetArrayPosition(kf)]<<" "<<trw->GetIndex()<<std::endl;
          trw->SetIndexMom(index_list[GetArrayPosition(pdg_mom)]);
          //std::cout<<"setindex: "<<index_list[GetArrayPosition(pdg_mom)]<<" "<<trw->GetIndexMom()<<std::endl;
          trw->SetPdgMother((kf==pdg_mom) ? 0:pdg_mom);
          //std::cout<<"setindex: "<<pdg_mom<<" "<<trw->GetPdgMother()<<std::endl;
          trw->SetPdg(kf);
          //std::cout<<"setindex: "<<kf<<" "<<trw->GetPdg()<<std::endl;
          new (aarrtr[icount]) KMCProbeFwd(*trw);
          icount++;
        }
      }
    }
    printf("Background in array = %i \n",icount);
    tree->Fill();

  }
  f->cd();
  tree->Write();
  f->Close();
  
  TFile *outfile = new TFile(Form("bkgdistributions%s.root",suffix.Data()), "recreate");
  hNPartRec->Divide(hNPartGen);

  outfile->cd();
  hNevents->Write();
  hNPartRec->Write();
  hPDGSaved->Write();
  hPDGSurv1->Write();
  hPDGSurv2->Write();
  hPDGGen->Write();
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
  if (Eint == 20)
    return 1.9;   // gaussian y mean - 20 GeV
  else if (Eint == 40)
    return 2.22; // gaussian y mean - 40 GeV
  else if (Eint == 60)
    return 2.42;  // gaussian y mean - 60 GeV
  else if (Eint == 80)
    return 2.57;  // gaussian y mean - 80 GeV
  else if (Eint == 158)
    return 2.9;   // gaussian y mean - 160 GeV
  else if (Eint == 400)
    return 3.37;  // gaussian y mean - 400 GeV
  return 2.22;
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


int GetArrayPosition(int pdgParticle){
  int counter = 3;

  if(TMath::Abs(pdgParticle) == 211) 
    return 0;
  if(TMath::Abs(pdgParticle) == 321) 
    return 1;
  if(TMath::Abs(pdgParticle) == 2212) 
    return 2;
  
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg)
      return counter;
    counter++;
  }

  return NParticles+3;


}