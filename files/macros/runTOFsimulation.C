#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
#include "../KMCProbeFwd.h"
#include "../KMCDetectorFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TROOT.h"
#include "../GenMUONLMR.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliAODRecoDecay.h"
#include "AliDecayer.h"
#include "AliDecayerEvtGen.h"
#include "TDatabasePDG.h"
#include "./HFUtils.C"
#include "TMVA/Reader.h"
#include "KFParticle.h"
#include "DCAFitterN.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "MeasurementUtils.h"
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;

double vX = 0, vY = 0, vZ = 0; // event vertex
TDatime dt;
/*
THnSparseF* CreateSparse(Double_t mass, Int_t nbody);
Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk, int pdgMother, int pdgDaughterPos, int pdgDaughterNeg);
//measurements for phi(pdg code 333) and K0s(pdg code 333) and lambda (pdg code 3312)
//number of energy of the colliding beam
const int NEnergy = 5;
//number of particles 
const int NParticles = 5;
//array of the energy
double Elab[NEnergy] = {20,30,40,80,158};

//T parameter in the exponential pT distribution [matter/antimatter][particle][beam energy] i
//                                                    phi                        K               Lambda                Omega                    Csi
double Tslope[2][NParticles][NEnergy] = {{{196.8,237.4,244.6,239.8,298.7},{0,0,229, 223.1, 229},{244,249,258,265,301},{0,0,218,0,267},{221,233,222,227,277}},//matter e rapidity distribution [matter/antimatter][particle][beam energy]
                                         {{196.8,237.4,244.6,239.8,298.7},{0,0,226,217,226},{339,284,301,292,303},{0,0,218,0,259},{311,277,255,321,0}}};//antimatter
//sigma parameter of the gaussians of th                    phi                        K                       Lambda                Omega                    Csi
double sigma_rapidity[2][NParticles][NEnergy] = {{{0.425,0.538,0.696,0.658,1.451},{0,0,0.674, 0.743, 0.84},{0.51,0.66,0.91,0.87,0},{0,0,0.6,0,1.2},{0.45,0.56,0.76,0.71,1.18}},//matter
                                                 {{0.425,0.538,0.696,0.658,1.451},{0,0,0,0,0},{0,0,0.71,0.85,0.95},{0,0,0.6,0,1.0},{0,0,0,0,0}}};//antimatter
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
Double_t KFVertexer(AliExternalTrackParam kTrack [], Double_t kMasses[], int n_dau, float bz);
Double_t O2Vertexer(AliExternalTrackParam kTrack[], Double_t kMasses[], int n_dau, float bz);
*/
void GenerateSignalCandidates(Int_t nevents = 100000, 
				double Eint = 40, 
        TString suffix = "_PHI_TEST",
				const char *setup = "../setups/setup-EHN1_BetheBloch.txt",//"../setups/setup_EHN1-H8_short_10pixel_1.5T_BB.txt",//
				const char *privateDecayTable = "../decaytables/USERTABPHI.DEC",
        int pdgParticle = 333,
				bool writeNtuple = kTRUE, 
				bool writeNtuplGen = kTRUE, 
				bool simulateBg=kFALSE,
        int nbody = 2,
        int pdg_unstable_dau = 0,
        bool matter = true,
        int minITSHits = 5){
  
  // Generate strange particle signals and simulate detector response for decay tracks
  //nevents *= GetMultiplicity(pdgParticle,Eint,matter);
  int refreshBg = 70;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);
  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");
  
  //define the pT and rapidity probability function
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();
  Double_t lifetime = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Lifetime();
  if(lifetime==0)
    lifetime = TMath::Power(8.59,-11);
  TF1 *fpt = new TF1("fpt", "x*exp(-TMath::Sqrt(x**2+[0]**2)/[1])", ptminSG, ptmaxSG);
  fpt->SetParameter(0,mass);
  fpt->SetParameter(1, GetTslope(pdgParticle,Eint,matter)/1000);
  TF1 *fy = new TF1("fy"," exp(-0.5*((x-[0]-[2])/[1])**2)+exp(-0.5*((x+[0]-[2])/[1])**2)", yminSG, ymaxSG);
  fy->SetParameter(0, GetY0Rapidity(pdgParticle, Eint, matter));
  fy->SetParameter(1, GetSigmaRapidity(pdgParticle, Eint, matter));
  fy->SetParameter(2, GetY0(Eint));
  TF1 *fm = new TF1("fm","1/((x-[0])**2+([1]/2)**2)",mass*0.8,mass*1.2);
  fm->SetParameter(0, mass);
  fm->SetParameter(1, TMath::Hbar()*6.242*TMath::Power(10,9)/lifetime);

  TFile *fout = new TFile(Form("Signal_histos%s.root",suffix.Data()), "recreate");
  
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

  double PrimVtxZ = det->GetLayer(0)->GetZ();
  // prepare decays
  TGenPhaseSpace decay;
  TLorentzVector parent, parentgen, daugen[3], daurec[3], parentrefl, daurecswapmass[2];
  KMCProbeFwd recProbe[3];  
  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); //read the default decay table DECAY.DEC and particle table
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
    printf("-- Use existing decay modes in aliroot\n");
    fDecayer->SetForceDecay(kHadronicD); 
  }
  TClonesArray *particles = new TClonesArray("TParticle", 1000);
  TLorentzVector *mom = new TLorentzVector();

  Int_t icount = 0;

  //read the pdg code of the daughters
  int pdg_dau[2] = {0,0};
  GetPDGDaughters(pdgParticle,pdg_dau,matter);
  //if decay daughter is unstable read the pdg code of its daughters
  int pdg_dau2[2] = {0,0};
  if(nbody == 3)
    GetPDGDaughters(pdg_unstable_dau,pdg_dau2,pdg_unstable_dau > 0);

  Double_t massPos = TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass();
  Double_t massNeg = TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass();

  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    if(iev%100==0) printf(" ***************  ev = %d \n", iev);
    int nrec = 0;
    int nfake = 0;
    bool IsReconstructable = true;
    if (simulateBg && iev%refreshBg==0) det->GenBgEvent(0.,0.,PrimVtxZ);

    Double_t ptGen = fpt->GetRandom(ptminSG,ptmaxSG); 
    Double_t yGen = fy->GetRandom(yminSG,ymaxSG);
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t pxGen = ptGen * TMath::Cos(phi);
    Double_t pyGen = ptGen * TMath::Sin(phi);
    Double_t massbw = (pdgParticle == 333) ? fm->GetRandom(mass*0.97,mass*1.03) : mass;
    Double_t mt = TMath::Sqrt(ptGen * ptGen + massbw * massbw);
    Double_t pzGen = mt * TMath::SinH(yGen);
    Double_t en = mt * TMath::CosH(yGen);

    mom->SetPxPyPzE(pxGen, pyGen, pzGen, en);
    Int_t np;
    do{
      fDecayer->Decay(pdgParticle, mom);
      np = fDecayer->ImportParticles(particles);
    } while (np < 0);

    Double_t secvertgen[3][3];
    // loop on decay products
    for (int i = 0; i < np; i++) {
      TParticle *iparticle = (TParticle *)particles->At(i);
      Int_t kf = iparticle->GetPdgCode();
      vX = iparticle->Vx();
      vY = iparticle->Vy();
      vZ = iparticle->Vz()+PrimVtxZ;
      if (kf == pdgParticle){
        // Mother particle
        hYGen->Fill(iparticle->Y());
        hPtGen->Fill(iparticle->Pt());
        hYPtGen->Fill(iparticle->Y(), iparticle->Pt());
      }
      //check if the particle is one of the final decay products
      bool IsDecayDaughter = pdg_dau[0] == kf || pdg_dau[1] == kf || pdg_dau2[0] == kf || pdg_dau2[1] == kf;
      bool IsStable = kf != pdg_unstable_dau;
      
      if (IsDecayDaughter && IsStable){
        Int_t crg = (iparticle->GetPdgCode() > 0) ? 1 : -1;
        Double_t ptdau = iparticle->Pt();
        Double_t ydau = iparticle->Y();
        Int_t crg_index = (crg == 1) ? 0 : 1;
        if(nbody == 2)
        {
          daurecswapmass[crg_index].SetXYZM(iparticle->Px(), iparticle->Py(), iparticle->Pz(),(iparticle->GetMass() == massPos)? massNeg : massPos);
          hptdau[crg_index]->Fill(ptdau,ptGen);
          hydau[crg_index]->Fill(ydau);
        }
        else{
          if((pdg_dau[0] == kf || pdg_dau[1] == kf) && nrec==0){
            hptdau[crg_index]->Fill(ptdau,ptGen);
            hydau[crg_index]->Fill(ydau);
          }
          else{
            hptdau[(kf == pdg_dau2[0]) ? 1 : 2]->Fill(ptdau,ptGen);
            hydau[(kf == pdg_dau2[0]) ? 1 : 2]->Fill(ydau);
          }
        }

        if(IsReconstructable){
          if (ntgen){
            arrntgen[0] = ptGen;
            arrntgen[1] = 1;
            arrntgen[2] = yGen;
            arrntgen[3] = TMath::Sqrt(vX*vX+vY*vY+vZ*vZ);
            arrntgen[4] = TMath::Sqrt(vX*vX+vY*vY+vZ*vZ)*mass/TMath::Sqrt(ptGen*ptGen+pzGen*pzGen);
            ntgen->Fill(arrntgen);
          }

          hDistx->Fill(vX);
          hDisty->Fill(vY);
          hDistz->Fill(vZ);
          hDistxy->Fill(TMath::Sqrt(vX*vX+vY*vY));
          hDistTot->Fill(TMath::Sqrt(vX*vX+vY*vY+vZ*vZ));
          hCt->Fill(TMath::Sqrt(vX*vX+vY*vY+vZ*vZ)*mass/TMath::Sqrt(ptGen*ptGen+pzGen*pzGen));
        }

        daugen[nrec].SetXYZM(iparticle->Px(), iparticle->Py(), iparticle->Pz(), iparticle->GetMass());

      }
    }
    
  } //event loop
  
  fout->cd();  
  hCentrality->Write();

  TH1F* hHitITS;
  TH2F* hHitITS2D;
  hHitITS = det->GetHhitITS();
  hHitITS->Write();
  hHitITS = det->GetHhitITSSum();
  hHitITS->Write();
  hHitITS2D = det->GetHhitITSvsX();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsXSum();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsY();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsYSum();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsZ();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsZSum();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsYRap();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsYRapSum();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsPt();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsPtSum();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsTheta();
  hHitITS2D->Write();
  hHitITS2D = det->GetHhitITSvsThetaSum();
  hHitITS2D->Write();

  hMassAll->Write();
  hMassTrue->Write();
  hMassFake->Write();
  hMassVsPt->Write();
  hMassVsY->Write();

  if(nbody == 2){
    hMassRefl->Write();
    hMassReflVsPt->Write();
    hMassReflVsY->Write();
    hDCA->Write();
    hDCAx->Write();
    hDCAy->Write();
    hDCAz->Write();
  }
  else{
    hDCA12->Write();
    hDCA13->Write();
    hDCA23->Write();
  }

  //efficiency vs pT /rapidity
  for(int i = 1; i <= hPtGen->GetNbinsX(); i++){
    double n = hPtGen->GetBinContent(i);
    if(n<0)
      hPtGen->SetBinContent(i,0);
  }
  
  for(int i = 1; i <= hYGen->GetNbinsX(); i++){
    double n = hYGen->GetBinContent(i);
    if(n<0)
      hYGen->SetBinContent(i,0);
  }
  
  hYEff->Divide(hYGen);
  hPtEff->Divide(hPtGen);

  for(int i = 1; i <= hPtEff->GetNbinsX(); i++){
    double eff = hPtEff->GetBinContent(i);
    double n = hPtGen->GetBinContent(i);
    if(n>0)
      hPtEff->SetBinError(i,TMath::Sqrt(eff*(1-eff)/n));
    else
      hPtEff->SetBinError(i,0);
  }

  for(int i = 1; i <= hYEff->GetNbinsX(); i++){
    double eff = hYEff->GetBinContent(i);
    double n = hYGen->GetBinContent(i);
    if(n>0)
      hYEff->SetBinError(i,TMath::Sqrt(eff*(1-eff)/n));
    else
      hYEff->SetBinError(i,0);
  }
  hMassVsOpen->Write();
  hYFake->Write();
  hPtFake->Write();
  hPxFake->Write();
  hPyFake->Write();
  hPzFake->Write();
  hYPtGen->Write();
  hPtGen->Write();
  hYGen->Write();
  hYEff->Write();
  hPtEff->Write();
  hYPtRecoAll->Write();  
  hYPtRecoFake->Write();
  hPtRecoAll->Write();
  hPtGenRecoAll->Write();
  hPtRecoVsGenAll->Write();
  hDiffPtRecoGenAll->Write();
  hYRecoAll->Write();
  hYGenRecoAll->Write();
  hPtRecoFake->Write();
  hArmPod->Write();
  hDistXY->Write();
  hDistx->Write();
  hDisty->Write();
  hDistxy->Write();
  hDistz->Write();
  hDistTot->Write();
  hDistxRec->Write();
  hDistyRec->Write();
  hDistxyRec->Write();
  hDistzRec->Write();
  hDistTotRec->Write();
  hCt->Write();
  hDist->Write();
  hDistgenXY->Write();
  hDistgen->Write();
  hCosp->Write();
  hResVxVsPt->Write();
  hResVyVsPt->Write();
  hResVzVsPt->Write();
  hResPtVsPt->Write();
  hResPtVsY->Write();
  hResVxVsY->Write();
  hResVyVsY->Write();
  hResVzVsY->Write();
  hResPxVsY->Write();
  hResPyVsY->Write();
  hResPzVsY->Write();
  hResDist->Write();
  hResDistXY->Write();
  hChi2True->Write();
  hChi2Fake->Write();
  hd0XYprod->Write();

  for(int i = 0; i < nbody; ++i){
    hd0XY[i]->Write();
    hResPxDauVsPt[i]->Write();
    hResPyDauVsPt[i]->Write();
    hResPzDauVsPt[i]->Write();
    hResPxDauVsY[i]->Write();
    hResPyDauVsY[i]->Write();
    hResPzDauVsY[i]->Write();
    hPxDauFake[i]->Write();
    hPyDauFake[i]->Write();
    hPzDauFake[i]->Write();
  }
  hNevents->Write();
  if (ntcand){
    fnt->cd();
    ntcand->Write();
    ntgen->Write();
    fnt->Close();
  }
  
  TFile fout2(Form("DecayHistos%s.root",suffix.Data()), "RECREATE");
  hPtGen->Write();
  hYGen->Write();
  hMGen->Write();
  if(nbody == 2){
    hptdau[0]->SetName("hptdauPos");
    hydau[0]->SetName("hydauPos");
    hptdau[1]->SetName("hptdauNeg");
    hydau[1]->SetName("hydauNeg");
    hydau2D->Write();
  }
  else{
    hptdau[0]->SetName(Form("hptdau_%i",(pdg_dau[0] == pdg_unstable_dau) ? pdg_dau[1] : pdg_dau[0]));
    hydau[0]->SetName(Form("hydau_%i",(pdg_dau[0] == pdg_unstable_dau) ? pdg_dau[1] : pdg_dau[0]));
    hptdau[1]->SetName(Form("hptdau2_%i",pdg_dau2[0]));
    hydau[1]->SetName(Form("hydau2_%i",pdg_dau2[0]));
    hptdau[2]->SetName(Form("hptdau2%i",pdg_dau2[1]));
    hydau[2]->SetName(Form("hydau2_%i",pdg_dau2[1]));
  }
  for(int i = 0; i < nbody; i++){
    hptdau[i]->Write();
    hydau[i]->Write();
  }
  fout2.Close();
  fout->Close();

  f->cd();
  tree->Write();
  f->Close();
  std::cout<<"rejected solve:"<<count_solve<<std::endl;
  std::cout<<"rejected mcwinner:"<<count_mcwin<<std::endl;
  std::cout<<"rejected chi2:"<<count_chi2<<std::endl;
  std::cout<<"tot rejected: "<<count_solve+count_mcwin+count_chi2<<std::endl;

}