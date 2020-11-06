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
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;

// settings for signal generation
double yminSG = -10.; // min y to generate
double ymaxSG = 10.;  //
double ptminSG = 0.;
double ptmaxSG = 10; //Elena's change, it was 3 GeV/c

double vX = 0, vY = 0, vZ = 0; // event vertex

THnSparseF* CreateSparse();
TDatime dt;
Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk, int pdgMother, int pdgDaughterPos, int pdgDaughterNeg);
//measurements for phi(pdg code 333) and K0s(pdg code 333)
//per il K0s ho usato le misure di NA49 del K+
const int NEnergy = 5;
const int NParticles =2
double Elab[NEnergy] = {20,30,40,80,158};
double Tslope[NParticles][NEnergy] = {{196.8,237.4,244.6,239.8,298.7},{0,0,244.6,239.8,298.7}};
double sigma_rapidity[NParticles][NEnergy] = {{0.425,0.538,0.696,0.658,1.2},{0,0,0.694,0.742,0.839}};
double y0_rapidity[NParticles][NEnergy] = {{0.425,0.538,0.487,0.682,1.451},{0,0,0.725,0.792,0.88}};
double multiplicity[NParticles][NEnergy] = {{1.89,1.84,2.55,4.04,8.46},{0,0,59.1,76.9,103.0}};
double pdg_code[NParticles] = {333,310};

double GetTslope(int pdgParticle, double Eint);
double GetSigmaRapidity(int pdgParticle, double Eint);
double GetY0Rapidity(int pdgParticle, double Eint);
double GetMultiplicity(int pdgParticle, double Eint);
double GetY0(double Eint);
const int nEv = 10000;

void GenerateSignalCandidates(Int_t nevents = 10000, 
				double Eint = 158, 
        TString suffix = "_k0s",
				const char *setup = "../setups/setup_short_5pixel_1.5T.txt",
				const char *privateDecayTable = "../decaytables/USERTABK0.DEC",
				bool writeNtuple = kTRUE, 
				bool simulateBg=kTRUE,
        int pdgParticle = 310){
  
  // Generate strange particle in two body decay signals and simulate detector response for decay tracks

  nevents *= GetMultiplicity(pdgParticle,Eint);
  int refreshBg = 1000;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);
  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");
  
  //define the pT and rapidity probability function
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();
  TF1 *fpt = new TF1("fpt","x*expo(-TMath::Sqrt(x**2+[0]**2)/[1])",ptminSG,ptmaxSG);
  fpt->SetParameter(0,mass);
  fpt->SetParameter(1,GetTslope(pdgParticle,Eint));
  TF1 *fy = new TF1("fy","expo(-0.5*((x-[0])/[1])**2)+expo(-0.5*((x+[0])/[1])**2) ",yminSG,ymaxSG);
  fy->SetParameter(0,GetY0Rapidity(pdgParticle,Eint));
  fy->SetParameter(1,GetSigmaRapidity(pdgParticle,Eint));

  TH2F *hptDauPlus = new TH2F("hptDauPlus", "positive daughter;#it{p}_{T+} (GeV/#it{c});#it{p}_{TM} (GeV/#it{c});counts", 50,0.,10.,50, 0., 10.);
  TH2F *hptDauMinus = new TH2F("hptDauMinus", "negative daughter;#it{p}_{T-} (GeV/#it{c});#it{p}_{TM} (GeV/#it{c});counts", 50, 0.,10.,50,0., 10.);
  TH1D *hyDauPlus = new TH1D("hyDauPlus", "y positive daughter;y_{+};counts", 50, 0., 5.);
  TH1D *hyDauMinus = new TH1D("hyDauMinus", "y negative daughter;y_{-};counts", 50, 0., 5.);
  TH2F *hyDau2D = new TH2F("hyDau2D", "y negative daughter vs y positive daughter;y_{-};y_{+};counts", 50, 0., 5., 50, 0., 5.);
  
  //Magnetic field and detector parameters
  MagField *mag = new MagField(1);
  int BNreg = mag->GetNReg();
  const double *BzMin = mag->GetZMin();
  const double *BzMax = mag->GetZMax();
  const double *BVal;
  printf("*************************************\n");
  printf("number of magnetic field regions = %d\n", BNreg);
  
  TFile *fout = new TFile(Form("Signal_histos%s.root",suffix.Data()), "recreate");
  
  for (int i = 0; i < BNreg; i++){
    BVal = mag->GetBVals(i);
    printf("*** Field region %d ***\n", i);
    if (i == 0){
      printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    }else if (i == 1){
      printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    }
  }
  
  //int outN = nev/10;
  //if (outN<1) outN=1;

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT

  det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); //NA60+
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
  det->SetMaxChi2Vtx(20);  // fiducial cut on chi2 of convergence to vtx  
  // IMPORTANT FOR NON-UNIFORM FIELDS
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); //NA60+
  //det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  det->BookControlHistos();
  //
  
  
  // prepare decays
  TGenPhaseSpace decay;
  TLorentzVector parentgen, daugen[2], parent, daurec[2], parentrefl, daurecswapmass[2]; 
  KMCProbeFwd recProbe[2];  
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
  //fDecayer->ForceDecay();
  TClonesArray *particles = new TClonesArray("TParticle", 1000);
  TLorentzVector *mom = new TLorentzVector();
  
  
  TH2F* hYPtGen = new TH2F("hYPtGen", "y-#it{p}_{T} corr match;y;#it{p}_{T};counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtGen = new TH1D("hPtGen", "#it{p}_{T} gen;#it{p}_{T} (GeV/#it{c});counts", 40, ptminSG, ptmaxSG);
  TH1D* hYGen = new TH1D("hYGen", "y full phase space;y;counts", 80., 1., 5.4);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "y-#it{p}_{T} all match;y;#it{p}_{T};counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Reconstructed #it{p}_{T} all match;#it{p}_{T};counts", 40, ptminSG, ptmaxSG);  
  TH1D* hPtGenRecoAll = new TH1D("hPtGenRecoAll", "Generated #it{p}_{T} all match;#it{p}_{T};counts", 40, ptminSG, ptmaxSG);
  TH2F* hPtRecoVsGenAll = new TH2F("hPtRecoVsGenAll"," ; Generated #it{p}_{T} ; Reconstructed #it{p}_{T}",40, ptminSG, ptmaxSG,40, ptminSG, ptmaxSG);
  TH2F* hDiffPtRecoGenAll = new TH2F("hDiffPtRecoGenAll"," ; Generated #it{p}_{T} ; Reco #it{p}_{T} - Gen #it{p}_{T}",40, ptminSG, ptmaxSG,100,-0.2,0.2);

  TH1D* hYRecoAll = new TH1D("hYRecoAll", "Reconstructed y all match;y;counts", 80., 1., 5.4);
  TH1D* hYGenRecoAll = new TH1D("hYGenRecoAll", "Generated y all match;y;counts", 80., 1., 5.4);
  TH2F* hYPtRecoFake = new TH2F("hYPtRecoFake", "y-#it{p}_{T} fake match;;counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoFake = new TH1D("hPtRecoFake", "#it{p}_{T} fake match;;counts", 40, ptminSG, ptmaxSG);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 1.5*mass);
  TH1D* hMassFake = new TH1D("hMassFake", "Mass fake match;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 1.5*mass);
  TH1D* hMassRefl = new TH1D("hMassRefl", "Mass reflections;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 1.5*mass);

  TH2F *hDistXY = new TH2F("hDistXY", ";;#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", ";;#it{p}_{T} (GeV/#it{c});counts", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", ";;#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", ";;#it{p}_{T} (GeV/#it{c});counts", 300, 0, 10, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", ";;#it{p}_{T} (GeV/#it{c});counts", 300, -1, 1, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", ";;#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.01, 0.01, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hMassVsPt = new TH2F("hMassVsPt", ";m (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", 200, 0.5*mass, 1.5*mass, 6, 0, 3);
  TH2F *hMassVsY = new TH2F("hMassVsY", ";m (GeV/#it{c}^{2});y;counts", 200, 0.5*mass, 1.5*mass, 10, 0, 5);
  TH2F *hMassReflVsPt = new TH2F("hMassReflVsPt", ";m (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", 200, 0.5*mass, 1.5*mass, 6, 0, 3);
  TH2F *hMassReflVsY = new TH2F("hMassReflVsY", ";m (GeV/#it{c}^{2});y;counts", 200, 0.5*mass, 1.5*mass, 10, 0, 5);
  
  TH2F *hResVx = new TH2F("hResVx", ";V_{xgen}-V_{xrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", ";V_{ygen}-V_{yrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", ";V_{zgen}-V_{zrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, 0, 3);
  TH2F *hResVxVsY = new TH2F("hResVxVsY", ";V_{xgen}-V_{xrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", ";V_{ygen}-V_{yrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", ";V_{zgen}-V_{zrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResPxVsY = new TH2F("hResPxVsY", ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5); //for Kaons
  TH2F *hResPyVsY = new TH2F("hResPyVsY", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5);
  TH2F *hResPzVsY = new TH2F("hResPzVsY", ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5);
  TH2F *hResDist = new TH2F("hResDist", ";d_{gen}-d_{rec} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", ";d_{xygen}-d_{xyrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH1D *hNevents = new TH1D("hNevents", ";;counts", 1, 0, 1);
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  if (writeNtuple)
  {
    fnt = new TFile(Form("fntSig%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "mass:pt:y:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax", 32000);
  }

  Float_t arrnt[11];
  double y0 = GetY0(Eint);
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    if(iev%100==0) printf(" ***************  ev = %d \n", iev);
    int nrec = 0;
    int nfake = 0;
    double pxyz[3];
    
    if (simulateBg && (iev%refreshBg)==0) det->GenBgEvent(0.,0.,0.);
    Double_t ptGenD = fpt->GetRandom(); // get D0 distribution from file
    Double_t yGenD = fy->GetRandom()+y0;
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t pxGenD = ptGenD * TMath::Cos(phi);
    Double_t pyGenD = ptGenD * TMath::Sin(phi);
    
    Double_t mt = TMath::Sqrt(ptGenD * ptGenD + mass * mass);
    Double_t pzGenD = mt * TMath::SinH(yGenD);
    Double_t en = mt * TMath::CosH(yGenD);

    Double_t massKaon = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    Double_t massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    mom->SetPxPyPzE(pxGenD, pyGenD, pzGenD, en);
    Int_t np;
    do{
      fDecayer->Decay(pdgParticle, mom);
      np = fDecayer->ImportParticles(particles);
    } while (np < 0);

    Int_t arrpdgdau[2];
    Double_t ptK = -999.;
    Double_t ptPi = -999.;
    Double_t yK=-999.;
    Double_t yPi = -999.;
    Int_t icount = 0;
    Double_t secvertgenK[3]={0.,0.,0.};
    Double_t secvertgenPi[3]={0.,0.,0.};
    Int_t pdgDaughterNeg = -312;
    Int_t pdgDaughterPos = 312;
    // loop on decay products
    for (int i = 0; i < np; i++) {
      TParticle *iparticle1 = (TParticle *)particles->At(i);
      Int_t kf = iparticle1->GetPdgCode();
      vX = iparticle1->Vx();
      vY = iparticle1->Vy();
      vZ = iparticle1->Vz();
      if (TMath::Abs(kf) == pdgParticle){
        // D0 particle
        hYGen->Fill(iparticle1->Y());
        hPtGen->Fill(iparticle1->Pt());
        hYPtGen->Fill(iparticle1->Y(), iparticle1->Pt());
        //printf("mother part = %d, code=%d, pt=%f, y=%f \n", i, kf, iparticle1->Pt(), iparticle1->Y());
      }

      bool IsDecayDaughter = (kf != 22) && (TMath::Abs(kf) != pdgParticle);
      if (IsDecayDaughter){
        // daughters that can be reconstructed: K and pi
        Double_t e = iparticle1->Energy();
        Double_t px = iparticle1->Px();
        Double_t py = iparticle1->Py();
        Double_t pz = iparticle1->Pz();	    
        TLorentzVector *pDecDau = new TLorentzVector(0., 0., 0., 0.);
        pDecDau->SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
        Int_t crg=1;
        if(iparticle1->GetPdgCode()<0) crg=-1;
        if (!det->SolveSingleTrack(pDecDau->Pt(), pDecDau->Rapidity(), pDecDau->Phi(), iparticle1->GetMass(), crg, vX, vY, vZ, 0, 1, 99)) continue;
        KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
        if (!trw) continue;
        if (trw->GetNormChi2(kTRUE) > ChiTot) continue;
        nrec++;
            
        nfake += trw->GetNFakeITSHits();
        trw->GetPXYZ(pxyz);
        
        if (kf > 0){
          // Positive daughter
          pdgDaughterPos = kf;
          ptK = iparticle1->Pt();
          yK = iparticle1->Y();
          hptDauPlus->Fill(ptGenD,ptK);
          hyDauPlus->Fill(iparticle1->Y());
          secvertgenK[0] = iparticle1->Vx();
          secvertgenK[1] = iparticle1->Vy();
          secvertgenK[2] = iparticle1->Vz();
          daugen[0].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
          daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], iparticle1->GetMass());
          daurecswapmass[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2],massPion);
          recProbe[0] = *trw;
        }else if (kf < 0){
          // Negative daughter
          pdgDaughterNeg = kf;
          ptPi = iparticle1->Pt();
          yPi = iparticle1->Y();
          hptDauMinus->Fill(ptGenD,ptPi);
          hyDauMinus->Fill(iparticle1->Y());
          secvertgenPi[0] = iparticle1->Vx();
          secvertgenPi[1] = iparticle1->Vy();
          secvertgenPi[2] = iparticle1->Vz();
          daugen[1].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
          daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2],  iparticle1->GetMass());
          daurecswapmass[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2],massKaon);
          recProbe[1] = *trw;
        }
      }
    }
    if (ptK > 0 && ptPi > 0) hyDau2D->Fill(yPi, yK);
    if (nrec < 2) continue;
    
    recProbe[0].PropagateToDCA(&recProbe[1]);
    
    parent = daurec[0];
    parent += daurec[1];
    parentgen = daugen[0];
    parentgen += daugen[1];
    parentrefl = daurecswapmass[0];
    parentrefl += daurecswapmass[1];
      
    Double_t  ptRecD=parent.Pt();
    Double_t  massRecD=parent.M();
    Double_t  massRecReflD=parentrefl.M();
    Double_t yRecD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
    hYPtRecoAll->Fill(yRecD, ptRecD);
    hPtRecoAll->Fill(ptRecD);
    hPtGenRecoAll->Fill(ptGenD);
    hPtRecoVsGenAll->Fill(ptGenD,ptRecD);
    hDiffPtRecoGenAll->Fill(ptGenD,(ptRecD-ptGenD));
    hYRecoAll->Fill(yRecD);
    hYGenRecoAll->Fill(yGenD);
    hMassAll->Fill(massRecD);
    hMassRefl->Fill(massRecReflD);
    if (nfake > 0){
      hYPtRecoFake->Fill(yRecD, ptRecD);
      hPtRecoFake->Fill(ptRecD);
      hMassFake->Fill(massRecD);
    }
    hMassVsPt->Fill(massRecD,ptRecD);
    hMassVsY->Fill(massRecD,yRecD);
    hMassReflVsPt->Fill(massRecReflD,ptRecD);
    hMassReflVsY->Fill(massRecReflD,yRecD);
   
    Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
    Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
    Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
    Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    // printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
    hDCA->Fill(dca, ptRecD);
    hDCAx->Fill(d1, ptRecD);
    hDCAy->Fill(d2, ptRecD);
    hDCAz->Fill(d3, ptRecD);
    
    //      Double_t xP = (recProbe[1].GetX() + recProbe[0].GetX()) / 2.;
    //      Double_t yP = (recProbe[1].GetY() + recProbe[0].GetY()) / 2.;
    //      Double_t zP = (recProbe[1].GetZ() + recProbe[0].GetZ()) / 2.;
    
    Double_t xP, yP, zP;
    ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
    Double_t residVx=10000.*(xP - secvertgenK[0]);
    Double_t residVy=10000.*(yP - secvertgenK[1]);
    Double_t residVz=10000.*(zP - secvertgenK[2]);
    hResVx->Fill(residVx, ptRecD);
    hResVy->Fill(residVy, ptRecD);
    hResVz->Fill(residVz, ptRecD);
    hResVxVsY->Fill(residVx, yRecD);
    hResVyVsY->Fill(residVy, yRecD);
    hResVzVsY->Fill(residVz, yRecD);
    
    hResPx->Fill(daurec[0].Px() - daugen[0].Px(), ptRecD);
    hResPy->Fill(daurec[0].Py() - daugen[0].Py(), ptRecD);
    hResPz->Fill(daurec[0].Pz() - daugen[0].Pz(), ptRecD);
    hResPx->Fill(daurec[1].Px() - daugen[1].Px(), ptRecD);
    hResPy->Fill(daurec[1].Py() - daugen[1].Py(), ptRecD);
    hResPz->Fill(daurec[1].Pz() - daugen[1].Pz(), ptRecD);

    hResPxVsY->Fill(daurec[0].Px() - daugen[0].Px(), yRecD);
    hResPyVsY->Fill(daurec[0].Py() - daugen[0].Py(), yRecD);
    hResPzVsY->Fill(daurec[0].Pz() - daugen[0].Pz(), yRecD);
    hResPxVsY->Fill(daurec[1].Px() - daugen[1].Px(), yRecD);
    hResPyVsY->Fill(daurec[1].Py() - daugen[1].Py(), yRecD);
    hResPzVsY->Fill(daurec[1].Pz() - daugen[1].Pz(), yRecD);
    
    // cout << "secvert generated Pion: " << secvertgenPi[0] << "  " << secvertgenPi[1] << "  " << secvertgenPi[2] << endl;
    // cout << "Reco Vert  Pion: " << xP << "  " << yP << "  " << zP << endl;
    
    Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
    Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
    Float_t distgen = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1] + secvertgenPi[2] * secvertgenPi[2]);
    Float_t distgenXY = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1]);
    // printf("dist = %f , distXY=%f , dx=%f, dy=%f, dz=%f z1=%f, z2=%f \n", dist, distXY, xP, yP, zP, recProbe[0].GetZ(), recProbe[1].GetZ());
    // printf("distgen = %f , distgenXY=%f \n", distgen, distgenXY);
    
    Double_t vsec[3] = {xP, yP, zP};
    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
    Double_t cts = CosThetaStar(parent,daurec[0],pdgParticle,pdgDaughterPos,pdgDaughterNeg);
    Double_t ipD = ImpParXY(vprim, vsec, parent);
    hCosp->Fill(cosp, ptRecD);
    // printf(" ***** ***** cos point = %f \n", cosp);
    //if (cosp < -0.98)
    //    printf("SMALL COSPOINT");
    
    hResDist->Fill(dist - distgen, ptRecD);
    hResDistXY->Fill(distXY - distgenXY, ptRecD);
    
    //recProbe[0].PropagateToDCA(&recProbe[1]);
    
    // hYPtAll->Fill(parent.Y(), ptRecD);
    // hPtAll->Fill(ptRecD);
    hDistXY->Fill(distXY, ptRecD);
    hDist->Fill(dist, ptRecD);
    hDistgenXY->Fill(distgenXY, ptRecD);
    hDistgen->Fill(distgen, ptRecD);
      
    //AliExternalTrackParam *track1 = (AliExternalTrackParam *)recProbe[0].GetTrack();
    //AliExternalTrackParam *track2 = (AliExternalTrackParam *)recProbe[1].GetTrack();
    recProbe[0].PropagateToZBxByBz(0);
    Double_t d0x1 = recProbe[0].GetX();
    Double_t d0y1 = recProbe[0].GetY();
    Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
    if (d0x1 < 0)
      d0xy1 *= -1;
    
    recProbe[1].PropagateToZBxByBz(0);
    Double_t d0x2 = recProbe[1].GetX();
    Double_t d0y2 = recProbe[1].GetY();
    Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
    if (d0x2 < 0)
      d0xy2 *= -1;
    
    // printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
    
    hd0XYprod->Fill(d0xy1 * d0xy2, ptRecD);
    hd0XY1->Fill(d0xy1, ptRecD);
    hd0XY2->Fill(d0xy2, ptRecD);
      
    arrsp[0] = massRecD;
    arrsp[1] = ptRecD;
    arrsp[2] = yRecD;
    arrsp[3] = dist;
    arrsp[4] = cosp;
    arrsp[5] = TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2));
    arrsp[6] = d0xy1 * d0xy2;
    arrsp[7] = dca;
    arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
    arrsp[9] = TMath::Abs(ipD);	    
    arrsp[10] = cts;      
    hsp->Fill(arrsp);
    
    if (ntcand){
      arrnt[0] = massRecD;
      arrnt[1] = ptRecD;
      arrnt[2] = yRecD;
      arrnt[3] = dist;
      arrnt[4] = cosp;
      arrnt[5] = d0xy1;
      arrnt[6] = d0xy2;
      arrnt[7] = d0xy1 * d0xy2;
      arrnt[8] = dca;
      arrnt[9] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      arrnt[10] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      ntcand->Fill(arrnt);
    }
  } //event loop
  
  hMassAll->SetLineColor(kBlue);
  hMassAll->Draw();
  hMassAll->SetMinimum(0.1);
  hMassFake->SetLineColor(kRed);
  hMassFake->Draw("same");
  
  fout->cd();  
  hMassAll->Write();
  hMassFake->Write();
  hMassRefl->Write();
  hMassVsPt->Write();
  hMassVsY->Write();
  hMassReflVsPt->Write();
  hMassReflVsY->Write();
  hYPtGen->Write();
  hPtGen->Write();
  hYGen->Write();
  hYPtRecoAll->Write();  
  hYPtRecoFake->Write();
  hPtRecoAll->Write();
  hPtGenRecoAll->Write();
  hPtRecoVsGenAll->Write();
  hDiffPtRecoGenAll->Write();
  hYRecoAll->Write();
  hYGenRecoAll->Write();
  hPtRecoFake->Write();
  hDistXY->Write();
  hDist->Write();
  hDistgenXY->Write();
  hDistgen->Write();
  hCosp->Write();
  hDCA->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAz->Write();
  hResVx->Write();
  hResVy->Write();
  hResVz->Write();
  hResPx->Write();
  hResPy->Write();
  hResPz->Write();
  hResVxVsY->Write();
  hResVyVsY->Write();
  hResVzVsY->Write();
  hResPxVsY->Write();
  hResPyVsY->Write();
  hResPzVsY->Write();
  hResDist->Write();
  hResDistXY->Write();
  hd0XYprod->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hNevents->Write();
  hsp->Write();
  if (ntcand){
    fnt->cd();
    ntcand->Write();
    fnt->Close();
  }
  // TCanvas *ccdau = new TCanvas();
  // ccdau->Divide(3, 2);
  // ccdau->cd(1)->SetLogy();
  // hPtGen->Draw();
  // ccdau->cd(2)->SetLogy();
  // hptDauPlus->Draw();
  // ccdau->cd(3)->SetLogy();
  // hptDauMinus->Draw();
  // ccdau->cd(4)->SetLogy();
  // hYGen->Draw();
  // ccdau->cd(5)->SetLogy();
  // hyDauPlus->Draw();
  // ccdau->cd(6)->SetLogy();
  // hyDauMinus->Draw();
  
  TFile fout2(Form("DecayHistos%s.root",suffix.Data()), "RECREATE");
  // TFile fout2("DecayHistostest.root", "RECREATE");
  hPtGen->Write();
  hptDauPlus->Write();
  hptDauMinus->Write();
  hyDauPlus->Write();
  hyDauMinus->Write();
  hYGen->Write();
  hyDau2D->Write();
  fout2.Close();

  fout->Close();
}



void MakeCombinBkgCandidates(const char* trackTreeFile="treeBkgEvents.root",
             TString suffix = "_phi",
			       Int_t nevents = 999999, 
			       Int_t writeNtuple = kTRUE,
             int pdgMother = 310,
             int pdgDaughterPos = 211,
             int pdgDaughterNeg = -211){

  // Read the TTree of tracks produced with runBkgVT.C
  // Create D0 combinatorial background candidates (= OS pairs of tracks)
  // Store in THnSparse and (optionally) TNtuple

  TFile *filetree = new TFile(trackTreeFile);
  TTree *tree = (TTree *)filetree->Get("tree");
  TClonesArray *arr = 0;
  tree->SetBranchAddress("tracks", &arr);
  Int_t entries = tree->GetEntries();
  printf("Number of events in tree = %d\n",entries);
  if(nevents>entries) nevents=entries;
  else printf(" --> Analysis performed on first %d events\n",nevents);

  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  
  TFile *fout = new TFile("Bkg-histos.root", "recreate");
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 0., 2.5);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", "", 100, -0.01, 0.01, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  
  TH2F *hResVx = new TH2F("hResVx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 100, -1, 1, 30, 0, 3);
  TH2F *hCosThStVsMass = new TH2F("hCosThStVsMass", "", 50, 1.5, 2.5, 40, -1, 1);
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  
  TH1D *hcand = new TH1D("hcand", "", 1000, 0, 500000);
  TH1D *hcandpeak = new TH1D("hcandpeak", "", 500, 0, 15000);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
    
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  Float_t arrnt[11];
  if (writeNtuple){
    fnt = new TFile(Form("fntBkg%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "mass:pt:y:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax", 32000);
  }

  KMCProbeFwd recProbe[2];
  TLorentzVector parent, daurec[2];
  Double_t massM = TDatabasePDG::Instance()->GetParticle(pdgMother)->Mass();
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    tree->GetEvent(iev);
    Int_t arrentr = arr->GetEntriesFast();
    
    for (Int_t itr = 0; itr < arrentr; itr++){
      KMCProbeFwd *tr1 = (KMCProbeFwd *)arr->At(itr);
      // cout << "tr P=" << tr1->GetP() << endl;
      Float_t ch1 = tr1->GetCharge();
      for (Int_t itr2 = itr; itr2 < arrentr; itr2++){
        KMCProbeFwd *tr2 = (KMCProbeFwd *)arr->At(itr2);
        Float_t ch2 = tr2->GetCharge();
        if (ch1 * ch2 > 0) continue;
        if (ch1 < 0){ //convention: first track negative
          recProbe[0] = *tr1;
          recProbe[1] = *tr2;
        }else if (ch2 < 0){
          recProbe[0] = *tr2;
          recProbe[1] = *tr1;
        }
        recProbe[0].PropagateToDCA(&recProbe[1]);
        Double_t pxyz[3];
        recProbe[0].GetPXYZ(pxyz);
        
        Double_t pxyz2[3];
        recProbe[1].GetPXYZ(pxyz2);
	
        for(Int_t iMassHyp=0; iMassHyp<2; iMassHyp++){
          // mass hypothesis: Kpi, piK
          Int_t iKaon=-1;
          if(iMassHyp==0){
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], TDatabasePDG::Instance()->GetParticle(pdgDaughterNeg)->Mass());
            daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdgDaughterPos)->Mass());
            iKaon=0;
          }else{
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], TDatabasePDG::Instance()->GetParticle(pdgDaughterPos)->Mass());
            daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdgDaughterNeg)->Mass());
            iKaon=1;
          }
          parent = daurec[0];
          parent += daurec[1];
          countCand++;
          Float_t ptD=parent.Pt();
          Float_t invMassD=parent.M();
          Float_t yD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
          hYPtRecoAll->Fill(yD, ptD);
          hPtRecoAll->Fill(ptD);
          hMassAll->Fill(invMassD);
          if(invMassD>0.5*massM  && invMassD<1.5*massM){
            // range to fill histos
            //TODO: controllare questi valori
            if(invMassD>0.8*massM && invMassD<1.2*massM) countCandInPeak++;
            Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
            Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
            Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
            Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
            
            //printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
            hDCA->Fill(dca, ptD);
            hDCAx->Fill(d1, ptD);
            hDCAy->Fill(d2, ptD);
            hDCAz->Fill(d3, ptD);
            // Float_t xP = (recProbe[1].GetX() + recProbe[0].GetX()) / 2.;
            // Float_t yP = (recProbe[1].GetY() + recProbe[0].GetY()) / 2.;
            // Float_t zP = (recProbe[1].GetZ() + recProbe[0].GetZ()) / 2.;
            Double_t xP, yP, zP;
            ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
            Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
            Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
            Double_t vsec[3] = {xP, yP, zP};
            Double_t cosp = CosPointingAngle(vprim, vsec, parent);
            Double_t cts = CosThetaStar(parent,daurec[iKaon],pdgMother,pdgDaughterPos,pdgDaughterNeg);
            Double_t ipD = ImpParXY(vprim, vsec, parent);
            hCosp->Fill(cosp, ptD);
            hCosThStVsMass->Fill(invMassD,cts);
            //printf(" ***** ***** cos point = %f \n", cosp);	    
            hDistXY->Fill(distXY, ptD);
            hDist->Fill(dist, ptD);
            
            recProbe[0].PropagateToZBxByBz(0);
            Double_t d0x1 = recProbe[0].GetX();
            Double_t d0y1 = recProbe[0].GetY();
            Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
            if (d0x1 < 0) d0xy1 *= -1;
            
            recProbe[1].PropagateToZBxByBz(0);
            Double_t d0x2 = recProbe[1].GetX();
            Double_t d0y2 = recProbe[1].GetY();
            Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
            if (d0x2 < 0) d0xy2 *= -1;
              
            //printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
            
            hd0XYprod->Fill(d0xy1 * d0xy2, ptD);
            hd0XY1->Fill(d0xy1, ptD);
            hd0XY2->Fill(d0xy2, ptD);
            arrsp[0] = invMassD;
            arrsp[1] = ptD;
            arrsp[2] = yD;
            arrsp[3] = dist;
            arrsp[4] = cosp;
            arrsp[5] = TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2));
            arrsp[6] = d0xy1 * d0xy2;
            arrsp[7] = dca;
            arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
            arrsp[9] = TMath::Abs(ipD);	    
            arrsp[10] = cts;
            hsp->Fill(arrsp);
            
            if (ntcand){
              arrnt[0] = invMassD;
              arrnt[1] = ptD;
              arrnt[2] = yD;
              arrnt[3] = dist;
              arrnt[4] = cosp;
              arrnt[5] = d0xy1;
              arrnt[6] = d0xy2;
              arrnt[7] = d0xy1 * d0xy2;
              arrnt[8] = dca;
              arrnt[9] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
              arrnt[10] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
              ntcand->Fill(arrnt);
            }
	    
          } // check on inv mass
        } // loop on mass hypothesis
      } // loop on first track
    } // loop on second track
    
    hcand->Fill(countCand);
    hcandpeak->Fill(countCandInPeak);
    printf(" --> Event %d, tot candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hNevents->Write();
  hcand->Write();
  hcandpeak->Write();  
  hMassAll->Write();
  hYPtRecoAll->Write();
  hPtRecoAll->Write();
  hDistXY->Write();
  hDist->Write();
  hDCA->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAz->Write();
  hCosp->Write();
  hCosThStVsMass->Write();
  hd0XYprod->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hsp->Write();
  fout->Close();
  if (ntcand){
    fnt->cd();
    ntcand->Write();
    fnt->Close();
  }
}


Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk, int pdgMother, int pdgDaughterPos, int pdgDaughterNeg) {

  Double_t massMoth = TDatabasePDG::Instance()->GetParticle(pdgMother)->Mass();
  Double_t massK = TDatabasePDG::Instance()->GetParticle(pdgDaughterPos)->Mass();
  Double_t massPi = TDatabasePDG::Instance()->GetParticle(pdgDaughterNeg)->Mass();

  Double_t pStar = TMath::Sqrt((massMoth*massMoth-massK*massK-massPi*massPi)*(massMoth*massMoth-massK*massK-massPi*massPi)-4.*massK*massK*massPi*massPi)/(2.*massMoth);

  Double_t pMoth=parent.P();
  Double_t e=TMath::Sqrt(massMoth*massMoth+pMoth*pMoth);
  Double_t beta = pMoth/e;
  Double_t gamma = e/massMoth;
  TVector3 momDau(dauk.Px(),dauk.Py(),dauk.Pz());
  TVector3 momMoth(parent.Px(),parent.Py(),parent.Pz());
  Double_t qlProng=momDau.Dot(momMoth)/momMoth.Mag();
  Double_t cts = (qlProng/gamma-beta*TMath::Sqrt(pStar*pStar+massK*massK))/pStar;

  return cts;
}

THnSparseF* CreateSparse(){
  const Int_t nAxes=11;
  TString axTit[nAxes]={"Inv. mass (GeV/c^{2})","#it{p}_{T} (GeV/c)","y",
			"Dec Len (cm)","cos(#vartheta_{p})",
			"d_0^{min} (cm)",
			"d_0*d_0 (cm^{2})","DCA",
			"#it{p}_{T}^{min} (GeV/c)",
			"d_0^{D} (cm)","cos(#theta*)"};
  Int_t bins[nAxes] =   {100,   5,  40, 30,  20,   10,   10,      12,      8,  16,  10}; 
  Double_t min[nAxes] = {1.65,  0., 1., 0., 0.98, 0.,   -0.0006, 0.0,   0.,  0.,  -1.};
  Double_t max[nAxes] = {2.15,  5., 5., 0.3, 1.,   0.05, 0.,      0.03,  4.,  0.04, 1.};  
  THnSparseF *hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
  for(Int_t iax=0; iax<nAxes; iax++) hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  return hsp;
}

double GetTslope(int pdgParticle, double Eint){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_code){
    if(pdgParticle == pdg)
      index_pdg = counter;
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E)
      index_E = counter;
    counter++;
  }
  return Tslope[index_pdg][index_E];
}

double GetSigmaRapidity(int pdgParticle, double Eint){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_code){
    if(pdgParticle == pdg)
      index_pdg = counter;
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E)
      index_E = counter;
    counter++;
  }
  return sigma_rapidity[index_pdg][index_E];
}

double GetY0Rapidity(int pdgParticle, double Eint){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_code){
    if(pdgParticle == pdg)
      index_pdg = counter;
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E)
      index_E = counter;
    counter++;
  }
  return y0_rapidity[index_pdg][index_E];
}

double GetMultiplicity(int pdgParticle, double Eint){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_code){
    if(pdgParticle == pdg)
      index_pdg = counter;
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E)
      index_E = counter;
    counter++;
  }
  return multiplicity[index_pdg][index_E];
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