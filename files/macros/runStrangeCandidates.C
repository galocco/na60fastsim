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
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "MeasurementUtils.h"
#include "TClonesArray.h"
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;
double zTOF = 0.40;
double vX = 0, vY = 0, vZ = 0; // event vertex
TDatime dt;

void GenerateSignalCandidates(Int_t nevents = 10000, 
				double Eint = 40, 
        TString suffix = "_PHI_TEST",
				const char *setup = "../setups/setup_EHN1-H8_short_10pixel_1.5T_BB.txt",//
				const char *privateDecayTable = "../decaytables/USERTABPHI.DEC",
        Int_t pdgParticle = 333,
				bool simulateBg=kFALSE,
				bool writeNtuple = kTRUE,
        int minITSHits = 5){
  
  // Generate strange particle signals and simulate detector response for decay tracks
  bool matter = pdgParticle > 0;
  pdgParticle = TMath::Abs(pdgParticle);
  int pdg_unstable_dau = 0;
  int nbody = GetNBody(pdgParticle, pdg_unstable_dau);
  int refreshBg = 70;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);
  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");
  
  //define the pT and rapidity probability function
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();
  Double_t lifetime = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Lifetime();
  std::cout<<"lifetime: "<<lifetime<<"\n";
  if(pdgParticle == 0)
    lifetime = TMath::Power(8.59, -11);
  else if(pdgParticle == 3334)
    lifetime = TMath::Power(8.21, -11);
  else if(pdgParticle == 3312)
    lifetime = TMath::Power(1.64, -10);
  TF1 *fpt = new TF1("fpt", "x*exp(-TMath::Sqrt(x**2+[0]**2)/[1])", ptminSG, ptmaxSG);
  fpt->SetParameter(0,mass);
  fpt->SetParameter(1, GetTslope(pdgParticle,Eint,matter)/1000);
  TF1 *fy = new TF1("fy"," exp(-0.5*((x-[0]-[2])/[1])**2)+exp(-0.5*((x+[0]-[2])/[1])**2)", yminSG, ymaxSG);
  fy->SetParameter(0, GetY0Rapidity(pdgParticle, Eint, matter));
  fy->SetParameter(1, GetSigmaRapidity(pdgParticle, Eint, matter));
  fy->SetParameter(2, GetY0(Eint));
  TF1 *fm = new TF1("fm","1/((x-[0])**2+([1]/2)**2)",mass*0.975,mass*1.025);
  fm->SetParameter(0, mass);
  fm->SetParameter(1, TMath::Hbar()*6.242*TMath::Power(10,9)/lifetime);
  TFile *fout = new TFile(Form("Signal_histos%s.root",suffix.Data()), "recreate");
  TH1D* hCentrality = new TH1D("hCentrality", ";centrality;counts", 100, 0, 100);

  int count_solve=0, count_mcwin=0, count_chi2=0;

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT
  det->SetMinITSHits(minITSHits); //NA60+
  //det->SetMinITSHits(det->GetNumberOfActiveLayersITS()-1); //NA60
  det->SetMinMSHits(0); //NA60+
  //det->SetMinMSHits(det->GetNumberOfActiveLayersMS()-1); //NA60
  det->SetMinTRHits(0);
  //
  // max number of seed on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000);
  //
  // set chi2 cuts
  det->SetMaxChi2Cl(10.);  // max track to cluster chi2
  det->SetMaxChi2NDF(3.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(-20);  // fiducial cut on chi2 of convergence to vtx  
  // IMPORTANT FOR NON-UNIFORM FIEL
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); //NA60+
  //det->SetMinP2Propagate(2); //NA60
  //
  //det->SetApplyBransonPCorrection(-1);
  //det->SetApplyBransonPCorrection(3.); // kind of syst error on vertex position precision
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  //det->UseTrackOriginAsVertex();
  det->BookControlHistos();
  Double_t Bz = 0;
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
	      printf("Bx = %f By = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
        Bz = BVal[0]/10.;
      }else if (i == 1){
	      printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
      }
    }
  }

  double PrimVtxZ = det->GetLayer(0)->GetZ();
  // prepare decays
  TGenPhaseSpace decay;
  TLorentzVector parent, daugen[3], daurec[3], parentrefl, daurecswapmass[2];
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
  
  
  TH2F* hYPtGen = new TH2F("hYPtGen", "y-#it{p}_{T} corr match;y;#it{p}_{T};counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtGen  = new TH1D("hPtGen", "#it{p}_{Tgen};#it{p}_{T} (GeV/#it{c});counts", 40, ptminSG, ptmaxSG);
  TH1D* hPtFake = new TH1D("hPtFake", "#it{p}_{T};#it{p}_{T} (GeV/#it{c});counts", 40, ptminSG, ptmaxSG);
  TH1D* hPxFake = new TH1D("hPxFake", "#it{p}_{x};#it{p}_{T} (GeV/#it{c});counts", 40, ptminSG, ptmaxSG);
  TH1D* hPyFake = new TH1D("hPyFake", "#it{p}_{y};#it{p}_{T} (GeV/#it{c});counts", 40, ptminSG, ptmaxSG);
  TH1D* hPzFake = new TH1D("hPzFake", "#it{p}_{z};#it{p}_{T} (GeV/#it{c});counts", 40, ptminSG, ptmaxSG);
  TH1D* hYFake = new TH1D("hYFake", ";y;counts", 160., yminSG,ymaxSG);  
  TH1D* hYGen = new TH1D("hYGen", "y full phase space;y;counts", 160., yminSG,ymaxSG);  
  TH1D* hMGen = new TH1D("hMGen", "Mass all match;m (GeV/#it{c}^{2});counts", 1000, 0.9*mass, 1.1*mass);
  TH1D* hPtEff = new TH1D("hPtEff", "#it{p}_{T} efficiency;#it{p}_{T} (GeV/#it{c});counts", 40., ptminSG, ptmaxSG);
  TH1D* hYEff = new TH1D("hYEff", "y efficiency;counts", 160., yminSG,ymaxSG);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "y-#it{p}_{T} all match;y;#it{p}_{T};counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Reconstructed #it{p}_{T} all match;#it{p}_{T};counts", 40, ptminSG, ptmaxSG);  
  TH1D* hPtGenRecoAll = new TH1D("hPtGenRecoAll", "Generated #it{p}_{T} all match;#it{p}_{T};counts", 40, ptminSG, ptmaxSG);
  TH2F* hPtRecoVsGenAll = new TH2F("hPtRecoVsGenAll"," ; Generated #it{p}_{T} ; Reconstructed #it{p}_{T}",40, ptminSG, ptmaxSG,40, ptminSG, ptmaxSG);
  TH2F* hDiffPtRecoGenAll = new TH2F("hDiffPtRecoGenAll"," ; Generated #it{p}_{T} ; Reco #it{p}_{T} - Gen #it{p}_{T}",40, ptminSG, ptmaxSG,100,-0.2,0.2);
  TH2F* hMassVsOpen = new TH2F("hMassVsOpen", "Mass vs opening angle", 40, mass*0.96, mass*1.04, 50, 0., TMath::Pi());
  TH1F *hnSigmaTOF = new TH1F("hnSigmaTOF", ";n#sigma_{#beta};counts", 100, -5., 5.);
  TH2F *hTOF = new TH2F("hTOF", ";#it{p}(GeV/#it{c});TOF #beta;counts", 200, 0, 8, 300, 0.8, 1.15);
  TH2F* hptdau[3];
  TH1D* hydau[3];
  TH2F* hResPxDauVsPt[3];
  TH2F* hResPyDauVsPt[3];
  TH2F* hResPzDauVsPt[3];
  TH2F* hResPxDauVsY[3];
  TH2F* hResPyDauVsY[3];
  TH2F* hResPzDauVsY[3];
  TH1F* hPxDauFake[3];
  TH1F* hPyDauFake[3];
  TH1F* hPzDauFake[3];

  for(int i = 0; i < nbody; i++){
    hptdau[i] = new TH2F(Form("hptdau%i",i), " ;#it{p}_{T} (GeV/#it{c});#it{p}_{TM} (GeV/#it{c});counts", 50,ptminSG, 3,50,ptminSG, ptmaxSG);
    hydau[i]  = new TH1D(Form("hydau%i",i), ";y_{};counts", 160, yminSG, ymaxSG);
    hResPxDauVsPt[i] = new TH2F(Form("hResPxDauVsPt%i",i), ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -0.01, 0.01, 30, ptminSG, ptmaxSG);
    hResPyDauVsPt[i] = new TH2F(Form("hResPyDauVsPt%i",i), ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -0.01, 0.01, 30, ptminSG, ptmaxSG);
    hResPzDauVsPt[i] = new TH2F(Form("hResPzDauVsPt%i",i), ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, ptminSG, ptmaxSG);
    hResPxDauVsY[i] = new TH2F(Form("hResPxDauVsY%i",i), ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});y;counts", 100, -0.01, 0.01, 50, 0, 5);
    hResPyDauVsY[i] = new TH2F(Form("hResPyDauVsY%i",i), ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});y;counts", 100, -0.01, 0.01, 50, 0, 5);
    hResPzDauVsY[i] = new TH2F(Form("hResPzDauVsY%i",i), ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5);
    hPxDauFake[i] = new TH1F(Form("hPxDauFake%i",i), ";#it{p}_{x} (GeV/#it{c});y;counts", 100, 0, 10);
    hPyDauFake[i] = new TH1F(Form("hPyDauFake%i",i), ";#it{p}_{y} (GeV/#it{c});y;counts", 100, 0, 10);
    hPzDauFake[i] = new TH1F(Form("hPzDauFake%i",i), ";#it{p}_{z} (GeV/#it{c});y;counts", 100, 0, 10);
  }
  TH2F *hydau2D = new TH2F("hydau2D", "y negative daughter vs y positive daughter;y_{-};y_{+};counts", 50, 0., 5., 50, 0., 5.);
  
  TH1D* hYRecoAll = new TH1D("hYRecoAll", "Reconstructed y all match;y;counts", 80., 1., 5.4);
  TH1D* hYGenRecoAll = new TH1D("hYGenRecoAll", "Generated y all match;y;counts", 80., 1., 5.4);
  TH2F* hYPtRecoFake = new TH2F("hYPtRecoFake", "y-#it{p}_{T} fake match;;counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoFake = new TH1D("hPtRecoFake", "#it{p}_{T} fake match;;counts", 40, ptminSG, ptmaxSG);
  TH1D* hMassTrue = new TH1D("hMassTrue", "Mass all match;m (GeV/#it{c}^{2});counts", 200, 0.96*mass, 1.04*mass);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 2*mass);
  TH1D* hMassFake = new TH1D("hMassFake", "Mass fake match;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 2*mass);
  TH2F *hMassVsPt = new TH2F("hMassVsPt", ";m (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", 200, 0.9*mass, 1.1*mass, 50, ptminSG, ptmaxSG);
  TH2F *hMassVsY = new TH2F("hMassVsY", ";m (GeV/#it{c}^{2});y;counts", 200, 0.9*mass, 1.1*mass, 50, yminSG, ymaxSG);

  TH2F *hArmPod = new TH2F("hArmPod", ";#alpha;#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 100, 0, 2);
  TH2F *hDistXY = new TH2F("hDistXY", ";d_{xy} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", ";d (cm);#it{p}_{T} (GeV/#it{c});counts", 300, 0, 10, 30, 0, 3);
  TH1D *hDistx = new TH1D("hDistx", "generated secondary vertex x;x (cm);counts", 400, -30, 30);
  TH1D *hDisty = new TH1D("hDisty", "generated secondary vertex y;y (cm);counts", 400, -30, 30);
  TH1D *hDistxy = new TH1D("hDistxy", "generated secondary vertex d_{xy};d_{xy} (cm);counts", 100, 0, 30);
  TH1D *hDistz = new TH1D("hDistz", "generated secondary vertex z;z (cm);counts", 400, 0, 200);
  TH1D *hDistTot = new TH1D("hDistTot", "generated secondary vertex distance;d (cm);counts", 400, 0, 200);

  TH1D *hDistxRec = new TH1D("hDistxRec", "reconstructed secondary vertex x;x (cm);counts", 400, -30, 30);
  TH1D *hDistyRec = new TH1D("hDistyRec", "reconstructed secondary vertex y;y (cm);counts", 400, -30, 30);
  TH1D *hDistxyRec = new TH1D("hDistxyRec", "reconstructed secondary vertex d_{xy};d_{xy} (cm);counts", 20, 0, 1);
  TH1D *hDistzRec = new TH1D("hDistzRec", "reconstructed secondary vertex z;z (cm);counts", 400, 0, 200);
  TH1D *hDistTotRec = new TH1D("hDistTotRec", "reconstructed secondary vertex distance;d (cm);counts", 400, 0, 200);

  TH1D *hCt = new TH1D("hCt", ";#it{c}t (cm);counts", 100, 0, 100);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", ";d_{genxy} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", ";d_{gen} (cm);#it{p}_{T} (GeV/#it{c});counts", 300, 0, 10, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", ";cos(#theta_{p});#it{p}_{T} (GeV/#it{c});counts", 300, -1, 1, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", ";DCA (cm);#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", ";DCA_{x} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", ";DCA_{y} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", ";DCA_{z} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY[3];
  for(int i = 0; i < nbody; ++i)
    hd0XY[i] = new TH2F(Form("hd0xy%i",i), ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.01, 0.01, 30, 0, 3);

  //histograms for the mass of candidates built from misidentified, only if the daughters are stables
  TH1D* hMassRefl = new TH1D("hMassRefl", "Mass reflections;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 2*mass);
  TH2F *hMassReflVsPt = new TH2F("hMassReflVsPt", ";m (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", 200, 0.5*mass, 2*mass, 6, 0, 3);
  TH2F *hMassReflVsY = new TH2F("hMassReflVsY", ";m (GeV/#it{c}^{2});y;counts", 200, 0.5*mass, 2*mass, 10, 0, 5);
  
  TH2F *hResVxVsPt = new TH2F("hResVxVsPt", ";V_{xgen}-V_{xrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResVyVsPt = new TH2F("hResVyVsPt", ";V_{ygen}-V_{yrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResVzVsPt = new TH2F("hResVzVsPt", ";V_{zgen}-V_{zrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResPxVsPt = new TH2F("hResPxVsPt", ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResPyVsPt = new TH2F("hResPyVsPt", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResPzVsPt = new TH2F("hResPzVsPt", ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResPtVsPt = new TH2F("hResPtVsPt", ";#it{p}_{Tgen}-#it{p}_{Trec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -0.02, 0.02, 30, ptminSG, ptmaxSG);

  TH2F *hResVxVsY = new TH2F("hResVxVsY", ";V_{xgen}-V_{xrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", ";V_{ygen}-V_{yrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", ";V_{zgen}-V_{zrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResPxVsY = new TH2F("hResPxVsY", ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});y;counts", 100, -0.01, 0.01, 50, 0, 5);
  TH2F *hResPyVsY = new TH2F("hResPyVsY", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});y;counts", 100, -0.01, 0.01, 50, 0, 5);
  TH2F *hResPzVsY = new TH2F("hResPzVsY", ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5);
  TH2F *hResPtVsY = new TH2F("hResPtVsY", ";#it{p}_{Tgen}-#it{p}_{Trec} (GeV/#it{c});y;counts", 100, -0.02, 0.02, 50, 0, 5);
  
  TH2F *hResDist = new TH2F("hResDist", ";d_{gen}-d_{rec} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", ";d_{xygen}-d_{xyrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH1D *hNevents = new TH1D("hNevents", ";;counts", 1, 0, 1);
  TH1D *hChi2True = new TH1D("hChi2True", ";#chi^{2};counts", 15, 0, 1.5);
  TH1D *hChi2Fake = new TH1D("hChi2Fake", ";#chi^{2};counts", 15, 0, 1.5);

  TFile *tracks_file = new TFile(Form("treeStrangeParticles%s.root",suffix.Data()), "RECREATE");
  TTree *tree = new TTree("tree", "tree Sig");
  TClonesArray *arrtr = new TClonesArray("KMCProbeFwd");
  TClonesArray &aarrtr = *arrtr;
  tree->Branch("tracks", &arrtr);
  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  TNtuple *ntgen = 0x0;
  if (writeNtuple)
  {
    fnt = new TFile(Form("fntSig%s.root",suffix.Data()), "recreate");
    ntgen = new TNtuple("ntgen", "ntgen", "pt:rapidity:dist", 32000); 
    if(nbody == 2)
      ntcand = new TNtuple("ntcand", "ntcand", "m:pt:rapidity:dist:cosp:d0prod:ptMin:ptMax:dca:thetad:arm:nsigma1:nsigma2", 32000);
    else
      ntcand = new TNtuple("ntcand", "ntcand", "m:pt:rapidity:dist:distD:cosp:cospD:bxy:bxyD:dca:dcaD", 32000);
  }

  Float_t* arrnt = new Float_t[(nbody == 2) ? 13 : 11];
  Float_t arrntgen[4];
  //read the pdg code of the daughters
  Int_t pdg_dau[2] = {0,0};
  GetPDGDaughters(pdgParticle, pdg_dau,matter);
  //if decay daughter is unstable read the pdg code of its daughters
  Int_t pdg_dau2[2] = {0, 0};
  Int_t pdgDaughter[3] = {0, 0, 0};
  Double_t massDau[3] = {0, 0, 0};
  Double_t massUnstable = 0;
  if(nbody == 3){
    GetPDGDaughters(pdg_unstable_dau,pdg_dau2, pdg_unstable_dau > 0);
    massUnstable = TDatabasePDG::Instance()->GetParticle(pdg_unstable_dau)->Mass();
    pdgDaughter[0] = (pdg_dau[0] != pdg_unstable_dau) ? pdg_dau[0] : pdg_dau[1];
    pdgDaughter[1] = (pdgDaughter[0]*pdg_dau2[0] < 0) ? pdg_dau2[0] : pdg_dau2[1];
    pdgDaughter[2] = (pdgDaughter[0]*pdg_dau2[0] < 0) ? pdg_dau2[1] : pdg_dau2[0];
  }
  
  if(nbody==2){
    massDau[0] = TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass();//negative
    massDau[1] = TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass();//positive
  }
  else{
    massDau[0] = TDatabasePDG::Instance()->GetParticle(pdgDaughter[0])->Mass();
    massDau[1] = TDatabasePDG::Instance()->GetParticle(pdgDaughter[1])->Mass();
    massDau[2] = TDatabasePDG::Instance()->GetParticle(pdgDaughter[2])->Mass();
  }

  Double_t multi_times_br = GetMultiplicity(pdgParticle, Eint, matter)*GetBRatio(pdgParticle);
  Int_t gen_part = gRandom->Poisson(multi_times_br);
  Int_t icount = 0;

  Double_t L = 0;
  Double_t pxyz[3] = {0,0,0};
  Double_t tHit = 0;
  Double_t nsigma1 = 0, nsigma2 = 0;
  Double_t T0 = smearT(0);
  Int_t iev = 0;
  Int_t istrange = 0; 
  while (SetEvent(nevents, iev, tree, aarrtr, T0, gen_part, multi_times_br, icount)){
    gen_part--;
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    int nrec = 0;
    int nfake = 0;
    bool IsReconstructable = true;
    if (simulateBg && iev%refreshBg==0) det->GenBgEvent(0. ,0. , PrimVtxZ);

    Double_t ptGen = fpt->GetRandom(ptminSG, ptmaxSG); 
    Double_t yGen = fy->GetRandom(yminSG, ymaxSG);
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t pxGen = ptGen * TMath::Cos(phi);
    Double_t pyGen = ptGen * TMath::Sin(phi);
    Double_t massbw = (pdgParticle == 333) ? fm->GetRandom(mass * 0.96, mass * 1.04) : mass;
    Double_t mt = TMath::Sqrt(ptGen * ptGen + massbw * massbw);
    Double_t pzGen = mt * TMath::SinH(yGen);
    Double_t en = mt * TMath::CosH(yGen);
    mom->SetPxPyPzE(pxGen, pyGen, pzGen, en);
    Int_t np;
    do{
      fDecayer->Decay(pdgParticle, mom);
      np = fDecayer->ImportParticles(particles);
    } while (np < 0);
    istrange++;
    Double_t secvertgen[3] = {0, 0, 0};
    Double_t trdvertgen[3] = {0, 0, 0};
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
      else if (kf == pdg_dau2[1] && nbody==3){
        trdvertgen[0] = vX;
        trdvertgen[1] = vY;
        trdvertgen[2] = vZ;
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
          daurecswapmass[crg_index].SetXYZM(iparticle->Px(), iparticle->Py(), iparticle->Pz(),(iparticle->GetMass() == massDau[1])? massDau[0] : massDau[1]);
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

        if(IsReconstructable && nrec==0){
          if (ntgen){
            arrntgen[0] = ptGen;
            arrntgen[1] = yGen;
            arrntgen[2] = TMath::Sqrt(vX*vX+vY*vY+vZ*vZ);
            ntgen->Fill(arrntgen);
          }

          secvertgen[0] = vX;
          secvertgen[1] = vY;
          secvertgen[2] = vZ;
          hDistx->Fill(vX);
          hDisty->Fill(vY);
          hDistz->Fill(vZ);
          hDistxy->Fill(TMath::Sqrt(vX*vX+vY*vY));
          hDistTot->Fill(TMath::Sqrt(vX*vX+vY*vY+vZ*vZ));
          hCt->Fill(TMath::Sqrt(vX*vX+vY*vY+vZ*vZ)*mass/TMath::Sqrt(ptGen*ptGen+pzGen*pzGen));
        }

        daugen[nrec].SetXYZM(iparticle->Px(), iparticle->Py(), iparticle->Pz(), iparticle->GetMass());
        
        if (!det->SolveSingleTrack(daugen[nrec].Pt(), daugen[nrec].Rapidity(), daugen[nrec].Phi(), iparticle->GetMass(), crg, vX, vY, vZ, 0, 1, 99)){
          if(IsReconstructable){ 
            IsReconstructable = false;
            count_solve++;
            }
          continue;
        }

        KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
        if (!trw){
          if(IsReconstructable){ 
            IsReconstructable = false;
            count_mcwin++;
          }
          continue;
        }
        
        if (trw->GetNormChi2(kTRUE) > ChiTot){
          if(IsReconstructable){ 
            IsReconstructable = false;
            count_chi2++;
          }
          continue;
        }

        trw->SetIndexMom(istrange);
        trw->SetIndex(icount);
        trw->SetPdg(kf);
        trw->SetPdgMother(pdgParticle);

        double test[3] = {iparticle->Px(), iparticle->Py(), iparticle->Pz()};
        tHit = smearT(hasTOF(pxyz, iparticle->GetMass(), L, zTOF));
        trw->SetTOF(tHit-T0);
        trw->SetL(L);

        new (aarrtr[icount]) KMCProbeFwd(*trw);
        icount++;
        nfake += trw->GetNFakeITSHits();
        if(trw->GetNFakeITSHits()==0)
          hChi2True->Fill(trw->GetNormChi2(kTRUE));
        else
          hChi2Fake->Fill(trw->GetNormChi2(kTRUE));

        if(nbody==2){
          recProbe[crg_index] = *trw;
          double betaMeas = betaTOF(L, tHit-T0);
          hTOF->Fill(TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]),betaMeas);
          double beta = betaGen(pxyz, iparticle->GetMass());
          if(crg_index == 0){
            nsigma1 = (betaMeas-beta)/sigmaBeta(zTOF);
            hnSigmaTOF->Fill(nsigma1);
          }
          else{
            nsigma2 = (betaMeas-beta)/sigmaBeta(zTOF);
            hnSigmaTOF->Fill(nsigma2);
          }
        }
        else{
          if(i == 2){
            recProbe[0] = *trw;
          }
          else if(i == 3){
            recProbe[1] = *trw;
          }
          else if(i == 4){
            recProbe[2] = *trw;
          }
        }
        nrec++;

      }
    }
    if(nrec < nbody)
      continue;//not alle the daughters are reconstructed

    hCentrality->Fill(5*gRandom->Rndm());
    
    double pos[3] = {0,0,0};
    KMCProbeFwd helper, fin_cand;
    if(nbody == 2){
      recProbe[0].PropagateToDCA(&recProbe[1]);  
      recProbe[0].GetPXYZ(pxyz);
      daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], recProbe[0].GetMass());
      recProbe[1].GetPXYZ(pxyz);
      daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2], recProbe[1].GetMass());
      parent = daurec[0];
      parent += daurec[1];
    }
    else{
      recProbe[1].PropagateToDCA(&recProbe[2]);
      ComputeVertex(recProbe[1],recProbe[2],pos[0],pos[1],pos[2]);
      recProbe[1].GetPXYZ(pxyz);
      daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], recProbe[1].GetMass());
      recProbe[2].GetPXYZ(pxyz);
      daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2], recProbe[2].GetMass());
      daurec[2] = daurec[0];
      daurec[2] += daurec[1];
  
      pxyz[0] = daurec[2].Px();
      pxyz[1] = daurec[2].Py();
      pxyz[2] = daurec[2].Pz();

      helper = KMCProbeFwd(pos, pxyz, recProbe[1].GetCharge()+recProbe[2].GetCharge());
      recProbe[0].PropagateToDCA(&helper);

      recProbe[0].GetPXYZ(pxyz);
      daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], recProbe[0].GetMass());
      helper.GetPXYZ(pxyz);
      daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2], massUnstable);

      parent = daurec[0];
      parent += daurec[1];
    }
 
    hMGen->Fill(massbw);

    Double_t  ptRec = parent.Pt();
    Double_t  massRec = parent.M();
    Double_t yRec = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));

    Float_t dca = 0, dcaD = 0;
    Double_t xP = 0, yP = 0, zP = 0, sigmaVert = 0, arm = 0;
    Double_t thetad = OpeningAngle(daurec[0],daurec[1]);
    hMassVsOpen->Fill(massRec, TMath::ACos(thetad));
    if(nbody == 2){
      hydau2D->Fill(daurec[0].Rapidity(), daurec[1].Rapidity());
      parentrefl = daurecswapmass[0];
      parentrefl += daurecswapmass[1];
      Double_t  massRecReflD=parentrefl.M();
      hMassRefl->Fill(massRecReflD);
      hMassReflVsPt->Fill(massRecReflD,ptRec);
      hMassReflVsY->Fill(massRecReflD,yRec);

      Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
      Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
      Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
      dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
      hDCA->Fill(dca, ptRec);
      hDCAx->Fill(d1, ptRec);
      hDCAy->Fill(d2, ptRec);
      hDCAz->Fill(d3, ptRec);
      ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
      arm = ArmenterosAlpha(daurec[1],daurec[0]);
      hArmPod->Fill(arm, ptRec);
    }
    else{
      ComputeVertex(recProbe[0], helper, xP, yP, zP);
      Float_t d1 = recProbe[0].GetX() - helper.GetX();
      Float_t d2 = recProbe[0].GetY() - helper.GetY();
      Float_t d3 = recProbe[0].GetZ() - helper.GetZ();
      dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
      d1 = recProbe[2].GetX() - recProbe[1].GetX();
      d2 = recProbe[2].GetY() - recProbe[1].GetY();
      d3 = recProbe[2].GetZ() - recProbe[1].GetZ();
      dcaD = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    }

    hYPtRecoAll->Fill(yRec, ptRec);
    hPtRecoAll->Fill(ptRec);
    hPtGenRecoAll->Fill(ptGen);
    hPtRecoVsGenAll->Fill(ptGen, ptRec);
    hDiffPtRecoGenAll->Fill(ptGen, ptRec-ptGen);
    hYRecoAll->Fill(yRec);
    hYGenRecoAll->Fill(yGen);
    hMassAll->Fill(massRec);

    if (nfake > 0){
      hYPtRecoFake->Fill(yRec, ptRec);
      hPtRecoFake->Fill(ptRec);
      hMassFake->Fill(massRec);
    }
    else{
      hMassTrue->Fill(massRec);
      hYEff->Fill(yRec);
      hPtEff->Fill(ptRec);
      hYFake->Fill(yRec);
      hPtFake->Fill(ptRec);
      hPxFake->Fill(parent.Px());
      hPyFake->Fill(parent.Py());
      hPzFake->Fill(parent.Pz());

      for(int i = 0; i < nbody; ++i){
        hPxDauFake[i]->Fill(daurec[i].Px());
        hPyDauFake[i]->Fill(daurec[i].Py());
        hPzDauFake[i]->Fill(daurec[i].Pz());
      }

    }

    hMassVsPt->Fill(massRec,ptRec);
    hMassVsY->Fill(massRec,yRec);
    
    Double_t residVx = 10000.*(xP - secvertgen[0]);
    Double_t residVy = 10000.*(yP - secvertgen[1]);
    Double_t residVz = 10000.*(zP - secvertgen[2]);

    hResVxVsPt->Fill(residVx, ptRec);
    hResVyVsPt->Fill(residVy, ptRec);
    hResVzVsPt->Fill(residVz, ptRec);

    hResPtVsPt->Fill(ptRec-ptGen,ptRec);
    hResPtVsY->Fill(ptRec-ptGen,yRec);
    
    hResVxVsY->Fill(residVx, yRec);
    hResVyVsY->Fill(residVy, yRec);
    hResVzVsY->Fill(residVz, yRec);
    
    hResPxVsPt->Fill(parent.Px() - pxGen, ptRec);
    hResPyVsPt->Fill(parent.Py() - pyGen, ptRec);
    hResPzVsPt->Fill(parent.Pz() - pzGen, ptRec);

    hResPxVsY->Fill(parent.Px() - pxGen, yRec);
    hResPyVsY->Fill(parent.Py() - pyGen, yRec);
    hResPzVsY->Fill(parent.Pz() - pzGen, yRec);

    for(int i = 0; i < nbody; ++i){
      hResPxDauVsPt[i]->Fill(daurec[i].Px()-daugen[i].Px(), ptRec);
      hResPyDauVsPt[i]->Fill(daurec[i].Py()-daugen[i].Py(), ptRec);
      hResPzDauVsPt[i]->Fill(daurec[i].Pz()-daugen[i].Pz(), ptRec);
      hResPxDauVsY[i]->Fill(daurec[i].Px()-daugen[i].Px(), yRec);
      hResPyDauVsY[i]->Fill(daurec[i].Py()-daugen[i].Py(), yRec);
      hResPzDauVsY[i]->Fill(daurec[i].Pz()-daugen[i].Pz(), yRec);
    }
    
    Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
    Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
    Float_t distgen = TMath::Sqrt(secvertgen[0] * secvertgen[0] + secvertgen[1] * secvertgen[1] + secvertgen[2] * secvertgen[2]);
    Float_t distgenXY = TMath::Sqrt(secvertgen[0] * secvertgen[0] + secvertgen[1] * secvertgen[1]);

    Float_t bxy = 0, bxyD = 0;
    Float_t distD = 0;
    Double_t vsec[3] = {xP, yP, zP};
    if(nbody==3){
      pxyz[0] = parent.Px();
      pxyz[1] = parent.Py();
      pxyz[2] = parent.Pz();
      fin_cand = KMCProbeFwd(pos, pxyz, recProbe[0].GetCharge()+helper.GetCharge());
      distD = TMath::Sqrt((pos[0] - vsec[0]) * (pos[0] - vsec[0]) + (pos[1] - vsec[1]) * (pos[1] - vsec[1]) + (pos[2] - vsec[2]) * (pos[2] - vsec[2]));

      fin_cand.PropagateToZBxByBz(0);
      Float_t d0x = fin_cand.GetX();
      Float_t d0y = fin_cand.GetY();
      bxy = TMath::Sqrt(d0x * d0x + d0y * d0y);
      if (d0x < 0)
        bxy *= -1;

      helper.PropagateToZBxByBz(0);
      d0x = helper.GetX();
      d0y = helper.GetY();
      bxyD = TMath::Sqrt(d0x * d0x + d0y * d0y);
      if (d0x < 0)
        bxyD *= -1;

    }

    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
    Double_t cospD = CosPointingAngle(vsec, pos, daurec[2]);
    Double_t ipD = ImpParXY(vprim, vsec, parent);
    hCosp->Fill(cosp, ptRec);
    
    hResDist->Fill(dist - distgen, ptRec);
    hResDistXY->Fill(distXY - distgenXY, ptRec);
    
    hDistXY->Fill(distXY, ptRec);
    hDist->Fill(dist, ptRec);
    hDistxRec->Fill(xP);
    hDistyRec->Fill(yP);
    hDistzRec->Fill(zP);
    hDistxyRec->Fill(distXY);
    hDistTotRec->Fill(dist);
    hDistgenXY->Fill(distgenXY, ptRec);
    hDistgen->Fill(distgen, ptRec);

    Double_t d0xy[3] = {-1,-1,-1};
    for(int i = 0; i < nbody; i++){
      recProbe[i].PropagateToZBxByBz(0);
      Double_t d0x = recProbe[i].GetX();
      Double_t d0y = recProbe[i].GetY();
      d0xy[i] = TMath::Sqrt(d0x * d0x + d0y * d0y);
      if (d0x < 0)
        d0xy[i] *= -1;
      hd0XY[i]->Fill(d0xy[i], ptRec);
    }
    
    hd0XYprod->Fill(d0xy[0] * d0xy[1], ptRec);

    if (writeNtuple){
      if(nbody == 2){
        arrnt[0] = massRec;
        arrnt[1] = ptRec;
        arrnt[2] = yRec;
        arrnt[3] = dist;
        arrnt[4] = cosp;
        arrnt[5] = d0xy[0] * d0xy[1];
        arrnt[6] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
        arrnt[7] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
        arrnt[8] = dca; 
        arrnt[9] = TMath::ACos(thetad); 
        arrnt[10] = arm; 
        arrnt[11] = nsigma1; 
        arrnt[12] = nsigma2; 
      }
      else{
        arrnt[0] = massRec;
        arrnt[1] = ptRec;
        arrnt[2] = yRec;
        arrnt[3] = dist;
        arrnt[4] = distD;
        arrnt[5] = cosp;
        arrnt[6] = cospD;
        arrnt[7] = bxy;
        arrnt[8] = bxyD;
        arrnt[9] = dca;
        arrnt[10] = dcaD;
      }
      ntcand->Fill(arrnt);
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

  //efficiency vs pT /rapidity
  for(int i = 1; i <= hPtGen->GetNbinsX(); i++){
    double n = hPtGen->GetBinContent(i);
    if(n<0)
      hPtGen->SetBinContent(i,0);
    double rec = hPtEff->GetBinContent(i);
    if(rec>n)
      hPtEff->SetBinContent(i,n);
      
  }
  
  for(int i = 1; i <= hYGen->GetNbinsX(); i++){
    double n = hYGen->GetBinContent(i);
    if(n<0)
      hYGen->SetBinContent(i,0);
    double rec = hYEff->GetBinContent(i);
    if(rec>n)
      hYEff->SetBinContent(i,n);
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
  hTOF->Write();
  hnSigmaTOF->Write();

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
  if (writeNtuple){
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

  tracks_file->cd();
  tree->Write();
  tracks_file->Close();
  std::cout<<"rejected solve:"<<count_solve<<std::endl;
  std::cout<<"rejected mcwinner:"<<count_mcwin<<std::endl;
  std::cout<<"rejected chi2:"<<count_chi2<<std::endl;
  std::cout<<"tot rejected: "<<count_solve+count_mcwin+count_chi2<<std::endl;

}


void MakeCombinBkgCandidates3Body(const char* trackTreeFileBkg="treeBkgEvents_layer5.root",
            const char* trackTreeFileSig="treeBkgEvents_layer5.root",
            TString suffix = "_Omega",
            const char *setup = "../setups/setup_EHN1-H8_short_10pixel_1.5T_BB.txt",
            Int_t pdgMother = 3334,
            Int_t fullsim = kTRUE,
			      Int_t writeNtuple = kTRUE){

  Int_t pdg_unstable_dau;
  Int_t chMother = TDatabasePDG::Instance()->GetParticle(pdgMother)->Charge()/3;//->Charge() gives the charge in units of |e|/3
  pdgMother = TMath::Abs(pdgMother);
  GetNBody(pdgMother, pdg_unstable_dau);
  // Read the TTree of tracks produced with runBkgVT.C
  // Create combinatorial background candidates (= OS pairs of tracks)
  // background tracks
  TFile *filetreeBkg = new TFile(trackTreeFileBkg);
  TTree *treeBkg = (TTree *)filetreeBkg->Get("tree");
  TClonesArray *arrBkg = 0;
  treeBkg->SetBranchAddress("tracks", &arrBkg);

  Int_t nevents = treeBkg->GetEntries();
  // signal tracks
  TFile *filetreeSig = new TFile(trackTreeFileSig);
  TTree *treeSig = (TTree *)filetreeSig->Get("tree");
  TClonesArray *arrSig = 0;
  treeSig->SetBranchAddress("tracks", &arrSig);
  if(treeSig->GetEntries()<nevents && fullsim)
    nevents = treeSig->GetEntries();
  printf("Number of events in tree = %d\n",nevents);
  
  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  
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


  TFile *fout = new TFile(Form("Bkg-histos%s.root",suffix.Data()), "recreate");
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);

  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 1., 3.5);
  TH1D* hMassKK = new TH1D("hMassKK", "KK Inv Mass ", 200, 0.9, 2.9);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY3 = new TH2F("hd0xy3", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 100, -1, 1, 30, 0, 3);
  
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
  
  // define mother particle mass
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgMother)->Mass();

  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  Float_t arrnt[12];
  if (writeNtuple){
    fnt = new TFile(Form("fntBkg%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "m:pt:rapidity:dist:distD:cosp:cospD:bxy:bxyD:dca:dcaD:true", 32000);
                                             
  }

  //read the pdg code of the daughters
  int pdg_dau[2] = {0,0};
  GetPDGDaughters(pdgMother, pdg_dau, pdgMother>0);
  //read the pdg code of the unstable daughter daughters
  int pdg_dau2[2] = {0,0};
  GetPDGDaughters(pdg_unstable_dau, pdg_dau2, pdg_unstable_dau > 0);
  int pdgDaughter[3];
  pdgDaughter[0] = (pdg_dau[0] != pdg_unstable_dau) ? pdg_dau[0] : pdg_dau[1];
  pdgDaughter[1] = (pdgDaughter[0]*pdg_dau2[0] < 0) ? pdg_dau2[0] : pdg_dau2[1];
  pdgDaughter[2] = (pdgDaughter[0]*pdg_dau2[0] < 0) ? pdg_dau2[1] : pdg_dau2[0];
  
  Double_t massUnstable = TDatabasePDG::Instance()->GetParticle(pdg_unstable_dau)->Mass();
  Double_t massDau[3];
  massDau[0] = TDatabasePDG::Instance()->GetParticle(pdgDaughter[0])->Mass();
  massDau[1] = TDatabasePDG::Instance()->GetParticle(pdgDaughter[1])->Mass();
  massDau[2] = TDatabasePDG::Instance()->GetParticle(pdgDaughter[2])->Mass();

  KMCProbeFwd recProbe[3];
  TLorentzVector parent, pair, daurec[3];
  Int_t trueCand = 0;
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    treeBkg->GetEvent(iev);
    Int_t arrentrBkg = arrBkg->GetEntriesFast();
    if(fullsim)
      treeSig->GetEvent(iev);
    Int_t arrentrSig = (fullsim) ? arrSig->GetEntriesFast() : 0;
    Double_t pxyz[3];
    for (Int_t itr = 0; itr < arrentrBkg+arrentrSig; itr++){
      KMCProbeFwd *tr1;
      if(itr<arrentrBkg)
        tr1 = (KMCProbeFwd *)arrBkg->At(itr);
      else
        tr1 = (KMCProbeFwd *)arrSig->At(itr-arrentrBkg);
      Float_t ch1 = tr1->GetCharge();

      for (Int_t itr2 = itr+1; itr2 < arrentrBkg+arrentrSig; itr2++){
        KMCProbeFwd *tr2;
        if(itr2<arrentrBkg)
          tr2 = (KMCProbeFwd *)arrBkg->At(itr2);
        else
          tr2 = (KMCProbeFwd *)arrSig->At(itr2-arrentrBkg);
        Float_t ch2 = tr2->GetCharge();
        // convention: charge signs are ordered as +-+ or -+-
        if (ch1 * ch2 > 0) continue;

        for (Int_t itr3 = itr2+1; itr3 < arrentrBkg+arrentrSig; itr3++){
          KMCProbeFwd *tr3;
          if(itr3<arrentrBkg)
            tr3 = (KMCProbeFwd *)arrBkg->At(itr3);
          else
            tr3 = (KMCProbeFwd *)arrSig->At(itr3-arrentrBkg);
          Float_t ch3 = tr3->GetCharge();
          // convention: charge signs are ordered as +-+ or -+-
          if (ch3 * ch2 > 0) continue;
          if ((ch1 + ch2 + ch3) != chMother) continue;

          ////////////////////////////////
          recProbe[0] = *tr1;
          recProbe[0].PropagateToZBxByBz(vprim[2]);
          Double_t d0x1 = recProbe[0].GetX();
          Double_t d0y1 = recProbe[0].GetY();
          Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
          if (d0x1 < 0) d0xy1 *= -1;
          ////////////////////////////////
          recProbe[1] = *tr2;
          recProbe[1].PropagateToZBxByBz(vprim[2]);
          Double_t d0x2 = recProbe[1].GetX();
          Double_t d0y2 = recProbe[1].GetY();
          Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
          if (d0x2 < 0) d0xy2 *= -1;
          ////////////////////////////////
          recProbe[2] = *tr3;
          recProbe[2].PropagateToZBxByBz(vprim[2]);
          Double_t d0x3 = recProbe[2].GetX();
          Double_t d0y3 = recProbe[2].GetY();
          Double_t d0xy3 = TMath::Sqrt(d0x3 * d0x3 + d0y3 * d0y3);
          if (d0x3 < 0) d0xy3 *= -1;
          ////////////////////////////////
          for(Int_t iMassHyp=0; iMassHyp<2; iMassHyp++){
            countCand++;
            Double_t VSec[3] = {0, 0, 0};
            Double_t VTrd[3] = {0, 0, 0};
            recProbe[1].PropagateToDCA(&recProbe[2*iMassHyp]);
            ComputeVertex(recProbe[1], recProbe[2*iMassHyp], VTrd[0], VTrd[1], VTrd[2]);
            recProbe[1].GetPXYZ(pxyz);
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], massDau[1]);
            recProbe[2*iMassHyp].GetPXYZ(pxyz);
            daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2], massDau[2]);
            pair = daurec[1];
            pair += daurec[0];

            Double_t  massPair = pair.M();
            //skip if the mass of the candidate unstable daughter is far from the expected mass
            if(massPair < massUnstable*0.97 || massPair > massUnstable*1.03)
              continue;

            pxyz[0] = pair.Px();
            pxyz[1] = pair.Py();
            pxyz[2] = pair.Pz();

            KMCProbeFwd helper(VTrd, pxyz,  recProbe[2*iMassHyp].GetCharge()+recProbe[1].GetCharge());
            recProbe[2*(1-iMassHyp)].PropagateToDCA(&helper);
            ComputeVertex(recProbe[2*(1-iMassHyp)], helper, VSec[0], VSec[1], VSec[2]);
            recProbe[2*(1-iMassHyp)].GetPXYZ(pxyz);
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], massDau[0]);
            helper.GetPXYZ(pxyz);
            daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2], massUnstable);

            
            Float_t distPair = TMath::Sqrt(VTrd[0] * VTrd[0] + VTrd[1] * VTrd[1] + VTrd[2] * VTrd[2]);
            Float_t dist = TMath::Sqrt(VSec[0] * VSec[0] + VSec[1] * VSec[1] + VSec[2] * VSec[2]);
            //skip if the tertiary decay occurs before the secondary decay
            if(dist < 1.)
              continue;
            if(dist > distPair)
              continue;

            Float_t distD = TMath::Sqrt((VSec[0]-VTrd[0]) * (VSec[0]-VTrd[0]) + (VSec[1]-VTrd[1]) * (VSec[1]-VTrd[1]) + (VSec[2]-VTrd[2]) * (VSec[2]-VTrd[2]));

            parent = daurec[0];
            parent += daurec[1];

            Float_t invMass = parent.M();
            if(invMass>mass*0.96  && invMass<mass*1.04){
              hMassKK->Fill(massPair);
              hMassAll->Fill(invMass);
              Float_t pt = parent.Pt();
              Float_t y = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
              hYPtRecoAll->Fill(y, pt);
              hPtRecoAll->Fill(pt);

              Double_t sigmaVert = 0;
              Float_t distXY = TMath::Sqrt(VSec[0] * VSec[0] + VSec[1] * VSec[1]);
              Double_t cosp = CosPointingAngle(vprim, VSec, parent);
              Double_t cospD = CosPointingAngle(VSec, VTrd, pair);
              //skip if the cosine of the pointing angle is too small
              if(cosp < 0.999)
                continue; 

              countCandInPeak++;
              Float_t dca = 0, dcaD = 0;
              Float_t d1 = recProbe[2*(1-iMassHyp)].GetX() - helper.GetX();
              Float_t d2 = recProbe[2*(1-iMassHyp)].GetY() - helper.GetY();
              Float_t d3 = recProbe[2*(1-iMassHyp)].GetZ() - helper.GetZ();
              dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
              d1 = recProbe[2*iMassHyp].GetX() - recProbe[1].GetX();
              d2 = recProbe[2*iMassHyp].GetY() - recProbe[1].GetY();
              d3 = recProbe[2*iMassHyp].GetZ() - recProbe[1].GetZ();
              dcaD = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
              Float_t bxy = 0, bxyD = 0;

              pxyz[0] = parent.Px();
              pxyz[1] = parent.Py();
              pxyz[2] = parent.Pz();
              KMCProbeFwd fin_cand(VSec, pxyz, recProbe[2*(1-iMassHyp)].GetCharge()+helper.GetCharge());
              fin_cand.PropagateToZBxByBz(0);
              Float_t d0x = fin_cand.GetX();
              Float_t d0y = fin_cand.GetY();
              bxy = TMath::Sqrt(d0x * d0x + d0y * d0y);
              if (d0x < 0)
                bxy *= -1;

              helper.PropagateToZBxByBz(0);
              d0x = helper.GetX();
              d0y = helper.GetY();
              bxyD = TMath::Sqrt(d0x * d0x + d0y * d0y);
              if (d0x < 0)
                bxyD *= -1;


              hCosp->Fill(cosp, pt);
              hDistXY->Fill(distXY, pt);
              hDist->Fill(dist, pt);
              hd0XY1->Fill(d0xy1, pt);
              hd0XY2->Fill(d0xy2, pt);
              hd0XY3->Fill(d0xy3, pt);

              trueCand = IsTrueCandidate(tr1, tr2, tr3, pdgMother, pdgDaughter, iMassHyp);
              if (writeNtuple){
                arrnt[0] = invMass;
                arrnt[1] = pt;
                arrnt[2] = y;
                arrnt[3] = dist;
                arrnt[4] = distD;
                arrnt[5] = cosp;
                arrnt[6] = cospD;
                arrnt[7] = bxy;
                arrnt[8] = bxyD;
                arrnt[9] = dca;
                arrnt[10] = dcaD;
                arrnt[11] = trueCand;
                
                ntcand->Fill(arrnt);
              }
            } // check on inv mass
          } // loop on mass hypothesis
	      } // loop on third track
      } // loop on second track
    }// loop on first track
    printf(" --> Event %d, tot candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hNevents->Write();
  hMassAll->Write();
  hMassKK->Write();
  hYPtRecoAll->Write();
  hPtRecoAll->Write();
  hDistXY->Write();
  hDist->Write();
  hDCA->Write();
  hCosp->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hd0XY3->Write();
  fout->Close();
  if (writeNtuple){
    fnt->cd();
    ntcand->Write();
    fnt->Close();
  }
}

void MakeCombinBkgCandidates2Body(const char* trackTreeFileBkg="treeBkgEvents_layer5.root",
            const char* trackTreeFileSig="treeBkgEvents_layer5.root",
            TString suffix = "_PHI_test",
            int pdgMother = 333,
            const char *setup = "../setups/setup-EHN1_BetheBloch.txt",
            Int_t fullsim = kTRUE,
			      Int_t writeNtuple = kTRUE,
            Int_t writeSparse = kTRUE){

  // Read the TTree of tracks produced with runBkgVT.C
  // Store in THnSparse and (optionally) TNtuple
  // background tracks
  TFile *filetreeBkg = new TFile(trackTreeFileBkg);
  TTree *treeBkg = (TTree *)filetreeBkg->Get("tree");
  TClonesArray *arrBkg = 0;
  treeBkg->SetBranchAddress("tracks", &arrBkg);
  Int_t nevents = treeBkg->GetEntries();
  // signal tracks
  TFile *filetreeSig = new TFile(trackTreeFileSig);
  TTree *treeSig = (TTree *)filetreeSig->Get("tree");
  TClonesArray *arrSig = 0;
  treeSig->SetBranchAddress("tracks", &arrSig);
  if(treeSig->GetEntries()<nevents && fullsim)
    nevents = treeSig->GetEntries();

  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  Double_t massM = TDatabasePDG::Instance()->GetParticle(pdgMother)->Mass();

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  
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
  
  TFile *fout = new TFile(Form("Bkg-histos%s.root",suffix.Data()), "recreate");
  TH1D* hCentrality = new TH1D("hCentrality", ";centrality;counts", 100, 0, 100);
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 0., 2.5);
  TH2F* hMassVsOpen = new TH2F("hMassVsOpen", "Mass vs opening angle", 40, massM*0.96, massM*1.04, 50, 0., TMath::Pi());

  TH2F *hDistXY = new TH2F("hDistXY", ";d_{xy} (cm); #it{p}_{T} (GeV/#it{c}) ;counts", 100, 0, 0.1, 100, 0, 3);
  TH2F *hDistXYPlane = new TH2F("hDistXYPlane", ";x (cm); y (cm);counts", 100, -0.004, 0.004, 100, -0.0000001, 0.00000001);
  TH2F *hDistVsPt = new TH2F("hDistVsPt", ";z (cm); #it{p}_{T} (GeV/#it{c});counts", 300, 0, 3, 30, 0, 3);
  TH2F *hDistzVsPt = new TH2F("hDistzVsPt", ";z (cm); #it{p}_{T} (GeV/#it{c});counts", 300, 0, 3, 30, 0, 3);
  TH1F *hDist = new TH1F("hDist", ";z (cm); counts", 100, 0, 3);
  TH1F *hDistz = new TH1F("hDistz", ";z (cm); counts", 100, 0, 0.0001);

  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", "", 100, -0.1, 0.1, 30, 0, 3);
  //TH2F *//hd0XYprod = new TH2F("//hd0xyprod", "", 100, -0.01, 0.01, 30, 0, 3);
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
  
  TH1D *hcand = new TH1D("hcand", "", 1000, 0, 500000);
  TH1D *hcandpeak = new TH1D("hcandpeak", "", 500, 0, 15000);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
    

  TFile *fnt = 0x0;
  TFile *fhsp = 0x0;
  TNtuple *ntcand = 0x0;
  const Int_t nAxes = 12;
  Double_t arrhsp[nAxes];
  Float_t  arrnt[nAxes+1];

  if (writeNtuple){
    fnt = new TFile(Form("fntBkg%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "m:pt:rapidity:dist:cosp:d0prod:dca:ptMin:ptMax:thetad:true:nsigma1:nsigma2", 32000);
  }

  THnSparse* hsp = 0x0;
  if (writeSparse){
    fhsp = new TFile(Form("fhspBkg%s.root",suffix.Data()), "recreate");
    TString axTit[nAxes]={"m",
                          "pt",
                          "rapidity",
                          "dist",
                          "cosp",
                          "d0prod",
                          "dca",
                          "ptMin",
                          "ptMax",
                          "thetad",
                          "nsigma1",
                          "nsigma2"};
    //                        m            pt    y    d  cosp     d0   dca   pti pta  the  true ns1 ns2
    Int_t bins[nAxes] =   {40,             30,   30,  30,   30,       30,   30,   30, 30,  30, 16, 16}; 
    Double_t min[nAxes] = {0.96*massM,  ptminSG,  1,  0.,-1.00, -0.00015,    0,    0,  0,   0,  0,  0};
    Double_t max[nAxes] = {1.04*massM,  ptmaxSG,  5,  1., 1.00,  0.00015,  0.1,   30, 30, 0.5,  8,  8};  
    hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
    for(Int_t iax=0; iax<nAxes; iax++)
      hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  }

  //read the pdg code of the daughters
  int pdg_dau[2] = {0,0};
  GetPDGDaughters(pdgMother,pdg_dau,true);
  int swap_mass = TMath::Abs(pdg_dau[0]) != TMath::Abs(pdg_dau[1]) ? 2 : 1;
  KMCProbeFwd recProbe[2],recProbeTo0[2];
  TLorentzVector parent, daurec[2];
  Int_t trueCand = 0;

  for (Int_t iev = 0; iev < nevents; iev++){
    hCentrality->Fill(1);
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    treeBkg->GetEvent(iev);
    Int_t arrentrBkg = arrBkg->GetEntriesFast();
    if(fullsim)
      treeSig->GetEvent(iev);
    Int_t arrentrSig = (fullsim) ? arrSig->GetEntriesFast() : 0;

    for (Int_t itr = 0; itr < arrentrBkg+arrentrSig; itr++){
       KMCProbeFwd *tr1;
      if(itr < arrentrBkg)
        tr1 = (KMCProbeFwd *)arrBkg->At(itr);
      else
        tr1 = (KMCProbeFwd *)arrSig->At(itr-arrentrBkg);
      Float_t ch1 = tr1->GetCharge();
      for (Int_t itr2 = itr+1; itr2 < arrentrBkg+arrentrSig; itr2++){ 
        KMCProbeFwd *tr2;
        if(itr2<arrentrBkg)
          tr2 = (KMCProbeFwd *)arrBkg->At(itr2);
        else
          tr2 = (KMCProbeFwd *)arrSig->At(itr2-arrentrBkg);

        Float_t ch2 = tr2->GetCharge();
        if (ch1 * ch2 > 0) continue;
        if (ch1 < 0){ //convention: first track negative
          recProbe[0] = *tr1;
          recProbe[1] = *tr2;
        }
        else if (ch2 < 0){
          recProbe[0] = *tr2;
          recProbe[1] = *tr1;
        }
        
        trueCand = IsTrueCandidate(tr1, tr2, pdgMother, pdg_dau);
        
        Double_t tTOF1 = tr1->GetTOF();
        Double_t L1 = tr1->GetL();
        Double_t beta1 = betaTOF(L1, tTOF1);

        Double_t tTOF2 = tr2->GetTOF();
        Double_t L2 = tr2->GetL();
        Double_t beta2 = betaTOF(L2, tTOF2); 

        Double_t pxyz[3] = {0, 0, 0};
        Double_t pxyz2[3] = {0, 0, 0};

        recProbe[0].PropagateToDCA(&recProbe[1]);
        recProbe[0].GetPXYZ(pxyz);
        recProbe[1].GetPXYZ(pxyz2);

        Double_t betaDau1 = 0; 
        Double_t betaDau2 = 0; 

        for(Int_t iMassHyp=0; iMassHyp < swap_mass; iMassHyp++){
          countCand++;
          Int_t iNeg=-1;
          if(iMassHyp==0){
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass());
            daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass());
            betaDau2 = betaGen(pxyz, TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass()); 
            betaDau1 = betaGen(pxyz2, TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass()); 
            iNeg=0;
          }else{
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass());
            daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass());
            betaDau1 = betaGen(pxyz, TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass()); 
            betaDau2 = betaGen(pxyz2, TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass()); 
            iNeg=1;
          }
          Double_t nsigma1 = TMath::Abs((beta1-betaDau1))/sigmaBeta(zTOF);
          Double_t nsigma2 = TMath::Abs((beta2-betaDau2))/sigmaBeta(zTOF);
          
          parent = daurec[0];
          parent += daurec[1];
          Float_t pt = parent.Pt();
          Float_t invMass = parent.M();
          if(invMass<0.96*massM  || invMass>1.04*massM) continue;
          Float_t y = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
          Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
          Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
          Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
          Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
          

          Double_t xP = 0, yP = 0, zP = 0;
          ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP); 
          Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
          if(dist < 0.5 && pdgMother!=333)
            continue;
          Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
          Double_t vsec[3] = {xP, yP, zP};
          Double_t thetad = OpeningAngle(daurec[0],daurec[1]);
          Double_t cosp = CosPointingAngle(vprim, vsec, parent);
          if(cosp < 0.999 && pdgMother!=333)
            continue;
          
          recProbeTo0[0] = recProbe[0];
          recProbeTo0[0].PropagateToZBxByBz(0);
          Double_t d0x1 = recProbeTo0[0].GetX();
          Double_t d0y1 = recProbeTo0[0].GetY();
          Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
          if (d0x1 < 0) 
            d0xy1 *= -1;
          
          recProbeTo0[1] = recProbe[1];
          recProbeTo0[1].PropagateToZBxByBz(0);
          Double_t d0x2 = recProbeTo0[1].GetX();
          Double_t d0y2 = recProbeTo0[1].GetY();
          Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
          if (d0x2 < 0)
            d0xy2 *= -1;
          
          hYPtRecoAll->Fill(y, pt);
          hPtRecoAll->Fill(pt);
          hMassAll->Fill(invMass);
          hDCA->Fill(dca, pt);
          hDCAx->Fill(d1, pt);
          hDCAy->Fill(d2, pt);
          hDCAz->Fill(d3, pt);
          hMassVsOpen->Fill(invMass, TMath::ACos(thetad));
          hCosp->Fill(cosp, pt);
          hDistXY->Fill(distXY, pt);
          hDistXYPlane->Fill(distXY, zP);
          hDistVsPt->Fill(dist, pt);
          hDistzVsPt->Fill(zP, pt);
          hDist->Fill(dist);
          hDistz->Fill(zP);
          hd0XY1->Fill(d0xy1, pt);
          hd0XY2->Fill(d0xy2, pt);

          countCandInPeak++;
          if (writeNtuple){
            arrnt[0] = invMass;
            arrnt[1] = pt;
            arrnt[2] = y;
            arrnt[3] = dist;
            arrnt[4] = cosp;
            arrnt[5] = d0xy1 * d0xy2;
            arrnt[6] = dca;
            arrnt[7] = TMath::Min(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
            arrnt[8] = TMath::Max(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
            arrnt[9] = TMath::ACos(thetad);
            arrnt[10] = trueCand;
            arrnt[11] = nsigma1;
            arrnt[12] = nsigma2;
            ntcand->Fill(arrnt);
          }

          if (writeSparse){
            arrhsp[0] = invMass;
            arrhsp[1] = pt;
            arrhsp[2] = y;
            arrhsp[3] = dist;
            arrhsp[4] = cosp;
            arrhsp[5] = d0xy1 * d0xy2;
            arrhsp[6] = dca;
            arrhsp[7] = TMath::Min(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
            arrhsp[8] = TMath::Max(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
            arrhsp[9] = TMath::ACos(thetad);
            arrhsp[10] = nsigma1;
            arrhsp[11] = nsigma2;
            hsp->Fill(arrhsp);
          }
	    
        } // loop on mass hypothesis
      } // loop on first track
    } // loop on second track
    
    hcand->Fill(countCand);
    hcandpeak->Fill(countCandInPeak);
    printf(" --> Event %d, tot candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hCentrality->Write();
  hNevents->Write();
  hcand->Write();
  hcandpeak->Write();  
  hMassAll->Write();
  hYPtRecoAll->Write();
  hMassVsOpen->Write();
  hPtRecoAll->Write();
  hDistXY->Write();
  hDistVsPt->Write();
  hDist->Write();
  hDistXYPlane->Write();
  hDistzVsPt->Write();
  hDistz->Write();
  hDCA->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAz->Write();
  hCosp->Write();
  hCosThStVsMass->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  fout->Close();
  if (writeNtuple){
    fnt->cd();
    ntcand->Write("", TObject::kWriteDelete);
    fnt->Close();
  }
  if(writeSparse){
    fhsp->cd();
    hsp->Write();
    fhsp->Close();
  }
}
