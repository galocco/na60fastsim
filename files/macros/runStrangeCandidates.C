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
  TH1D* hMassTrue = new TH1D("hMassTrue", "Mass all match;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 2*mass);
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
  TH2F *hDCA12 = new TH2F("hDCA12", ";DCA_{12} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCA13 = new TH2F("hDCA13", ";DCA_{13} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCA23 = new TH2F("hDCA23", ";DCA_{23} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
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

  TFile *f = new TFile(Form("treeStrangeParticles%s.root",suffix.Data()), "RECREATE");
  TTree *tree = new TTree("tree", "tree Strange");
  TClonesArray *arrtr = new TClonesArray("KMCProbeFwd");
  TClonesArray &aarrtr = *arrtr;
  tree->Branch("tracks", &arrtr);
  Int_t icount = 0;

  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  TNtuple *ntgen = 0x0;
  if (writeNtuple)
  {
    fnt = new TFile(Form("fntSig%s.root",suffix.Data()), "recreate");
    ntgen = new TNtuple("ntgen", "ntgen", "pt:centrality:rapidity:dist:ct", 32000); 
    if(nbody == 2)
      ntcand = new TNtuple("ntcand", "ntcand", "m:pt:centrality:rapidity:dist:ct:cosp:d01:d02:d0prod:ipD:ptMin:ptMax:dca:thetad:arm", 32000);
    else
      ntcand = new TNtuple("ntcand", "ntcand", "m:pt:centrality:rapidity:dist:ct:cosp:d01:d02:sigmavert:ipD:ptMin:ptMax:dca12:dca13:dca23", 32000);
  }

  Float_t* arrnt = new Float_t[(nbody == 2) ? 16 : 16];
  Float_t arrntgen[5];
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

        if(IsReconstructable && nrec==0){
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


        trw->SetIndexMom(iev);
        trw->SetIndex(icount);
        trw->SetPdg(kf);
        trw->SetPdgMother(pdgParticle);
        new (aarrtr[icount]) KMCProbeFwd(*trw);
        icount++;
        nfake += trw->GetNFakeITSHits();
        if(trw->GetNFakeITSHits()==0)
          hChi2True->Fill(trw->GetNormChi2(kTRUE));
        else
          hChi2Fake->Fill(trw->GetNormChi2(kTRUE));

        if(nbody==2)
          recProbe[crg_index] = *trw;
        else
          recProbe[nrec] = *trw;
        nrec++;
      }
    }
    if(nrec < nbody)
      continue;//not both daughters are reconstructed

    hCentrality->Fill(5*gRandom->Rndm());
    if(nbody == 2){
      recProbe[0].PropagateToDCA(&recProbe[1]);
    }
    if(nbody == 3){
      recProbe[0].PropagateToDCA(&recProbe[2]);
      recProbe[1].PropagateToDCA(&recProbe[2]);
    }


    Double_t pxyz[3]={0,0,0};
    for(int i = 0; i < nbody; ++i){
      recProbe[i].GetPXYZ(pxyz);
      daurec[i].SetXYZM(pxyz[0], pxyz[1], pxyz[2], recProbe[i].GetMass());
      if(i==0){
        parent = daurec[0];
        parentgen = daugen[0];
      }
      else{
        parent += daurec[i];
        parentgen += daugen[i];
      }

      secvertgen[0][i] = recProbe[i].GetX();
      secvertgen[1][i] = recProbe[i].GetY();
      secvertgen[2][i] = recProbe[i].GetZ();
    }
    hMGen->Fill(parentgen.M());

    Double_t  ptRec=parent.Pt();
    Double_t  massRec=parent.M();
    Double_t yRec = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));

    Float_t dca = 0, dca12 = 0, dca13 = 0, dca23 = 0;
    Double_t xP = 0, yP = 0, zP = 0, sigmaVert, arm = 0;
    Double_t thetad = OpeningAngle(daurec[0],daurec[1]);;
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
      ComputeVertex(recProbe[0],recProbe[1],recProbe[2],vprim[2],xP,yP,zP,sigmaVert);

      Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
      Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
      Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
      dca12 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
      d1 = recProbe[2].GetX() - recProbe[0].GetX();
      d2 = recProbe[2].GetY() - recProbe[0].GetY();
      d3 = recProbe[2].GetZ() - recProbe[0].GetZ();
      dca13 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
      d1 = recProbe[2].GetX() - recProbe[1].GetX();
      d2 = recProbe[2].GetY() - recProbe[1].GetY();
      d3 = recProbe[2].GetZ() - recProbe[1].GetZ();
      dca23 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
      hDCA12->Fill(dca12, ptRec);
      hDCA13->Fill(dca13, ptRec);
      hDCA23->Fill(dca23, ptRec);
    }

    hYPtRecoAll->Fill(yRec, ptRec);
    hPtRecoAll->Fill(ptRec);
    hPtGenRecoAll->Fill(ptGen);
    hPtRecoVsGenAll->Fill(ptGen,ptRec);
    hDiffPtRecoGenAll->Fill(ptGen,(ptRec-ptGen));
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
    
    Double_t residVx = 10000.*(xP - secvertgen[0][0]);
    Double_t residVy = 10000.*(yP - secvertgen[1][0]);
    Double_t residVz = 10000.*(zP - secvertgen[2][0]);

    hResVxVsPt->Fill(residVx, ptRec);
    hResVyVsPt->Fill(residVy, ptRec);
    hResVzVsPt->Fill(residVz, ptRec);

    hResPtVsPt->Fill(ptRec-ptGen,ptRec);
    hResPtVsY->Fill(ptRec-ptGen,yRec);
    
    hResVxVsY->Fill(residVx, yRec);
    hResVyVsY->Fill(residVy, yRec);
    hResVzVsY->Fill(residVz, yRec);
    
    hResPxVsPt->Fill(parent.Px() - parentgen.Px(), ptRec);
    hResPyVsPt->Fill(parent.Py() - parentgen.Py(), ptRec);
    hResPzVsPt->Fill(parent.Pz() - parentgen.Pz(), ptRec);

    hResPxVsY->Fill(parent.Px() - parentgen.Px(), yRec);
    hResPyVsY->Fill(parent.Py() - parentgen.Py(), yRec);
    hResPzVsY->Fill(parent.Pz() - parentgen.Pz(), yRec);

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
    Float_t distgen = TMath::Sqrt(secvertgen[0][0] * secvertgen[0][0] + secvertgen[1][0] * secvertgen[1][0] + secvertgen[2][0] * secvertgen[2][0]);
    Float_t distgenXY = TMath::Sqrt(secvertgen[0][0] * secvertgen[0][0] + secvertgen[1][0] * secvertgen[1][0]);

    Double_t vsec[3] = {xP, yP, zP};
    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
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

    if (ntcand){
      arrnt[0] = massRec;
      arrnt[1] = ptRec;
      arrnt[2] = 1;
      arrnt[3] = yRec;
      arrnt[4] = (zP > 0) ? dist : -dist;
      arrnt[5] = mass*dist/TMath::Sqrt(parent.Pt()*parent.Pt()+parent.Pz()*parent.Pz());
      arrnt[6] = cosp;
      arrnt[7] = d0xy[0];
      arrnt[8] = d0xy[1];
      arrnt[9] = (nbody == 2) ? d0xy[0] * d0xy[1] : sigmaVert;
      arrnt[10] = ipD;
      if(nbody == 2){
        arrnt[11] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
        arrnt[12] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
        arrnt[13] = dca; 
        arrnt[14] = thetad; 
        arrnt[15] = arm; 
      }
      else{
        arrnt[11] = TMath::Min(TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
        arrnt[12] = TMath::Max(TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
        arrnt[13] = dca12;
        arrnt[14] = dca13;
        arrnt[15] = dca23;
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




void MakeCombinBkgCandidates3Body(const char* trackTreeFile="treeBkgEvents.root",
             TString suffix = "_Omega",
             const char *setup = "../setups/setup_short_5pixel_1.5T.txt",
			       Int_t nevents = 999999, 
			       Int_t writeNtuple = kFALSE,
			       Bool_t usePID=kFALSE,
             int pdgMother = 3334,
             int pdg_unstable_dau = 3122,
             int minITSHits = 5){

  // Read the TTree of tracks produced with runBkgVT.C
  // Create combinatorial background candidates (= OS pairs of tracks)
  // Store in THnSparse and (optionally) TNtuple

  TFile *filetree = new TFile(trackTreeFile);
  TTree *tree = (TTree *)filetree->Get("tree");
  TClonesArray *arr = 0;
  tree->SetBranchAddress("tracks", &arr);
  Int_t entries = tree->GetEntries();
  printf("Number of events in tree = %d\n",entries);
  if(nevents>entries) nevents=entries;
  else printf(" --> Analysis performed on first %d events\n",nevents);
  if(usePID) printf("Rough PID cuts will be used\n");
  
  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

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
  det->ReadSetup(setup, setup);
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


  TFile *fout = new TFile(Form("Bkg-histos%s.root",suffix.Data()), "recreate");
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);

  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 1., 3.5);
  TH1D* hMassKK = new TH1D("hMassKK", "KK Inv Mass ", 200, 0.9, 2.9);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistZ = new TH2F("hDistZ", "", 100, 0, 0.2, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY3 = new TH2F("hd0xy3", "", 100, -0.1, 0.1, 30, 0, 3);

  TH1D* hMomPion = new TH1D("hMomPion","",200,0.,10.);
  TH1D* hMomKaon = new TH1D("hMomKaon","",200,0.,10.);
  TH1D* hMomProton = new TH1D("hMomProton","",200,0.,10.);
  
  TH2F *hResVx = new TH2F("hResVx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 100, -1, 1, 30, 0, 3);
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  
  TH1D *hcand = new TH1D("hcand", "", 1000, 0, 500000000);
  TH1D *hcandpeak = new TH1D("hcandpeak", "", 500, 0, 15000000);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
  
  // define mother particle mass
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgMother)->Mass();

  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  Float_t arrnt[14];
  if (writeNtuple){
    fnt = new TFile(Form("fntBkg%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "m:pt:centrality:rapidity:dist:cosp:d01:d02:sigmavert:ptMin:ptMax:dca12:dca13:dca23", 32000);
  }

  //read the pdg code of the daughters
  int pdg_dau[2] = {0,0};
  GetPDGDaughters(pdgMother, pdg_dau, pdgMother>0);
  //read the pdg code of the unstable daughter daughters
  int pdg_dau2[2] = {0,0};
  GetPDGDaughters(pdg_unstable_dau, pdg_dau2, pdg_unstable_dau > 0);
  int pdgDaughter1 = (pdg_dau[0] == pdg_unstable_dau) ? pdg_dau[0] : pdg_dau[1];
  int pdgDaughter2 = pdg_dau2[0];
  int pdgDaughter3 = pdg_dau2[1];
  KMCProbeFwd recProbe[3];
  TLorentzVector parent, pair13, daurec[3];

  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    tree->GetEvent(iev);
    Int_t arrentr = arr->GetEntriesFast();

    Double_t pxyz0[3],pxyz1[3],pxyz2[3];
    for (Int_t itr = 0; itr < arrentr; itr++){
      KMCProbeFwd *tr1 = (KMCProbeFwd *)arr->At(itr);
      // cout << "tr P=" << tr1->GetP() << endl;
      Float_t ch1 = tr1->GetCharge();
      recProbe[0] = *tr1;
      recProbe[0].PropagateToZBxByBz(vprim[2]);
      recProbe[0].GetPXYZ(pxyz0);
      Double_t d0x1 = recProbe[0].GetX();
      Double_t d0y1 = recProbe[0].GetY();
      Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
      if (d0x1 < 0) d0xy1 *= -1;

      for (Int_t itr2 = 0; itr2 < arrentr; itr2++){
        if(itr2==itr) continue;
        KMCProbeFwd *tr2 = (KMCProbeFwd *)arr->At(itr2);
        Float_t ch2 = tr2->GetCharge();
        // convention: charge signs are ordered as +-+ or -+-
        if (ch1 * ch2 > 0) continue;
        recProbe[1] = *tr2;
        recProbe[1].PropagateToZBxByBz(vprim[2]);
        recProbe[1].GetPXYZ(pxyz1);
        Double_t d0x2 = recProbe[1].GetX();
        Double_t d0y2 = recProbe[1].GetY();
        Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
        if (d0x2 < 0) d0xy2 *= -1;

        for (Int_t itr3 = itr2+1; itr3 < arrentr; itr3++){
          if(itr3==itr) continue;
          KMCProbeFwd *tr3 = (KMCProbeFwd *)arr->At(itr3);
          Float_t ch3 = tr3->GetCharge();
          // convention: charge signs are ordered as +-+ or -+-
          if (ch3 * ch2 > 0) continue;
          recProbe[2] = *tr3;
          recProbe[2].PropagateToZBxByBz(vprim[2]);
          recProbe[2].GetPXYZ(pxyz2);
          Double_t d0x3 = recProbe[2].GetX();
          Double_t d0y3 = recProbe[2].GetY();
          Double_t d0xy3 = TMath::Sqrt(d0x3 * d0x3 + d0y3 * d0y3);
          if (d0x3 < 0) d0xy3 *= -1;
          //printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);

          for(Int_t iMassHyp=0; iMassHyp<2; iMassHyp++){
            // mass hypothesis: 123, 321, 
            Double_t mom1=0;
            Double_t mom2=0;
            Double_t mom3=0;
            if(iMassHyp==0){
              daurec[0].SetXYZM(pxyz0[0], pxyz0[1], pxyz0[2], TDatabasePDG::Instance()->GetParticle(pdgDaughter1)->Mass());
              daurec[1].SetXYZM(pxyz1[0], pxyz1[1], pxyz1[2], TDatabasePDG::Instance()->GetParticle(pdgDaughter2)->Mass());
              daurec[2].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdgDaughter3)->Mass());
              mom1=recProbe[2].GetTrack()->P();
              mom2=recProbe[1].GetTrack()->P();
              mom3=recProbe[0].GetTrack()->P();
              pair13 = daurec[1];
              pair13 += daurec[0];
            }else{
              daurec[0].SetXYZM(pxyz0[0], pxyz0[1], pxyz0[2], TDatabasePDG::Instance()->GetParticle(pdgDaughter3)->Mass());
              daurec[1].SetXYZM(pxyz1[0], pxyz1[1], pxyz1[2], TDatabasePDG::Instance()->GetParticle(pdgDaughter2)->Mass());
              daurec[2].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdgDaughter1)->Mass());
              mom1=recProbe[0].GetTrack()->P();
              mom2=recProbe[1].GetTrack()->P();
              mom3=recProbe[2].GetTrack()->P();
              pair13 = daurec[1];
              pair13 += daurec[2];
            }
            parent = daurec[0];
            parent += daurec[1];
            parent += daurec[2];
          
            countCand++;
            Float_t pt=parent.Pt();
            Float_t invMass=parent.M();
            Float_t y = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
            Double_t  massRec13=pair13.M();
            hYPtRecoAll->Fill(y, pt);
            hPtRecoAll->Fill(pt);
            hMassAll->Fill(invMass);
            hMassKK->Fill(massRec13);
            if(invMass>mass*0.5  && invMass<mass*1.5){
                    // range to fill histos
              if(TMath::Abs(invMass-mass)<0.06) countCandInPeak++;
              hMomPion->Fill(mom1);
              hMomKaon->Fill(mom2);
              hMomKaon->Fill(mom3);

              Double_t xV, yV, zV;
              Double_t sigmaVert;
              ComputeVertex(recProbe[0],recProbe[1],recProbe[2],vprim[2],xV,yV,zV,sigmaVert);
              Float_t dist = TMath::Sqrt(xV * xV + yV * yV + zV * zV);
              Float_t distXY = TMath::Sqrt(xV * xV + yV * yV);
              Float_t distZ = zV;
              Double_t vsec[3] = {xV, yV, zV};
              Double_t cosp = CosPointingAngle(vprim, vsec, parent);
              Double_t ipD = ImpParXY(vprim, vsec, parent);
              hCosp->Fill(cosp, pt);
              //printf(" ***** ***** cos point = %f \n", cosp);	    
              hDistXY->Fill(distXY, pt);
              hDistZ->Fill(zV, pt);
              hDist->Fill(dist, pt);
              hd0XY1->Fill(d0xy1, pt);
              hd0XY2->Fill(d0xy2, pt);
              hd0XY3->Fill(d0xy3, pt);

              recProbe[0].PropagateToDCA(&recProbe[1]);
              recProbe[0].PropagateToDCA(&recProbe[2]);
              recProbe[1].PropagateToDCA(&recProbe[2]);

              Float_t dca12,dca13,dca23;
              Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
              Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
              Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
              dca12 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
              d1 = recProbe[2].GetX() - recProbe[0].GetX();
              d2 = recProbe[2].GetY() - recProbe[0].GetY();
              d3 = recProbe[2].GetZ() - recProbe[0].GetZ();
              dca13 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
              d1 = recProbe[2].GetX() - recProbe[1].GetX();
              d2 = recProbe[2].GetY() - recProbe[1].GetY();
              d3 = recProbe[2].GetZ() - recProbe[1].GetZ();
              dca23 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);

              if (ntcand){
                arrnt[0] = invMass;
                arrnt[1] = pt;
                arrnt[2] = 1;
                arrnt[3] = y;
                arrnt[4] = dist;
                arrnt[5] = cosp;
                arrnt[6] = d0xy1;
                arrnt[7] = d0xy2;
                arrnt[8] = sigmaVert;
                arrnt[9] = TMath::Min(TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
                arrnt[10] = TMath::Max(TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
                arrnt[11] = dca12;
                arrnt[12] = dca13;
                arrnt[13] = dca23;
                
                ntcand->Fill(arrnt);
              }
            } // check on inv mass
          } // loop on mass hypothesis
	      } // loop on third track
      } // loop on second track
    }// loop on first track
    hcand->Fill(countCand);
    hcandpeak->Fill(countCandInPeak);
    printf(" --> Event %d, tot  candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hNevents->Write();
  hcand->Write();
  hcandpeak->Write();  
  hMassAll->Write();
  hMassKK->Write();
  hYPtRecoAll->Write();
  hPtRecoAll->Write();
  hDistXY->Write();
  hDistZ->Write();
  hDist->Write();
  hDCA->Write();
  hCosp->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hd0XY3->Write();
  hMomPion->Write();
  hMomKaon->Write();
  hMomProton->Write();
  fout->Close();
  if (ntcand){
    fnt->cd();
    ntcand->Write();
    fnt->Close();
  }
}
 void MakeCombinBkgCandidates2Body(const char* trackTreeFileBkg="treeBkgEvents_layer5.root",
            const char* trackTreeFileSig="treeBkgEvents_layer5.root",
            TString suffix = "_PHI",
            int pdgMother = 333,
            const char *setup = "../setups/setup-EHN1_BetheBloch.txt",
			      Int_t writeNtuple = kTRUE,
            Int_t fullsim = kFALSE){

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
  TNtuple *ntcand = 0x0;
  Float_t arrnt[17];
  if (writeNtuple){
    fnt = new TFile(Form("fntBkg%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "m:pt:centrality:rapidity:dist:ct:cosp:d01:d02:d0prod:ipD:dca:ptMin:ptMax:thetad:arm:true", 32000);
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
    Int_t arrentrSig = (fullsim) ? arrBkg->GetEntriesFast() : 0; 
    
    for (Int_t itr = 0; itr < arrentrBkg+arrentrSig; itr++){
       KMCProbeFwd *tr1;
      if(itr<arrentrBkg)
        tr1 = (KMCProbeFwd *)arrBkg->At(itr);
      else
        tr1 = (KMCProbeFwd *)arrSig->At(itr);
      Bool_t pdgMom1 = tr1->GetPdgMother() == pdgMother;
      Bool_t pdgDa1 = tr1->GetPdg() == pdg_dau[0] || tr1->GetPdg() == pdg_dau[1];
      Float_t ch1 = tr1->GetCharge();
      for (Int_t itr2 = itr; itr2 < arrentrBkg+arrentrSig; itr2++){ 
        KMCProbeFwd *tr2;
        if(itr2<arrentrBkg)
          tr2 = (KMCProbeFwd *)arrBkg->At(itr2);
        else
          tr2 = (KMCProbeFwd *)arrSig->At(itr2);

        Float_t ch2 = tr2->GetCharge();
        if (ch1 * ch2 > 0) continue;
        if (ch1 < 0){ //convention: first track negative
          recProbe[0] = *tr1;
          recProbe[1] = *tr2;
        }else if (ch2 < 0){
          recProbe[0] = *tr2;
          recProbe[1] = *tr1;
        }
        
        Bool_t pdgMom2 = tr2->GetPdgMother() == pdgMother;
        Bool_t index = tr2->GetIndexMom() == tr1->GetIndexMom();
        Bool_t pdgDa2 = (tr2->GetPdg() == pdg_dau[0] || tr2->GetPdg() == pdg_dau[1]) && tr2->GetPdg()!=tr1->GetPdg() ;
        if (pdgMom1 && pdgMom2 && index && pdgDa1 && pdgDa2)
          trueCand = 1;
        else
          trueCand = 0;

        Double_t pxyz[3] = {0, 0, 0};
        Double_t pxyz2[3] = {0, 0, 0};

        recProbe[0].PropagateToDCA(&recProbe[1]);
        recProbe[0].GetPXYZ(pxyz);
        recProbe[1].GetPXYZ(pxyz2);
	
        for(Int_t iMassHyp=0; iMassHyp < swap_mass; iMassHyp++){
          Int_t iNeg=-1;
          if(iMassHyp==0){
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass());
            daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass());
            iNeg=0;
          }else{
            daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass());
            daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass());
            iNeg=1;
          }
          parent = daurec[0];
          parent += daurec[1];
          countCand++;
          Float_t pt = parent.Pt();
          Float_t invMass = parent.M();
          Float_t y = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
          hYPtRecoAll->Fill(y, pt);
          hPtRecoAll->Fill(pt);
          hMassAll->Fill(invMass);
          if(invMass>0.8*massM  && invMass<1.2*massM){
            // range to fill histos
            if(invMass>0.8*massM && invMass<1.2*massM) countCandInPeak++;
            Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
            Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
            Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
            Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
            
            hDCA->Fill(dca, pt);
            hDCAx->Fill(d1, pt);
            hDCAy->Fill(d2, pt);
            hDCAz->Fill(d3, pt);

            Double_t xP = 0, yP = 0, zP = 0;
            ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP); 
            Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
            Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
            Double_t vsec[3] = {xP, yP, zP};
            Double_t cosp = CosPointingAngle(vprim, vsec, parent);
            Double_t thetad = OpeningAngle(daurec[0],daurec[1]);
            Double_t cts = CosThetaStar(parent,daurec[iNeg],pdgMother,pdg_dau[0],pdg_dau[1]);
            Double_t arm = ArmenterosAlpha(daurec[0],daurec[1]);
            Double_t ipD = ImpParXY(vprim, vsec, parent);

            hCosp->Fill(cosp, pt);
            hCosThStVsMass->Fill(invMass,cts);

            hDistXY->Fill(distXY, pt);
            hDistXYPlane->Fill(xP, zP);
            hDistVsPt->Fill(dist, pt);
            hDistzVsPt->Fill(zP, pt);
            hDist->Fill(dist);
            hDistz->Fill(zP);
            
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
            
            hd0XY1->Fill(d0xy1, pt);
            hd0XY2->Fill(d0xy2, pt);

            if (ntcand){
              arrnt[0] = invMass;
              arrnt[1] = pt;
              arrnt[2] = 1;
              arrnt[3] = y;
              arrnt[4] = dist;
              arrnt[5] = dist*massM/TMath::Sqrt(parent.Pt() * parent.Pt() + parent.Pz() * parent.Pz());
              arrnt[6] = cosp;
              arrnt[7] = d0xy1;
              arrnt[8] = d0xy2;
              arrnt[9] = d0xy1 * d0xy2;
              arrnt[10] = ipD;
              arrnt[11] = dca;
              arrnt[12] = TMath::Min(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
              arrnt[13] = TMath::Max(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
              arrnt[14] = thetad;
              arrnt[15] = arm;
              arrnt[16] = trueCand;
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
  hCentrality->Write();
  hNevents->Write();
  hcand->Write();
  hcandpeak->Write();  
  hMassAll->Write();
  hYPtRecoAll->Write();
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
  if (ntcand){
    fnt->cd();
    ntcand->Write();
    fnt->Close();
  }
}
/*
Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk, int pdgMother, int pdgDaughterPos, int pdgDaughterNeg) {

  Double_t massMoth = TDatabasePDG::Instance()->GetParticle(pdgMother)->Mass();
  Double_t massPos = TDatabasePDG::Instance()->GetParticle(pdgDaughterPos)->Mass();
  Double_t massNeg = TDatabasePDG::Instance()->GetParticle(pdgDaughterNeg)->Mass();

  Double_t pStar = TMath::Sqrt((massMoth*massMoth-massPos*massPos-massNeg*massNeg)*(massMoth*massMoth-massPos*massPos-massNeg*massNeg)-4.*massPos*massPos*massNeg*massNeg)/(2.*massMoth);

  Double_t pMoth=parent.P();
  Double_t e=TMath::Sqrt(massMoth*massMoth+pMoth*pMoth);
  Double_t beta = pMoth/e;
  Double_t gamma = e/massMoth;
  TVector3 momDau(dauk.Px(),dauk.Py(),dauk.Pz());
  TVector3 momMoth(parent.Px(),parent.Py(),parent.Pz());
  Double_t qlProng=momDau.Dot(momMoth)/momMoth.Mag();
  Double_t cts = (qlProng/gamma-beta*TMath::Sqrt(pStar*pStar+massPos*massPos))/pStar;

  return cts;
}

THnSparseF* CreateSparse(Double_t mass, Int_t nbody = 2){
  const Int_t nAxes=10;
  THnSparseF *hsp;
  if(nbody == 3){
    TString axTit[nAxes]={"Inv. mass (GeV/c^{2})",
                          "#it{p}_{T} (GeV/c)",
                          "y",
                          "ct (cm)",
                          "cos(#vartheta_{p})",
                          "#it{p}_{T}^{min} (GeV/c)",
                          "d_0^{min} (cm)",
                          "#bar{DCA} (cm)",
                          "b (cm)",
                          "sigmaVert (cm)"};
    Int_t bins[nAxes] =   {100,            10,     40,  30,   15,  10,    8,   12,   16,  10}; 
    Double_t min[nAxes] = {0.5*mass,  ptminSG, yminSG,  0., 0.98,  0., 0.00, 0.00, 0.00,  0.};
    Double_t max[nAxes] = {1.5*mass,  ptmaxSG, ymaxSG, 30., 1.00,  10., 0.05, 0.03, 0.04,  2.};  
    hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
    for(Int_t iax=0; iax<nAxes; iax++)
      hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  }
  else{
    TString axTit[nAxes]={"Inv. mass (GeV/c^{2})",
                          "#it{p}_{T} (GeV/c)",
                          "y",
                          "ct (cm)",
                          "cos(#vartheta_{p})",
                          "#it{p}_{T}^{min} (GeV/c)",
                          "d_0^{min} (cm)",
                          "DCA (dca)",
                          "b (cm)",
                          "cos(#theta*)"};
    Int_t bins[nAxes] =   {100,            10,     40,  30,   15,  10,    8,   12,   16,  20}; 
    Double_t min[nAxes] = {0.5*mass,  ptminSG, yminSG,  0., 0.98,  0., 0.00, 0.00, 0.00, -1.};
    Double_t max[nAxes] = {1.5*mass,  ptmaxSG, ymaxSG, 30., 1.00,  10., 0.05, 0.03, 0.04,  1.};  
    hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
    for(Int_t iax=0; iax<nAxes; iax++)
      hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  }
  return hsp;
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
  int counter = 0;
  int sign = matter ? 1 : -1;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      for(int i = 0; i < 2; i++)
        pdgDaughters[i] = pdg_daughter[index_pdg][i]*sign;
      return;
    }
    counter++;
  }
}

Double_t KFVertexer(AliExternalTrackParam kTrack [], Double_t kMasses[], int n_dau, float bz){
  double posmom[6],cov[21];
  KFParticle helper[3];

  for (int iT=0; iT < n_dau; iT++) {
    kTrack[iT].GetXYZ(posmom);
    kTrack[iT].GetPxPyPz(posmom+3);
    kTrack[iT].GetCovarianceXYZPxPyPz(cov);
    helper[iT].Create(posmom,cov,kTrack[iT].Charge(),kMasses[iT]);
    helper[iT].Chi2() = kTrack[iT].GetTPCchi2();
    helper[iT].NDF() = kTrack[iT].GetNumberOfTPCClusters() * 2;
    //helper[iT].SetFieldCoeff(bz,0);
  }
  
  KFParticle oneCandidate;
  oneCandidate.AddDaughter(helper[0]);
  KFParticle twoCandidate{oneCandidate};
  twoCandidate.AddDaughter(helper[1]);
  if(n_dau == 2)
    return twoCandidate.GetMass();
  else{
    KFParticle threeCandidate{twoCandidate};
    threeCandidate.AddDaughter(helper[2]);
    return threeCandidate.GetMass();

  }
}

Double_t O2Vertexer(AliExternalTrackParam kTrack[], Double_t kMasses[], int n_dau, float bz){  
  o2::vertexing::DCAFitter3 fVertexer3;
  o2::vertexing::DCAFitter2 fVertexer2;
  fVertexer3.setBz(bz);
  fVertexer2.setBz(bz);
  fVertexer3.setMaxChi2(10e9);
  fVertexer2.setMaxChi2(10e9);
  o2::track::TrackParCov helper[3];
  for (int iT=0; iT < n_dau; iT++)
    helper[iT] = *(o2::track::TrackParCov*)((AliExternalTrackParam*)&kTrack[iT]);
  int nVert{0};
  try {
    if(n_dau == 2)
      nVert = fVertexer2.process(helper[0], helper[1]);
    else if (n_dau == 3)
      nVert = fVertexer3.process(helper[0], helper[1], helper[2]);
    else{
      std::cout << "MyException caught" << std::endl;
      return -1;
    }
  }
  catch (std::runtime_error& e) {
    std::cout << "MyException caught" << std::endl;
    std::cout << e.what() << std::endl;
    return -1;
  }
  if (nVert) {
    if(n_dau == 2){
      auto vert = fVertexer2.getPCACandidate();
      auto& dau1 = fVertexer2.getTrack(0);
      auto& dau2 = fVertexer2.getTrack(1);
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> vector1{dau1.Pt(), dau1.Eta(), dau1.Phi(), kMasses[0]};
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> vector2{dau2.Pt(), dau2.Eta(), dau2.Phi(), kMasses[1]};
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> mother{vector1 + vector2};
      return mother.mass();
    }
    else{
      auto vert = fVertexer3.getPCACandidate();
      auto& dau1 = fVertexer3.getTrack(0);
      auto& dau2 = fVertexer3.getTrack(1);
      auto& dau3 = fVertexer3.getTrack(2);
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> vector1{dau1.Pt(), dau1.Eta(), dau1.Phi(), kMasses[0]};
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> vector2{dau2.Pt(), dau2.Eta(), dau2.Phi(), kMasses[1]};
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> vector3{dau3.Pt(), dau3.Eta(), dau3.Phi(), kMasses[2]};
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> mother{vector1 + vector2 + vector3};
      return mother.mass();
    }
  }
  else 
    return -1;
}
*/