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
// rapidity range
double yminSG = -10.; 
double ymaxSG = 10.;
// pT range
double ptminSG = 0.;
double ptmaxSG = 3.; 

double vX = 0, vY = 0, vZ = 0; // event vertex

THnSparseF* CreateSparse(Double_t mass, Int_t nbody);
TDatime dt;
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
double Tslope[2][NParticles][NEnergy] = {{{196.8,237.4,244.6,239.8,298.7},{0,0,228.9, 223.1, 228.9},{244,249,258,265,301},{0,0,218,0,267},{221,233,222,227,277}},//matter
                                         {{196.8,237.4,244.6,239.8,298.7},{0,0,226,217,226},{339,284,301,292,303},{0,0,218,0,259},{311,277,255,321,0}}};//antimatter
//sigma parameter of the gaussians of the rapidity distribution [matter/antimatter][particle][beam energy]
//                                                          phi                        K                       Lambda                Omega                    Csi
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

void GenerateSignalCandidates(Int_t nevents = 10000, 
				double Eint = 158, 
        TString suffix = "_phi",
				const char *setup = "../setups/setup-EHN1_BetheBloch.txt",
				const char *privateDecayTable = "../decaytables/USERTABPHI.DEC",
				bool writeNtuple = kTRUE, 
				bool simulateBg=kTRUE,
        int pdgParticle = 333,
        int nbody = 2,
        int pdg_unstable_dau = 0,
        bool matter = true){
  
  // Generate strange particle signals and simulate detector response for decay tracks
  //nevents *= GetMultiplicity(pdgParticle,Eint,matter);
  int refreshBg = 1000;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);
  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");
  
  //define the pT and rapidity probability function
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();
  TF1 *fpt = new TF1("fpt","x*exp(-TMath::Sqrt(x**2+[0]**2)/[1])",ptminSG,ptmaxSG);
  fpt->SetParameter(0,mass);
  fpt->SetParameter(1,GetTslope(pdgParticle,Eint,matter)/1000);
  TF1 *fy = new TF1("fy","exp(-0.5*((x-[0])/[1])**2)+exp(-0.5*((x+[0])/[1])**2) ",yminSG,ymaxSG);
  fy->SetParameter(0,GetY0Rapidity(pdgParticle,Eint,matter));
  fy->SetParameter(1,GetSigmaRapidity(pdgParticle,Eint,matter));
  int count1=0,count2=0,count3=0;
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
  // max number of see on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000);
  //
  // set chi2 cuts
  det->SetMaxChi2Cl(10.);  // max track to cluster chi2
  det->SetMaxChi2NDF(3.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(20);  // fiducial cut on chi2 of convergence to vtx  
  // IMPORTANT FOR NON-UNIFORM FIEL
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); //NA60+
  //det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  det->BookControlHistos();
  
  // prepare decays
  TGenPhaseSpace decay;
  TLorentzVector parentgen, daugen[3], parent, daurec[3], parentrefl, daurecswapmass[2]; 
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
  //fDecayer->ForceDecay();
  TClonesArray *particles = new TClonesArray("TParticle", 1000);
  TLorentzVector *mom = new TLorentzVector();
  
  
  TH2F* hYPtGen = new TH2F("hYPtGen", "y-#it{p}_{T} corr match;y;#it{p}_{T};counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtGen = new TH1D("hPtGen", "#it{p}_{Tgen};#it{p}_{T} (GeV/#it{c});counts", 40, ptminSG, ptmaxSG);
  TH1D* hYGen = new TH1D("hYGen", "y full phase space;y;counts", 160., yminSG,ymaxSG);
  TH1D* hPtEff = new TH1D("hPtEff", "#it{p}_{T} efficiency;#it{p}_{T} (GeV/#it{c});counts", 40., ptminSG, ptmaxSG);
  TH1D* hYEff = new TH1D("hYEff", "y efficiency;counts", 160., yminSG,ymaxSG);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "y-#it{p}_{T} all match;y;#it{p}_{T};counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Reconstructed #it{p}_{T} all match;#it{p}_{T};counts", 40, ptminSG, ptmaxSG);  
  TH1D* hPtGenRecoAll = new TH1D("hPtGenRecoAll", "Generated #it{p}_{T} all match;#it{p}_{T};counts", 40, ptminSG, ptmaxSG);
  TH2F* hPtRecoVsGenAll = new TH2F("hPtRecoVsGenAll"," ; Generated #it{p}_{T} ; Reconstructed #it{p}_{T}",40, ptminSG, ptmaxSG,40, ptminSG, ptmaxSG);
  TH2F* hDiffPtRecoGenAll = new TH2F("hDiffPtRecoGenAll"," ; Generated #it{p}_{T} ; Reco #it{p}_{T} - Gen #it{p}_{T}",40, ptminSG, ptmaxSG,100,-0.2,0.2);

  TH2F* hptau[3];
  TH1D* hyau[3];
  for(int i = 0; i < nbody; i++){
    hptau[i] = new TH2F(Form("hptdau%i",i), " ;#it{p}_{T} (GeV/#it{c});#it{p}_{TM} (GeV/#it{c});counts", 50,ptminSG, 3,50,ptminSG, ptmaxSG);
    hyau[i] = new TH1D(Form("hydau%i",i), ";y_{};counts", 50, 0., 5.);
  }
  TH2F *hyau2D = new TH2F("hydau2D", "y negative daughter vs y positive daughter;y_{-};y_{+};counts", 50, 0., 5., 50, 0., 5.);
  
  TH1D* hYRecoAll = new TH1D("hYRecoAll", "Reconstructed y all match;y;counts", 80., 1., 5.4);
  TH1D* hYGenRecoAll = new TH1D("hYGenRecoAll", "Generated y all match;y;counts", 80., 1., 5.4);
  TH2F* hYPtRecoFake = new TH2F("hYPtRecoFake", "y-#it{p}_{T} fake match;;counts", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoFake = new TH1D("hPtRecoFake", "#it{p}_{T} fake match;;counts", 40, ptminSG, ptmaxSG);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 1.5*mass);
  TH1D* hMassFake = new TH1D("hMassFake", "Mass fake match;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 1.5*mass);
  TH2F *hMassVsPt = new TH2F("hMassVsPt", ";m (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", 200, 0.5*mass, 1.5*mass, 6, 0, 3);
  TH2F *hMassVsY = new TH2F("hMassVsY", ";m (GeV/#it{c}^{2});y;counts", 200, 0.5*mass, 1.5*mass, 10, 0, 5);

  TH2F *hDistXY = new TH2F("hDistXY", ";d_{xy} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", ";d (cm);#it{p}_{T} (GeV/#it{c});counts", 300, 0, 10, 30, 0, 3);
  TH1D *hDistx = new TH1D("hDistx", "generated secondary vertex x;x (cm);counts", 100, 0, 30);
  TH1D *hDisty = new TH1D("hDisty", "generated secondary vertex y;y (cm);counts", 100, 0, 30);
  TH1D *hDistxy = new TH1D("hDistxy", "generated secondary vertex d_{xy};d_{xy} (cm);counts", 100, 0, 30);
  TH1D *hDistz = new TH1D("hDistz", "generated secondary vertex z;z (cm);counts", 200, 0, 200);
  TH1D *hDistTot = new TH1D("hDistTot", "generated secondary vertex distance;d (cm);counts", 200, 0, 200);

  TH1D *hDistxRec = new TH1D("hDistxRec", "reconstructed secondary vertex x;x (cm);counts", 30, 0, 1.5);
  TH1D *hDistyRec = new TH1D("hDistyRec", "reconstructed secondary vertex y;y (cm);counts", 30, 0, 1.5);
  TH1D *hDistxyRec = new TH1D("hDistxyRec", "reconstructed secondary vertex d_{xy};d_{xy} (cm);counts", 20, 0, 1);
  TH1D *hDistzRec = new TH1D("hDistzRec", "reconstructed secondary vertex z;z (cm);counts", 30, 0, 15);
  TH1D *hDistTotRec = new TH1D("hDistTotRec", "reconstructed secondary vertex distance;d (cm);counts", 30, 0, 15);

  TH1D *hCtz = new TH1D("hctz", ";ctz (cm);counts", 100, 0, 100);
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
  //TH2F *hd0XYprod = new TH2F("hd0xyprod", ";;#it{p}_{T} (GeV/#it{c});counts", 100, -0.01, 0.01, 30, 0, 3);
  //histograms for the mass of candidates built from misidentified, only if the daughters are stables
  TH1D* hMassRefl = new TH1D("hMassRefl", "Mass reflections;m (GeV/#it{c}^{2});counts", 200, 0.5*mass, 1.5*mass);
  TH2F *hMassReflVsPt = new TH2F("hMassReflVsPt", ";m (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", 200, 0.5*mass, 1.5*mass, 6, 0, 3);
  TH2F *hMassReflVsY = new TH2F("hMassReflVsY", ";m (GeV/#it{c}^{2});y;counts", 200, 0.5*mass, 1.5*mass, 10, 0, 5);
  
  TH2F *hResVx = new TH2F("hResVx", ";V_{xgen}-V_{xrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResVy = new TH2F("hResVy", ";V_{ygen}-V_{yrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResVz = new TH2F("hResVz", ";V_{zgen}-V_{zrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 200, -1000., 1000., 30, ptminSG, ptmaxSG);
  TH2F *hResPx = new TH2F("hResPx", ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, ptminSG, ptmaxSG);
  TH2F *hResPy = new TH2F("hResPy", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, ptminSG, ptmaxSG);
  TH2F *hResPz = new TH2F("hResPz", ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, ptminSG, ptmaxSG);
  TH2F *hResPt = new TH2F("hResPt", ";#it{p}_{Tgen}-#it{p}_{Trec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});counts", 100, -1, 1, 30, ptminSG, ptmaxSG);
  TH2F *hResVxVsY = new TH2F("hResVxVsY", ";V_{xgen}-V_{xrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", ";V_{ygen}-V_{yrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", ";V_{zgen}-V_{zrec} (cm);y;counts", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResPxVsY = new TH2F("hResPxVsY", ";#it{p}_{xgen}-#it{p}_{xrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5);
  TH2F *hResPyVsY = new TH2F("hResPyVsY", ";#it{p}_{ygen}-#it{p}_{yrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5);
  TH2F *hResPzVsY = new TH2F("hResPzVsY", ";#it{p}_{zgen}-#it{p}_{zrec} (GeV/#it{c});y;counts", 100, -1, 1, 50, 0, 5);
  TH2F *hResDist = new TH2F("hResDist", ";d_{gen}-d_{rec} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", ";d_{xygen}-d_{xyrec} (cm);#it{p}_{T} (GeV/#it{c});counts", 100, -0.1, 0.1, 30, 0, 3);
  TH1D *hNevents = new TH1D("hNevents", ";;counts", 1, 0, 1);
  TH1D *hChi2 = new TH1D("hChi2", ";#chi^{2};counts", 40, 0, 20);
  
  THnSparseF *hsp = CreateSparse(mass, nbody);
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  TNtuple *ntgen = 0x0;
  if (writeNtuple)
  {
    fnt = new TFile(Form("fntSig%s.root",suffix.Data()), "recreate");
    ntgen = new TNtuple("ntgen", "ntgen", "pt:y", 32000); 
    if(nbody == 2)
      ntcand = new TNtuple("ntcand", "ntcand", "mass:pt:y:dist:cosp:d01:d02:d0prod:ptMin:ptMax:dca", 32000);
    else
      ntcand = new TNtuple("ntcand", "ntcand", "mass:pt:y:dist:cosp:d01:d02:sigmavert:ptMin:ptMax:dca12:dca13:dca23", 32000);
  }

  Float_t* arrnt = new Float_t[(nbody == 2) ? 11 : 13];
  Float_t arrntgen[2];
  double y0 = GetY0(Eint);
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
    double pxyz[3];
    bool cond1 = true,cond2 = true,cond3 = true;
    if (simulateBg && (iev%refreshBg)==0) det->GenBgEvent(0.,0.,0.);
    Double_t ptGen = fpt->GetRandom(); 
    Double_t yGen = fy->GetRandom(yminSG,ymaxSG)+y0;
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t pxGen = ptGen * TMath::Cos(phi);
    Double_t pyGen = ptGen * TMath::Sin(phi);

    if (ntgen){
      arrntgen[0] = ptGen;
      arrntgen[1] = yGen;
      ntgen->Fill(arrntgen);
    }
    
    Double_t mt = TMath::Sqrt(ptGen * ptGen + mass * mass);
    Double_t pzGen = mt * TMath::SinH(yGen);
    Double_t en = mt * TMath::CosH(yGen);

    mom->SetPxPyPzE(pxGen, pyGen, pzGen, en);
    Int_t np;
    do{
      fDecayer->Decay(pdgParticle, mom);
      np = fDecayer->ImportParticles(particles);
    } while (np < 0);
    Double_t ptau[3];
    Double_t yau[3];
    Double_t secvertgen[3][3];
    Int_t pdgDaughter[3];
    // loop on decay products
    for (int i = 0; i < np; i++) {
      TParticle *iparticle1 = (TParticle *)particles->At(i);
      Int_t kf = iparticle1->GetPdgCode();
      vX = iparticle1->Vx();
      vY = iparticle1->Vy();
      vZ = iparticle1->Vz();
      if (kf == pdgParticle){
        // Mother particle
        hYGen->Fill(iparticle1->Y());
        hPtGen->Fill(iparticle1->Pt());
        hYPtGen->Fill(iparticle1->Y(), iparticle1->Pt());
        //printf("mother part = %d, code=%d, pt=%f, y=%f \n", i, kf, iparticle1->Pt(), iparticle1->Y());
      }
      //check if the particle is one of the final decay products
      bool IsDecayaughter = pdg_dau[0] == kf || pdg_dau[1] == kf;
      if(nbody == 3)
        IsDecayaughter = IsDecayaughter || pdg_dau2[0] == kf || pdg_dau2[1] == kf;
      bool IsStable = kf != pdg_unstable_dau;
      if (IsDecayaughter && IsStable){
        // daughters that can be reconstructed: K and pi
        Double_t e = iparticle1->Energy();
        Double_t px = iparticle1->Px();
        Double_t py = iparticle1->Py();
        Double_t pz = iparticle1->Pz();
        if(cond1 && nrec==0){
          hDistx->Fill(vX);
          hDisty->Fill(vY);
          hDistz->Fill(vZ);
          hDistxy->Fill(TMath::Sqrt(vX*vX+vY*vY));
          hDistTot->Fill(TMath::Sqrt(vX*vX+vY*vY+vZ*vZ));
          hCtz->Fill(vZ*mass/iparticle1->P());
        }
        TLorentzVector *pDecDau = new TLorentzVector(0., 0., 0., 0.);
        pDecDau->SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
        Int_t crg=1;
        if(iparticle1->GetPdgCode()<0) crg=-1;
        if (!det->SolveSingleTrack(pDecDau->Pt(), pDecDau->Rapidity(), pDecDau->Phi(), iparticle1->GetMass(), crg, vX, vY, vZ, 0, 1, 99)){
          if(cond1){ cond1 = false;count1++;}
          continue;
        }
        KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
        if (!trw){
          if(cond1){ cond1 = false;count2++;}
          continue;
        }
        hChi2->Fill(trw->GetNormChi2(kTRUE));
        if (trw->GetNormChi2(kTRUE) > ChiTot){
          if(cond1){ cond1 = false;count3++;}
          continue;
        }
        
        nfake += trw->GetNFakeITSHits();
        trw->GetPXYZ(pxyz);

        ptau[nrec] = iparticle1->Pt();
        yau[nrec] = iparticle1->Y();
        secvertgen[0][nrec] = iparticle1->Vx();
        secvertgen[1][nrec] = iparticle1->Vy();
        secvertgen[2][nrec] = iparticle1->Vz();
        daugen[nrec].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
        daurec[nrec].SetXYZM(pxyz[0], pxyz[1], pxyz[2], iparticle1->GetMass());
        if(nbody == 2)
        {
          daurecswapmass[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2],(iparticle1->GetMass() == massPos)? massNeg : massPos);
          hptau[(crg == 1)? 0 : 1]->Fill(ptau[nrec],ptGen);
          hyau[(crg == 1)? 0 : 1]->Fill(yau[nrec]);
        }
        else{
          if((pdg_dau[0] == kf || pdg_dau[1] == kf) && nrec==0){
            hptau[(crg == 1)? 0 : 1]->Fill(ptau[nrec],ptGen);
            hyau[(crg == 1)? 0 : 1]->Fill(yau[nrec]);
          }
          else{
            hptau[(kf == pdg_dau2[0])? 1 : 2]->Fill(ptau[nrec],ptGen);
            hyau[(kf == pdg_dau2[0])? 1 : 2]->Fill(yau[nrec]);
          }
        }
        recProbe[nrec] = *trw;
        nrec++;
      }
    }
    if(nrec < nbody)
      continue;

    if(nbody == 2)
      recProbe[0].PropagateToDCA(&recProbe[1]);
    if(nbody == 5){
      recProbe[0].PropagateToDCA(&recProbe[2]);
      recProbe[1].PropagateToDCA(&recProbe[2]);
    }

    parent = daurec[0];
    parentgen = daugen[0];
    for(int i = 1; i < nbody; ++i){
      parent += daurec[i];
      parentgen += daugen[i];
    }

    Double_t  ptRec=parent.Pt();
    Double_t  massRec=parent.M();
    Double_t yRec = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));

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
      hYEff->Fill(yRec);
      hPtEff->Fill(ptRec);
    }
    hMassVsPt->Fill(massRec,ptRec);
    hMassVsY->Fill(massRec,yRec);

    Float_t dca = 0, dca12 = 0, dca13 = 0, dca23 = 0;
    Double_t xP, yP, zP, sigmaVert, cts = 0;
    if(nbody == 2){
      if (ptau[0] > 0 && ptau[1] > 0) hyau2D->Fill(yau[0], yau[1]);
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
      // printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
      hDCA->Fill(dca, ptRec);
      hDCAx->Fill(d1, ptRec);
      hDCAy->Fill(d2, ptRec);
      hDCAz->Fill(d3, ptRec);

      ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
      cts = CosThetaStar(parent,daurec[0],pdgParticle,pdg_dau[0],pdg_dau[1]);
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
      // printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
      hDCA12->Fill(dca12, ptRec);
      hDCA13->Fill(dca13, ptRec);
      hDCA23->Fill(dca23, ptRec);
    }

    Double_t residVx=10000.*(xP - secvertgen[0][0]);
    Double_t residVy=10000.*(yP - secvertgen[1][0]);
    Double_t residVz=10000.*(zP - secvertgen[2][0]);
    hResVx->Fill(residVx, ptRec);
    hResVy->Fill(residVy, ptRec);
    hResVz->Fill(residVz, ptRec);
    hResVxVsY->Fill(residVx, yRec);
    hResVyVsY->Fill(residVy, yRec);
    hResVzVsY->Fill(residVz, yRec);
    
    hResPx->Fill(daurec[0].Px() - daugen[0].Px(), ptRec);
    hResPy->Fill(daurec[0].Py() - daugen[0].Py(), ptRec);
    hResPz->Fill(daurec[0].Pz() - daugen[0].Pz(), ptRec);
    hResPx->Fill(daurec[1].Px() - daugen[1].Px(), ptRec);
    hResPy->Fill(daurec[1].Py() - daugen[1].Py(), ptRec);
    hResPz->Fill(daurec[1].Pz() - daugen[1].Pz(), ptRec);
    hResPt->Fill(ptRec-ptGen,ptRec);

    hResPxVsY->Fill(daurec[0].Px() - daugen[0].Px(), yRec);
    hResPyVsY->Fill(daurec[0].Py() - daugen[0].Py(), yRec);
    hResPzVsY->Fill(daurec[0].Pz() - daugen[0].Pz(), yRec);
    hResPxVsY->Fill(daurec[1].Px() - daugen[1].Px(), yRec);
    hResPyVsY->Fill(daurec[1].Py() - daugen[1].Py(), yRec);
    hResPzVsY->Fill(daurec[1].Pz() - daugen[1].Pz(), yRec);
    
    // cout << "secvert generated Pion: " << secvertgenPi[0] << "  " << secvertgenPi[1] << "  " << secvertgenPi[2] << endl;
    // cout << "Reco Vert  Pion: " << xP << "  " << yP << "  " << zP << endl;
    
    Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
    Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
    Float_t distgen = TMath::Sqrt(secvertgen[0][0] * secvertgen[0][0] + secvertgen[1][0] * secvertgen[1][0] + secvertgen[2][0] * secvertgen[2][0]);
    Float_t distgenXY = TMath::Sqrt(secvertgen[0][0] * secvertgen[0][0] + secvertgen[1][0] * secvertgen[1][0]);
    // printf("dist = %f , distXY=%f , dx=%f, dy=%f, dz=%f z1=%f, z2=%f \n", dist, distXY, xP, yP, zP, recProbe[0].GetZ(), recProbe[1].GetZ());
    // printf("distgen = %f , distgenXY=%f \n", distgen, distgenXY);
    
    Double_t vsec[3] = {xP, yP, zP};
    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
    Double_t ipD = ImpParXY(vprim, vsec, parent);
    hCosp->Fill(cosp, ptRec);
    // printf(" ***** ***** cos point = %f \n", cosp);
    //if (cosp < -0.98)
    //    printf("SMALL COSPOINT");
    
    hResDist->Fill(dist - distgen, ptRec);
    hResDistXY->Fill(distXY - distgenXY, ptRec);
    
    hDistXY->Fill(distXY, ptRec);
    hDist->Fill(dist, ptRec);
    hDistxRec->Fill(xP);
    hDistyRec->Fill(yP);
    hDistxyRec->Fill(TMath::Sqrt(xP*xP+yP*yP));
    hDistzRec->Fill(zP);
    hDistTotRec->Fill(dist);
    hDistgenXY->Fill(distgenXY, ptRec);
    hDistgen->Fill(distgen, ptRec);
      
    //AliExternalTrackParam *track1 = (AliExternalTrackParam *)recProbe[0].GetTrack();
    //AliExternalTrackParam *track2 = (AliExternalTrackParam *)recProbe[1].GetTrack();

    Double_t d0x[3];
    Double_t d0y[3];
    Double_t d0xy[3] = {10000,10000,10000};
    Double_t ptdau[3] = {10000,10000,10000};

    for(int i = 0; i < nbody; i++){
      recProbe[0].PropagateToZBxByBz(0);
      d0x[i] = recProbe[i].GetX();
      d0y[i] = recProbe[i].GetY();
      d0xy[i] = TMath::Sqrt(d0x[i] * d0x[i] + d0y[i] * d0y[i]);
      if (d0x[i] < 0)
        d0xy[i] *= -1;
      for(int i = 0; i < nbody; ++i)
        hd0XY[i]->Fill(d0xy[i], ptRec);
      // printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
      //hd0XYprod->Fill(d0xy1 * d0xy2, ptRec);
      ptdau[i] = recProbe[i].GetTrack()->Pt();

    }
      
    arrsp[0] = massRec;
    arrsp[1] = ptRec;
    arrsp[2] = yRec;
    arrsp[3] = mass*dist/parent.P();
    arrsp[4] = cosp;
    arrsp[5] = TMath::Min(TMath::Min(ptdau[0],ptdau[1]),ptdau[2]);
    arrsp[6] = TMath::Min(TMath::Min(TMath::Abs(d0xy[0]),TMath::Abs(d0xy[1])),TMath::Abs(d0xy[2]));
    arrsp[7] = (nbody == 2) ? dca : (dca12+dca13+dca23)/3.;
    arrsp[8] = TMath::Abs(ipD);	    
    arrsp[9] = (nbody == 2) ? cts : sigmaVert;      
    hsp->Fill(arrsp);

    if (ntcand){
      arrnt[0] = massRec;
      arrnt[1] = ptRec;
      arrnt[2] = yRec;
      arrnt[3] = dist;
      arrnt[4] = cosp;
      arrnt[5] = d0xy[0];
      arrnt[6] = d0xy[1];
      arrnt[7] = (nbody == 2) ? d0xy[0] * d0xy[1] : sigmaVert;
      if(nbody == 2){
        arrnt[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
        arrnt[9] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
        arrnt[10] = dca; 
      }
      else{
        arrnt[8] = TMath::Min(TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
        arrnt[9] = TMath::Max(TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
        arrnt[10] = dca12;
        arrnt[11] = dca13;
        arrnt[12] = dca23;
      }
      ntcand->Fill(arrnt);
    }
  } //event loop
  
  fout->cd();  
  hMassAll->Write();
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

  hYEff->Divide(hYGen);
  hPtEff->Divide(hPtGen);
  for(int i = 1; i <= hPtEff->GetNbinsX(); i++){
    double eff = hPtEff->GetBinContent(i);
    double n = hPtGen->GetBinContent(i);
    hPtEff->SetBinError(i,TMath::Sqrt(eff*(1-eff)/n));
  }
  for(int i = 1; i <= hYEff->GetNbinsX(); i++){
    double eff = hYEff->GetBinContent(i);
    double n = hYGen->GetBinContent(i);
    hYEff->SetBinError(i,TMath::Sqrt(eff*(1-eff)/n));
  }
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
  hCtz->Write();
  hDist->Write();
  hDistgenXY->Write();
  hDistgen->Write();
  hCosp->Write();
  hResVx->Write();
  hResVy->Write();
  hResVz->Write();
  hResPx->Write();
  hResPy->Write();
  hResPz->Write();
  hResPt->Write();
  hResVxVsY->Write();
  hResVyVsY->Write();
  hResVzVsY->Write();
  hResPxVsY->Write();
  hResPyVsY->Write();
  hResPzVsY->Write();
  hResDist->Write();
  hResDistXY->Write();
  hChi2->Write();
  //hd0XYprod->Write();
  for(int i = 0; i < nbody; ++i)
    hd0XY[i]->Write();
  hNevents->Write();
  hsp->Write();
  if (ntcand){
    fnt->cd();
    ntcand->Write();
    ntgen->Write();
    fnt->Close();
  }
  
  TFile fout2(Form("DecayHistos%s.root",suffix.Data()), "RECREATE");
  hPtGen->Write();
  hYGen->Write();
  if(nbody == 2){
    hptau[0]->SetName("hptdauPos");
    hyau[0]->SetName("hydauPos");
    hptau[1]->SetName("hptdauNeg");
    hyau[1]->SetName("hydauNeg");
    hyau2D->Write();
  }
  else{
    hptau[0]->SetName(Form("hptdau_%i",(pdg_dau[0] == pdg_unstable_dau) ? pdg_dau[1] : pdg_dau[0]));
    hyau[0]->SetName(Form("hydau_%i",(pdg_dau[0] == pdg_unstable_dau) ? pdg_dau[1] : pdg_dau[0]));
    hptau[1]->SetName(Form("hptdau2_%i",pdg_dau2[0]));
    hyau[1]->SetName(Form("hydau2_%i",pdg_dau2[0]));
    hptau[2]->SetName(Form("hptdau2%i",pdg_dau2[1]));
    hyau[2]->SetName(Form("hydau2_%i",pdg_dau2[1]));
  }
  for(int i = 0; i < nbody; i++){
    hptau[i]->Write();
    hyau[i]->Write();
  }
  fout2.Close();
  fout->Close();
  std::cout<<count1<<std::endl;
  std::cout<<count2<<std::endl;
  std::cout<<count3<<std::endl;
  std::cout<<count1+count2+count3<<std::endl;

}




void MakeCombinBkgCandidates3Body(const char* trackTreeFile="treeBkgEvents.root",
             TString suffix = "_Omega",
			       Int_t nevents = 999999, 
			       Int_t writeNtuple = kFALSE,
			       Bool_t usePID=kFALSE,
             int pdgMother = 3334,
             int pdg_unstable_dau = 3122){

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
  THnSparseF *hsp = CreateSparse(mass, 3);
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntcand = 0x0;
  Float_t arrnt[13];
  if (writeNtuple){
    fnt = new TFile(Form("fntBkg%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "mass:pt:y:dist:cosp:d01:d02:sigmavert:ptMin:ptMax:dca12:dca13:dca23", 32000);
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

              arrsp[0] = invMass;
              arrsp[1] = pt;
              arrsp[2] = y;
              arrsp[3] = mass*dist/parent.P();
              arrsp[4] = cosp;
              arrsp[5] = TMath::Min(TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
              arrsp[6] = TMath::Min(TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2)),TMath::Abs(d0xy3));
              arrsp[7] = (dca12+dca13+dca23)/3.;
              arrsp[8] = TMath::Abs(ipD);	    
              arrsp[9] = sigmaVert;      
              hsp->Fill(arrsp);

              if (ntcand){
                arrnt[0] = invMass;
                arrnt[1] = pt;
                arrnt[2] = y;
                arrnt[3] = dist;
                arrnt[4] = cosp;
                arrnt[5] = d0xy1;
                arrnt[6] = d0xy2;
                arrnt[7] = sigmaVert;
                arrnt[8] = TMath::Min(TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
                arrnt[9] = TMath::Max(TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt()),recProbe[2].GetTrack()->Pt());
                arrnt[10] = dca12;
                arrnt[11] = dca13;
                arrnt[12] = dca23;
                
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
  hsp->Write();
  fout->Close();
  if (ntcand){
    fnt->cd();
    ntcand->Write();
    fnt->Close();
  }
}

void MakeCombinBkgCandidates2Body(const char* trackTreeFile="treeBkgEvents.root",
             TString suffix = "_phi",
			       Int_t nevents = 999999, 
			       Int_t writeNtuple = kTRUE,
             int pdgMother = 333){

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
  int pdg_dau[2] = {0,0};
  GetPDGDaughters(pdgMother,pdg_dau,pdgMother > 0);

  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  
  TFile *fout = new TFile(Form("Bkg-histos%s.root",suffix.Data()), "recreate");
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
    
  Double_t massM = TDatabasePDG::Instance()->GetParticle(pdgMother)->Mass();
  THnSparseF *hsp = CreateSparse(massM, 2);
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
          Float_t pt=parent.Pt();
          Float_t invMass=parent.M();
          Float_t y = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
          hYPtRecoAll->Fill(y, pt);
          hPtRecoAll->Fill(pt);
          hMassAll->Fill(invMass);
          if(invMass>0.5*massM  && invMass<1.5*massM){
            // range to fill histos
            if(invMass>0.8*massM && invMass<1.2*massM) countCandInPeak++;
            Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
            Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
            Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
            Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
            
            //printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
            hDCA->Fill(dca, pt);
            hDCAx->Fill(d1, pt);
            hDCAy->Fill(d2, pt);
            hDCAz->Fill(d3, pt);
            Double_t xP, yP, zP;
            ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
            Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
            Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
            Double_t vsec[3] = {xP, yP, zP};
            Double_t cosp = CosPointingAngle(vprim, vsec, parent);
            Double_t cts = CosThetaStar(parent,daurec[iNeg],pdgMother,pdg_dau[0],pdg_dau[1]);
            Double_t ipD = ImpParXY(vprim, vsec, parent);
            hCosp->Fill(cosp, pt);
            hCosThStVsMass->Fill(invMass,cts);
            //printf(" ***** ***** cos point = %f \n", cosp);	    
            hDistXY->Fill(distXY, pt);
            hDist->Fill(dist, pt);
            
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
            
            ////hd0XYprod->Fill(d0xy1 * d0xy2, pt);
            hd0XY1->Fill(d0xy1, pt);
            hd0XY2->Fill(d0xy2, pt);
            
            arrsp[0] = invMass;
            arrsp[1] = pt;
            arrsp[2] = y;
            arrsp[3] = dist*massM*parent.P(),
            arrsp[4] = cosp;
            arrsp[5] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
            arrsp[6] = TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2));
            arrsp[7] = dca;
            arrsp[8] = TMath::Abs(ipD);	    
            arrsp[9] = cts;
            hsp->Fill(arrsp);
            
            if (ntcand){
              arrnt[0] = invMass;
              arrnt[1] = pt;
              arrnt[2] = y;
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
  //hd0XYprod->Write();
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