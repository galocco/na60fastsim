#ifndef MEASUREMENT_UTILS_H
#define MEASUREMENT_UTILS_H
#endif

#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
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
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliAODRecoDecay.h"
#include "AliDecayer.h"
#include "AliDecayerEvtGen.h"
#include "TDatabasePDG.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

// settings for signal generation
// rapidity range
const double yminSG = -2.; 
const double ymaxSG = 8.;
// pT range
const double ptminSG = 0;
const double ptmaxSG = 10; 

//number of energy of the colliding beam
const int NEnergy = 5;
//number of particles 
const int NParticles = 5;
//array of the energy
double Elab[NEnergy] = {20,30,40,80,158};

//T parameter in the exponential pT distribution [matter/antimatter][particle][beam energy] i
//                                                    phi                        K               Lambda                Omega                    Xi
double Tslope[2][NParticles][NEnergy] = {{{196.8,237.4,244.6,239.8,298.7},{0,0,229, 223.1, 229},{244,249,258,265,301},{0,0,218,0,267},{221,233,222,227,277}},//matter e rapidity distribution [matter/antimatter][particle][beam energy]
                                         {{196.8,237.4,244.6,239.8,298.7},{0,0,226,217,226},{339,284,301,292,303},{0,0,218,0,259},{311,277,255,321,0}}};//antimatter
//sigma parameter of the gaussians of th                    phi                        K                       Lambda                Omega                    Xi
double sigma_rapidity[2][NParticles][NEnergy] = {{{0.425,0.538,0.696,0.658,1.451},{0,0,0.725,0.792,0.88},{0.51,0.66,0.91,0.87,0},{0,0,0.6,0,1.2},{0.45,0.56,0.76,0.71,1.18}},//matter
                                                 {{0.425,0.538,0.696,0.658,1.451},{0,0,0.635,0.705,0.81},{0.62,0.69,0.77,0.83,1.00},{0,0,0.6,0,1.0},{0,0.76,0.65,0.87,0.73}}};//antimatter
//mu paramter of the gaussian of the rapidity distribution [matter/antimatter][particle][beam energy]
//                                                       phi                        K                   Lambda                   Omega                    Xi
double y0_rapidity[2][NParticles][NEnergy] = {{{0.425,0.538,0.487,0.682,0.},{0,0,0.694,0.742,0.839},{0.49,0.59,0.65,0.94,0},{0,0,0,0,0},{0.45,0.47,0.54,0.68,0}},//matter
                                              {{0.425,0.538,0.487,0.682,0.},{0,0,0.569,0.668,0.727},{0,0,0,0,0},{0,0,0.0,0},{0,0,0,0,0}}};//antimatter
//multiplicity for event [matter/antimatter][particle][beam energy]
//                                                       phi                        K                   Lambda                   Omega                    Xi
double multiplicity[2][NParticles][NEnergy] = {{{1.89,1.84,2.55,4.04,8.46},{0,0,59.1,76.9,103.0},{27.1,36.9,43.1,50.1,44.9},{0,0,0.14,0,0.43},{1.50,2.42,2.96,3.80,4.04}},//matter
                                               {{1.89,1.84,2.55,4.04,8.46},{0,0,19.2,32.4,51.9},{0.16,0.39,0.68,1.82,3.07},{0,0,0.14,0,0.19},{0,0.12,0.13,0.58,0.66}}};//antimatter
//branching ratio [particle]
//                            phi     K   Lambda  Omega   Xi
double bratio[NParticles] = {0.489, 0.692, 0.639, 0.433, 0.638};//for Omega (Xi) it's given by 0.678*0.639 (0.999*0.639)
//name of the particles
TString particle_name[NParticles] = {"phi","K0s","Lambda","Omega","Xi"};
TString all_particle_name[NParticles+3] = {"pion","kaon","proton","phi","K0s","Lambda","Omega","Xi"};
//pdg code of the particles
double pdg_mother[NParticles] = {333,310,3122,3334,3312};
//pdg code of the daughters
//                                     K+    K-   pi+  pi-    p    pi-  Lambda  K-  Lambda pi-
double pdg_daughter[NParticles][2] = {{321,-321},{211,-211},{2212,-211},{3122,-321},{3122,-211}};


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

Double_t ArmPod(TVector3 plus, TVector3 minus, TVector3 mom, float &qT){
  float qlp = TMath::Abs(mom.Dot(plus));
  float qlm = TMath::Abs(mom.Dot(minus));
  double alpha = (qlp-qlm)/(qlp+qlm);
  qT = TMath::Abs(plus.Cross(mom).Mag())/mom.Mag();
  return alpha;
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
  if(TMath::Abs(pdgParticle)==310)
    return (Tslope[0][index_pdg][index_E]+Tslope[1][index_pdg][index_E])/2.;
  else
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
  
  if(TMath::Abs(pdgParticle)==310)
    return (sigma_rapidity[0][index_pdg][index_E]+sigma_rapidity[1][index_pdg][index_E])/2.;
  else
    return sigma_rapidity[matter ? 0 : 1][index_pdg][index_E];
}

double GetBRatio(int pdgParticle){
  int index_pdg = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  return bratio[index_pdg];
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
  if(TMath::Abs(pdgParticle)==310)
    return (y0_rapidity[0][index_pdg][index_E]+y0_rapidity[1][index_pdg][index_E])/2.;
  else
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
  
  if(TMath::Abs(pdgParticle)==310)
    return (multiplicity[0][index_pdg][index_E]+multiplicity[1][index_pdg][index_E])/2.;
  else
    return multiplicity[matter ? 0 : 1][index_pdg][index_E];
}

double GetY0(double Eint){
  double y0BG = 2.9;
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

int GetNBody(int pdgParticle, int& pdg_unstable_dau, bool matter = true){
  pdgParticle = TMath::Abs(pdgParticle);
  int pdgDaughters[2];
  GetPDGDaughters(pdgParticle, pdgDaughters, true);
  int sign = matter ? 1 : -1;
  for(auto& pdg: pdg_mother){
    if(pdgDaughters[0] == pdg){
      pdg_unstable_dau = pdg*sign;
      return 3;
    }
    else if(pdgDaughters[1] == pdg){
      pdg_unstable_dau = pdg*sign;
      return 3;
    }
  }

  return 2;
}

const double sigmaTOF = 20;//ps (10^(-12))
const double lTOF = 0.15;//m

double sigmaBeta(double zTOF){
  double sigmas[7] = {0.020, 0.017, 0.014, 0.012, 0.010, 0.009, 0.008};
  if(zTOF == 0.4)
    return sigmas[0];
  if(zTOF == 0.5)
    return sigmas[1];
  if(zTOF == 0.6)
    return sigmas[2];
  if(zTOF == 0.7)
    return sigmas[3];
  if(zTOF == 0.8)
    return sigmas[4];
  if(zTOF == 0.9)
    return sigmas[5];
  if(zTOF == 1.0)
    return sigmas[6];
  return sigmas[0];
}

double minT(double zTOF){
  double sigmas[7] = {1200, 1200, 1200, 1200, 1200, 1200, 3200};
  if(zTOF == 0.4)
    return sigmas[0];
  if(zTOF == 0.5)
    return sigmas[1];
  if(zTOF == 0.6)
    return sigmas[2];
  if(zTOF == 0.7)
    return sigmas[3];
  if(zTOF == 0.8)
    return sigmas[4];
  if(zTOF == 0.9)
    return sigmas[5];
  if(zTOF == 1.0)
    return sigmas[6];
  return sigmas[0];
}

double maxT(double zTOF){
  double sigmas[7] = {1600, 3800, 3800, 3800, 3800, 3800, 3800};
  if(zTOF == 0.4)
    return sigmas[0];
  if(zTOF == 0.5)
    return sigmas[1];
  if(zTOF == 0.6)
    return sigmas[2];
  if(zTOF == 0.7)
    return sigmas[3];
  if(zTOF == 0.8)
    return sigmas[4];
  if(zTOF == 0.9)
    return sigmas[5];
  if(zTOF == 1.0)
    return sigmas[6];
  return sigmas[0];
}

double betaGen(double p[], double mass){
  double ptot2 = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
  return TMath::Sqrt(ptot2/(ptot2+mass*mass));
}

double hasTOF(double p[], double mass, double &L, double zTOF){
  double speed = betaGen(p, mass)*TMath::C();
  double ptot = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  double v[3] = {p[0]/ptot*speed, p[1]/ptot*speed, p[2]/ptot*speed};
  double t = zTOF/v[2];
  L = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*t;
  return t*TMath::Power(10,12);

  if(TMath::Abs(v[0]*t) > lTOF && TMath::Abs(v[1]*t) > lTOF){
    L = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*t;
    return t*TMath::Power(10,12);
  }
  else{
    L = -1;
    return -1;
  }
}

double smearT(double dt){
  return dt+gRandom->Gaus(0,sigmaTOF);
}

double betaTOF(double L, double dt){
  return L/(TMath::C()*dt*TMath::Power(10,-12));
}

int IsTrueCandidate2Body(KMCProbeFwd tr[], int pdgMother){
  if (tr[0].GetIndexMom() != tr[1].GetIndexMom()) return 0;
  for(int i=0; i<2; i++){
    if (TMath::Abs(tr[i].GetPdgMother()) != TMath::Abs(pdgMother)) return 0;
  }
  return 1;
}

int IsTrueCandidate3Body(KMCProbeFwd tr[], int pdgMother, int pdgguess[]){
  if (tr[2].GetPdg() != pdgguess[2]) return 0;
  if (tr[1].GetPdg() != pdgguess[0]) return 0;
  if (tr[0].GetPdg() != pdgguess[1]) return 0;

  for(int i=0; i<3; i++){
    if (TMath::Abs(tr[i].GetPdgMother()) != TMath::Abs(pdgMother)) return 0;
    if (i < 2)
      if (tr[i].GetIndexMom() != tr[i+1].GetIndexMom()) return 0;
  }
  return 1;
}

bool SetEvent(Int_t nevents, Int_t &iev, TTree* tree, TClonesArray &aarrtr, Double_t &T0, Int_t &gen_part, Double_t multi_times_br, Int_t &icount){
  while(gen_part == 0)
  {
    iev++;
    if(iev%100==0) printf(" ***************  ev = %d of %d \n", iev, nevents);
    if(iev>nevents) return false;
    icount = 0;
    tree->Fill();
    aarrtr.Clear();
    gen_part  = gRandom->Poisson(multi_times_br);
    T0 = smearT(0);
  }
  return true;
}