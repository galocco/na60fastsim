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
#include "TStopwatch.h"
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;

double vX = 0, vY = 0, vZ = 0; // event vertex

void MakeEventMixing(Int_t n_next = 1,
                    const char* trackTreeFileBkg="treeBkgEvents_layer5.root",
                    const char* trackTreeFileSig="treeBkgEvents_layer5.root",
                    TString suffix = "_PHI_test",
                    int pdgMother = 333,
                    const char *setup = "../setups/setup-EHN1_BetheBloch.txt",
                    Int_t fullsim = kFALSE,
                    Int_t writeNtuple = kTRUE,
                    Int_t writeSparse = kTRUE){
  TStopwatch timer;
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
  TFile *fhsp = 0x0;
  TNtuple *ntcand = 0x0;
  const Int_t nAxes = 10;
  Double_t arrhsp[nAxes];
  Float_t  arrnt[nAxes];

  if (writeNtuple){
    fnt = new TFile(Form("fntMixed%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "m:pt:rapidity:dist:cosp:d0prod:dca:ptMin:ptMax:thetad", 32000);
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
                          "thetad"};
    //                        m            pt    y    d  cosp     d0   dca       pti pta  the  ns1 ns2
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
  bool verbose = false;
  for (Int_t iev = 0; iev < nevents-n_next; iev++){
    hCentrality->Fill(1);
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    treeBkg->GetEvent(iev);
    TClonesArray arrBkgCopy = *arrBkg;
    Int_t arrentrBkg = arrBkg->GetEntriesFast();
    TClonesArray arrSigCopy;
    if(fullsim){
      treeSig->GetEvent(iev);
      arrSigCopy = *arrSig;
    }
    Int_t arrentrSig = (fullsim) ? arrSig->GetEntriesFast() : 0;

    
    for (Int_t itr = 0; itr < arrentrBkg+arrentrSig; itr++){
      KMCProbeFwd *tr1;
      if(itr < arrentrBkg)
        tr1 = (KMCProbeFwd *)arrBkgCopy->At(itr);
      else
        tr1 = (KMCProbeFwd *)arrSigCopy->At(itr-arrentrBkg);
      Float_t ch1 = tr1->GetCharge();
      for(int inext=1; inext <= n_next; inext++){
        treeBkg->GetEvent(iev+inext);
        Int_t arrentrBkgNext = arrBkg->GetEntriesFast();
        if(fullsim){
          treeSig->GetEvent(iev+inext);
          arrNextSigCopy = *arrSig;
        }
        Int_t arrentrSigNext = (fullsim) ? arrSig->GetEntriesFast() : 0;

        for (Int_t itr2 = itr; itr2 < arrentrBkgNext+arrentrSigNext; itr2++){ 
          KMCProbeFwd *tr2;
          if(itr2 < arrentrBkgNext)
            tr2 = (KMCProbeFwd *)arrBkg->At(itr2);
          else
            tr2 = (KMCProbeFwd *)arrSig->At(itr2-arrentrBkgNext);
          Float_t ch2 = tr2->GetCharge(); 
          if (ch1*ch2 > 0) continue;
          if (ch1 < 0){ //convention: first track negative
            recProbe[0] = *tr1;
            recProbe[1] = *tr2;
          }else if (ch2 < 0){
            recProbe[0] = *tr2;
            recProbe[1] = *tr1;
          }

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
            if(invMass>0.96*massM  && invMass<1.04*massM){
              // range to fill histos
              if(invMass>0.96*massM && invMass<1.04*massM) countCandInPeak++;
              Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
              Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
              Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
              Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
              
              //if(dca>0.018) continue;
              hDCA->Fill(dca, pt);
              hDCAx->Fill(d1, pt);
              hDCAy->Fill(d2, pt);
              hDCAz->Fill(d3, pt);

              Double_t xP = 0, yP = 0, zP = 0;
              ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP); 
              Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
              //if(dist > 1) continue;
              Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
              Double_t vsec[3] = {xP, yP, zP};
              Double_t thetad = OpeningAngle(daurec[0],daurec[1]);
              //PRESELECTION: dist<1 and  and thetad > 0.97 and arm < 0.8 and arm > -0.8
              //if(TMath::ACos(thetad)>0.15 || TMath::ACos(thetad)<0.005) continue;
              Double_t arm = ArmenterosAlpha(daurec[0],daurec[1]);
              //if(arm < -0.8 || arm > 0.8) continue;
              //Double_t ipD = ImpParXY(vprim, vsec, parent);
              //Double_t cts = CosThetaStar(parent,daurec[iNeg],pdgMother,pdg_dau[0],pdg_dau[1]);
              Double_t cosp = CosPointingAngle(vprim, vsec, parent);

              hCosp->Fill(cosp, pt);
              //hCosThStVsMass->Fill(invMass,cts);

              hDistXY->Fill(distXY, pt);
              hDistXYPlane->Fill(distXY, zP);
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
              
              //if(d0xy1 * d0xy2 > 0.1 || d0xy1 * d0xy2 < -0.1 ) continue;
              hd0XY1->Fill(d0xy1, pt);
              hd0XY2->Fill(d0xy2, pt);

              if (writeNtuple){
                arrnt[0] = invMass;
                arrnt[1] = pt;
                arrnt[2] = y;
                arrnt[3] = dist;
                //arrnt[5] = dist*massM/TMath::Sqrt(parent.Pt() * parent.Pt() + parent.Pz() * parent.Pz());
                arrnt[4] = cosp;
                arrnt[5] = d0xy1 * d0xy2;
                arrnt[6] = dca;
                arrnt[7] = TMath::Min(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
                arrnt[8] = TMath::Max(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
                arrnt[9] = TMath::ACos(thetad);
                ntcand->Fill(arrnt);
              }
              if (writeSparse){
                arrhsp[0] = invMass;
                arrhsp[1] = pt;
                arrhsp[2] = y;
                arrhsp[3] = dist;
                //arrhsp[5] = dist*massM/TMath::Sqrt(parent.Pt() * parent.Pt() + parent.Pz() * parent.Pz());
                arrhsp[4] = cosp;
                arrhsp[5] = d0xy1 * d0xy2;
                arrhsp[6] = dca;
                arrhsp[7] = TMath::Min(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
                arrhsp[8] = TMath::Max(recProbe[0].GetTrack()->Pt(), recProbe[1].GetTrack()->Pt());
                arrhsp[9] = TMath::ACos(thetad);
                hsp->Fill(arrhsp);
              }
        
            } // check on inv mass
            if(verbose) std::cout<<"ends mass loop\n";
          } // loop on mass hypothesis
        } // loop on first track 
      
      }
      if(verbose) std::cout<<"inner-most loop ended, ev "<<iev<<"\n";

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
  if (writeNtuple){
    fnt->cd();
    ntcand->Write();
    fnt->Close();
  }
  if (writeSparse){
    fhsp->cd();
    hsp->Write();
    fhsp->Close();
  }
  
}
