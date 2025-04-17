#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
#include "KMCProbeFwd.h"
#include "KMCDetectorFwd.h"
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
#include "GenMUONLMR.h"
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
                    const char *setup = "../setups/setup_proposal_40GeV.txt",
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
  
  TFile *filetreeBkgNext = new TFile(trackTreeFileBkg);
  TTree *treeBkgNext = (TTree *)filetreeBkgNext->Get("tree");
  TClonesArray *arrBkgNext = 0;
  treeBkgNext->SetBranchAddress("tracks", &arrBkgNext);
  // signal tracks
  TFile *filetreeSigNext = new TFile(trackTreeFileSig);
  TTree *treeSigNext = (TTree *)filetreeSigNext->Get("tree");
  TClonesArray *arrSigNext = 0;
  treeSigNext->SetBranchAddress("tracks", &arrSigNext);

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
  TH2F *hYPtReco = new TH2F("hYPtReco", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);
  TH1D *hMass = new TH1D("hMass", "Mass all match", 250, 0., 2.5);
  TH2F *hDistXY = new TH2F("hDistXY", ";d_{xy} (cm); #it{p}_{T} (GeV/#it{c}) ;counts", 100, 0, 0.1, 100, 0, 3);
  TH2F *hDistXYPlane = new TH2F("hDistXYPlane", ";x (cm); y (cm);counts", 100, -0.004, 0.004, 100, -0.0000001, 0.00000001);
  TH2F *hDistVsPt = new TH2F("hDistVsPt", ";z (cm); #it{p}_{T} (GeV/#it{c});counts", 300, 0, 3, 30, 0, 3);
  TH2F *hDistzVsPt = new TH2F("hDistzVsPt", ";z (cm); #it{p}_{T} (GeV/#it{c});counts", 300, 0, 3, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", "", 100, -0.01, 0.01, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 100, -1, 1, 30, 0, 3);
  TH2F *hCosThStVsMass = new TH2F("hCosThStVsMass", "", 50, 1.5, 2.5, 40, -1, 1);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
    
  TFile *fnt = 0x0;
  TFile *fhsp = 0x0;
  TNtuple *ntcand = 0x0;
  const Int_t nAxes = 7;
  Double_t arrhsp[nAxes];
  Float_t  arrnt[12];

  if (writeNtuple){
    fnt = new TFile(Form("fntMixed%s.root",suffix.Data()), "recreate");
    ntcand = new TNtuple("ntcand", "ntcand", "m:pt:rapidity:dist:cosp:d0prod:dca:ptMin:ptMax:thetad:histp:hitsm", 32000);
  }

  THnSparse* hsp = 0x0;
  if (writeSparse){
    fhsp = new TFile(Form("fhspBkg%s.root",suffix.Data()), "recreate");
    TString axTit[nAxes]={"m",
                      "pt",
                      "dist",
                      "d0prod",
                      "dca"
                      "hitsp",
                      "hitsm"};
    Int_t binpt = (ptmaxSG-ptminSG)*10;
    //                        m            pt    y    d         d0     dca    pti pta  the 
    Int_t bins[nAxes] =   {120 ,    binpt,  20,       30,   10, 7, 7}; 
    Double_t min[nAxes] = {0.98,  ptminSG,  0., -0.00015,    0, 3.5, 10.5};
    Double_t max[nAxes] = {1.10,  ptmaxSG,  1.,  0.00015,  0.1, 3.5, 10.5};  
    hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
    for(Int_t iax=0; iax<nAxes; iax++)
      hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  }

  //read the pdg code of the daughters
  int pdg_dau[2] = {0,0};
  GetPDGDaughters(pdgMother,pdg_dau,true);
  KMCProbeFwd recProbe[2];
  TLorentzVector parent, daurec[2];
  KMCProbeFwd *tr1;
  KMCProbeFwd *tr2;
  for (Int_t iev = 0; iev < nevents-n_next; ++iev){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    treeBkg->GetEvent(iev);
    Int_t arrentrBkg = arrBkg->GetEntriesFast();
    if(fullsim)
      treeSig->GetEvent(iev);
    Int_t arrentrSig = (fullsim) ? arrSig->GetEntriesFast() : 0;

    for(int inext=1; inext <= n_next; ++inext){
      treeBkgNext->GetEvent(iev+inext);
      Int_t arrentrBkgNext = arrBkgNext->GetEntriesFast();
      if(fullsim)
        treeSigNext->GetEvent(iev+inext);

      Int_t arrentrSigNext = (fullsim) ? arrSigNext->GetEntriesFast() : 0;
      for (Int_t itr = 0; itr < arrentrBkg+arrentrSig; ++itr){
        if(itr < arrentrBkg)
          tr1 = (KMCProbeFwd *)arrBkg->At(itr);
        else
          tr1 = (KMCProbeFwd *)arrSig->At(itr-arrentrBkg);
        Float_t ch1 = tr1->GetCharge();
        for (Int_t itr2 = itr; itr2 < arrentrBkgNext+arrentrSigNext; ++itr2){
          if(itr2 < arrentrBkgNext)
            tr2 = (KMCProbeFwd *)arrBkgNext->At(itr2);
          else
            tr2 = (KMCProbeFwd *)arrSigNext->At(itr2-arrentrBkgNext);
          Float_t ch2 = tr2->GetCharge(); 
          if (ch1*ch2 > 0) continue;
          if (ch1 < 0){ //convention: first track negative
            recProbe[0] = *tr1;
            recProbe[1] = *tr2;
          }else if (ch2 < 0){
            recProbe[0] = *tr2;
            recProbe[1] = *tr1;
          }

          countCand++;
          Double_t pxyz_pos[3] = {0, 0, 0};
          Double_t pxyz_neg[3] = {0, 0, 0};

          recProbe[0].PropagateToDCA(&recProbe[1]);
          recProbe[0].GetPXYZ(pxyz_pos);
          recProbe[1].GetPXYZ(pxyz_neg);

          daurec[0].SetXYZM(pxyz_pos[0], pxyz_pos[1], pxyz_pos[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[1])->Mass());
          daurec[1].SetXYZM(pxyz_neg[0], pxyz_neg[1], pxyz_neg[2], TDatabasePDG::Instance()->GetParticle(pdg_dau[0])->Mass());


          parent = daurec[0];
          parent += daurec[1];
          Float_t invMass = parent.M();
          if(invMass<0.98 && invMass>1.1) continue;

          Double_t vsec[3] = {0,0,0};
          ComputeVertex(recProbe[0], recProbe[1], vsec[0], vsec[1], vsec[2]); 
          Float_t dist = TMath::Sqrt(vsec[0] * vsec[0] + vsec[1] * vsec[1] + vsec[2] * vsec[2]);
          if(dist > 0.1) continue;
          
          Float_t pt = parent.Pt();
          Float_t y = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
          Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
          Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
          Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
          Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
          
          hDCA->Fill(dca, pt);
          hDCAx->Fill(d1, pt);
          hDCAy->Fill(d2, pt);
          hDCAz->Fill(d3, pt);

          Float_t distXY = TMath::Sqrt(vsec[0] * vsec[0] + vsec[1] * vsec[1]);
          Double_t thetad = OpeningAngle(daurec[0],daurec[1]);
          //Double_t arm = ArmenterosAlpha(daurec[0],daurec[1], parent);
          Double_t cosp = CosPointingAngle(vprim, vsec, parent);

          
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
          
          hd0XY1->Fill(d0xy1, pt);
          hd0XY2->Fill(d0xy2, pt);
          hCosp->Fill(cosp, pt);
          hCosThStVsMass->Fill(invMass,thetad);
          hDistXY->Fill(distXY, pt);
          hDistXYPlane->Fill(distXY, vsec[2]);
          hDistVsPt->Fill(dist, pt);
          hDistzVsPt->Fill(vsec[2], pt);
          hYPtReco->Fill(y, pt);
          hMass->Fill(invMass);
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
            arrnt[10] = recProbe[1].GetNITSHits();
            arrnt[11] = recProbe[0].GetNITSHits();
            ntcand->Fill(arrnt);
          }

          if (writeSparse){
            arrhsp[0] = invMass;
            arrhsp[1] = pt;
            arrhsp[2] = dist;
            arrhsp[3] = d0xy1 * d0xy2;
            arrhsp[4] = dca;
            arrhsp[5] = recProbe[1].GetNITSHits();
            arrhsp[6] = recProbe[0].GetNITSHits();
            hsp->Fill(arrhsp);
          }

        } // loop on first track 
      
      }// loop on second track
    } // loop on the next events
    
    printf(" --> Event %d, tot candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hNevents->Write();
  hMass->Write();
  hYPtReco->Write();
  hDistXY->Write();
  hDistVsPt->Write();
  hDistXYPlane->Write();
  hDistzVsPt->Write();
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
