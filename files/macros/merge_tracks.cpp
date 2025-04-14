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
#include "TClonesArray.h"
#endif

void MergeTreeTracks(TString trackTreeFile1="treeStrangeParticles.root",
                     TString trackTreeFile2="treeStrangeParticles.root",
                     TString suffix = "_PHI_test",
                     Int_t clean = kFALSE){

  TFile *tracks_file = new TFile(Form("treeStrangeParticles%s.root",suffix.Data()), "RECREATE");
  TTree *tree = new TTree("tree", "tree Sig");
  TClonesArray *arrtr = new TClonesArray("KMCProbeFwd");
  TClonesArray &aarrtr = *arrtr;
  tree->Branch("tracks", &arrtr);
  TFile *filetree1 = new TFile(trackTreeFile1.Data());
  TTree *tree1 = (TTree *)filetree1->Get("tree");
  TClonesArray *arr1 = 0;
  tree1->SetBranchAddress("tracks", &arr1);
  // signal tracks
  TFile *filetree2 = new TFile(trackTreeFile2.Data());
  TTree *tree2 = (TTree *)filetree2->Get("tree");
  TClonesArray *arr2 = 0;
  tree2->SetBranchAddress("tracks", &arr2);

  Int_t nevents = (tree1->GetEntries() < tree2->GetEntries()) ? tree1->GetEntries(): tree2->GetEntries();
  std::cout<<"nevents: "<<nevents<<"\n";
  for (Int_t iev = 0; iev < nevents; iev++){
    tree1->GetEvent(iev);
    Int_t arrentr1 = arr1->GetEntriesFast();
    tree2->GetEvent(iev);
    Int_t arrentr2 = arr2->GetEntriesFast();
    for (Int_t itr = 0; itr < arrentr1+arrentr2; itr++){
        KMCProbeFwd *trw;
        if(itr<arrentr1)
            trw = (KMCProbeFwd *)arr1->At(itr);
        else
            trw = (KMCProbeFwd *)arr2->At(itr-arrentr1);
        new (aarrtr[itr]) KMCProbeFwd(*trw);
    } // loop on second track
    tree->Fill();
    aarrtr.Clear();
    printf(" --> Event %d\n",iev);
  }
  filetree1->Close();
  filetree2->Close();
  tracks_file->cd();
  tree->Write();
  tracks_file->Close();
  if(clean){
    std::cout<<"deleting the file "<<trackTreeFile1.Data()<<"\n";
	  gSystem->Exec(Form("rm %s", trackTreeFile1.Data()));
  }
}
