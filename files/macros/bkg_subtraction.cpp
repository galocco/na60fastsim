#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TF1.h>
#include <TNtuple.h>
#include "TRandom.h"
#include "THnSparse.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TDatabasePDG.h"
#include "TMVA/HyperParameterOptimisation.h"
#include "TMVA/MethodBase.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TStopwatch.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TDirectory.h"
#endif
float eff_presel = 0.528;//pt > 0.9 0.191
float multiplicity = 2.55;
float bratio = 0.489;

TH2D* fill_hist(TString filename, Int_t pdg_code, Int_t nevents, const TString branch);
Int_t normalize_evt_mix(TH2D* hist_mix, TH2D* hist_bkg);
void save_projections(TH2D* hist, TDirectory* dir);
void save_and_fit_projections(TH2D* hist, TDirectory* dir);
int evt_mix_study(TString cmb_bkg_file="fntBkg_PHI_L5_E40_Data1.root",
                  TString evt_mix_file="fntMixed_PHI_evt_mix_L5_E40.root",
                  TString events_file="Bkg-histos_PHI_L5_E40_Data1.root",
                  TString sig_file="fntSig_PHI_L5_E40_Data.root",
                  Int_t   pdg_code=333,
                  TString suffix="",
                  Int_t fullsim = kFALSE
                  ){
    TStopwatch timer;
    TH2D* hist_bkg;
    TH2D* hist_mix;
    TH2D* hist_sig;
    TH2D* hist_gen;
    //load the background and the event mixing
    printf("loading the combinatorial bgk:\n");
    hist_bkg = fill_hist(cmb_bkg_file, pdg_code, 0, "ntcand");
    hist_bkg->SetName("hist_bkg");
    printf("loading the event-mixing bgk:\n");
    hist_mix = fill_hist(evt_mix_file, pdg_code, 0, "ntcand");
    hist_mix->SetName("hist_mix");
    normalize_evt_mix(hist_mix, hist_bkg);
    //read the number of events
    TFile *events = new TFile(events_file.Data(),"read");
    TH1D* hist_ev = (TH1D*) events->Get("hNevents");
    //inputFile->GetObject("hscore_bdt", hscore_data);
    Int_t n_ev = hist_ev->GetBinContent(1);
    Int_t nrec = eff_presel*multiplicity*n_ev*bratio;
    printf("loading the signal:\n");
    hist_sig = fill_hist(sig_file, pdg_code, nrec, "ntcand");
    hist_sig->SetName("hist_sig");

    hist_gen = fill_hist(sig_file, pdg_code, nrec, "ntgen");
    hist_gen->SetName("hist_gen");
    TH1D* hist_eff = hist_sig->ProjectionY("hist_eff");
    TH1D* hist_gen_1d = hist_gen->ProjectionY("hist_gen_1d");
    for(int i = 1; i <= hist_eff->GetNbinsX(); i++){
        double rec = hist_eff->GetBinContent(i);
        double gen = hist_gen_1d->GetBinContent(i);
        if(gen<rec){
            hist_eff->SetBinContent(i, 1);
            hist_eff->SetBinError(i, 0);
        }
        else if(gen>0){
            double eff = rec/gen;
            hist_eff->SetBinContent(i,eff);
            hist_eff->SetBinError(i,TMath::Sqrt(eff*(1-eff)/gen));
        }
        else{
            hist_eff->SetBinContent(i, 0);
            hist_eff->SetBinError(i, 0);
        }
        
    }
    TH2D* hist_all = (TH2D*) hist_bkg->Clone();
    hist_all->SetName("hist_all");
    
    if(!fullsim) hist_all->Add(hist_sig);
    TH2D* hist_sub = (TH2D*) hist_all->Clone();
    hist_sub->SetName("hist_sub");
    hist_sub->Add(hist_mix, -1);

    TH2D* hist_ratio = (TH2D*) hist_all->Clone();
    hist_ratio->SetName("hist_ratio");
    hist_ratio->Divide(hist_mix);

    TFile *output = new TFile(Form("event_mixing_phi_%s.root",suffix.Data()), "recreate");
    hist_ev->Write();
    hist_mix->Write();
    if(!fullsim){
        hist_bkg->Write();
        hist_sig->Write();
    }
    hist_all->Write();
    hist_sub->Write();
    hist_ratio->Write();
    hist_eff->Write();
    output->Close();
    std::cout<<"analysis compleated\n";
    std::cout<<"real time: "<<timer.RealTime()<<"\n";
    std::cout<<"cpu time: "<<timer.CpuTime()<<"\n";
    return 0;
}

TH2D* fill_hist(TString filename, Int_t pdg_code, Int_t nevents, const TString branch="ntcand"){
    Int_t test = nevents;
    TFile *input = TFile::Open(filename.Data());
    TTree* theTree = (TTree*)input->Get(branch.Data());
    Float_t pt, rapidity, dist, cosp, d0prod, dca, mass, thetad;
    theTree->SetBranchAddress("m",&mass);    
    theTree->SetBranchAddress("pt",&pt);
    theTree->SetBranchAddress("rapidity",&rapidity);
    //theTree->SetBranchAddress("dist",&cdist);
    theTree->SetBranchAddress("cosp",&cosp);
    theTree->SetBranchAddress("d0prod",&d0prod);
    theTree->SetBranchAddress("dca",&dca);
    theTree->SetBranchAddress("thetad",&thetad);
    Double_t mass_ref = TDatabasePDG::Instance()->GetParticle(pdg_code)->Mass();
    Int_t nbins_mass = 40;
    Int_t nbins_pt = 10;
    TH2D* hist = new TH2D("hist","; m (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c}); counts", nbins_mass, mass_ref*0.96, mass_ref*1.04, nbins_pt, 0., 3);
    if(nevents==0)
        nevents = theTree->GetEntries();
    
    for(Long64_t itr=0; itr<nevents; itr++){
        if(itr%10000==0)
            std::cout<<"candidate: "<<itr<<"\n";

        theTree->GetEntry(itr);
        hist->Fill(mass, pt);
        
    }
    return hist;
}

Int_t normalize_evt_mix(TH2D* hist_mix, TH2D* hist_cmb){
    for(int bin_y=1; bin_y<=hist_mix->GetNbinsY(); bin_y++){
        Float_t evt_mix_counts = 0;
        Float_t evt_cmb_counts = 0;
        for(int bin_x=1; bin_x<=hist_mix->GetNbinsX(); bin_x++){//il range su cui si normalizza va pensato più accuratamente
            evt_mix_counts += hist_mix->GetBinContent(bin_x, bin_y);
            evt_cmb_counts += hist_cmb->GetBinContent(bin_x, bin_y);
        }
        Float_t normalization = evt_cmb_counts/evt_mix_counts;
        for(int bin_x=1; bin_x<=hist_mix->GetNbinsX(); bin_x++){
            hist_mix->SetBinContent(bin_x, bin_y, int(normalization*hist_mix->GetBinContent(bin_x, bin_y)));
            hist_mix->SetBinContent(bin_x, bin_y, TMath::Sqrt(hist_mix->GetBinContent(bin_x, bin_y)));
        }
    }
    return 0;
}

std::vector<TString> files_list(TString dirname="C:/root/folder/", TString sub_dir="BkgSub", TString ext=".root"){
    TSystemDirectory dir("main_dir", dirname.Data());
    TList *files = dir.GetListOfFiles();
    std::vector<TString> namelist;
    if(files){
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) { 
            fname = file->GetName();
            if (file->IsDirectory() && fname.BeginsWith(sub_dir.Data())){
                TSystemDirectory dirSub("sub_dir", Form("%s/%s",dirname.Data(), fname.Data()));
                TList *filesSub = dirSub.GetListOfFiles();
                bool broken=false;
                if(filesSub){
                    TSystemFile *fileSub;
                    TString fnameSub;
                    TIter nextSub(filesSub);
                    while ((fileSub=(TSystemFile*)nextSub())) { 
                        fnameSub = fileSub->GetName();
                        if (!fileSub->IsDirectory() && fnameSub.EndsWith(ext.Data())){
                            namelist.push_back(Form("%s/%s/%s", dirname.Data(), fname.Data(), fnameSub.Data()));
                            broken=true;
                            break;
                        }
                    }
                    if(!broken)
                        std::cout<<fname.Data()<<"\n";
                }
            }
        }
    }
    return namelist;
}

std::vector<TString> dir_list(TString dirname="C:/root/folder/", TString sub_dir="BkgSub"){
    TSystemDirectory dir("main_dir", dirname.Data());
    TList *files = dir.GetListOfFiles();
    std::vector<TString> namelist;
    if(files){
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) { 
            fname = file->GetName();
            if (file->IsDirectory() && fname.BeginsWith(sub_dir.Data())){
                TSystemDirectory dirSub("sub_dir", Form("%s/%s",dirname.Data(), fname.Data()));
                TList *filesSub = dirSub.GetListOfFiles();
                namelist.push_back(Form("%s/%s", dirname.Data(), fname.Data()));
            }
        }
    }
    return namelist;
}

int sum_results(TString data_path = "/gr3-str/NA60plus/phi_analysis/results",
                TString sub_dir = "BkgSub",
                TString suffix = "",
                Int_t fullsim = kFALSE){
    std::vector<TString> namelist = files_list(data_path.Data(),sub_dir.Data(),".root");
    TH2D* hist_bkg;
    TH2D* hist_mix;
    TH2D* hist_sig;
    TH2D* hist_all;
    TH2D* hist_sub;
    TH1D* hist_ev;
    int count=0;
    for (auto it = begin(namelist); it != end(namelist); ++it){
        TFile* inputFile = new TFile(it->Data(),"read");
        TH1D* hist_ev_it = (TH1D*) inputFile->Get("hNevents");
        TH2D* hist_mix_it = (TH2D*) inputFile->Get("hist_mix");
        TH2D* hist_all_it = (TH2D*) inputFile->Get("hist_all");
        TH2D* hist_sub_it = (TH2D*) inputFile->Get("hist_sub");
        TH2D* hist_sig_it;
        TH2D* hist_bkg_it;
        if(fullsim){
            hist_sig_it = (TH2D*) inputFile->Get("hist_sig");
            hist_bkg_it = (TH2D*) inputFile->Get("hist_bkg");
        }
        if(count==0){
            hist_ev = (TH1D*) hist_ev_it->Clone();
            hist_mix = (TH2D*) hist_mix_it->Clone();
            hist_all = (TH2D*) hist_all_it->Clone();
            hist_sub = (TH2D*) hist_sub_it->Clone();
            if(!fullsim){
                hist_bkg = (TH2D*) hist_bkg_it->Clone();
                hist_sig = (TH2D*) hist_sig_it->Clone();
            }
        }
        else{
            hist_ev->Add(hist_ev_it);
            hist_mix->Add(hist_mix_it);
            hist_all->Add(hist_all_it);
            hist_sub->Add(hist_sub_it);
            if(!fullsim){
                hist_sig->Add(hist_sig_it);
                hist_bkg->Add(hist_bkg_it);
            }
        }
        count++;
    }
    std::cout<<"nfile: "<<count<<"\n";
    TFile *outputFile = new TFile(Form("sum_event_mixing_phi%s.root", suffix.Data()), "recreate");

    TDirectory* dir_mix = outputFile->mkdir("event mixing");
    TDirectory* dir_sub = outputFile->mkdir("subtraction");
    TDirectory* dir_all = outputFile->mkdir("all");
    save_projections(hist_mix, dir_mix);
    save_projections(hist_sub, dir_sub);
    save_projections(hist_all, dir_all);
    if(!fullsim){
        TDirectory* dir_bkg = outputFile->mkdir("background");
        TDirectory* dir_sig = outputFile->mkdir("signal");
        save_projections(hist_sig, dir_sig);
        save_projections(hist_bkg, dir_bkg);
    }
    outputFile->cd();
    hist_mix->Write();
    hist_all->Write();
    hist_sub->Write();
    hist_ev->Write();
    if(!fullsim){
        hist_sig->Write();
        hist_bkg->Write();
    }
    outputFile->Close();
    return 0;
}


int merge_sparse(TString data_path = "/gr3-str/NA60plus/phi_analysis/combinatorial",
                 TString mixing_path = "/gr3-str/NA60plus/phi_analysis/event_mixing",
                 TString sub_dir_suffix = "PHI_L10_E40",
                 TString suffix = "PHI_L10_E40"){
    std::vector<TString> datalist = dir_list(data_path.Data(),  Form("CmbBkg_%s",sub_dir_suffix.Data()));
    std::vector<TString> mixlist  = dir_list(mixing_path.Data(),Form("EvtMix_%s",sub_dir_suffix.Data()));
    THnSparse* hsp_mix;
    THnSparse* hsp_data;
    TH1D* hist_ev;
    int count_data=0;
    for (auto it = begin(datalist); it != end(datalist); ++it){
        TFile* inputFile = new TFile(Form("%s/fhspBkg_%s",it->Data(),sub_dir_suffix.Data()),"read");
        TFile* eventFile = new TFile(Form("%s/Bkg-histos_%s",it->Data(),sub_dir_suffix.Data()),"read");
        THnSparse* hist_ev_it = (THnSparse*) inputFile->Get("hsp");
        TH1D* hsp_it = (TH1D*) eventFile->Get("hNevents");
        if(count_data==0){
            hist_ev = (TH1D*) hist_ev_it->Clone();
            hsp_mix = (TH2D*) hsp_it->Clone();
        }
        else{
            hist_ev->Add(hist_ev_it);
            hsp_mix->Add(hsp_it);
        }
        count_data++;
        eventFile->Close();
        inputFile->Close();
    }
    int count_mix=0;
    for (auto it = begin(mixlist); it != end(mixlist); ++it){
        TFile* inputFile = new TFile(Form("%s/fhspBkg_%s",it->Data(),suffix.Data()),"read");
        if(count_mix==0)
            hsp_mix = (TH2D*) hsp_it->Clone();
        else
            hsp_mix->Add(hsp_it);
        count_mix++;
        inputFile->Close();
    }
    std::cout<<"nfile data: "<<count_data<<"\n";
    std::cout<<"nfile event-mixing: "<<count_mixing<<"\n";
    TFile *outputFile = new TFile(Form("Results%s.root", suffix.Data()), "recreate");
    hsp_data->Write();
    hsp_mix->Write();
    hist_ev->Write();
    outputFile->Close();
    return 0;
}


void save_projections(TH2D* hist, TDirectory* dir){
    dir->cd();
    for(int bin=1; bin<=hist->GetNbinsY(); bin++){
        Float_t ptbin_lw = hist->GetYaxis()->GetBinLowEdge(bin);
        Float_t ptbin_up = hist->GetYaxis()->GetBinUpEdge(bin);
        TH1D* proj = hist->ProjectionX(Form("%s_pt_%.1f_%.1f",hist->GetName(), ptbin_lw, ptbin_up), bin, bin);
        proj->Write();
    }
    TH1D* proj = hist->ProjectionX(Form("%s_all_pt",hist->GetName()), 1, hist->GetNbinsY());
    proj->Write();

}
/*
void save_and_fit_projections(TH2D* hist, TDirectory* dir){
    dir->cd();
    Int_t nbins = hist->GetNbinsY();
    Float_t min_bin = hist->GetYaxis()->GetBinLowEdge(1);
    Float_t max_bin = hist->GetYaxis()->GetBinUpEdge(nbins);
    
    for(int bin=1; bin<=nbins; bin++){
        Float_t ptbin_lw = hist->GetYaxis()->GetBinLowEdge(bin);
        Float_t ptbin_up = hist->GetYaxis()->GetBinUpEdge(bin);
        TH1D* proj = hist->ProjectionX(, bin, bin);
        //TF1* fit_func = new TF1("fit_func","gausn(0)+pol2(3)");
        //proj->Fit("fit_func","M+");
        proj->Write();
    }
    TH1D* proj = hist->ProjectionX(Form("%s_all_pt",hist->GetName()), 1, nbins);
    proj->Write();

}
*/


int fit_results(int pdg_code, TString dir, TString suffix =""){
    
    TH1D* pt_spectra_raw = new TH1D("pt_spectra_raw","; #it{p}_{T} (GeV/#it{c}); dN/d#it{p}_{T} (GeV/#it{c})^{-1}",nbins, min_bin, max_bin);
    TH1D* pt_spectra_cor = new TH1D("pt_spectra_cor","; #it{p}_{T} (GeV/#it{c}); dN/d#it{p}_{T} (GeV/#it{c})^{-1}",nbins, min_bin, max_bin);
    TH1D* sigma_value = new TH1D("sigma_value","; #it{p}_{T} (GeV/#it{c}); #sigma (GeV/#it{c}^{2})",nbins, min_bin, max_bin);
    TH1D* mass_value = new TH1D("mass_value","; #it{p}_{T} (GeV/#it{c}); mass (GeV/#it{c}^{2})",nbins, min_bin, max_bin);
    Double_t mass_ref = TDatabasePDG::Instance()->GetParticle(pdg_code)->Mass();
    TFile *inputFile = new TFile(Form("%/sum_event_mixing_phi%s.root", dir.Data(), suffix.Data()), "recreate");
    TFile *outputFile = new TFile("analysed_data.root","recreate");
    inputFile->cd();
    TH2D* hist_sub =  (TH2D*) inputFile->Get("hist_sub"););
    TH1D* hist_eff =  (TH1D*) inputFile->Get("hist_eff");
    int nbins  = hist_sub->GetNbinsY();
    TF1* fit_func = new TF1("fit_func","gausn(0)+pol2(3)");
    TF1* fit_pt = new TF1("fit_pt","[0]*x*TMath::exp(-TMath::Sqrt(x**2+[1]**2)/[2])");
    fit_pt->FitParameter(1, mass);
    TH1D* hist_proj;
    TDirectory* dir_sub = inputFile->mkdir("");
    for(int pt_bin = 1; pt_bin<nbins; pt_bin++){
        inputFile->cd();
        Float_t ptbin_lw = hist->GetYaxis()->GetBinLowEdge(bin);
        Float_t ptbin_up = hist->GetYaxis()->GetBinUpEdge(bin);
        hist_proj = (TH1D*) inputFile->Get(Form("subtraction/hist_sub_pt_%.1f_%.1f",hist->GetName(), ptbin_lw, ptbin_up));
        hist_proj->Fit("fit_func","M+");
        inputFile->cd();
        hist_proj->Write();
        mass_value->SetBinContent(pt_bin, fit_func->GetParameter(1));
        mass_value->SetBinError(pt_bin, fit_func->GetParError(1));
        sigma_value->SetBinContent(pt_bin, fit_func->GetParameter(2));
        sigma_value->SetBinError(pt_bin, fit_func->GetParError(2));
        pt_spectra_raw->SetBinContent(pt_bin, fit_func->GetParameter(0));
        pt_spectra_raw->SetBinError(pt_bin, fit_func->GetParError(0));
        pt_spectra_cor->SetBinContent(pt_bin, fit_func->GetParameter(0)/hist_eff->GetBinContent(pt_bin));
        pt_spectra_cor->SetBinError(pt_bin, fit_func->GetParError(0)/hist_eff->GetBinContent(pt_bin));
        
    }
    pt_spectra_corr->Fit("fit_pt","M+");
    pt_spectra_raw->Write();
    pt_spectra_cor->Write();
    sigma_value->Write();
    mass_value->Write();
    inputFile->Close();
    outputFile->Close();
    return 0;
}