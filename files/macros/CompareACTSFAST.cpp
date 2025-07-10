

#include "../KMCProbeFwd.h"

void setCanvasForTH2D()
{
    gStyle->SetTitleSize(0.07, "XYZ"); // Sets the size for X, Y, and Z axis titles
    gStyle->SetTitleFontSize(0.06);    // Sets the size for the histogram title
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadLeftMargin(0.1);   // Set left margin
    gStyle->SetPadRightMargin(0.1);  // Set right margin
    gStyle->SetPadTopMargin(0.2);    // Set top margin
    gStyle->SetPadBottomMargin(0.1); // Set bottom margin
    gStyle->SetHistLineWidth(2);
}

void setCanvasForTH1D()
{
    gStyle->SetTitleSize(0.07, "XYZ"); // Sets the size for X, Y, and Z axis titles
    gStyle->SetTitleFontSize(0.07);    // Sets the size for the histogram title
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadLeftMargin(0.15);   // Set left margin
    gStyle->SetPadRightMargin(0.01);  // Set right margin
    gStyle->SetPadTopMargin(0.03);    // Set top margin
    gStyle->SetPadBottomMargin(0.15); // Set bottom margin
    gStyle->SetHistLineWidth(2);
}
void SaveImpParFast(TString fileName = "/home/giacomo/na60fastsim_galocco/files/macros/treeBkgEvents.root", TString output = "res.root")
{
    TFile *filetreeFast = new TFile(fileName.Data());
    TTree *treeFast = (TTree *)filetreeFast->Get("tree");
    TClonesArray *arrFast = 0;
    treeFast->SetBranchAddress("tracks", &arrFast);
    TH2F *hIPVsPFast = new TH2F("hIPVsPFast", "", 20, 0, 10, 200, -0.03, 0.03);
    TH1F *hIPFast_GetX = new TH1F("hIPFast_GetX", "", 200, -0.03, 0.03);
    TH1F *hIPFast_GetY = new TH1F("hIPFast_GetY", "", 200, -0.03, 0.03);
    TH1F *hIPFast_GetZ = new TH1F("hIPFast_GetZ", "", 200, -0.03, 0.03);
    TH2F *hIPVsEtaFast = new TH2F("hIPVsEtaFast", "", 30, 2, 5, 200, -0.03, 0.03);
    for (Int_t iev = 0; iev < treeFast->GetEntries(); iev++)
    {
        Double_t vprim[3] = {0, 0, 0};
        treeFast->GetEvent(iev);
        Int_t arrentr = arrFast->GetEntriesFast();
        for (Int_t itr = 0; itr < arrentr; itr++)
        {
            KMCProbeFwd *tr1 = (KMCProbeFwd *)arrFast->At(itr);
            Double_t pxyz[3];
            tr1->GetPXYZ(pxyz);
            Double_t impParXY = tr1->GetY();
            Double_t pt = TMath::Sqrt(pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1]);
            Double_t p = TMath::Sqrt(pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1] + pxyz[2] * pxyz[2]);
            hIPVsPFast->Fill(p, impParXY);

            hIPFast_GetX->Fill(tr1->GetX());
            hIPFast_GetY->Fill(tr1->GetY());
            hIPFast_GetZ->Fill(tr1->GetZ());
            double e = TMath::Sqrt(p * p + tr1->GetMass() * tr1->GetMass());
            Double_t eta = 0.5 * TMath::Log((e + pxyz[2]) / (e - pxyz[2]));
            hIPVsEtaFast->Fill(eta, impParXY);
        }
    }
    TFile *file = new TFile(output.Data(), "RECREATE");

    TCanvas *c1Fast = new TCanvas("c1Fast", "", 1600, 800);
    c1Fast->Divide(4, 4);
    c1Fast->cd(1);
    hIPVsPFast->Draw("colz");
    TGraphErrors *gGausWFast = new TGraphErrors(0);
    for (int ib = 1; ib <= hIPVsPFast->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hIPVsPFast->ProjectionY(Form("ImpFastVsP_%d", ib), ib, ib);
        c1Fast->cd(1 + ib);
        htmp->Draw();
        htmp->Write();
        htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
        TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
        if (fg)
        {
            gGausWFast->SetPoint(ib - 1, hIPVsPFast->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) * 10000.);
            gGausWFast->SetPointError(ib - 1, 0.5 * hIPVsPFast->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) * 10000.);
        }
    }
    std::cout << "FIT FAST eta" << std::endl;

    TCanvas *c1FastEta = new TCanvas("c1FastEta", "", 1500, 800);
    c1FastEta->Divide(4, 4);
    c1FastEta->cd(1);
    hIPVsEtaFast->Draw("colz");
    TGraphErrors *gGausWFastEta = new TGraphErrors(0);
    for (int ib = 1; ib <= hIPVsEtaFast->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hIPVsEtaFast->ProjectionY(Form("ImpFastVsEta_%d", ib), ib, ib);
        c1FastEta->cd(1 + ib);
        htmp->Draw();
        htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
        TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
        if (fg)
        {
            gGausWFastEta->SetPoint(ib - 1, hIPVsEtaFast->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) * 10000.);
            gGausWFastEta->SetPointError(ib - 1, 0.5 * hIPVsEtaFast->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) * 10000.);
        }
    }

    hIPVsPFast->Write();
    hIPVsEtaFast->Write();

    gGausWFast->SetName("ImpVsP");
    gGausWFastEta->SetName("ImpVsEta");
    gGausWFast->Write();
    gGausWFastEta->Write();
    c1Fast->Write();
    c1FastEta->Write();

    hIPFast_GetX->Write();
    hIPFast_GetY->Write();
    hIPFast_GetZ->Write();
    file->Close();
}

void SaveImpParACTS(TString fileName = "output_events_40GeV", TString output = "res.root")
{

    gStyle->SetPadRightMargin(0.02);  // Set right margin
    gStyle->SetPadRightMargin(0.02);  // Set right margin
    gStyle->SetPadLeftMargin(0.1);  // Set right margin
    gStyle->SetPadLeftMargin(0.1);  // Set right margin
    TFile *filetreeACTS = new TFile((fileName + "/tracksummary_ambi.root").Data());
    TTree *treeACTS = (TTree *)filetreeACTS->Get("tracksummary");

    TFile *filetreeACTSll = new TFile((fileName + "/tracksummary_ambill.root").Data());
    TTree *treeACTSll = (TTree *)filetreeACTSll->Get("tracksummary");

    TFile *filetreeACTSllll = new TFile((fileName + "/tracksummary_ambillll.root").Data());
    TTree *treeACTSllll = (TTree *)filetreeACTSllll->Get("tracksummary");

    TFile *filetreeACTSllllll = new TFile((fileName + "/tracksummary_ambillllll.root").Data());
    TTree *treeACTSllllll = (TTree *)filetreeACTSllllll->Get("tracksummary");
    std::vector<float> *t_pT = nullptr;
    std::vector<float> *t_pTll = nullptr;
    std::vector<float> *t_pTllll = nullptr;
    std::vector<float> *t_pTllllll = nullptr;

    std::vector<float> *d0 = nullptr;
    std::vector<float> *d0ll = nullptr;
    std::vector<float> *d0llll = nullptr;
    std::vector<float> *d0llllll = nullptr;

    std::vector<float> *eta = nullptr;
    std::vector<float> *etall = nullptr;
    std::vector<float> *etallll = nullptr;
    std::vector<float> *etallllll = nullptr;

    TH2F *hIPVsPACTS = new TH2F("hIPVsPACTS", ";p [GeV/c];d_{0} [cm]", 20, 0, 10, 200, -0.03, 0.03);
    TH2F *hIPVsEtaACTS = new TH2F("hIPVsEtaACTS", ";#eta [GeV/c];d_{0} [cm]", 30, 2, 5, 200, -0.03, 0.03);

    treeACTS->SetBranchAddress("t_p", &t_pT);
    treeACTSll->SetBranchAddress("t_p", &t_pTll);
    treeACTSllll->SetBranchAddress("t_p", &t_pTllll);
    treeACTSllllll->SetBranchAddress("t_p", &t_pTllllll);

    treeACTS->SetBranchAddress("res_eLOC0_fit", &d0);
    treeACTSll->SetBranchAddress("res_eLOC0_fit", &d0ll);
    treeACTSllll->SetBranchAddress("res_eLOC0_fit", &d0llll);
    treeACTSllllll->SetBranchAddress("res_eLOC0_fit", &d0llllll);

    treeACTS->SetBranchAddress("t_eta", &eta);
    treeACTSll->SetBranchAddress("t_eta", &etall);
    treeACTSllll->SetBranchAddress("t_eta", &etallll);
    treeACTSllllll->SetBranchAddress("t_eta", &etallllll);

    for (Int_t iev = 0; iev < treeACTS->GetEntries(); iev++)
    {
        treeACTS->GetEntry(iev);
        treeACTSll->GetEntry(iev);
        treeACTSllll->GetEntry(iev);
        treeACTSllllll->GetEntry(iev);

        //  loop over tracks in this event
        Int_t nTracks = t_pT->size(); // t_pT, d0, etc. all have same size
        // std::cout<<t_pT->size()<<"  "<<d0->size()<<"  "<<eta->size()<<std::endl;
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            hIPVsPACTS->Fill(t_pT->at(itr), d0->at(itr) / 10.);
            hIPVsEtaACTS->Fill(eta->at(itr), d0->at(itr) / 10.);
        }
        // std::cout<<"ACTS pTll"<<std::endl;
        //  loop over tracks in this event
        nTracks = t_pTll->size(); // t_pT, d0, etc. all have same size
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            hIPVsPACTS->Fill(t_pTll->at(itr), d0ll->at(itr) / 10.);
            hIPVsEtaACTS->Fill(etall->at(itr), d0ll->at(itr) / 10.);
        } // loop over tracks in this event
        nTracks = t_pTllll->size(); // t_pT, d0, etc. all have same size
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            hIPVsPACTS->Fill(t_pTllll->at(itr), d0llll->at(itr) / 10.);
            hIPVsEtaACTS->Fill(etallll->at(itr), d0llll->at(itr) / 10.);
        } // loop over tracks in this event

        nTracks = t_pTllllll->size(); // t_pT, d0, etc. all have same size
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            hIPVsPACTS->Fill(t_pTllllll->at(itr), d0llllll->at(itr) / 10.);
            hIPVsEtaACTS->Fill(etallllll->at(itr), d0llllll->at(itr) / 10.);
        }
    }
    TFile *file = new TFile(output.Data(), "RECREATE");

    TCanvas *c1ACTS = new TCanvas("c1ACTS", "", 1500, 800);
    c1ACTS->Divide(4, 4);
    c1ACTS->cd(1);
    hIPVsPACTS->Draw("colz");
    TGraphErrors *gGausWACTS = new TGraphErrors(0);
    for (int ib = 1; ib <= hIPVsPACTS->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hIPVsPACTS->ProjectionY(Form("ImpACTSVsP_%d", ib), ib, ib);

        c1ACTS->cd(1 + ib);
        htmp->Draw();
        htmp->Write();
        htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
        TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
        if (fg)
        {
            gGausWACTS->SetPoint(ib - 1, hIPVsPACTS->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) * 10000.);
            gGausWACTS->SetPointError(ib - 1, 0.5 * hIPVsPACTS->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) * 10000.);
        }
    }
    std::cout << "FIT ACTS eta" << std::endl;
    TCanvas *c1ACTSEta = new TCanvas("c1ACTSEta", "", 1500, 800);
    c1ACTSEta->Divide(4, 4);
    c1ACTSEta->cd(1);
    hIPVsEtaACTS->Draw("colz");
    TGraphErrors *gGausWACTSEta = new TGraphErrors(0);
    for (int ib = 1; ib <= hIPVsEtaACTS->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hIPVsEtaACTS->ProjectionY(Form("ImpFastVsEta_%d", ib), ib, ib);
        c1ACTSEta->cd(1 + ib);
        htmp->Draw();
        htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
        TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
        if (fg)
        {
            gGausWACTSEta->SetPoint(ib - 1, hIPVsEtaACTS->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) * 10000.);
            gGausWACTSEta->SetPointError(ib - 1, 0.5 * hIPVsEtaACTS->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) * 10000.);
        }
    }

    hIPVsPACTS->Write();
    hIPVsEtaACTS->Write();

    gGausWACTS->SetName("ImpVsP");
    gGausWACTSEta->SetName("ImpVsEta");
    gGausWACTS->Write();
    gGausWACTSEta->Write();
    c1ACTS->Write();
    c1ACTSEta->Write();
    file->Close();
}

void SavePResFast(TString fileName = "/home/giacomo/na60fastsim_galocco/files/macros/treeBkgEvents.root", TString output = "res.root")

{

    TFile *fileh = new TFile(fileName.Data());
    TH2F *hResidPVsP = (TH2F *)fileh->Get("hResidPVsP");
    hResidPVsP->RebinX(2);
    TH2F *hResidPVsEta = (TH2F *)fileh->Get("hResidPVsEta");
    TH2F* hRecPVsEta = (TH2F *)fileh->Get("hRecPVsEta");


    TFile *file = new TFile(output.Data(), "RECREATE");

    TH1D *hResidPVsEtaProj = hResidPVsEta->ProjectionX("hResidPVsEtaProj");
    hResidPVsEtaProj->Write();
    TH1D *hResidPVsPProj = hResidPVsP->ProjectionX("hResidPVsPProj");
    hResidPVsPProj->Write();
    hRecPVsEta->Write();

    TCanvas *c1 = new TCanvas("c1", "", 1500, 800);
    c1->Divide(8, 8);
    c1->cd(1);
    hResidPVsP->Draw("colz");
    TGraphErrors *gRelGaus_resPVsP = new TGraphErrors(0);
    TGraphErrors *gRelGaus_meanP = new TGraphErrors(0);

    TDirectory *dirVsP = file->mkdir("ResPACTSVsP");
    dirVsP->cd();
    for (int ib = 1; ib <= hResidPVsP->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hResidPVsP->ProjectionY(Form("ResPFastVsP_%d", ib), ib, ib);
        c1->cd(1 + ib);
        htmp->Draw();
        htmp->Write();
        if (htmp->GetEntries() > 100)
        {
            int n = gRelGaus_resPVsP->GetN();
            htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
            TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
            if (fg)
            {
                gRelGaus_resPVsP->SetPoint(n, hResidPVsP->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) / hResidPVsP->GetXaxis()->GetBinCenter(ib));
                gRelGaus_resPVsP->SetPointError(n, 0.5 * hResidPVsP->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) / hResidPVsP->GetXaxis()->GetBinCenter(ib));
                gRelGaus_meanP->SetPoint(n, hResidPVsP->GetXaxis()->GetBinCenter(ib), fg->GetParameter(1));
                gRelGaus_meanP->SetPointError(n, 0.5 * hResidPVsP->GetXaxis()->GetBinWidth(ib), fg->GetParError(1));
            }
        }
    }
    c1->Write();
    TCanvas *c2 = new TCanvas("c2", "", 1500, 800);
    c2->Divide(8, 8);
    c2->cd(1);
    TGraphErrors *gRelGaus_resPVsEta = new TGraphErrors(0);
    TGraphErrors *gRelGaus_meanPVsEta = new TGraphErrors(0);
    TDirectory *dirVsEta = file->mkdir("ResPACTSVsEta");
    dirVsEta->cd();
    for (int ib = 1; ib <= hResidPVsEta->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hResidPVsEta->ProjectionY(Form("ResPFastVsEta_%d", ib), ib, ib);
        c2->cd(1 + ib);
        htmp->Draw();
        htmp->Write();
        if (htmp->GetEntries() > 100)
        {
            int n = gRelGaus_resPVsEta->GetN();
            htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
            TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
            if (fg)
            {
                gRelGaus_resPVsEta->SetPoint(n, hResidPVsEta->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) / hResidPVsEta->GetXaxis()->GetBinCenter(ib));
                gRelGaus_resPVsEta->SetPointError(n, 0.5 * hResidPVsEta->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) / hResidPVsEta->GetXaxis()->GetBinCenter(ib));
                gRelGaus_meanPVsEta->SetPoint(n, hResidPVsEta->GetXaxis()->GetBinCenter(ib), fg->GetParameter(1));
                gRelGaus_meanPVsEta->SetPointError(n, 0.5 * hResidPVsEta->GetXaxis()->GetBinWidth(ib), fg->GetParError(1));
            }
        }
    }
    file->cd();
    c2->Write();
    hResidPVsP->Write();
    hResidPVsEta->Write();

    gRelGaus_resPVsP->SetName("ResPVsP");
    gRelGaus_resPVsEta->SetName("ResPVsEta");
    gRelGaus_resPVsP->Write();
    gRelGaus_resPVsEta->Write();
    gRelGaus_meanP->SetName("MeanPVsP");
    gRelGaus_meanPVsEta->SetName("MeanPVsEta");
    gRelGaus_meanP->Write();
    gRelGaus_meanPVsEta->Write();
    file->Close();
}

void SavePResACTS(TString fileName = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_VTOnly", TString output = "res.root", bool fivehits = false)
{

    TFile *filetreeACTS = new TFile((fileName + "/tracksummary_ambi.root").Data());
    TTree *treeACTS = (TTree *)filetreeACTS->Get("tracksummary");

    TFile *filetreeACTSll = new TFile((fileName + "/tracksummary_ambill.root").Data());
    TTree *treeACTSll = (TTree *)filetreeACTSll->Get("tracksummary");

    TFile *filetreeACTSllll = new TFile((fileName + "/tracksummary_ambillll.root").Data());
    TTree *treeACTSllll = (TTree *)filetreeACTSllll->Get("tracksummary");

    TFile *filetreeACTSllllll = new TFile((fileName + "/tracksummary_ambillllll.root").Data());
    TTree *treeACTSllllll = (TTree *)filetreeACTSllllll->Get("tracksummary");
    std::vector<float> *t_p = nullptr;
    std::vector<float> *t_pll = nullptr;
    std::vector<float> *t_pllll = nullptr;
    std::vector<float> *t_pllllll = nullptr;
    std::vector<float> *t_eta = nullptr;
    std::vector<float> *t_etall = nullptr;
    std::vector<float> *t_etallll = nullptr;
    std::vector<float> *t_etallllll = nullptr;

    std::vector<float> *t_meas = nullptr;
    std::vector<float> *t_measll = nullptr;
    std::vector<float> *t_measllll = nullptr;
    std::vector<float> *t_measllllll = nullptr;

    std::vector<float> *t_maj = nullptr;
    std::vector<float> *t_majll = nullptr;
    std::vector<float> *t_majllll = nullptr;
    std::vector<float> *t_majllllll = nullptr;

    std::vector<float> *resp = nullptr;
    std::vector<float> *respll = nullptr;
    std::vector<float> *respllll = nullptr;
    std::vector<float> *respllllll = nullptr;

    treeACTS->SetBranchAddress("t_p", &t_p);
    treeACTSll->SetBranchAddress("t_p", &t_pll);
    treeACTSllll->SetBranchAddress("t_p", &t_pllll);
    treeACTSllllll->SetBranchAddress("t_p", &t_pllllll);

    treeACTS->SetBranchAddress("t_eta", &t_eta);
    treeACTSll->SetBranchAddress("t_eta", &t_etall);
    treeACTSllll->SetBranchAddress("t_eta", &t_etallll);
    treeACTSllllll->SetBranchAddress("t_eta", &t_etallllll);

    treeACTS->SetBranchAddress("nMeasurements", &t_meas);
    treeACTSll->SetBranchAddress("nMeasurements", &t_measll);
    treeACTSllll->SetBranchAddress("nMeasurements", &t_measllll);
    treeACTSllllll->SetBranchAddress("nMeasurements", &t_measllllll);

    treeACTS->SetBranchAddress("nMajorityHits", &t_maj);
    treeACTSll->SetBranchAddress("nMajorityHits", &t_majll);
    treeACTSllll->SetBranchAddress("nMajorityHits", &t_majllll);
    treeACTSllllll->SetBranchAddress("nMajorityHits", &t_majllllll);

    treeACTS->SetBranchAddress("eQOP_fit", &resp);
    treeACTSll->SetBranchAddress("eQOP_fit", &respll);
    treeACTSllll->SetBranchAddress("eQOP_fit", &respllll);
    treeACTSllllll->SetBranchAddress("eQOP_fit", &respllllll);

    TFile *fileh = new TFile("bkgdistributions_40GeV_Carbon_5hits.root", "READ");
    TH2F *hResidPVsP = (TH2F *)fileh->Get("hResidPVsP");
    Int_t nbinsX = hResidPVsP->GetNbinsX() / 2.;
    Int_t nbinsY = hResidPVsP->GetNbinsY();
    double xmin = hResidPVsP->GetXaxis()->GetXmin();
    double xmax = hResidPVsP->GetXaxis()->GetXmax();
    double ymin = hResidPVsP->GetYaxis()->GetXmin();
    double ymax = hResidPVsP->GetYaxis()->GetXmax();

    TH2F *hResidPVsPACTS = new TH2F("hResidPVsPACTS", ";p [GeV/c];#sigma_p",
                                    nbinsX, xmin, xmax,
                                    nbinsY, ymin, ymax);

    TH2F *hResidPVsEta = (TH2F *)fileh->Get("hResidPVsEta");
    TH2F *hResidPVsEtaACTS = new TH2F("hResidPVsEtaACTS", ";eta;#sigma_p", hResidPVsEta->GetNbinsX(), hResidPVsEta->GetXaxis()->GetXmin(), hResidPVsEta->GetXaxis()->GetXmax(), hResidPVsEta->GetNbinsY(), hResidPVsEta->GetYaxis()->GetXmin(), hResidPVsEta->GetYaxis()->GetXmax());
    TH2F* hRecPVsEta = new TH2F("hRecPVsEta","",40,1.,5.,200,0.,20.);

    TFile *file = new TFile(output.Data(), "RECREATE");

    TCanvas *c1ACTS = new TCanvas("c1ACTS", "", 1500, 800);
    c1ACTS->Divide(8, 8);
    c1ACTS->cd(1);
    hResidPVsPACTS->Draw("colz");
    
    for (Int_t iev = 0; iev < treeACTS->GetEntries(); iev++)
    {
        treeACTS->GetEntry(iev);
        treeACTSll->GetEntry(iev);
        treeACTSllll->GetEntry(iev);
        treeACTSllllll->GetEntry(iev);

        Int_t nTracks = t_p->size(); // t_p, d0, etc. all have same size
        // std::cout<<t_p->size()<<"  "<<d0->size()<<"  "<<eta->size()<<std::endl;
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            if (fivehits && t_meas->at(itr) < 5)//&& t_maj->at(itr) != 5)
                continue; // skip tracks with less than 5 hits
            Float_t prec = TMath::Abs(1. / (resp->at(itr)));
            
            hResidPVsPACTS->Fill(t_p->at(itr), prec - t_p->at(itr));
            hResidPVsEtaACTS->Fill(t_eta->at(itr), prec - t_p->at(itr));
            hRecPVsEta->Fill(t_eta->at(itr), t_p->at(itr));
        }
        nTracks = t_pll->size(); // t_p, resp, etc. all have same size
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            if (fivehits && t_measll->at(itr) < 5)//&& t_majll->at(itr) != 5)
                continue; // skip tracks with less than 5 hits
            Float_t prec = TMath::Abs(1. / (respll->at(itr)));
            hResidPVsPACTS->Fill(t_pll->at(itr), prec - t_pll->at(itr));
            hResidPVsEtaACTS->Fill(t_etall->at(itr), prec - t_pll->at(itr));
            hRecPVsEta->Fill(t_etall->at(itr), t_pll->at(itr));
        } // loop over tracks in this event

        nTracks = t_pllll->size(); // t_p, resp, etc. all have same size
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            if (fivehits && t_measllll->at(itr) < 5)//&& t_majllll->at(itr) != 5)
                continue; // skip tracks with less than 5 hits
            Float_t prec = TMath::Abs(1. / (respllll->at(itr)));
            hResidPVsPACTS->Fill(t_pllll->at(itr), prec - t_pllll->at(itr));            
            hResidPVsEtaACTS->Fill(t_etallll->at(itr), prec - t_pllll->at(itr));
            hRecPVsEta->Fill(t_etallll->at(itr), t_pllll->at(itr));

        } // loop over tracks in this event
        nTracks = t_pllllll->size(); // t_p, resp, etc. all have same size
        for (Int_t itr = 0; itr < nTracks; itr++)
        {
            if (fivehits && t_measllllll->at(itr) < 5)//&& t_majllllll->at(itr) != 5)
                continue; // skip tracks with less than 5 hits
            Float_t prec = TMath::Abs(1. / (respllllll->at(itr)));
            hResidPVsPACTS->Fill(t_pllllll->at(itr), prec - t_pllllll->at(itr));
            hResidPVsEtaACTS->Fill(t_etallllll->at(itr), prec - t_pllllll->at(itr));
            hRecPVsEta->Fill(t_etallllll->at(itr), t_pllllll->at(itr));
        }
    }

    TGraphErrors *gRelGausACTS = new TGraphErrors(0);
    TGraphErrors *gMeanGausACTS = new TGraphErrors(0);
    TDirectory *dirVsP = file->mkdir("ResPACTSVsP");
    dirVsP->cd();

    for (int ib = 1; ib <= hResidPVsPACTS->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hResidPVsPACTS->ProjectionY(Form("ResPACTSVsP_%d", ib), ib, ib);
        htmp->Write();
        c1ACTS->cd(1 + ib);
        htmp->Draw();
        if (htmp->GetEntries() > 100)
        {
            int n = gRelGausACTS->GetN();
            htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
            TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
            if (fg)
            {
                gRelGausACTS->SetPoint(n, hResidPVsPACTS->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) / hResidPVsPACTS->GetXaxis()->GetBinCenter(ib));
                gRelGausACTS->SetPointError(n, 0.5 * hResidPVsPACTS->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) / hResidPVsPACTS->GetXaxis()->GetBinCenter(ib));
                gMeanGausACTS->SetPoint(n, hResidPVsPACTS->GetXaxis()->GetBinCenter(ib), fg->GetParameter(1));
                gMeanGausACTS->SetPointError(n, 0.5 * hResidPVsPACTS->GetXaxis()->GetBinWidth(ib), fg->GetParError(1));
            }
        }
    }

    TGraphErrors *gRelGausVsEtaACTS = new TGraphErrors(0);
    TGraphErrors *gMeanGausVsEtaACTS = new TGraphErrors(0);
    TDirectory *dirVsEta = file->mkdir("ResPACTSVsEta");
    dirVsEta->cd();
    for (int ib = 1; ib <= hResidPVsEtaACTS->GetXaxis()->GetNbins(); ib++)
    {
        TH1D *htmp = hResidPVsEtaACTS->ProjectionY(Form("ResPACTSVsEta_%d", ib), ib, ib);
        htmp->Write();
        c1ACTS->cd(1 + ib);
        htmp->Draw();
        if (htmp->GetEntries() > 100)
        {
            int n = gRelGausVsEtaACTS->GetN();
            htmp->Fit("gaus", "Q", "", -2. * htmp->GetRMS(), 2. * htmp->GetRMS());
            TF1 *fg = (TF1 *)htmp->GetListOfFunctions()->FindObject("gaus");
            if (fg)
            {
                gRelGausVsEtaACTS->SetPoint(n, hResidPVsEtaACTS->GetXaxis()->GetBinCenter(ib), fg->GetParameter(2) / hResidPVsEtaACTS->GetXaxis()->GetBinCenter(ib));
                gRelGausVsEtaACTS->SetPointError(n, 0.5 * hResidPVsEtaACTS->GetXaxis()->GetBinWidth(ib), fg->GetParError(2) / hResidPVsEtaACTS->GetXaxis()->GetBinCenter(ib));
                gMeanGausVsEtaACTS->SetPoint(n, hResidPVsEtaACTS->GetXaxis()->GetBinCenter(ib), fg->GetParameter(1));
                gMeanGausVsEtaACTS->SetPointError(n, 0.5 * hResidPVsEtaACTS->GetXaxis()->GetBinWidth(ib), fg->GetParError(1));
            }
        }
    }
    file->cd();
    gRelGausACTS->SetName("ResPVsP");
    gRelGausVsEtaACTS->SetName("ResPVsEta");
    gRelGausACTS->Write();
    gRelGausVsEtaACTS->Write();
    hResidPVsPACTS->Write();
    
    gMeanGausACTS->SetName("MeanPVsP");
    gMeanGausVsEtaACTS->SetName("MeanPVsEta");
    gMeanGausACTS->Write();
    gMeanGausVsEtaACTS->Write();

    TH1D *hResidPVsEtaProj = hResidPVsEtaACTS->ProjectionX("hResidPVsEtaProj");
    hResidPVsEtaProj->Write();
    TH1D *hResidPVsPProj = hResidPVsPACTS->ProjectionX("hResidPVsPProj");
    hResidPVsPProj->Write();

    hResidPVsEtaACTS->Write();
    c1ACTS->Write();
    hRecPVsEta->Write();
    file->Close();
    std::cout << "ACTS PRes saved in " << output.Data() << std::endl;
}

void CompareACTSFAST(bool analyze = true)
{
    gROOT->SetBatch(true);

    if (analyze)
    {
        SaveImpParFast("/home/giacomo/na60fastsim_galocco/files/macros/treeBkgEvents_40GeV_Void_5hits.root", "output_fast/ImpPar_FAST_Void_5hits.root");
        SaveImpParFast("/home/giacomo/na60fastsim_galocco/files/macros/treeBkgEvents_40GeV_SiOnly_5hits.root", "output_fast/ImpPar_FAST_SiOnly_5hits.root");
        SaveImpParFast("/home/giacomo/na60fastsim_galocco/files/macros/treeBkgEvents_40GeV_Carbon_5hits.root", "output_fast/ImpPar_FAST_FullSetup_5hits.root");

        SavePResFast("/home/giacomo/na60fastsim_galocco/files/macros/bkgdistributions_40GeV_Void_5hits.root", "output_fast/PRes_FAST_Void_5hits.root");
        SavePResFast("/home/giacomo/na60fastsim_galocco/files/macros/bkgdistributions_40GeV_SiOnly_5hits.root", "output_fast/PRes_FAST_SiOnly_5hits.root");
        SavePResFast("/home/giacomo/na60fastsim_galocco/files/macros/bkgdistributions_40GeV_Carbon_5hits_175um.root", "output_fast/PRes_FAST_FullSetup_5hits.root");
        SavePResFast("/home/giacomo/na60fastsim_galocco/files/macros/bkgdistributions_40GeV_Carbon_5hits_compact.root", "output_fast/PRes_FAST_FullSetup_5hits_compact.root");

        SaveImpParACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_FullSetup", "output_acts/ImpPar_ACTS_FullSetup.root");
        SaveImpParACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_VTOnlySiOnly", "output_acts/ImpPar_ACTS_VTOnlySiOnly.root");
        SaveImpParACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_VTOnly", "output_acts/ImpPar_ACTS_VTOnly.root");
        SaveImpParACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_NoMaterial", "output_acts/ImpPar_ACTS_NoMaterial.root");

        SavePResACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_FullSetup", "output_acts/PRes_ACTS_FullSetup.root", true);
        SavePResACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_VTOnly", "output_acts/PRes_ACTS_VTOnly.root", true);
        SavePResACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_VTOnlyTest", "output_acts/PRes_ACTS_VTOnly.root", true);
        SavePResACTS("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_events_40GeV_muons_NoMaterial", "output_acts/PRes_ACTS_NoMaterial.root", true);
    }

    ////////////////////////////////////////////

    TFile *PRes_FAST_Void = new TFile("output_fast/PRes_FAST_Void_5hits.root", "read");
    TFile *PRes_FAST_SiOnly = new TFile("output_fast/PRes_FAST_SiOnly_5hits.root", "read");
    TFile *PRes_FAST_FullSetup = new TFile("output_fast/PRes_FAST_FullSetup_5hits.root", "read");

    TFile *PRes_ACTS_FullSetup = new TFile("output_acts/PRes_ACTS_VTOnly.root", "read");
    TFile *PRes_ACTS_VTOnlySiOnly = new TFile("output_acts/PRes_ACTS_VTOnlySiOnly.root", "read");
    TFile *PRes_ACTS_VTOnly = new TFile("output_acts/PRes_ACTS_VTOnly.root", "read");
    TFile *PRes_ACTS_NoMaterial = new TFile("output_acts/PRes_ACTS_NoMaterial.root", "read");

    float resMaxvsP = 0.02;
    float resMaxvsEta = 0.05;
    float maxP = 10.0; // GeV/c

    TGraphErrors *PRes_FAST_Void_graph_vs_p = (TGraphErrors *)PRes_FAST_Void->Get("ResPVsP");
    PRes_FAST_Void_graph_vs_p->SetName("PRes_FAST_Void_graph_vs_p");
    PRes_FAST_Void_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PRes_FAST_Void_graph_vs_p->GetYaxis()->SetRangeUser(0., resMaxvsP);
    PRes_FAST_Void_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PRes_FAST_SiOnly_graph_vs_p = (TGraphErrors *)PRes_FAST_SiOnly->Get("ResPVsP");
    PRes_FAST_SiOnly_graph_vs_p->SetName("PRes_FAST_SiOnly_graph_vs_p");
    PRes_FAST_SiOnly_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PRes_FAST_SiOnly_graph_vs_p->GetYaxis()->SetRangeUser(0., resMaxvsP);
    PRes_FAST_SiOnly_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PRes_FAST_FullSetup_graph_vs_p = (TGraphErrors *)PRes_FAST_FullSetup->Get("ResPVsP");
    PRes_FAST_FullSetup_graph_vs_p->SetName("PRes_FAST_FullSetup_graph_vs_p");
    PRes_FAST_FullSetup_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PRes_FAST_FullSetup_graph_vs_p->GetYaxis()->SetRangeUser(0., resMaxvsP);
    PRes_FAST_FullSetup_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);


    TGraphErrors *PRes_ACTS_FullSetup_graph_vs_p = (TGraphErrors *)PRes_ACTS_FullSetup->Get("ResPVsP");
    PRes_ACTS_FullSetup_graph_vs_p->SetName("PRes_ACTS_FullSetup_graph_vs_p");
    PRes_ACTS_FullSetup_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_FullSetup_graph_vs_p->GetYaxis()->SetRangeUser(0., resMaxvsP);
    PRes_ACTS_FullSetup_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PRes_ACTS_SiOnly_graph_vs_p = (TGraphErrors *)PRes_ACTS_VTOnlySiOnly->Get("ResPVsP");
    PRes_ACTS_SiOnly_graph_vs_p->SetName("PRes_ACTS_SiOnly_graph_vs_p");
    PRes_ACTS_SiOnly_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_SiOnly_graph_vs_p->GetYaxis()->SetRangeUser(0., resMaxvsP);
    PRes_ACTS_SiOnly_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PRes_ACTS_VTOnly_graph_vs_p = (TGraphErrors *)PRes_ACTS_VTOnly->Get("ResPVsP");
    PRes_ACTS_VTOnly_graph_vs_p->SetName("PRes_ACTS_VTOnly_graph_vs_p");
    PRes_ACTS_VTOnly_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_VTOnly_graph_vs_p->GetYaxis()->SetRangeUser(0., resMaxvsP);
    PRes_ACTS_VTOnly_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PRes_ACTS_NoMaterial_graph_vs_p = (TGraphErrors *)PRes_ACTS_NoMaterial->Get("ResPVsP");
    PRes_ACTS_NoMaterial_graph_vs_p->SetName("PRes_ACTS_NoMaterial_graph_vs_p");
    PRes_ACTS_NoMaterial_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_NoMaterial_graph_vs_p->GetYaxis()->SetRangeUser(0., resMaxvsP);
    PRes_ACTS_NoMaterial_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PRes_FAST_Void_graph_vs_eta = (TGraphErrors *)PRes_FAST_Void->Get("ResPVsEta");
    PRes_FAST_Void_graph_vs_eta->SetName("PRes_FAST_Void_graph_vs_eta");
    PRes_FAST_Void_graph_vs_eta->SetTitle(";#eta;#sigma_{#it{p}}/#it{p}");
    PRes_FAST_Void_graph_vs_eta->GetYaxis()->SetRangeUser(0., resMaxvsEta);


    TGraphErrors *PRes_FAST_SiOnly_graph_vs_eta = (TGraphErrors *)PRes_FAST_SiOnly->Get("ResPVsEta");
    PRes_FAST_SiOnly_graph_vs_eta->SetName("PRes_FAST_SiOnly_graph_vs_eta");
    PRes_FAST_SiOnly_graph_vs_eta->SetTitle(";#eta;#sigma_{#it{p}}/#it{p}");
    PRes_FAST_SiOnly_graph_vs_eta->GetYaxis()->SetRangeUser(0., resMaxvsEta);


    TGraphErrors *PRes_FAST_FullSetup_graph_vs_eta = (TGraphErrors *)PRes_FAST_FullSetup->Get("ResPVsEta");
    PRes_FAST_FullSetup_graph_vs_eta->SetName("PRes_FAST_FullSetup_graph_vs_eta");
    PRes_FAST_FullSetup_graph_vs_eta->SetTitle(";#eta;#sigma_{#it{p}}/#it{p}");
    PRes_FAST_FullSetup_graph_vs_eta->GetYaxis()->SetRangeUser(0., resMaxvsEta);


    TGraphErrors *PRes_ACTS_FullSetup_graph_vs_eta = (TGraphErrors *)PRes_ACTS_FullSetup->Get("ResPVsEta");
    PRes_ACTS_FullSetup_graph_vs_eta->SetName("PRes_ACTS_FullSetup_graph_vs_eta");
    PRes_ACTS_FullSetup_graph_vs_eta->SetTitle(";#eta;#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_FullSetup_graph_vs_eta->GetYaxis()->SetRangeUser(0., resMaxvsEta);


    TGraphErrors *PRes_ACTS_SiOnly_graph_vs_eta = (TGraphErrors *)PRes_ACTS_VTOnlySiOnly->Get("ResPVsEta");
    PRes_ACTS_SiOnly_graph_vs_eta->SetName("PRes_ACTS_SiOnly_graph_vs_eta");
    PRes_ACTS_SiOnly_graph_vs_eta->SetTitle(";#eta;#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_SiOnly_graph_vs_eta->GetYaxis()->SetRangeUser(0., resMaxvsEta);


    TGraphErrors *PRes_ACTS_VTOnly_graph_vs_eta = (TGraphErrors *)PRes_ACTS_VTOnly->Get("ResPVsEta");
    PRes_ACTS_VTOnly_graph_vs_eta->SetName("PRes_ACTS_VTOnly_graph_vs_eta");
    PRes_ACTS_VTOnly_graph_vs_eta->SetTitle(";#eta;#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_VTOnly_graph_vs_eta->GetYaxis()->SetRangeUser(0., resMaxvsEta);


    TGraphErrors *PRes_ACTS_NoMaterial_graph_vs_eta = (TGraphErrors *)PRes_ACTS_NoMaterial->Get("ResPVsEta");
    PRes_ACTS_NoMaterial_graph_vs_eta->SetName("PRes_ACTS_NoMaterial_graph_vs_eta");
    PRes_ACTS_NoMaterial_graph_vs_eta->SetTitle(";#eta;#sigma_{#it{p}}/#it{p}");
    PRes_ACTS_NoMaterial_graph_vs_eta->GetYaxis()->SetRangeUser(0., resMaxvsEta);


    TCanvas *c1 = new TCanvas("c1", "");

    c1->Update();
    c1->Modified();

    PRes_FAST_SiOnly_graph_vs_eta->SetMarkerColor(kBlack);
    PRes_ACTS_SiOnly_graph_vs_eta->SetMarkerColor(kBlack);

    PRes_FAST_FullSetup_graph_vs_eta->SetMarkerColor(kBlue);
    PRes_ACTS_VTOnly_graph_vs_eta->SetMarkerColor(kBlue);

    PRes_FAST_Void_graph_vs_eta->SetMarkerColor(kRed);
    PRes_ACTS_NoMaterial_graph_vs_eta->SetMarkerColor(kRed);

    PRes_FAST_SiOnly_graph_vs_eta->SetMarkerStyle(20);
    PRes_ACTS_SiOnly_graph_vs_eta->SetMarkerStyle(24);

    PRes_FAST_FullSetup_graph_vs_eta->SetMarkerStyle(21);
    PRes_ACTS_VTOnly_graph_vs_eta->SetMarkerStyle(25);

    PRes_FAST_Void_graph_vs_eta->SetMarkerStyle(22);
    PRes_ACTS_NoMaterial_graph_vs_eta->SetMarkerStyle(26);

    //////////////////////////////

    PRes_FAST_SiOnly_graph_vs_eta->Draw("AP");
    PRes_ACTS_SiOnly_graph_vs_eta->Draw("P SAME");

    PRes_FAST_FullSetup_graph_vs_eta->Draw("P SAME");
    PRes_ACTS_VTOnly_graph_vs_eta->Draw("P SAME");

    PRes_FAST_Void_graph_vs_eta->Draw("P SAME");
    PRes_ACTS_NoMaterial_graph_vs_eta->Draw("P SAME");

    TLegend *leg = new TLegend(0.35, 0.65, 0.7, 0.95);

    leg->AddEntry(PRes_FAST_SiOnly_graph_vs_eta, "FAST Si only", "p");
    leg->AddEntry(PRes_ACTS_SiOnly_graph_vs_eta, "ACTS Si only", "p");
    leg->AddEntry(PRes_FAST_FullSetup_graph_vs_eta, "FAST full setup", "p");
    leg->AddEntry(PRes_ACTS_VTOnly_graph_vs_eta, "ACTS full setup", "p");
    leg->AddEntry(PRes_FAST_Void_graph_vs_eta, "FAST No material budget", "p");
    leg->AddEntry(PRes_ACTS_NoMaterial_graph_vs_eta, "ACTS No material budget", "p");

    leg->Draw("SAME");

    c1->Update();
    c1->SaveAs("PResVsEta.png");
    c1->SaveAs("PResVsEta.pdf");

    PRes_FAST_SiOnly_graph_vs_p->SetMarkerColor(kBlack);
    PRes_ACTS_SiOnly_graph_vs_p->SetMarkerColor(kBlack);

    PRes_FAST_FullSetup_graph_vs_p->SetMarkerColor(kBlue);
    PRes_ACTS_VTOnly_graph_vs_p->SetMarkerColor(kBlue);

    PRes_FAST_Void_graph_vs_p->SetMarkerColor(kRed);
    PRes_ACTS_NoMaterial_graph_vs_p->SetMarkerColor(kRed);

    PRes_FAST_SiOnly_graph_vs_p->SetMarkerStyle(20);
    PRes_ACTS_SiOnly_graph_vs_p->SetMarkerStyle(24);

    PRes_FAST_FullSetup_graph_vs_p->SetMarkerStyle(21);
    PRes_ACTS_VTOnly_graph_vs_p->SetMarkerStyle(25);

    PRes_FAST_Void_graph_vs_p->SetMarkerStyle(22);
    PRes_ACTS_NoMaterial_graph_vs_p->SetMarkerStyle(26);

    //////////////////////////////

    PRes_FAST_SiOnly_graph_vs_p->Draw("AP");
    PRes_ACTS_SiOnly_graph_vs_p->Draw("P SAME");

    PRes_FAST_FullSetup_graph_vs_p->Draw("P SAME");
    PRes_ACTS_VTOnly_graph_vs_p->Draw("P SAME");

    PRes_FAST_Void_graph_vs_p->Draw("P SAME");
    PRes_ACTS_NoMaterial_graph_vs_p->Draw("P SAME");

    leg->Draw("SAME");

    c1->Update();
    c1->SaveAs("PResVsP.png");
    c1->SaveAs("PResVsP.pdf");

    TH1D* hResidPVsEtaProjACTS = (TH1D*) PRes_ACTS_VTOnlySiOnly->Get("hResidPVsEtaProj");
    TH1D* hResidPVsPProjACTS = (TH1D*) PRes_ACTS_VTOnlySiOnly->Get("hResidPVsPProj");
    TH1D* hResidPVsEtaProjFAST = (TH1D*) PRes_FAST_SiOnly->Get("hResidPVsEtaProj");
    TH1D* hResidPVsPProjFAST = (TH1D*) PRes_FAST_SiOnly->Get("hResidPVsPProj");

    hResidPVsPProjFAST->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hResidPVsPProjFAST->GetYaxis()->SetTitle("Counts");
    hResidPVsPProjFAST->Draw();
    hResidPVsPProjACTS->SetLineColor(kRed);
    hResidPVsPProjACTS->Draw("SAME");
    hResidPVsPProjFAST->GetYaxis()->SetRangeUser(0., std::max(hResidPVsPProjFAST->GetMaximum(), hResidPVsPProjACTS->GetMaximum()) * 1.2);
    TLegend *leg2 = new TLegend(0.65, 0.65, 0.95, 0.95);
    leg2->AddEntry(hResidPVsPProjFAST, "FAST", "l");
    leg2->AddEntry(hResidPVsPProjACTS, "ACTS", "l");

    leg2->Draw("SAME");

    c1->Update();
    c1->SaveAs("RecPdistr.png");
    c1->SaveAs("RecPdistr.pdf");

    hResidPVsEtaProjFAST->GetXaxis()->SetTitle("#eta");
    hResidPVsEtaProjFAST->GetYaxis()->SetTitle("Counts");
    hResidPVsEtaProjFAST->Draw();
    hResidPVsEtaProjACTS->SetLineColor(kRed);
    hResidPVsEtaProjACTS->Draw("SAME");
    hResidPVsEtaProjFAST->GetYaxis()->SetRangeUser(0., std::max(hResidPVsEtaProjFAST->GetMaximum(), hResidPVsEtaProjACTS->GetMaximum()) * 1.2);

    leg2->Draw("SAME");
    c1->Update();
    c1->SaveAs("RecEtadistr.png");
    c1->SaveAs("RecEtadistr.pdf");


    //setCanvasForTH2D();
    TH1D* hRecPVsEtaACTS = (TH1D*) PRes_ACTS_VTOnlySiOnly->Get("hRecPVsEta");
    hRecPVsEtaACTS->SetTitle(";#eta;#it{p} (GeV/c)");
    hRecPVsEtaACTS->Draw("colz");
    c1->Update();
    //setCanvasForTH2D();
    c1->SaveAs("RecPVsEtaACTS_SiOnly.png");
    TH1D* hRecPVsEtaFAST = (TH1D*) PRes_FAST_SiOnly->Get("hRecPVsEta");
    hRecPVsEtaFAST->SetTitle(";#eta;#it{p} (GeV/c)");
    hRecPVsEtaFAST->Draw("colz");
    //setCanvasForTH2D();
    c1->Update();
    c1->SaveAs("RecPVsEtaFAST_SiOnly.png");

    //setCanvasForTH2D();
    TH1D* hRecPVsEtaACTS_Void = (TH1D*) PRes_ACTS_NoMaterial->Get("hRecPVsEta");
    hRecPVsEtaACTS_Void->SetTitle(";#eta;#it{p} (GeV/c)");
    hRecPVsEtaACTS_Void->Draw("colz");
    c1->Update();
    //setCanvasForTH2D();
    c1->SaveAs("RecPVsEtaACTS_Void.png");
    TH1D* hRecPVsEtaFAST_Void = (TH1D*) PRes_FAST_Void->Get("hRecPVsEta");
    hRecPVsEtaFAST_Void->SetTitle(";#eta;#it{p} (GeV/c)");
    hRecPVsEtaFAST_Void->Draw("colz");
    //setCanvasForTH2D();
    c1->Update();
    c1->SaveAs("RecPVsEtaFAST_Void.png");


    //setCanvasForTH2D();
    TH1D* hRecPVsEtaACTS_FullSetup = (TH1D*) PRes_ACTS_FullSetup->Get("hRecPVsEta");
    hRecPVsEtaACTS_FullSetup->SetTitle(";#eta;#it{p} (GeV/c)");
    hRecPVsEtaACTS_FullSetup->Draw("colz");
    c1->Update();
    //setCanvasForTH2D();
    c1->SaveAs("RecPVsEtaACTS_FullSetup.png");
    TH1D* hRecPVsEtaFAST_FullSetup = (TH1D*) PRes_FAST_FullSetup->Get("hRecPVsEta");
    hRecPVsEtaFAST_FullSetup->SetTitle(";#eta;#it{p} (GeV/c)");
    hRecPVsEtaFAST_FullSetup->Draw("colz");
    //setCanvasForTH2D();
    c1->Update();
    c1->SaveAs("RecPVsEtaFAST_FullSetup.png");
    gROOT->SetBatch(false);

    //////////////////////////////
    /*
    TGraphErrors *PMean_FAST_Void_graph_vs_p = (TGraphErrors *)PRes_FAST_Void->Get("MeanPVsP");
    PMean_FAST_Void_graph_vs_p->SetName("PMean_FAST_Void_graph_vs_p");
    PMean_FAST_Void_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PMean_FAST_Void_graph_vs_p->GetYaxis()->SetRangeUser(-resMaxvsP, resMaxvsP);
    PMean_FAST_Void_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PMean_FAST_SiOnly_graph_vs_p = (TGraphErrors *)PRes_FAST_SiOnly->Get("MeanPVsP");
    PMean_FAST_SiOnly_graph_vs_p->SetName("PMean_FAST_SiOnly_graph_vs_p");
    PMean_FAST_SiOnly_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PMean_FAST_SiOnly_graph_vs_p->GetYaxis()->SetRangeUser(-resMaxvsP, resMaxvsP);
    PMean_FAST_SiOnly_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PMean_FAST_FullSetup_graph_vs_p = (TGraphErrors *)PRes_FAST_FullSetup->Get("MeanPVsP");
    PMean_FAST_FullSetup_graph_vs_p->SetName("PMean_FAST_FullSetup_graph_vs_p");
    PMean_FAST_FullSetup_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PMean_FAST_FullSetup_graph_vs_p->GetYaxis()->SetRangeUser(-resMaxvsP, resMaxvsP);
    PMean_FAST_FullSetup_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PMean_ACTS_SiOnly_graph_vs_p = (TGraphErrors *)PRes_ACTS_VTOnlySiOnly->Get("MeanPVsP");
    PMean_ACTS_SiOnly_graph_vs_p->SetName("PMean_ACTS_SiOnly_graph_vs_p");
    PMean_ACTS_SiOnly_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PMean_ACTS_SiOnly_graph_vs_p->GetYaxis()->SetRangeUser(-resMaxvsP, resMaxvsP);
    PMean_ACTS_SiOnly_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PMean_ACTS_VTOnly_graph_vs_p = (TGraphErrors *)PRes_ACTS_VTOnly->Get("MeanPVsP");
    PMean_ACTS_VTOnly_graph_vs_p->SetName("PMean_ACTS_VTOnly_graph_vs_p");
    PMean_ACTS_VTOnly_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PMean_ACTS_VTOnly_graph_vs_p->GetYaxis()->SetRangeUser(-resMaxvsP, resMaxvsP);
    PMean_ACTS_VTOnly_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);

    TGraphErrors *PMean_ACTS_NoMaterial_graph_vs_p = (TGraphErrors *)PRes_ACTS_NoMaterial->Get("MeanPVsP");
    PMean_ACTS_NoMaterial_graph_vs_p->SetName("PMean_ACTS_NoMaterial_graph_vs_p");
    PMean_ACTS_NoMaterial_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
    PMean_ACTS_NoMaterial_graph_vs_p->GetYaxis()->SetRangeUser(-resMaxvsP, resMaxvsP);
    PMean_ACTS_NoMaterial_graph_vs_p->GetXaxis()->SetRangeUser(0., maxP);


    TLegend *leg3 = new TLegend(0.35, 0.65, 0.7, 0.95);

    leg3->AddEntry(PMean_FAST_SiOnly_graph_vs_p, "FAST Si only", "p");
    leg3->AddEntry(PMean_ACTS_SiOnly_graph_vs_p, "ACTS Si only", "p");
    leg3->AddEntry(PMean_FAST_FullSetup_graph_vs_p, "FAST full setup", "p");
    leg3->AddEntry(PMean_ACTS_VTOnly_graph_vs_p, "ACTS full setup", "p");
    leg3->AddEntry(PMean_FAST_Void_graph_vs_p, "FAST No material budget", "p");
    leg3->AddEntry(PMean_ACTS_NoMaterial_graph_vs_p, "ACTS No material budget", "p");

    leg3->Draw("SAME");


    PMean_FAST_SiOnly_graph_vs_p->SetMarkerColor(kBlack);
    PMean_ACTS_SiOnly_graph_vs_p->SetMarkerColor(kBlack);

    PMean_FAST_FullSetup_graph_vs_p->SetMarkerColor(kBlue);
    PMean_ACTS_VTOnly_graph_vs_p->SetMarkerColor(kBlue);

    PMean_FAST_Void_graph_vs_p->SetMarkerColor(kRed);
    PMean_ACTS_NoMaterial_graph_vs_p->SetMarkerColor(kRed);

    PMean_FAST_SiOnly_graph_vs_p->SetMarkerStyle(20);
    PMean_ACTS_SiOnly_graph_vs_p->SetMarkerStyle(24);

    PMean_FAST_FullSetup_graph_vs_p->SetMarkerStyle(21);
    PMean_ACTS_VTOnly_graph_vs_p->SetMarkerStyle(25);

    PMean_FAST_Void_graph_vs_p->SetMarkerStyle(22);
    PMean_ACTS_NoMaterial_graph_vs_p->SetMarkerStyle(26);

    //////////////////////////////
    PMean_FAST_SiOnly_graph_vs_p->SetTitle(";#it{p} (GeV/#it{c});#mu_{#it{p}} (GeV/#it{c})");
    PMean_FAST_SiOnly_graph_vs_p->Draw("AP");
    PMean_ACTS_SiOnly_graph_vs_p->Draw("P SAME");

    PMean_FAST_FullSetup_graph_vs_p->Draw("P SAME");
    PMean_ACTS_VTOnly_graph_vs_p->Draw("P SAME");

    PMean_FAST_Void_graph_vs_p->Draw("P SAME");
    PMean_ACTS_NoMaterial_graph_vs_p->Draw("P SAME");

    leg->Draw("SAME");


    gStyle->SetTitleSize(0.07, "XYZ"); // Sets the size for X, Y, and Z axis titles
    gStyle->SetTitleFontSize(0.07);    // Sets the size for the histogram title
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadLeftMargin(0.2);   // Set left margin
    gStyle->SetPadRightMargin(0.01);  // Set right margin
    gStyle->SetPadTopMargin(0.03);    // Set top margin
    gStyle->SetPadBottomMargin(0.15); // Set bottom margin
    gStyle->SetHistLineWidth(2);
    c1->Update();
    c1->SaveAs("MeanPVsP.png");
    c1->SaveAs("MeanPVsP.pdf");
    */

}