import uproot
import pandas as pd
import ROOT
from ROOT import TF1, TH1, TH1D, TH2D, TFile, gDirectory
import matplotlib.pyplot as plt
import math
import numpy as np
from array import array

particle_list = ['K0s','phi']
for particle in particle_list:
    ROOT.gSystem.Exec('mkdir '+particle+'_folder')
    
#read the reconstruction efficiency
signal_file = ROOT.TFile('Signal_histos_'+particle+'.root')
hist_rec_eff = signal_file.Get('hPtEff')
#read the number of events for the background
background_file = ROOT.TFile('Bkg-histos_'+particle+'.root')
hist_ev = background_file.Get('hNevents')
n_ev = hist_ev.GetBinContent(1)

#test for the phi with Ebeam = 158 GeV
selection_list = ['cosp > 0.99','cosp > 0.999','cosp > 0.9999']
multiplicity = 8.46
PT_BINS = [0.5,1.0,1.5,2.0,2.5,3.0]
BINS = np.asarray(PT_BINS, dtype=np.float64) 
pdg_code = 333
Tslope = 0.200
mass = 1.019445 #TDatabasePDG::Instance()->GetParticle(pdg_code)->Mass()
fpt = ROOT.TF1("fpt","x*exp(-TMath::Sqrt(x**2+[0]**2)/[1])",0,10)
fpt.SetParameter(0,mass)
fpt.SetParameter(1,Tslope)

fPeak = ROOT.TF1("fPeak","[2]*exp(-0.5*((x-[0])/[1])**2)",1.05,1.05)
fPeak.SetParameter(0,mass)
fPeak.SetParLimits(0,1.015,1.025)
fPeak.SetParLimits(1,0.001,0.003)
fPeak.SetParLimits(2,100,3000)
fBkg = ROOT.TF1("fBkg","expo(0)")

normalization = fpt.Integral(0,10)
density_prob = []
particle_list = ['phi']
for particle in particle_list:
    signal_file = uproot.open('fntSig_'+particle+'.root')
    df_sig = signal_file["ntcand"].pandas.df()
    df_gen = signal_file["ntgen"].pandas.df()
    df_bkg = uproot.open('fntBkg_'+particle+'.root')["ntcand"].pandas.df()
    results_file = ROOT.TFile(particle+'_folder/'+particle+"_sign.root","recreate")
    
    counts_rec, _ = np.histogram(df_sig['pt'], len(PT_BINS)-1, range=[PT_BINS[0],PT_BINS[len(PT_BINS)-1]])
    hist_rec = ROOT.TH1D('hist_rec', ';#it{p}_{T} (GeV/c);counts', len(PT_BINS)-1, PT_BINS[0], PT_BINS[len(PT_BINS)-1])
    counts_gen, _ = np.histogram(df_gen['pt'], len(PT_BINS)-1, range=[PT_BINS[0],PT_BINS[len(PT_BINS)-1]])
    hist_gen = ROOT.TH1D('hist_gen', ';#it{p}_{T} (GeV/c);counts', len(PT_BINS)-1, PT_BINS[0], PT_BINS[len(PT_BINS)-1])
    hist_acceff = ROOT.TH1D('hist_acceff', ';#it{p}_{T} (GeV/c);efficiency', len(PT_BINS)-1 ,PT_BINS[0], PT_BINS[len(PT_BINS)-1])
    
    for index in range(0, len(PT_BINS)-1):
        hist_rec.SetBinContent(index + 1, counts_rec[index])
        hist_rec.SetBinError(index + 1, math.sqrt(counts_rec[index]))
        hist_gen.SetBinContent(index + 1, counts_gen[index])
        hist_gen.SetBinError(index + 1, math.sqrt(counts_gen[index]))
        eff = counts_rec[index]/counts_gen[index]
        hist_acceff.SetBinContent(index + 1, eff)
        hist_acceff.SetBinError(index + 1, math.sqrt(eff*(1-eff)/counts_gen[index]))
    hist_acceff.Write()
    for selection in selection_list:
        sel_dir = results_file.mkdir(selection)
        sel_dir.cd()       
            
        counts_sel, _ = np.histogram(df_sig.query(selection)['pt'], len(PT_BINS)-1, range=[PT_BINS[0],PT_BINS[len(PT_BINS)-1]])
        hist_eff = ROOT.TH1D('hist_eff', selection+';#it{p}_{T} (GeV/c);efficiency', len(PT_BINS)-1 ,PT_BINS[0], PT_BINS[len(PT_BINS)-1])
        hist_sigma = ROOT.TH1D('hist_sigma', selection+';#it{p}_{T} (GeV/c);#sigma (GeV/c^{2})', len(PT_BINS)-1 ,PT_BINS[0], PT_BINS[len(PT_BINS)-1])
        hist_significance = ROOT.TH1D('hist_significance', selection+';#it{p}_{T} (GeV/c);significance/#sqrt{n_{ev}}', len(PT_BINS)-1 ,PT_BINS[0], PT_BINS[len(PT_BINS)-1])
        print(counts_rec)
        for index in range(0, len(PT_BINS)-1):
            sub_dir = sel_dir.mkdir(f'pt_{PT_BINS[index]}_{PT_BINS[index + 1]}')
            sub_dir.cd()       
            
            eff = counts_sel[index]/counts_rec[index]
            hist_eff.SetBinContent(index + 1, eff)
            hist_eff.SetBinError(index + 1, math.sqrt(eff*(1-eff)/counts_rec[index]))
            density_prob.append(fpt.Integral(PT_BINS[index],PT_BINS[index + 1])/normalization)
            
            ptmin = PT_BINS[index]
            ptmax = PT_BINS[index+1]
            counts_sig, _ = np.histogram(df_sig.query(selection+'and pt > @ptmin and pt < @ptmax')['mass'], 40, range=[1, 1.05])
            counts_bkg, _ = np.histogram(df_bkg.query(selection+'and pt > @ptmin and pt < @ptmax')['mass'], 40, range=[1., 1.05])
            print(counts_sig)
            hist_sig = ROOT.TH1D('hist_sig', ';m (GeV/c^{2});counts', 40, 1.0, 1.05)
            hist_bkg = ROOT.TH1D('hist_bkg', ';m (GeV/c^{2});counts', 40, 1.0, 1.05)
            
            for index_mass in range(0, 40):
                hist_sig.SetBinContent(index_mass + 1, counts_sig[index_mass])
                hist_sig.SetBinError(index_mass + 1, math.sqrt(counts_sig[index_mass]))
                hist_bkg.SetBinContent(index_mass + 1, counts_bkg[index_mass])
                hist_bkg.SetBinError(index_mass + 1, math.sqrt(counts_bkg[index_mass]))
            
            hist_sig.Fit(fPeak,'MR+')
            mu = fPeak.GetParameter(0)
            sigma = fPeak.GetParameter(1)
            sigma_error = fPeak.GetParError(1)
            hist_sigma.SetBinContent(index + 1, sigma)
            hist_sigma.SetBinError(index + 1, sigma_error)
            bin_width = 0.05/40
            hist_bkg.Fit(fBkg,'M+')
            nbkg = fBkg.Integral(mu-3*sigma,mu+3*sigma)/bin_width
            nsig = multiplicity*density_prob[index]*eff*n_ev*hist_acceff.GetBinContent(index+1)
            print(mu-3*sigma,"--",mu+3*sigma)
            print(nsig)
            print(nbkg)
            print(density_prob)
            significance = nsig/math.sqrt(nsig+nbkg)
            sig_over_ev = significance/math.sqrt(n_ev)
            hist_significance.SetBinContent(index+1,sig_over_ev)
            
            hist_sig.Write()
            hist_bkg.Write()
            
        sel_dir.cd()
        hist_significance.Write()
        hist_rec.Write()
        hist_eff.Write()
        hist_sigma.Write()
    results_file.Close()