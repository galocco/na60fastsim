import uproot
import pandas as pd
import ROOT
from ROOT import TF1, TH1, TH1D, TH2D, TFile, gDirectory
import matplotlib.pyplot as plt
import math
import numpy as np

particle_list = ['K0s','phi']
for particle in particle_list:
    ROOT.gSystem.Exec('mkdir '+particle+'_folder')

ROOT.gROOT.ForceStyle()
ROOT.gStyle.SetOptStat(0)
hist_list = ['hPtGen','hYGen','hptdauPos','hptdauNeg','hydauPos','hydauNeg']
cv = ROOT.TCanvas("cv","cv")
#draw y,pT distributions of mother particles its daughters
for particle in particle_list:
    input_file = ROOT.TFile('DecayHistos_'+particle+'.root')
    for hist in hist_list:
        histogram = input_file.Get(hist)
        if 'hptdau' in hist:
            histogram = histogram.ProjectionX()
        if 'y' in hist:
            histogram.GetXaxis().SetRangeUser(3-5,3+5)
        if 'Y' in hist:
            histogram.GetXaxis().SetRangeUser(3-5,3+5)
        histogram.SetMarkerColor(ROOT.kBlue)
        histogram.SetLineColor(ROOT.kBlue)
        histogram.SetMarkerStyle(23)
        histogram.GetYaxis().SetTitle('counts')
        histogram.SetTitle('')
        histogram.SetMarkerStyle(23)
        histogram.Draw("e")
        cv.SetLeftMargin(10)
        cv.SaveAs(particle+'_folder/'+particle+'_'+hist+'.png')

#compute efficiencies and resolutions
hist_list = ['hResVx','hResVy','hResVz','hResPt']
title = ['#it{x}_{rec}-#it{x}_{gen}','#it{y}_{rec}-#it{y}_{gen}','#it{z}_{rec}-#it{z}_{gen}','#it{p}_{Trec}-#it{p}_{Tgen}']
unit = ['(cm)','(cm)','(cm)','(GeV/#it{c})']
for particle in particle_list:
    input_file = ROOT.TFile('Signal_histos_'+particle+'.root')
    
    histogram = input_file.Get('hYEff')
    histogram.SetMarkerColor(ROOT.kBlue)
    histogram.SetLineColor(ROOT.kBlue)
    histogram.SetMarkerStyle(23)
    histogram.GetXaxis().SetRangeUser(0,5)
    histogram.SetTitle(';y;efficiency x acceptance')
    histogram.Draw('e')
    cv.SaveAs(particle+'_folder/'+particle+'_YEff.png')
    
    histogram = input_file.Get('hPtEff')
    histogram.SetMarkerColor(ROOT.kBlue)
    histogram.SetLineColor(ROOT.kBlue)
    histogram.SetMarkerStyle(23)
    histogram.GetYaxis().SetTitle('efficiency x acceptance')
    histogram.SetTitle('')
    histogram.Draw('e')
    cv.SaveAs(particle+'_folder/'+particle+'_PtEff.png')
    output_file = ROOT.TFile(particle+'_folder/'+particle+'_resolution.root','recreate')
    counter = 0
    for hist in hist_list:
        histogram = input_file.Get(hist)
        biny = 12#histogram.GetNbinsY()
        low_edge = 0#histogram.GetYaxis().GetXmin();
        upper_edge = 6#histogram.GetYaxis().GetXmax();
        
        hist_mean = ROOT.TH1D(hist+'_mean', title[counter]+' mean;#it{p}_{T} (GeV/#it{c});mean'+unit[counter],biny,low_edge,upper_edge)
        hist_dev = ROOT.TH1D(hist+'_dev', title[counter]+' std dev;#it{p}_{T} (GeV/#it{c});std dev'+unit[counter],biny,low_edge,upper_edge)
        for ptbin in range(1,biny+1): 
            proj = histogram.ProjectionX("projection",ptbin,ptbin+1)
            if counter == 3:
                const = 1
            else:
                const = 1./10000
            mean = proj.GetMean()*const
            std = proj.GetStdDev()*const
            mean_err = proj.GetMeanError()*const
            std_err = proj.GetStdDevError()*const
            hist_mean.SetBinContent(ptbin,mean)
            hist_mean.SetBinError(ptbin,mean_err)
            hist_dev.SetBinContent(ptbin,std)
            hist_dev.SetBinError(ptbin,std_err)
        output_file.cd()
        hist_mean.Write()
        hist_mean.Draw()
        cv.SetLeftMargin(0.15)
        cv.SaveAs(particle+'_folder/'+particle+'_'+hist+'_mean.png')
        hist_dev.Write()
        hist_dev.Draw()
        cv.SetLeftMargin(0.15)
        cv.SaveAs(particle+'_folder/'+particle+'_'+hist+'_std.png')
        counter = counter + 1
    input_file.Close()

#draw the comparison between the distributions of  the topological variables of signal and background
columns = ['pt','y','dist','cosp','d01','d02','d0prod','dca','ptMin','ptMax']
xlabel = ['#it{p}_{T} (GeV/c)','y','d (cm)','cosp','d_{xy+} (cm)','d_{xy-} (cm)','d_{xy+}d_{xy-} (cm^{2})','dca (cm)',
           '#it{p}_{Tmin} (GeV/c)','#it{p}_{Tmax} (GeV/c)']
min_value = [0.,  3-5, 0., 0.9, -0.02, -0.02, -0.00004, 0,  0.,  0.]
max_value = [10., 3+5, 30, 1.0,  0.02,  0.02,  0.00004, 0.01, 30., 30.]
for particle in particle_list:
    df_sig = uproot.open('fntSig_'+particle+'.root')["ntcand"].pandas.df()
    df_bkg = uproot.open('fntBkg_'+particle+'.root')["ntcand"].pandas.df()
    results_file = ROOT.TFile(particle+'_folder/'+particle+"_variables.root","recreate")
    counter = 0
    for variable in columns:
        min_value.append(min(df_sig[variable].min(),df_bkg[variable].min()))
        max_value.append(max(df_sig[variable].max(),df_bkg[variable].max()))
        counts_sig, _ = np.histogram(df_sig[variable], bins=80, range=[min_value[counter],max_value[counter]])
        counts_bkg, _ = np.histogram(df_bkg[variable], bins=80, range=[min_value[counter],max_value[counter]])
        
        
        hist_sig = ROOT.TH1D('hist_'+variable+'_sig', ';'+xlabel[counter]+';pdf',80,min_value[counter],max_value[counter])
        hist_sig.SetMarkerColor(ROOT.kRed)
        hist_sig.SetLineColor(ROOT.kRed)
        hist_sig.SetMarkerStyle(23)
        for index in range(0, 80):
            hist_sig.SetBinContent(index+1, counts_sig[index])
            hist_sig.SetBinError(index + 1, math.sqrt(counts_sig[index]))
        hist_bkg = ROOT.TH1D('hist_'+variable+'_bkg', ';'+xlabel[counter]+';pdf',80,min_value[counter],max_value[counter])
        hist_bkg.SetMarkerColor(ROOT.kBlue)
        hist_bkg.SetLineColor(ROOT.kBlue)
        hist_bkg.SetMarkerStyle(23)
        for index in range(0, 80):
            hist_bkg.SetBinContent(index+1, counts_bkg[index])
            hist_bkg.SetBinError(index + 1, math.sqrt(counts_bkg[index]))
        hist_sig.Write()
        hist_bkg.Write()
        cv = ROOT.TCanvas("cv_"+variable,"cv_"+variable)    
        maximum = max(hist_sig.GetBinContent(hist_sig.GetMaximumBin())/sum(counts_sig),hist_bkg.GetBinContent(hist_bkg.GetMaximumBin())/sum(counts_bkg))
        hist_sig.Scale(1./sum(counts_sig))
        hist_bkg.Scale(1./sum(counts_bkg))
        hist_sig.SetAxisRange(0., maximum*1.5,"Y")
        hist_bkg.SetAxisRange(0., maximum*1.5,"Y")
        hist_sig.Draw()
        hist_bkg.Draw("same")
        legend = ROOT.TLegend(0.4,0.7,0.7,0.9)
        legend.SetHeader("Candidates "+particle,"C")
        legend.AddEntry('hist_'+variable+'_bkg',"background","p")
        legend.AddEntry('hist_'+variable+'_sig',"signal","p")
        legend.Draw()
        cv.SaveAs(particle+'_folder/'+particle+'_'+variable+'.png')
        cv.Write()
        counter = counter + 1
    results_file.Close()