//To compare MaCh3 and P-theta sensitivities for t2k-nova joint analysis
#include <iostream>
#include <string.h>
#include "draw1Dcreds_3comp.h"

bool nsteps = false;

struct Binning {
    Int_t nbins;
    Double_t bin_min;
    Double_t bin_max;
    const char *xtitle;
    const char *date;
};

void SetHistStyle(TH1* h) {
    h->GetXaxis()->SetTitleSize(0);
    h->GetXaxis()->SetLabelSize(0);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.40);
}

void SetRatioStyle(TH1* h, const char* title) {
    h->SetTitle("");
    h->GetXaxis()->SetTitle(title);
    h->GetXaxis()->SetTitleSize(0.11);
    h->GetXaxis()->SetLabelSize(0.10);
    h->GetXaxis()->SetTitleOffset(1.00);
    
    h->GetYaxis()->SetTitle("Ratio");
    h->GetYaxis()->SetTitleSize(0.10);
    h->GetYaxis()->SetLabelSize(0.10);
    h->GetYaxis()->SetTitleOffset(0.70);
    h->GetYaxis()->SetNdivisions(2000510); //309
    std::cout << 0.8*h->GetBinContent(h->GetMinimumBin()) << " " << 1.2*h->GetBinContent(h->GetMaximumBin()) << std::endl;
    h->GetYaxis()->SetRangeUser(std::max(0.0, 0.7*h->GetBinContent(h->GetMinimumBin())), std::min(2.5, 1.3*h->GetBinContent(h->GetMaximumBin())));
}

void SetTopPadStyle(TPad* pad) {
    pad->SetBottomMargin(0.038);
    //pad->SetRightMargin(0.05);
    pad->SetLeftMargin(0.14);
    //pad->SetRightMargin(0.005730659);
    pad->SetTickx();
    pad->SetTicky();
    pad->Draw();             // Draw the upper pad: pad1
    pad->cd();
}

void SetBottomPadStyle(TPad* pad) {
    pad->SetBottomMargin(0.258);
    pad->SetLeftMargin(0.14);
    pad->SetTopMargin(0.045);
    //pad->SetLeftMargin(0.005730659);
    //pad->SetRightMargin(0.1919771);
    pad->SetTickx(1);
    //pad->SetTopMargin(0.05649718);
    //pad->SetLeftMargin(0.1638418);
    pad->SetTicky();
    pad->SetGridy();
    pad->Draw();
    pad->cd();
}

void GetRatioHist(TH1D* hp, TH1D* hm, TH1D* ha, TH1D* hmb, TLegend *leg, Binning binning, bool isDchi2, const char* output_name) {
    TCanvas *c=new TCanvas("","",700,700);
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.346, 1, 1.0);
    SetTopPadStyle(pad1);
    
    Double_t bm = hm->GetBinContent(hm->GetMaximumBin());
    Double_t bp = hp->GetBinContent(hp->GetMaximumBin());
    Double_t ba = ha->GetBinContent(ha->GetMaximumBin());
    Double_t bmb = hmb->GetBinContent(hmb->GetMaximumBin());
    
    std::vector<Double_t> cont_vec={bm, bp, ba, bmb};
    Double_t bcmax = *std::max_element(cont_vec.begin(), cont_vec.end());
    std::cout << bm << " " << bp << " " << ba << " " << bmb << " " << bcmax << std::endl;
    
    if(isDchi2)  hm->GetYaxis()->SetRangeUser(0,10.0);
    else hm->GetYaxis()->SetRangeUser(0.0, 1.2*bcmax);
    //else hm->GetYaxis()->SetRangeUser(std::max(0.0, 0.8*hm->GetBinContent(hm->GetMinimumBin())), 1.2*bcmax);
    
    
    SetHistStyle(hm);
    
    hm->SetLineColor(kBlue);
    ha->SetLineColor(kGreen);
    hp->SetLineColor(kRed);
    hmb->SetLineColor(kBlack);
    
    hm->SetLineWidth(2);
    ha->SetLineWidth(2);
    hp->SetLineWidth(2);
    hmb->SetLineWidth(2);
    
    hm->Draw("HIST");
    hmb->Draw("HIST same");
    ha->Draw("HIST same");
    hp->Draw("HIST same");
    
    leg->Draw("same");
    
    std::vector<std::array<double,2> > mdm2_cred1sig = get1Dcred_disjoint(hm,0.6826);
    std::vector<std::array<double,2> > mdm2_cred2sig = get1Dcred_disjoint(hm,0.9544);
    std::vector<std::array<double,2> > mbdm2_cred1sig = get1Dcred_disjoint(hmb,0.6826);
    std::vector<std::array<double,2> > mbdm2_cred2sig = get1Dcred_disjoint(hmb,0.9544);
    std::vector<std::array<double,2> > adm2_cred1sig = get1Dcred_disjoint(ha,0.6826);
    std::vector<std::array<double,2> > adm2_cred2sig = get1Dcred_disjoint(ha,0.9544);
    std::vector<std::array<double,2> > pdm2_cred1sig = get1Dcred_disjoint(hp,0.6826);
    std::vector<std::array<double,2> > pdm2_cred2sig = get1Dcred_disjoint(hp,0.9544);
    
    std::cout << " dm2_cred1sig.size() " << mdm2_cred1sig.size() << " " << " dm2_cred2sig.size() " << mdm2_cred2sig.size() << std::endl;
    for(int i=0;i<mdm2_cred1sig.size();++i) {
        std::cout << " 1 sig " << i << " " << mdm2_cred1sig[i][0] << " " << mdm2_cred1sig[i][1] << std::endl;
    }
    
    for(int i=0;i<mdm2_cred2sig.size();++i) {
        std::cout << " 2 sig " << i << " " << mdm2_cred2sig[i][0] << " " << mdm2_cred2sig[i][1] << std::endl;
    }
    
    //Double_t mh_min = hm->GetBinContent(hm->GetMaximumBin());
    Double_t mh_min = std::max(0.0, 0.8*hm->GetBinContent(hm->GetMinimumBin()));
    Double_t mh_max = hm->GetBinContent(hm->GetMaximumBin());
    cout << " mh_min " << mh_min << " mh_max " << mh_max << endl;
    
    if(!isDchi2) {
        for(unsigned int a=0;a<mdm2_cred1sig.size();a++) {
            drawCredInterval(pad1,mdm2_cred1sig[a][0],mdm2_cred1sig[a][1],0.9*mh_max,1.2*mh_max,hm->GetLineColor(),0.25,9);
        }
        
        for(unsigned int a=0;a<mdm2_cred2sig.size();a++) {
            drawCredInterval(pad1,mdm2_cred2sig[a][0],mdm2_cred2sig[a][1],0.9*mh_max,1.2*mh_max,hm->GetLineColor(),0.25,2);
        }
        
        for(unsigned int a=0;a<mbdm2_cred1sig.size();a++) {
            drawCredInterval(pad1,mbdm2_cred1sig[a][0],mbdm2_cred1sig[a][1],0.3*mh_max,0.6*mh_max,hmb->GetLineColor(),0.25,9);
        }
        
        for(unsigned int a=0;a<mbdm2_cred2sig.size();a++) {
            drawCredInterval(pad1,mbdm2_cred2sig[a][0],mbdm2_cred2sig[a][1],0.3*mh_max,0.6*mh_max,hmb->GetLineColor(),0.25,2);
        }
        
        for(unsigned int a=0;a<adm2_cred1sig.size();a++) {
            drawCredInterval(pad1,adm2_cred1sig[a][0],adm2_cred1sig[a][1],0.6*mh_max,0.9*mh_max,ha->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<adm2_cred2sig.size();a++) {
            drawCredInterval(pad1,adm2_cred2sig[a][0],adm2_cred2sig[a][1],0.6*mh_max,0.9*mh_max,ha->GetLineColor(),0.25,2);
        }
        
        for(unsigned int a=0;a<pdm2_cred1sig.size();a++) {
            drawCredInterval(pad1,pdm2_cred1sig[a][0],pdm2_cred1sig[a][1],0,0.3*mh_max,hp->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<pdm2_cred2sig.size();a++) {
            drawCredInterval(pad1,pdm2_cred2sig[a][0],pdm2_cred2sig[a][1],0,0.3*mh_max,hp->GetLineColor(),0.25,2);
        }
    }
    
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.358);
    SetBottomPadStyle(pad2);
    
    TH1D* hratio_mmb = (TH1D*)hm->Clone();
    TH1D* hratio_pmb = (TH1D*)hp->Clone();
    TH1D* hratio_amb = (TH1D*)ha->Clone();
    
    hratio_mmb->SetMarkerStyle(20);
    hratio_mmb->SetMarkerColor(kBlue);
    hratio_pmb->SetMarkerStyle(20);
    hratio_pmb->SetMarkerColor(kRed);
    hratio_amb->SetMarkerStyle(20);
    hratio_amb->SetMarkerColor(kGreen);
    
    hratio_mmb->Divide(hmb);
    hratio_pmb->Divide(hmb);
    hratio_amb->Divide(hmb);
    
    //    for(int i=0;i<hratio_mmb->GetNbinsX();++i) {
    //        std::cout << hm->GetBinCenter(i+1) << " " << hm->GetBinContent(i+1) << " " << hmb->GetBinContent(i+1) << " ratio " << (hm->GetBinContent(i+1)/hmb->GetBinContent(i+1)) << " " << hratio_mmb->GetBinContent(i+1)
    //        << " " << hp->GetBinContent(i+1) << " " << hmb->GetBinContent(i+1) << " ratio " << (hp->GetBinContent(i+1)/hmb->GetBinContent(i+1)) << " " << hratio_pmb->GetBinContent(i+1)
    //        << " " << ha->GetBinContent(i+1) << " " << hmb->GetBinContent(i+1) << " ratio " << (ha->GetBinContent(i+1)/hmb->GetBinContent(i+1)) << " " << hratio_amb->GetBinContent(i+1)
    //        << std::endl;
    //    }
    
    SetRatioStyle(hratio_mmb, binning.xtitle);
    
    hratio_mmb->Draw("phist");
    hratio_pmb->Draw("phistsame");
    hratio_amb->Draw("phistsame");
    
    TLine *line = new TLine(binning.bin_min,1,binning.bin_max,1);
    line->Draw("same");
    
    c->SaveAs(output_name);
}

const char* oscparname;

Binning dcp_binning = {51,-3.2044245,3.2044245};
Binning dm2_binning = {81,2.25e-3,2.75e-3};
Binning dm2_io_binning = {81,-2.75e-3,-2.25e-3};
Binning dm2_both_binning = {972, -3.0e-3, 3.0e-3};
Binning th23_binning = {81,0.38,0.62};
Binning binning;
Binning binning_io;
Binning binning_both;

void Compare_MaCh3_Ptheta_Aria_MaCh3postBANFF_t2knova_dm2_21Aug2022(const char* expt, const char* osc_par, const char* outstring) {
    
    gStyle->SetOptStat(0);
    
    if(strstr(osc_par,"dCP") != NULL) {
        std::cout << " dcp " << std::endl;
        binning = dcp_binning;
        binning_io = dcp_binning;
        oscparname = "dcp";
        binning.xtitle = "#delta_{CP}";
    }
    else if(strstr(osc_par,"dm2") != NULL) {
        std::cout << " dm2 " << std::endl;
        binning = dm2_binning;
        binning_io = dm2_io_binning;
        binning_both = dm2_both_binning;
        oscparname = "dm2";
        binning.date = "30June2022";
        binning.xtitle = "#Delta m^{2}_{32}";
    }
    else if(strstr(osc_par,"th23") != NULL) {
        std::cout << " th23 " << std::endl;
        binning = th23_binning;
        binning_io = th23_binning;
        oscparname = "sin223";
        binning.xtitle = "sin^{2}#theta_{23}";
    }
    else {
        std::cout << "Please provide valid osc par.!" << std::endl;
        exit(1);
    }
    //TFile *fm=new TFile("MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
    //TFile *fm=new TFile("MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
    
    TFile* fm;
    if(strstr(expt,"T2KNOvA") != NULL) {
        std::cout << " MaCh3 T2KNOvA input file: " << std::endl;
        fm=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
    }
    else {
        std::cout << " MaCh3 T2K only input file: " << std::endl;
        //fm=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
        //fm=new TFile("MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
        fm=new TFile("~/Downloads/MaCh3_t2konly_Asimov0_03Sep2021_job3_ch0-999_red_reweighted.root"); //Asimov0
    }
    
    
    TTree *trm = (TTree*)fm->Get("osc_posteriors");
    Double_t dcp;
    Double_t dm2;
    Double_t th23;
    Double_t RCreweight;
    trm->SetBranchAddress("dcp",&dcp);
    trm->SetBranchAddress("dm23",&dm2);
    trm->SetBranchAddress("theta23",&th23);
    trm->SetBranchAddress("RCreweight",&RCreweight);
    
    Int_t nentries = trm->GetEntries();
    
    TH1D* hm=new TH1D("hm",Form("Both Hierarchies;%s;Posterior probability", binning.xtitle),           binning_both.nbins, binning_both.bin_min, binning_both.bin_max);
    TH1D* hm_both_no =new TH1D("hm_both_no",Form("Both Hierarchies;%s;Posterior probability", binning.xtitle), binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hm_both_io =new TH1D("hm_both_io",Form("Both Hierarchies;%s;Posterior probability", binning.xtitle), binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    TH1D* hm_no =new TH1D("hm_no",Form("Normal Hierarchy;%s;Posterior probability", binning.xtitle), binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hm_io =new TH1D("hm_io",Form("Inverted Hierarchy;%s;Posterior probability", binning.xtitle), binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    
    for(int i=0;i<((nsteps)?nentries:100000);++i) { //nentries
        trm->GetEntry(i);
        if(strstr(osc_par,"dCP") != NULL) {
            hm->Fill(dcp,RCreweight);
            if(dm2>0) hm_no->Fill(dcp,RCreweight); //Normal hierarchy
            else if (dm2<0) hm_io->Fill(dcp,RCreweight);
        }
        else if(strstr(osc_par,"dm2") != NULL) {
            hm->Fill(dm2,RCreweight);
            if(dm2>0) hm_no->Fill(dm2,RCreweight); //Normal hierarchy
            else if (dm2<0) hm_io->Fill(dm2,RCreweight);
        }
        else if(strstr(osc_par,"th23") != NULL) {
            hm->Fill(th23,RCreweight);
            if(dm2>0) hm_no->Fill(th23,RCreweight); //Normal hierarchy
            else if (dm2<0) hm_io->Fill(th23,RCreweight);
        }
    }
    //hm->Scale(1/hm->Integral()); //normalize
    hm_no->Scale(1/hm_no->Integral()); //normalize
    hm_io->Scale(1/hm_io->Integral()); //normalize
    
    TSpline3* sm = new TSpline3(hm);
    for(int i=0;i<hm_no->GetNbinsX();++i) {
        hm_both_no->SetBinContent(i+1, sm->Eval(hm_both_no->GetBinCenter(i+1)));
        hm_both_io->SetBinContent(i+1, sm->Eval(hm_both_io->GetBinCenter(i+1)));
        //std::cout << i+1 << " hm_no->GetBinCenter(i+1) " << hm_no->GetBinCenter(i+1) << " Eval " << sm->Eval(hm_no->GetBinCenter(i+1)) << " hm_io->GetBinCenter(i+1) " << hm_io->GetBinCenter(i+1) << " Eval " << sm->Eval(hm_io->GetBinCenter(i+1)) << std::endl;
    }
    
    hm_both_no->Scale(1/hm->Integral()); //normalize
    hm_both_io->Scale(1/hm->Integral()); //normalize
    
    for(int i=0;i<hm_no->GetNbinsX();++i) {
        //std::cout << i+1 << " hm_no->GetBinCenter(i+1) " << hm_no->GetBinCenter(i+1) << " Eval " << sm->Eval(hm_no->GetBinCenter(i+1)) << " hm_io->GetBinCenter(i+1) " << hm_io->GetBinCenter(i+1) << " Eval " << sm->Eval(hm_io->GetBinCenter(i+1)) << std::endl;
    }
    
    //    TH1D* hm_iozoom = (TH1D*)hm->Clone();
    //    hm_iozoom->SetLineColor(kBlue);
    //    hm_iozoom->SetLineWidth(2);
    //    hm_iozoom->GetXaxis()->SetRangeUser(-0.00275,-0.00225);
    //    hm->GetXaxis()->SetRangeUser(0.00225, 0.00275);
    
    hm->SetLineWidth(2);
    hm_no->SetLineWidth(2);
    hm_io->SetLineWidth(2);
    
    hm->SetLineColor(kBlue);
    hm_no->SetLineColor(kBlue);
    hm_io->SetLineColor(kBlue);
    
    //************************//
    //== Aria ==//
    
    TFile* fa;
    if(strstr(expt,"T2KNOvA") != NULL) {
        std::cout << " Aria T2KNOvA input file: " << std::endl;
        fa=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root");
    }
    else {
        std::cout << " Aria T2K only input file: " << std::endl;
        //fa=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root");
        //fa=new TFile("~/Downloads/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root"); //Asimov1
        fa=new TFile("~/Downloads/mcmcsamples_t2kcafana_joint_both_systs_noreactor_nddata_thinned1_postburnin20000.root"); //Asimov0
    }
    
    TTree *tra = (TTree*)fa->Get("samples");
    Double_t adcp;
    Double_t adm2;
    Double_t ath23;
    Double_t ath13;
    Double_t amh;
    tra->SetBranchAddress("delta(pi)",&adcp);
    tra->SetBranchAddress("dmsq32",&adm2);
    tra->SetBranchAddress("th23",&ath23);
    tra->SetBranchAddress("th13",&ath13);
    tra->SetBranchAddress("MH",&amh);
    
    Int_t anentries = tra->GetEntries();
    
    TH1D* ha=new TH1D("ha","ha",           binning_both.nbins, binning_both.bin_min, binning_both.bin_max);
    TH1D* ha_both_no =new TH1D("ha_both_no",Form("Both Hierarchies;%s;Posterior probability", binning.xtitle), binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* ha_both_io =new TH1D("ha_both_io",Form("Both Hierarchies;%s;Posterior probability", binning.xtitle), binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    TH1D* ha_no =new TH1D("ha_no",Form("Normal Hierarchy;%s;Posterior probability", binning.xtitle), binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* ha_io =new TH1D("ha_io",Form("Inverted Hierarchy;%s;Posterior probability", binning.xtitle), binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    
    for(int i=0;i<((nsteps)?anentries:100000);++i) { //anentries
        tra->GetEntry(i);
        
        Double_t ssth13 = TMath::Power(TMath::Sin(ath13),2);
        Double_t ssth23 = TMath::Power(TMath::Sin(ath23),2);
        double ssth23_weight = 2.*TMath::Sqrt((1-ssth23)*ssth23); // To go from flat prior in th23 -> flat prior in sin^2(th23)
        double ssth13_weight = 2.*TMath::Sqrt((1-ssth13)*ssth13); // To go from flat prior in th13 -> flat prior in sin^2(th13)
        Double_t rc_weight = TMath::Gaus(ssth13, 2.18e-2, 0.07e-2);
        
        if(strstr(osc_par,"dCP") != NULL) {
            //NOvA dCP is in NOvA units from 0 to 2, so need to do the conversion.
            if(adcp > 1.0) adcp -= 2.0;
            adcp *= TMath::Pi();
            
            ha->Fill(adcp, rc_weight*ssth13_weight*ssth23_weight);
            if(amh>0) ha_no->Fill(adcp, rc_weight*ssth13_weight*ssth23_weight); //Normal hierarchy
            else if (amh<0) ha_io->Fill(adcp, rc_weight*ssth13_weight*ssth23_weight);
        }
        else if(strstr(osc_par,"dm2") != NULL) {
            ha->Fill(adm2, rc_weight*ssth13_weight*ssth23_weight);
            if(amh>0) ha_no->Fill(adm2, rc_weight*ssth13_weight*ssth23_weight); //Normal hierarchy
            else if (amh<0) ha_io->Fill(adm2, rc_weight*ssth13_weight*ssth23_weight);
        }
        else if(strstr(osc_par,"th23") != NULL) {
            ha->Fill(ssth23, rc_weight*ssth13_weight*ssth23_weight);
            if(amh>0) ha_no->Fill(ssth23, rc_weight*ssth13_weight*ssth23_weight); //Normal hierarchy
            else if (amh<0) ha_io->Fill(ssth23, rc_weight*ssth13_weight*ssth23_weight);
        }
    }
    //ha->Scale(1/ha->Integral()); //normalize
    ha_no->Scale(1/ha_no->Integral()); //normalize
    ha_io->Scale(1/ha_io->Integral()); //normalize
    
    TSpline3* sa = new TSpline3(ha);
    for(int i=0;i<ha_no->GetNbinsX();++i) {
        ha_both_no->SetBinContent(i+1, sa->Eval(ha_both_no->GetBinCenter(i+1)));
        ha_both_io->SetBinContent(i+1, sa->Eval(ha_both_io->GetBinCenter(i+1)));
    }
    
    ha_both_no->Scale(1/ha->Integral()); //normalize
    ha_both_io->Scale(1/ha->Integral()); //normalize
    
    //    TH1D* ha_iozoom = (TH1D*)ha->Clone();
    //    ha_iozoom->SetLineColor(kGreen);
    //    ha_iozoom->SetLineWidth(2);
    //    ha_iozoom->GetXaxis()->SetRangeUser(-0.00275,-0.00225);
    //    ha->GetXaxis()->SetRangeUser(0.00225, 0.00275);
    
    ha->SetLineWidth(2);
    ha_no->SetLineWidth(2);
    ha_io->SetLineWidth(2);
    
    ha->SetLineColor(kBlue);
    ha_no->SetLineColor(kGreen);
    ha_io->SetLineColor(kGreen);
    
    //************************//
    //== MaCh3 post-BANFF ==//
    
    TFile* fmbanff;
    if(strstr(expt,"T2KNOvA") != NULL) {
        std::cout << " MaCh3 T2KNOvA input file: " << std::endl;
        fmbanff=new TFile("MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
    }
    else {
        std::cout << " MaCh3 T2K only input file: " << std::endl;
        //fmbanff=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/t2konly_banff_fit_job0-8_reweighted.root");
        //fmbanff=new TFile("~/Downloads/t2konly_banff_fit_job0-8_reweighted.root"); // RC 2018
        //fmbanff=new TFile("~/Downloads/t2konly_banff_job0-8_reweighted.root"); // RC 2019 Asimov1
        fmbanff=new TFile("t2konly_postbanff_Asimov0_BURNIN_CUT.root"); //Asimov0
    }
    
    TTree *trmbanff = (TTree*)fmbanff->Get("posteriors");
    Double_t mbanff_dcp;
    Double_t mbanff_dm2;
    Double_t mbanff_th23;
    Double_t mbanff_RCreweight;
    trmbanff->SetBranchAddress("delta_cp",&mbanff_dcp);
    trmbanff->SetBranchAddress("delm2_23",&mbanff_dm2);
    trmbanff->SetBranchAddress("sin2th_23",&mbanff_th23);
    trmbanff->SetBranchAddress("RCreweight",&mbanff_RCreweight);
    
    Int_t nentries_mbanff = trmbanff->GetEntries();
    
    TH1D* hmbanff=new TH1D("hmbanff","hmbanff",           binning_both.nbins, binning_both.bin_min, binning_both.bin_max);
    TH1D* hmbanff_both_no =new TH1D("hmbanff_both_no",Form("Both Hierarchies;%s;Posterior probability", binning.xtitle), binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hmbanff_both_io =new TH1D("hmbanff_both_io",Form("Both Hierarchies;%s;Posterior probability", binning.xtitle), binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    TH1D* hmbanff_no =new TH1D("hmbanff_no",Form("Normal Hierarchy;%s;Posterior probability", binning.xtitle), binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hmbanff_io =new TH1D("hmbanff_io",Form("Inverted Hierarchy;%s;Posterior probability", binning.xtitle), binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    
    for(int i=0;i<((nsteps)?nentries_mbanff:100000);++i) { //nentries_mbanff
        trmbanff->GetEntry(i);
        if(strstr(osc_par,"dCP") != NULL) {
            hmbanff->Fill(mbanff_dcp,mbanff_RCreweight);
            if(mbanff_dm2>0) hmbanff_no->Fill(mbanff_dcp,mbanff_RCreweight); //Normal hierarchy
            else if (mbanff_dm2<0) hmbanff_io->Fill(mbanff_dcp,mbanff_RCreweight);
        }
        else if(strstr(osc_par,"dm2") != NULL) {
            hmbanff->Fill(mbanff_dm2,mbanff_RCreweight);
            if(mbanff_dm2>0) hmbanff_no->Fill(mbanff_dm2,mbanff_RCreweight); //Normal hierarchy
            else if (mbanff_dm2<0) hmbanff_io->Fill(mbanff_dm2,mbanff_RCreweight);
        }
        else if(strstr(osc_par,"th23") != NULL) {
            hmbanff->Fill(mbanff_th23,mbanff_RCreweight);
            if(mbanff_dm2>0) hmbanff_no->Fill(mbanff_th23,mbanff_RCreweight); //Normal hierarchy
            else if (mbanff_dm2<0) hmbanff_io->Fill(mbanff_th23,mbanff_RCreweight);
        }
    }
    //hmbanff->Scale(1/hmbanff->Integral()); //normalize
    hmbanff_no->Scale(1/hmbanff_no->Integral()); //normalize
    hmbanff_io->Scale(1/hmbanff_io->Integral()); //normalize
    
    TSpline3* smbanff = new TSpline3(hmbanff);
    for(int i=0;i<hmbanff_no->GetNbinsX();++i) {
        hmbanff_both_no->SetBinContent(i+1, smbanff->Eval(hmbanff_both_no->GetBinCenter(i+1)));
        hmbanff_both_io->SetBinContent(i+1, smbanff->Eval(hmbanff_both_io->GetBinCenter(i+1)));
    }
    
    hmbanff_both_no->Scale(1/hmbanff->Integral()); //normalize
    hmbanff_both_io->Scale(1/hmbanff->Integral()); //normalize
    
    //    TH1D* hmbanff_iozoom = (TH1D*)ha->Clone();
    //    hmbanff_iozoom->SetLineColor(kBlack);
    //    hmbanff_iozoom->SetLineWidth(2);
    //    hmbanff_iozoom->GetXaxis()->SetRangeUser(-0.00275,-0.00225);
    //    hmbanff->GetXaxis()->SetRangeUser(0.00225, 0.00275);
    
    hmbanff->SetLineWidth(2);
    hmbanff_no->SetLineWidth(2);
    hmbanff_io->SetLineWidth(2);
    
    hmbanff->SetLineColor(kBlack);
    hmbanff_no->SetLineColor(kBlack);
    hmbanff_io->SetLineColor(kBlack);
    
    //TCanvas *c=new TCanvas("c","c",700,400);
    //c->SetLeftMargin(0.14);
    //c->Divide(2,1);
    //hm->Draw("HIST");
    
    //TFile *fp = TFile::Open("credible_T2KNOvA_both.root");
    //TFile *fp_no = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_%s/credible_%s.root",osc_par,binning.date,expt));
    TFile *fp_no = TFile::Open(Form("credible_%s_%s_asimov0_both.root",osc_par,expt));
    TH1D* hp_no=(TH1D*)fp_no->Get("LL");
    hp_no->SetLineColor(kRed);
    
    if(strstr(osc_par,"dCP") != NULL) hp_no->Rebin(51);
    else hp_no->Rebin(81);
    
    TH1D* htemp_both_no=new TH1D("htemp_both_no", "hemp_both_no",binning.nbins,binning.bin_min,binning.bin_max);
    TH1D* htemp_no=new TH1D("htemp_no", "hemp_no",binning.nbins,binning.bin_min,binning.bin_max);
    
    for(int i=0;i<binning.nbins;++i) {
        htemp_both_no->SetBinContent(i+1, hp_no->GetBinContent(i+1));
        htemp_no->SetBinContent(i+1, hp_no->GetBinContent(i+1));
    }
    
    //TFile *fp_io = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_%s/credible_%s_IH.root",osc_par, binning.date,expt));
    //TFile *fp_io = TFile::Open(Form("credible_%s_%s_IH.root",osc_par,expt));
    TFile *fp_io = TFile::Open(Form("credible_%s_%s_asimov0_both.root",osc_par,expt));
    TH1D* hp_io=(TH1D*)fp_io->Get("LL_inv");
    
    if(strstr(osc_par,"dCP") != NULL) hp_io->Rebin(51);
    else hp_io->Rebin(81);
    
    //hp_no->Draw("hist");
    //hp_io->Draw("hist same");
    
    TSpline3 *sp_io=new TSpline3(hp_io);
    
    TH1D* htemp_both_io=new TH1D("htemp_both_io", "hemp_both_io",binning_io.nbins,binning_io.bin_min,binning_io.bin_max);
    TH1D* htemp_io=new TH1D("htemp_io", "hemp_io",binning_io.nbins,binning_io.bin_min,binning_io.bin_max);
    
    for(int i=0;i<binning.nbins;++i) {
        //htemp_io->SetBinContent(i+1, hp_io->GetBinContent(i+1));
        double dm32 = hm_io->GetBinCenter(i+1);
        double dm21 = 7.53e-5;
        double dm31 = dm32 + dm21;
        double abs_dm31 = fabs(dm31);
        //std::cout << " abs_dm31 " << abs_dm31 << " sp_io->Eval(abs_dm31) " << sp_io->Eval(abs_dm31) << std::endl;
        htemp_both_io->SetBinContent(i+1, sp_io->Eval(abs_dm31));
        htemp_io->SetBinContent(i+1, sp_io->Eval(abs_dm31));
    }
    
    Double_t totalbc=0;
    for(int i=0;i<htemp_no->GetNbinsX();++i) {
        totalbc += htemp_no->GetBinContent(i+1) + htemp_io->GetBinContent(i+1);
    }
    
    //htemp_both_no->Draw("hist");
    std::cout << " totlabc " << totalbc << " integral " << (htemp_no->Integral() + htemp_io->Integral()) << std::endl;
    
//    htemp_both_no->Scale(1/totalbc);
//    htemp_both_io->Scale(1/totalbc);
//
//    htemp_no->Scale(1/htemp_no->Integral());
//    htemp_io->Scale(1/htemp_io->Integral());
    
    //std::cout << " totlabc " << totalbc << " integral " << htemp_no->Integral() << " " << htemp_io->Integral() << " " << " " << htemp_both_no->Integral() << " " << htemp_both_io->Integral() << std::endl;

    //htemp_both_no->Draw("hist");
    //htemp_no->Draw("hist same");
    
    hm_both_no->GetXaxis()->SetMaxDigits(2);
    hm_both_io->GetXaxis()->SetMaxDigits(2);
    hm_no->GetXaxis()->SetMaxDigits(2);
    hm_io->GetXaxis()->SetMaxDigits(2);
    
    TLegend* leg_no;
    //leg=new TLegend(0.6,0.67,0.75,0.85);
    leg_no=new TLegend(0.18,0.65,0.28,0.85);
    leg_no->SetBorderSize(0);
    leg_no->SetTextSize(0.04);

    TLegend* leg_io;
    //leg=new TLegend(0.6,0.67,0.75,0.85);
    leg_io=new TLegend(0.6181948,0.6534149,0.7686246,0.852758);
    leg_io->SetBorderSize(0);
    leg_io->SetTextSize(0.04);
    
    TLegend* leg_dchi2;
    //leg=new TLegend(0.6,0.67,0.75,0.85);
    if(strstr(outstring,"asimov1") != NULL)
        leg_dchi2=new TLegend(0.36,0.65,0.46,0.85);
    else
        leg_dchi2=new TLegend(0.26,0.65,0.36,0.85);
    leg_dchi2->SetBorderSize(0);
    leg_dchi2->SetTextSize(0.04);
    
    leg_no->AddEntry(hm,"MaCh3","l");
    leg_no->AddEntry(hp_no,"Ptheta","l");
    leg_no->AddEntry(ha_no,"Aria","l");
    leg_no->AddEntry(hmbanff,"MaCh3 post-BANFF","l");
    
    leg_io->AddEntry(hm,"MaCh3","l");
    leg_io->AddEntry(hp_no,"Ptheta","l");
    leg_io->AddEntry(ha_no,"Aria","l");
    leg_io->AddEntry(hmbanff,"MaCh3 post-BANFF","l");
    
    leg_dchi2->AddEntry(hm,"MaCh3","l");
    leg_dchi2->AddEntry(hp_no,"Ptheta","l");
    leg_dchi2->AddEntry(ha_no,"Aria","l");
    leg_dchi2->AddEntry(hmbanff,"MaCh3 post-BANFF","l");
    
    std::cout << " #bins hm " << hm->GetNbinsX() << " hmbanff " << hmbanff->GetNbinsX() << " ha " << ha->GetNbinsX() << " hp " << htemp_no->GetNbinsX() << std::endl; //(strstr(osc_par,"asimov1") != NULL)?leg_no:leg_io
    
    GetRatioHist(htemp_both_no, hm_both_no, ha_both_no, hmbanff_both_no, (strstr(outstring,"asimov1") != NULL)?leg_no:leg_io, binning, false,  Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_Both_NO_%s_%s.pdf",osc_par,expt,outstring)); //Both ierarchies meaning NO/IO normalized with both.

    GetRatioHist(htemp_both_io, hm_both_io, ha_both_io, hmbanff_both_io, (strstr(outstring,"asimov1") != NULL)?leg_io:leg_no, binning, false,  Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_Both_IO_%s_%s.pdf",osc_par,expt,outstring));

    GetRatioHist(htemp_no, hm_no, ha_no, hmbanff_no, (strstr(outstring,"asimov1") != NULL)?leg_no:leg_io, binning, false,  Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_NO_%s_%s.pdf",osc_par,expt,outstring));

    GetRatioHist(htemp_io, hm_io, ha_io, hmbanff_io, (strstr(outstring,"asimov1") != NULL)?leg_io:leg_no, binning, false,  Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_IO_%s_%s.pdf",osc_par,expt,outstring));
    
//    GetRatioHist(htemp_no, hm_both_no, ha_both_no, hmbanff_both_no, (strstr(outstring,"asimov1") != NULL)?leg_no:leg_io, binning, false,  Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_Both_NO_%s_%s.pdf",osc_par,expt,outstring)); //Both ierarchies meaning NO/IO normalized with both.
//
//    GetRatioHist(htemp_io, hm_both_io, ha_both_io, hmbanff_both_io, (strstr(outstring,"asimov1") != NULL)?leg_io:leg_no, binning, false,  Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_Both_IO_%s_%s.pdf",osc_par,expt,outstring));
    
    //(strstr(osc_par,"th23") != NULL)?leg:leg1
    
 /*
    //Plot Dchi2
    
    Double_t bincontmax_no=hm_no->GetBinContent(hm_no->GetMaximumBin());
    Double_t bincontmax_io=hm_io->GetBinContent(hm_io->GetMaximumBin());
    
    Double_t abincontmax_no=ha_no->GetBinContent(ha_no->GetMaximumBin());
    Double_t abincontmax_io=ha_io->GetBinContent(ha_io->GetMaximumBin());
    
    Double_t mbanffbincontmax_no=hmbanff_no->GetBinContent(hmbanff_no->GetMaximumBin());
    Double_t mbanffbincontmax_io=hmbanff_io->GetBinContent(hmbanff_io->GetMaximumBin());
    
    std::cout << bincontmax_no << " " << bincontmax_io << abincontmax_no << " " << abincontmax_io << " " << mbanffbincontmax_no << " " << mbanffbincontmax_io << hm->GetNbinsX() << std::endl;
    
    TH1D* hmdchi2_no=(TH1D*)hm_no->Clone();
    TH1D* hmdchi2_io=(TH1D*)hm_io->Clone();
    
    TH1D* hadchi2_no=(TH1D*)ha_no->Clone();
    TH1D* hadchi2_io=(TH1D*)ha_io->Clone();
    
    TH1D* hmbanffdchi2_no=(TH1D*)hmbanff_no->Clone();
    TH1D* hmbanffdchi2_io=(TH1D*)hmbanff_io->Clone();
    
    //    hmbanffdchi2_no->SetLineColor(kViolet);
    //    hmbanffdchi2_io->SetLineColor(kViolet);
    //
    //    hmbanffdchi2_no->SetLineWidth(2);
    //    hmbanffdchi2_io->SetLineWidth(2);
    
    
    for(int i=0;i<binning.nbins;++i) {
        std::cout << i+1 << " " << hm_no->GetBinContent(i+1) << " " << hm_io->GetBinContent(i+1) << " " << -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no) << " " << -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io) << std::endl;
        if(hm_no->GetBinContent(i+1) > 0) hmdchi2_no->SetBinContent(i+1, -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no));
        if(hm_io->GetBinContent(i+1) > 0) hmdchi2_io->SetBinContent(i+1, -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io));
        
        if(ha_no->GetBinContent(i+1) > 0) hadchi2_no->SetBinContent(i+1, -2*TMath::Log(ha_no->GetBinContent(i+1)/abincontmax_no));
        if(ha_io->GetBinContent(i+1) > 0) hadchi2_io->SetBinContent(i+1, -2*TMath::Log(ha_io->GetBinContent(i+1)/abincontmax_io));
        
        if(hmbanff_no->GetBinContent(i+1) > 0) hmbanffdchi2_no->SetBinContent(i+1, -2*TMath::Log(hmbanff_no->GetBinContent(i+1)/mbanffbincontmax_no));
        if(hmbanff_io->GetBinContent(i+1) > 0) hmbanffdchi2_io->SetBinContent(i+1, -2*TMath::Log(hmbanff_io->GetBinContent(i+1)/mbanffbincontmax_io));
    }
    
    // NO
    //TCanvas *c1=new TCanvas("c1","c1",700,700);
    
    //TFile *fpdchi2 = new TFile(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_%s/contour_%s_100k_%s.root",osc_par,binning.date,osc_par,expt));
    TFile *fpdchi2 = new TFile(Form("contour_%s_100k_%s.root",osc_par,expt));
    TH1D *hpdchi2_no=(TH1D*)fpdchi2->Get(Form("SA_%s",oscparname));
    TH1D *hpdchi2_io=(TH1D*)fpdchi2->Get(Form("SA_%s_IH",oscparname));
    TH1D *hpdchi2_temp=(TH1D*)hpdchi2_io->Clone();
    
    TH1D *htemp1=new TH1D("htemp1", "htemp1",binning_io.nbins,binning_io.bin_min,binning_io.bin_max);
    
    std::cout << hmdchi2_no->GetNbinsX() << " " << hmdchi2_no->GetXaxis()->GetXmin() << " " << hmdchi2_no->GetXaxis()->GetXmax()
    << " " << hadchi2_no->GetNbinsX() << " " << hadchi2_no->GetXaxis()->GetXmin() << " " << hadchi2_no->GetXaxis()->GetXmax()
    << " " << hpdchi2_no->GetNbinsX() << " " << hpdchi2_no->GetBinLowEdge(1) << " " << hpdchi2_no->GetBinCenter(81)
    << " " << hmbanff_no->GetNbinsX() << " " << hmbanff_no->GetBinLowEdge(1) << " " << hmbanff_no->GetBinCenter(81)
    << std::endl;
    
    //    hpdchi2_no->SetLineColor(kRed);
    //    hpdchi2_io->SetLineColor(kRed);
    //    hpdchi2_no->SetLineWidth(2);
    //    hpdchi2_io->SetLineWidth(2);
    
    hmdchi2_no->GetYaxis()->SetTitle("#Delta#chi^{2}");
    hmdchi2_io->GetYaxis()->SetTitle("#Delta#chi^{2}");
    
    //(strstr(osc_par,"th23") != NULL)?hmdchi2_no->GetYaxis()->SetRangeUser(0,30):hmdchi2_no->GetYaxis()->SetRangeUser(0,15);
    hmdchi2_no->GetXaxis()->SetMaxDigits(2);
    hmdchi2_io->GetXaxis()->SetMaxDigits(2);
    
    GetRatioHist(hpdchi2_no, hmdchi2_no, hadchi2_no, hmbanffdchi2_no, leg_dchi2, binning, true, Form("compare_dchi2_mach3_ptheta_aria_mach3postbanff_%s_NO_%s_%s.pdf",osc_par,expt,outstring));
    
    TSpline3 *sp=new TSpline3(hpdchi2_io);
    //for(int i=0;i<binning.nbins;++i) hpdchi2_io->SetBinContent(i+1, hpdchi2_temp->GetBinContent(i+1)-locmin);
    
    for(int i=0;i<binning.nbins;++i) {
        double dm32 = htemp1->GetBinCenter(i+1);
        double dm21 = 7.53e-5;
        double dm31 = dm32 + dm21;
        double abs_dm31 = fabs(dm31);
        //std::cout << " abs_dm31 " << abs_dm31 << std::endl;
        htemp1->SetBinContent(i+1, sp->Eval(abs_dm31));
    }
    
    GetRatioHist(htemp1, hmdchi2_io, hadchi2_io, hmbanffdchi2_io, leg_dchi2, binning, true, Form("compare_dchi2_mach3_ptheta_aria_mach3postbanff_%s_IO_%s_%s.pdf",osc_par,expt,outstring));
 */
}
