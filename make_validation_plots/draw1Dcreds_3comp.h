//This is to draw credible intervals
#include "TPaveText.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TMath.h"
#include <iostream>

void drawCredLine(TCanvas* canv, double val, double ymin, double ymax, int color, int style, int width) {
    
    TLine* line = new TLine(val,ymin,val,ymax);
    line->SetLineWidth(width);
    line->SetLineColorAlpha(color,1.);
    line->SetLineStyle(style);
    canv->cd();
    line->Draw();
    
}

void drawCredLine(TPad* canv, double val, double ymin, double ymax, int color, int style, int width) {
    
    TLine* line = new TLine(val,ymin,val,ymax);
    line->SetLineWidth(width);
    line->SetLineColorAlpha(color,1.);
    line->SetLineStyle(style);
    canv->cd();
    line->Draw();
    
}

void drawCredInterval(TCanvas* canv, double xmin, double xmax, double ymin, double ymax, int color, double alpha, int style) {
    
    TBox* cred = new TBox(xmin,ymin,xmax,ymax);
    cred->SetLineWidth(0.);
    cred->SetLineColorAlpha(0.,0.);
    cred->SetFillColorAlpha(color,alpha);
    drawCredLine(canv,xmin,ymin,ymax,color,style,3);
    drawCredLine(canv,xmax,ymin,ymax,color,style,3);
    canv->cd();
    cred->Draw();
}

void drawCredInterval(TPad* canv, double xmin, double xmax, double ymin, double ymax, int color, double alpha, int style) {
    
    TBox* cred = new TBox(xmin,ymin,xmax,ymax);
    cred->SetLineWidth(0.);
    cred->SetLineColorAlpha(0.,0.);
    cred->SetFillColorAlpha(color,alpha);
    drawCredLine(canv,xmin,ymin,ymax,color,style,3);
    drawCredLine(canv,xmax,ymin,ymax,color,style,3);
    canv->cd();
    cred->Draw();
}

TH1D* get1Dcred(TH1D* hist, double level, double &bf, double &up, double &low, int &nbins) {
    
    if(level>1.) {
        std::cout << "credible interval must be < 1." << std::endl;
        throw;
    }
    
    TH1D* hist_copy = (TH1D*)hist->Clone("hist_copy");
    TH1D* hist_ret = (TH1D*)hist->Clone("hist_ret");
    double integral, tsum=0.;
    integral = hist_copy->Integral();
    double xval;
    double xwidth;
    
    // get best fit
    bf = hist_copy->GetXaxis()->GetBinCenter(hist_copy->GetMaximumBin());
    up = -9999999.;
    low = 9999999.;
    nbins = 0;
    
    while((tsum/integral)<level) {
        double tmax = hist_copy->GetMaximum();
        int bin = hist_copy->GetMaximumBin();
        xval = hist_copy->GetXaxis()->GetBinCenter(bin);
        xwidth = hist_copy->GetXaxis()->GetBinWidth(bin);
        if((tsum/integral)<level) {
            hist_copy->SetBinContent(bin,-1.0);
            hist_ret->SetBinContent(bin,0.);
            if(xval<low && xval<bf) low = xval - xwidth/2.;
            if(xval>up && xval>bf) up = xval + xwidth/2.;
            nbins++;
        }
        tsum+=tmax;
    }
    
    return hist_ret;
    
}

std::vector<std::array<double,2> > get1Dcred_disjoint(TH1D* hist, double level) {
    
    if(level>1.) {
        std::cout << "credible interval must be < 1." << std::endl;
        throw;
    }
    
    // a vector of lower and upper bounds for each disjoint credible intervals
    // (can be a single interval)
    std::vector<std::array<double,2> > cred;
    
    TH1D* hist_copy = (TH1D*)hist->Clone("hist_copy");
    double integral, tsum=0.;
    integral = hist_copy->Integral();
    std::vector<int> binsIn;
    
    while((tsum/integral)<level) {
        double tmax = hist_copy->GetBinContent(hist_copy->GetMaximumBin()); //->GetMaximum();
        int bin = hist_copy->GetMaximumBin();
        //std::cout << " integral " << integral << " tmax " << tmax << " bin " << bin << " tsum " << tsum << " (tsum/integral) " << (tsum/integral) << std::endl;
        if((tsum/integral)<level) {
            hist_copy->SetBinContent(bin,-1.0);
            binsIn.push_back(bin);
        }
        tsum+=tmax;
    }
    
    std::sort(binsIn.begin(),binsIn.end());
    
    std::array<double,2> firstlb = {hist->GetXaxis()->GetBinCenter(binsIn[0])-hist->GetXaxis()->GetBinWidth(binsIn[0])/2.,-999.};
    cred.push_back(firstlb);
    int N_disjoint=0;
    
    for(unsigned int i=1; i<binsIn.size();i++) {
        if(binsIn[i]-binsIn[i-1]==1) continue;
        else {
            cred[N_disjoint][1]=hist->GetXaxis()->GetBinCenter(binsIn[i-1])+hist->GetXaxis()->GetBinWidth(binsIn[i-1])/2.;
            std::array<double,2> temp_arr = {hist->GetXaxis()->GetBinCenter(binsIn[i])-hist->GetXaxis()->GetBinWidth(binsIn[i])/2.,-999.};
            cred.push_back(temp_arr);
            N_disjoint++;
        }
    }
    
    cred[N_disjoint][1]=hist->GetXaxis()->GetBinCenter(binsIn.back())+hist->GetXaxis()->GetBinWidth(binsIn.back())/2.;
    
    for(int i=0;i<cred.size();++i) {
        std::cout << " sig " << i << " " << cred[i][0] << " " << cred[i][1] << std::endl;
    }
    
    return cred;
    
}


// hieararchy = 1-->NH, 0-->both, -1-->IH
// app = true for appearance contours, false for disappearance
void makePlot(int hierarchy,string postNoCorr, string postCorr, string postAnticorr, string outTag = "nom", string asimov = "AsimovNA") {
    
    bool drawAsimovPoint = true;
    bool draw2sig = true;
    
    double asim_dcp, asim_th13, asim_th23, asim_dm2;
    
    if(asimov=="Asimov0") {
        asim_dcp = 2.576;
        asim_th13 = 0.0218;
        asim_th23 = 0.57;
        asim_dm2 = 2.41e-3;
    }
    else if(asimov=="Asimov1") {
        asim_dcp = -1.601;
        asim_th13 = 0.0218;
        asim_th23 = 0.528;
        asim_dm2 = 2.509e-3;
    }
    else if(asimov=="Asimov4") {
        asim_dcp = -1.58;
        asim_th13 = 0.0218;
        asim_th23 = 0.55;
        asim_dm2 = -2.45e-3;
    }
    
    // acor == antiCorr
    // corr == corr
    // ncor == noCorr
    
    TFile* f_acor = new TFile(postAnticorr.c_str());
    TTree* t_acor = (TTree*)f_acor->Get("osc_posteriors");
    TFile* f_corr = new TFile(postCorr.c_str());
    TTree* t_corr = (TTree*)f_corr->Get("osc_posteriors");
    TFile* f_ncor = new TFile(postNoCorr.c_str());
    TTree* t_ncor = (TTree*)f_ncor->Get("osc_posteriors");
    
    TH1D* h_dcp_acor = new TH1D("h_dcp_acor",";#delta_{CP};posterior probability",100,-1.*TMath::Pi(),TMath::Pi());
    TH1D* h_dcp_corr = new TH1D("h_dcp_corr",";#delta_{CP};posterior probability",100,-1.*TMath::Pi(),TMath::Pi());
    TH1D* h_dcp_ncor = new TH1D("h_dcp_ncor",";#delta_{CP};posterior probability",100,-1.*TMath::Pi(),TMath::Pi());
    
    TH1D* h_th23_acor = new TH1D("h_th23_acor",";sin^{2}#theta_{3};posterior probability",125,0.35,0.65);
    TH1D* h_th23_corr = new TH1D("h_th23_corr",";sin^{2}#theta_{23};posterior probability",125,0.35,0.65);
    TH1D* h_th23_ncor = new TH1D("h_th23_ncor",";sin^{2}#theta_{23};posterior probability",125,0.35,0.65);
    
    TH1D* h_dm32_acor = new TH1D("h_dm32_acor",";#Delta m^{2}_{32};posterior probability",1000,-0.003,0.003);
    TH1D* h_dm32_corr = new TH1D("h_dm32_corr",";#Delta m^{2}_{32};posterior probability",1000,-0.003,0.003);
    TH1D* h_dm32_ncor = new TH1D("h_dm32_ncor",";#Delta m^{2}_{32};posterior probability",1000,-0.003,0.003);
    
    TH1D* h_th13_acor = new TH1D("h_th13_acor",";sin^{2}#theta_{13};posterior probability",100,0.019,0.025);
    TH1D* h_th13_corr = new TH1D("h_th13_corr",";sin^{2}#theta_{13};posterior probability",100,0.019,0.025);
    TH1D* h_th13_ncor = new TH1D("h_th13_ncor",";sin^{2}#theta_{13};posterior probability",100,0.019,0.025);
    
    TPaveText* hlab = new TPaveText(0.83,0.84,0.88,0.9,"NDC");
    
    if(hierarchy==1) { // normal hiearchy --> dm32 > 0.
        ///*
        t_ncor->Draw("dcp>>h_dcp_ncor","(dm23>0.)*(step>200000)*RCreweight");
        t_ncor->Draw("theta13>>h_th13_ncor","(dm23>0.)*(step>200000)*RCreweight");
        t_ncor->Draw("dm23>>h_dm32_ncor","(dm23>0.)*(step>200000)*RCreweight");
        t_ncor->Draw("theta23>>h_th23_ncor","(dm23>0.)*(step>200000)*RCreweight");
        //*/
        /*
         t_ncor->Draw("dcp>>h_dcp_ncor","(dm23>0.)*(step>200000)");
         t_ncor->Draw("theta13>>h_th13_ncor","(dm23>0.)*(step>200000)");
         t_ncor->Draw("dm23>>h_dm32_ncor","(dm23>0.)*(step>200000)");
         t_ncor->Draw("theta23>>h_th23_ncor","(dm23>0.)*(step>200000)");
         */
        ///*
        t_corr->Draw("dcp>>h_dcp_corr","(dm23>0.)*(step>200000)*RCreweight");
        t_corr->Draw("theta13>>h_th13_corr","(dm23>0.)*(step>200000)*RCreweight");
        t_corr->Draw("dm23>>h_dm32_corr","(dm23>0.)*(step>200000)*RCreweight");
        t_corr->Draw("theta23>>h_th23_corr","(dm23>0.)*(step>200000)*RCreweight");
        //*/
        /*
         t_corr->Draw("dcp>>h_dcp_corr","(dm23>0.)*(step>200000)");
         t_corr->Draw("theta13>>h_th13_corr","(dm23>0.)*(step>200000)");
         t_corr->Draw("dm23>>h_dm32_corr","(dm23>0.)*(step>200000)");
         t_corr->Draw("theta23>>h_th23_corr","(dm23>0.)*(step>200000)");
         */
        ///*
        t_acor->Draw("dcp>>h_dcp_acor","(dm23>0.)*(step>200000)*RCreweight");
        t_acor->Draw("theta13>>h_th13_acor","(dm23>0.)*(step>200000)*RCreweight");
        t_acor->Draw("dm23>>h_dm32_acor","(dm23>0.)*(step>200000)*RCreweight");
        t_acor->Draw("theta23>>h_th23_acor","(dm23>0.)*(step>200000)*RCreweight");
        //*/
        /*
         t_acor->Draw("dcp>>h_dcp_acor","(dm23>0.)*(step>200000)");
         t_acor->Draw("theta13>>h_th13_acor","(dm23>0.)*(step>200000)");
         t_acor->Draw("dm23>>h_dm32_acor","(dm23>0.)*(step>200000)");
         t_acor->Draw("theta23>>h_th23_acor","(dm23>0.)*(step>200000)");
         */
        
        
        hlab->AddText("NH");
    }
    else if(hierarchy==0) {
        
        /*
         t_ncor->Draw("dcp>>h_dcp_ncor","(step>200000)");
         t_ncor->Draw("theta13>>h_th13_ncor","(step>200000)");
         t_ncor->Draw("dm23>>h_dm32_ncor","(step>200000)");
         t_ncor->Draw("theta23>>h_th23_ncor","(step>200000)");
         */
        ///*
        t_ncor->Draw("dcp>>h_dcp_ncor","(step>200000)*RCreweight");
        t_ncor->Draw("theta13>>h_th13_ncor","(step>200000)*RCreweight");
        t_ncor->Draw("dm23>>h_dm32_ncor","(step>200000)*RCreweight");
        t_ncor->Draw("theta23>>h_th23_ncor","(step>200000)*RCreweight");
        //*/
        ///*
        t_corr->Draw("dcp>>h_dcp_corr","(step>200000)*RCreweight");
        t_corr->Draw("theta13>>h_th13_corr","(step>200000)*RCreweight");
        t_corr->Draw("dm23>>h_dm32_corr","(step>200000)*RCreweight");
        t_corr->Draw("theta23>>h_th23_corr","(step>200000)*RCreweight");
        //*/
        /*
         t_corr->Draw("dcp>>h_dcp_corr","(step>200000)");
         t_corr->Draw("theta13>>h_th13_corr","(step>200000)");
         t_corr->Draw("dm23>>h_dm32_corr","(step>200000)");
         t_corr->Draw("theta23>>h_th23_corr","(step>200000)");
         */
        ///*
        t_acor->Draw("dcp>>h_dcp_acor","(step>200000)*RCreweight");
        t_acor->Draw("theta13>>h_th13_acor","(step>200000)*RCreweight");
        t_acor->Draw("dm23>>h_dm32_acor","(step>200000)*RCreweight");
        t_acor->Draw("theta23>>h_th23_acor","(step>200000)*RCreweight");
        //*/
        /*
         t_acor->Draw("dcp>>h_dcp_acor","(step>200000)");
         t_acor->Draw("theta13>>h_th13_acor","(step>200000)");
         t_acor->Draw("dm23>>h_dm32_acor","(step>200000)");
         t_acor->Draw("theta23>>h_th23_acor","(step>200000)");
         */
        
        hlab->AddText("NH+IH");
    }
    else if(hierarchy==-1) { // inverted hiearchy --> dm32 < 0.
        ///*
        t_ncor->Draw("dcp>>h_dcp_ncor","(dm23<0.)*(step>200000)*RCreweight");
        t_ncor->Draw("theta13>>h_th13_ncor","(dm23<0.)*(step>200000)*RCreweight");
        t_ncor->Draw("dm23>>h_dm32_ncor","(dm23<0.)*(step>200000)*RCreweight");
        t_ncor->Draw("theta23>>h_th23_ncor","(dm23<0.)*(step>200000)*RCreweight");
        //*/
        /*
         t_ncor->Draw("dcp>>h_dcp_ncor","(dm23<0.)*(step>200000)");
         t_ncor->Draw("theta13>>h_th13_ncor","(dm23<0.)*(step>200000)");
         t_ncor->Draw("dm23>>h_dm32_ncor","(dm23<0.)*(step>200000)");
         t_ncor->Draw("theta23>>h_th23_ncor","(dm23<0.)*(step>200000)");
         */
        ///*
        t_corr->Draw("dcp>>h_dcp_corr","(dm23<0.)*(step>200000)*RCreweight");
        t_corr->Draw("theta13>>h_th13_corr","(dm23<0.)*(step>200000)*RCreweight");
        t_corr->Draw("dm23>>h_dm32_corr","(dm23<0.)*(step>200000)*RCreweight");
        t_corr->Draw("theta23>>h_th23_corr","(dm23<0.)*(step>200000)*RCreweight");
        //*/
        /*
         t_corr->Draw("dcp>>h_dcp_corr","(dm23<0.)*(step>200000)");
         t_corr->Draw("theta13>>h_th13_corr","(dm23<0.)*(step>200000)");
         t_corr->Draw("dm23>>h_dm32_corr","(dm23<0.)*(step>200000)");
         t_corr->Draw("theta23>>h_th23_corr","(dm23<0.)*(step>200000)");
         */
        ///*
        t_acor->Draw("dcp>>h_dcp_acor","(dm23<0.)*(step>200000)*RCreweight");
        t_acor->Draw("theta13>>h_th13_acor","(dm23<0.)*(step>200000)*RCreweight");
        t_acor->Draw("dm23>>h_dm32_acor","(dm23<0.)*(step>200000)*RCreweight");
        t_acor->Draw("theta23>>h_th23_acor","(dm23<0.)*(step>200000)*RCreweight");
        //*/
        /*
         t_acor->Draw("dcp>>h_dcp_acor","(dm23<0.)*(step>200000)");
         t_acor->Draw("theta13>>h_th13_acor","(dm23<0.)*(step>200000)");
         t_acor->Draw("dm23>>h_dm32_acor","(dm23<0.)*(step>200000)");
         t_acor->Draw("theta23>>h_th23_acor","(dm23<0.)*(step>200000)");
         */
        hlab->AddText("IH");
    }
    
    
    h_dm32_ncor->GetXaxis()->SetMaxDigits(2);
    h_dm32_corr->GetXaxis()->SetMaxDigits(2);
    h_dm32_acor->GetXaxis()->SetMaxDigits(2);
    
    h_dcp_ncor->Scale(1./h_dcp_ncor->Integral(),"width");
    h_dcp_corr->Scale(1./h_dcp_corr->Integral(),"width");
    h_dcp_acor->Scale(1./h_dcp_acor->Integral(),"width");
    
    h_dm32_ncor->Scale(1./h_dm32_ncor->Integral(),"width");
    h_dm32_corr->Scale(1./h_dm32_corr->Integral(),"width");
    h_dm32_acor->Scale(1./h_dm32_acor->Integral(),"width");
    
    h_th23_ncor->Scale(1./h_th23_ncor->Integral(),"width");
    h_th23_corr->Scale(1./h_th23_corr->Integral(),"width");
    h_th23_acor->Scale(1./h_th23_acor->Integral(),"width");
    
    h_th13_ncor->Scale(1./h_th13_ncor->Integral(),"width");
    h_th13_corr->Scale(1./h_th13_corr->Integral(),"width");
    h_th13_acor->Scale(1./h_th13_acor->Integral(),"width");
    
    std::cout << "4.25" << std::endl;
    
    // get credible intervals:
    std::vector<std::array<double,2> > dcp_cred1sig_regions_corr = get1Dcred_disjoint(h_dcp_corr,0.6826);
    std::vector<std::array<double,2> > dcp_cred2sig_regions_corr = get1Dcred_disjoint(h_dcp_corr,0.9544);
    std::vector<std::array<double,2> > th23_cred1sig_regions_corr = get1Dcred_disjoint(h_th23_corr,0.6826);
    std::vector<std::array<double,2> > th23_cred2sig_regions_corr = get1Dcred_disjoint(h_th23_corr,0.9544);
    std::vector<std::array<double,2> > th13_cred1sig_regions_corr = get1Dcred_disjoint(h_th13_corr,0.6826);
    std::vector<std::array<double,2> > th13_cred2sig_regions_corr = get1Dcred_disjoint(h_th13_corr,0.9544);
    std::vector<std::array<double,2> > dm32_cred1sig_regions_corr = get1Dcred_disjoint(h_dm32_corr,0.6826);
    std::vector<std::array<double,2> > dm32_cred2sig_regions_corr = get1Dcred_disjoint(h_dm32_corr,0.9544);
    
    std::cout << "4.30" << std::endl;
    std::vector<std::array<double,2> > dcp_cred1sig_regions_ncor = get1Dcred_disjoint(h_dcp_ncor,0.6826);
    std::vector<std::array<double,2> > dcp_cred2sig_regions_ncor = get1Dcred_disjoint(h_dcp_ncor,0.9544);
    std::vector<std::array<double,2> > th23_cred1sig_regions_ncor = get1Dcred_disjoint(h_th23_ncor,0.6826);
    std::vector<std::array<double,2> > th23_cred2sig_regions_ncor = get1Dcred_disjoint(h_th23_ncor,0.9544);
    std::vector<std::array<double,2> > th13_cred1sig_regions_ncor = get1Dcred_disjoint(h_th13_ncor,0.6826);
    std::vector<std::array<double,2> > th13_cred2sig_regions_ncor = get1Dcred_disjoint(h_th13_ncor,0.9544);
    std::vector<std::array<double,2> > dm32_cred1sig_regions_ncor = get1Dcred_disjoint(h_dm32_ncor,0.6826);
    std::vector<std::array<double,2> > dm32_cred2sig_regions_ncor = get1Dcred_disjoint(h_dm32_ncor,0.9544);
    
    std::cout << "4.35" << std::endl;
    std::vector<std::array<double,2> > dcp_cred1sig_regions_acor = get1Dcred_disjoint(h_dcp_acor,0.6826);
    std::vector<std::array<double,2> > dcp_cred2sig_regions_acor = get1Dcred_disjoint(h_dcp_acor,0.9544);
    std::vector<std::array<double,2> > th23_cred1sig_regions_acor = get1Dcred_disjoint(h_th23_acor,0.6826);
    std::vector<std::array<double,2> > th23_cred2sig_regions_acor = get1Dcred_disjoint(h_th23_acor,0.9544);
    std::vector<std::array<double,2> > th13_cred1sig_regions_acor = get1Dcred_disjoint(h_th13_acor,0.6826);
    std::vector<std::array<double,2> > th13_cred2sig_regions_acor = get1Dcred_disjoint(h_th13_acor,0.9544);
    std::vector<std::array<double,2> > dm32_cred1sig_regions_acor = get1Dcred_disjoint(h_dm32_acor,0.6826);
    std::vector<std::array<double,2> > dm32_cred2sig_regions_acor = get1Dcred_disjoint(h_dm32_acor,0.9544);
    
    // TODO: write credible intervals and HPD to screen
    
    TCanvas* c_dcp = new TCanvas("c_dcp","c_dcp",600,600);
    TCanvas* c_dm2 = new TCanvas("c_dm2","c_dm2",600,600);
    TCanvas* c_dm2_logy = new TCanvas("c_dm2_logy","c_dm2_logy",600,600);
    TCanvas* c_th23 = new TCanvas("c_th23","c_th23",600,600);
    TCanvas* c_th13 = new TCanvas("c_th13","c_th13",600,600);
    //c->Print("jarl.pdf]");
    
    //TLegend* leg1 = new TLegend(0.35,0.7,0.7,0.85);
    TH1D* shade1sig = new TH1D("shade1sig","shade1sig",1,0.,1.);
    TH1D* shade2sig = new TH1D("shade2sig","shade2sig",1,0.,1.);
    TH1D* asim_l = new TH1D("asim_l","asim_l",1,0.,1.);
    shade1sig->SetFillColorAlpha(kGray,0.50);
    shade2sig->SetFillColorAlpha(kGray,0.25);
    shade1sig->SetLineStyle(9);
    shade2sig->SetLineStyle(2);
    shade1sig->SetLineWidth(3);
    shade2sig->SetLineWidth(3);
    asim_l->SetLineWidth(3);
    asim_l->SetLineColor(kBlack);
    TLegend* leg1 = new TLegend(0.12,0.72,0.85,0.87);
    leg1->SetNColumns(2);
    leg1->SetTextSize(0.035);
    
    if(outTag=="sens") {
        leg1->AddEntry(h_dcp_ncor,"T2K+NOvA","l");
        leg1->AddEntry(shade1sig,"1#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_corr,"T2K-only","l");
        leg1->AddEntry(shade2sig,"2#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_acor,"NOvA-only","l");
        if(drawAsimovPoint && asim_dm2*hierarchy>=0.) leg1->AddEntry(asim_l,"True Simulated Value");
    }
    else if(outTag=="stat") {
        leg1->AddEntry(h_dcp_ncor,"All systematics","l");
        leg1->AddEntry(shade1sig,"1#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_corr,"No systematics","l");
        leg1->AddEntry(shade2sig,"2#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_acor,"Only SK detector systematics","l");
        if(drawAsimovPoint && asim_dm2*hierarchy>=0.) leg1->AddEntry(asim_l,"True Simulated Value");
    }
    else if(outTag=="knob") {
        leg1->AddEntry(h_dcp_ncor,"Old 3-knob","l");
        leg1->AddEntry(shade1sig,"1#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_corr,"New 2-knob","l");
        leg1->AddEntry(shade2sig,"2#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_acor,"New 3-knob","l");
        if(drawAsimovPoint && asim_dm2*hierarchy>=0.) leg1->AddEntry(asim_l,"True Simulated Value");
    }
    else {
        leg1->AddEntry(h_dcp_ncor,"No Correlation","l");
        leg1->AddEntry(shade1sig,"1#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_corr,"100\% Correlation","l");
        leg1->AddEntry(shade2sig,"2#sigma Credible Intervals");
        leg1->AddEntry(h_dcp_acor,"100\% Anticorrelation","l");
        if(drawAsimovPoint && asim_dm2*hierarchy>=0.) leg1->AddEntry(asim_l,"True Simulated Value");
    }
    
    leg1->SetFillStyle(0);
    leg1->SetLineColorAlpha(0,0);
    leg1->SetTextFont(132);
    
    TLegend* leg2 = new TLegend(0.1,0.7,0.7,0.85);
    if(outTag=="sens") {
        leg2->AddEntry(h_dcp_ncor,"T2K+NOvA","l");
        leg2->AddEntry(h_dcp_corr,"T2K-only","l");
        leg2->AddEntry(h_dcp_acor,"NOvA-only","l");
    }
    else if(outTag=="knob") {
        leg2->AddEntry(h_dcp_ncor,"Old 3-knob","l");
        leg2->AddEntry(h_dcp_corr,"New 2-knob","l");
        leg2->AddEntry(h_dcp_acor,"New 3-knob","l");
    }
    
    else {
        leg2->AddEntry(h_dcp_ncor,"No Correlation","l");
        leg2->AddEntry(h_dcp_corr,"100\% Correlation","l");
        leg2->AddEntry(h_dcp_acor,"100\% Anticorrelation","l");
    }
    
    leg2->SetFillStyle(0);
    leg2->SetLineColorAlpha(0,0);
    leg2->SetTextFont(132);
    leg2->SetTextSize(0.07);
    
    h_dcp_ncor->GetYaxis()->SetRangeUser(0.,h_dcp_ncor->GetMaximum()*1.5);
    h_dcp_corr->GetYaxis()->SetRangeUser(0.,h_dcp_corr->GetMaximum()*1.5);
    h_dcp_acor->GetYaxis()->SetRangeUser(0.,h_dcp_acor->GetMaximum()*1.5);
    
    h_th23_ncor->GetYaxis()->SetRangeUser(0.,h_th23_ncor->GetMaximum()*1.5);
    h_th23_corr->GetYaxis()->SetRangeUser(0.,h_th23_corr->GetMaximum()*1.5);
    h_th23_acor->GetYaxis()->SetRangeUser(0.,h_th23_acor->GetMaximum()*1.5);
    
    h_dm32_ncor->GetYaxis()->SetRangeUser(0.000001,h_dm32_ncor->GetMaximum()*1.5);
    h_dm32_corr->GetYaxis()->SetRangeUser(0.000001,h_dm32_corr->GetMaximum()*1.5);
    h_dm32_acor->GetYaxis()->SetRangeUser(0.000001,h_dm32_acor->GetMaximum()*1.5);
    
    /*
     h_dm32_ncor->GetYaxis()->SetRangeUser(0.,h_dm32_ncor->GetMaximum()*1.5);
     h_dm32_corr->GetYaxis()->SetRangeUser(0.,h_dm32_corr->GetMaximum()*1.5);
     h_dm32_acor->GetYaxis()->SetRangeUser(0.,h_dm32_acor->GetMaximum()*1.5);
     */
    
    h_th13_ncor->GetYaxis()->SetRangeUser(0.,h_th13_ncor->GetMaximum()*1.5);
    h_th13_corr->GetYaxis()->SetRangeUser(0.,h_th13_corr->GetMaximum()*1.5);
    h_th13_acor->GetYaxis()->SetRangeUser(0.,h_th13_acor->GetMaximum()*1.5);
    
    h_dcp_ncor->GetYaxis()->SetLabelSize(0.);
    h_dm32_ncor->GetYaxis()->SetLabelSize(0.);
    h_th23_ncor->GetYaxis()->SetLabelSize(0.);
    h_th13_ncor->GetYaxis()->SetLabelSize(0.);
    h_dcp_ncor->GetXaxis()->SetLabelSize(0.04);
    h_dm32_ncor->GetXaxis()->SetLabelSize(0.04);
    h_th23_ncor->GetXaxis()->SetLabelSize(0.04);
    h_th13_ncor->GetXaxis()->SetLabelSize(0.04);
    
    h_dcp_ncor->GetYaxis()->SetTitleSize(0.05);
    h_dm32_ncor->GetYaxis()->SetTitleSize(0.05);
    h_th23_ncor->GetYaxis()->SetTitleSize(0.05);
    h_th13_ncor->GetYaxis()->SetTitleSize(0.05);
    h_dcp_ncor->GetXaxis()->SetTitleSize(0.05);
    h_dm32_ncor->GetXaxis()->SetTitleSize(0.05);
    h_th23_ncor->GetXaxis()->SetTitleSize(0.05);
    h_th13_ncor->GetXaxis()->SetTitleSize(0.05);
    
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    
    c_dcp->Draw();
    c_dcp->cd();
    h_dcp_ncor->SetLineColor(kGreen+3);
    h_dcp_ncor->SetLineWidth(2);
    h_dcp_corr->SetLineColor(kBlue);
    h_dcp_corr->SetLineWidth(2);
    h_dcp_acor->SetLineColor(kRed);
    h_dcp_acor->SetLineWidth(2);
    h_dcp_ncor->Draw("hist");
    h_dcp_corr->Draw("hist same");
    h_dcp_acor->Draw("hist same");
    for(unsigned int a=0;a<dcp_cred1sig_regions_corr.size();a++) {
        drawCredInterval(c_dcp,dcp_cred1sig_regions_corr[a][0],dcp_cred1sig_regions_corr[a][1],h_dcp_ncor->GetMaximum()*0.25,h_dcp_ncor->GetMaximum()*0.50,h_dcp_corr->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<dcp_cred2sig_regions_corr.size();a++) {
        drawCredInterval(c_dcp,dcp_cred2sig_regions_corr[a][0],dcp_cred2sig_regions_corr[a][1],h_dcp_ncor->GetMaximum()*0.25,h_dcp_ncor->GetMaximum()*0.50,h_dcp_corr->GetLineColor(),0.25,2);
    }
    for(unsigned int a=0;a<dcp_cred1sig_regions_ncor.size();a++) {
        drawCredInterval(c_dcp,dcp_cred1sig_regions_ncor[a][0],dcp_cred1sig_regions_ncor[a][1],h_dcp_ncor->GetMaximum()*0.50,h_dcp_ncor->GetMaximum()*0.75,h_dcp_ncor->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<dcp_cred2sig_regions_ncor.size();a++) {
        drawCredInterval(c_dcp,dcp_cred2sig_regions_ncor[a][0],dcp_cred2sig_regions_ncor[a][1],h_dcp_ncor->GetMaximum()*0.50,h_dcp_ncor->GetMaximum()*0.75,h_dcp_ncor->GetLineColor(),0.25,2);
    }
    for(unsigned int a=0;a<dcp_cred1sig_regions_acor.size();a++) {
        drawCredInterval(c_dcp,dcp_cred1sig_regions_acor[a][0],dcp_cred1sig_regions_acor[a][1],h_dcp_ncor->GetMaximum()*0.00,h_dcp_ncor->GetMaximum()*0.25,h_dcp_acor->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<dcp_cred2sig_regions_acor.size();a++) {
        drawCredInterval(c_dcp,dcp_cred2sig_regions_acor[a][0],dcp_cred2sig_regions_acor[a][1],h_dcp_ncor->GetMaximum()*0.00,h_dcp_ncor->GetMaximum()*0.25,h_dcp_acor->GetLineColor(),0.25,2);
    }
    TLine* l_asim_dcp = new TLine(asim_dcp,0.,asim_dcp,0.77*h_dcp_ncor->GetMaximum());
    l_asim_dcp->SetLineWidth(3);
    if(drawAsimovPoint && asim_dm2*hierarchy>=0.) l_asim_dcp->Draw();
    leg1->Draw();
    gPad->SetTickx();
    gPad->SetTicky();
    gStyle->SetLegendBorderSize(0);
    gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.13);
    //gPad->SetTopMargin(0.08);
    gPad->RedrawAxis();
    
    c_th23->Draw();
    c_th23->cd();
    h_th23_ncor->SetLineColor(kGreen+3);
    h_th23_ncor->SetLineWidth(2);
    h_th23_corr->SetLineColor(kBlue);
    h_th23_corr->SetLineWidth(2);
    h_th23_acor->SetLineColor(kRed);
    h_th23_acor->SetLineWidth(2);
    h_th23_ncor->Draw("hist");
    h_th23_corr->Draw("hist same");
    h_th23_acor->Draw("hist same");
    for(unsigned int a=0;a<th23_cred1sig_regions_corr.size();a++) {
        drawCredInterval(c_th23,th23_cred1sig_regions_corr[a][0],th23_cred1sig_regions_corr[a][1],h_th23_ncor->GetMaximum()*0.25,h_th23_ncor->GetMaximum()*0.50,h_th23_corr->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<th23_cred2sig_regions_corr.size();a++) {
        drawCredInterval(c_th23,th23_cred2sig_regions_corr[a][0],th23_cred2sig_regions_corr[a][1],h_th23_ncor->GetMaximum()*0.25,h_th23_ncor->GetMaximum()*0.50,h_th23_corr->GetLineColor(),0.25,2);
    }
    for(unsigned int a=0;a<th23_cred1sig_regions_ncor.size();a++) {
        drawCredInterval(c_th23,th23_cred1sig_regions_ncor[a][0],th23_cred1sig_regions_ncor[a][1],h_th23_ncor->GetMaximum()*0.50,h_th23_ncor->GetMaximum()*0.75,h_th23_ncor->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<th23_cred2sig_regions_ncor.size();a++) {
        drawCredInterval(c_th23,th23_cred2sig_regions_ncor[a][0],th23_cred2sig_regions_ncor[a][1],h_th23_ncor->GetMaximum()*0.50,h_th23_ncor->GetMaximum()*0.75,h_th23_ncor->GetLineColor(),0.25,2);
    }
    for(unsigned int a=0;a<th23_cred1sig_regions_acor.size();a++) {
        drawCredInterval(c_th23,th23_cred1sig_regions_acor[a][0],th23_cred1sig_regions_acor[a][1],h_th23_ncor->GetMaximum()*0.00,h_th23_ncor->GetMaximum()*0.25,h_th23_acor->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<th23_cred2sig_regions_acor.size();a++) {
        drawCredInterval(c_th23,th23_cred2sig_regions_acor[a][0],th23_cred2sig_regions_acor[a][1],h_th23_ncor->GetMaximum()*0.00,h_th23_ncor->GetMaximum()*0.25,h_th23_acor->GetLineColor(),0.25,2);
    }
    TLine* l_asim_th23 = new TLine(asim_th23,0.,asim_th23,0.77*h_th23_ncor->GetMaximum());
    l_asim_th23->SetLineWidth(3);
    if(drawAsimovPoint && asim_dm2*hierarchy>=0.) l_asim_th23->Draw();
    leg1->Draw();
    gPad->SetTickx();
    gPad->SetTicky();
    gStyle->SetLegendBorderSize(0);
    gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.13);
    //gPad->SetTopMargin(0.08);
    gPad->RedrawAxis();
    
    c_th13->Draw();
    c_th13->cd();
    h_th13_ncor->SetLineColor(kGreen+3);
    h_th13_ncor->SetLineWidth(2);
    h_th13_corr->SetLineColor(kBlue);
    h_th13_corr->SetLineWidth(2);
    h_th13_acor->SetLineColor(kRed);
    h_th13_acor->SetLineWidth(2);
    h_th13_ncor->Draw("hist");
    h_th13_corr->Draw("hist same");
    h_th13_acor->Draw("hist same");
    
    for(unsigned int a=0;a<th13_cred1sig_regions_corr.size();a++) {
        drawCredInterval(c_th13,th13_cred1sig_regions_corr[a][0],th13_cred1sig_regions_corr[a][1],h_th13_ncor->GetMaximum()*0.25,h_th13_ncor->GetMaximum()*0.50,h_th13_corr->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<th13_cred2sig_regions_corr.size();a++) {
        drawCredInterval(c_th13,th13_cred2sig_regions_corr[a][0],th13_cred2sig_regions_corr[a][1],h_th13_ncor->GetMaximum()*0.25,h_th13_ncor->GetMaximum()*0.50,h_th13_corr->GetLineColor(),0.25,2);
    }
    for(unsigned int a=0;a<th13_cred1sig_regions_ncor.size();a++) {
        drawCredInterval(c_th13,th13_cred1sig_regions_ncor[a][0],th13_cred1sig_regions_ncor[a][1],h_th13_ncor->GetMaximum()*0.50,h_th13_ncor->GetMaximum()*0.75,h_th13_ncor->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<th13_cred2sig_regions_ncor.size();a++) {
        drawCredInterval(c_th13,th13_cred2sig_regions_ncor[a][0],th13_cred2sig_regions_ncor[a][1],h_th13_ncor->GetMaximum()*0.50,h_th13_ncor->GetMaximum()*0.75,h_th13_ncor->GetLineColor(),0.25,2);
    }
    for(unsigned int a=0;a<th13_cred1sig_regions_acor.size();a++) {
        drawCredInterval(c_th13,th13_cred1sig_regions_acor[a][0],th13_cred1sig_regions_acor[a][1],h_th13_ncor->GetMaximum()*0.00,h_th13_ncor->GetMaximum()*0.25,h_th13_acor->GetLineColor(),0.25,9);
    }
    for(unsigned int a=0;a<th13_cred2sig_regions_acor.size();a++) {
        drawCredInterval(c_th13,th13_cred2sig_regions_acor[a][0],th13_cred2sig_regions_acor[a][1],h_th13_ncor->GetMaximum()*0.00,h_th13_ncor->GetMaximum()*0.25,h_th13_acor->GetLineColor(),0.25,2);
    }
    TLine* l_asim_th13 = new TLine(asim_th13,0.,asim_th13,0.77*h_th13_ncor->GetMaximum());
    l_asim_th13->SetLineWidth(3);
    if(drawAsimovPoint && asim_dm2*hierarchy>=0.) l_asim_th13->Draw();
    leg1->Draw();
    gPad->SetTickx();
    gPad->SetTicky();
    gStyle->SetLegendBorderSize(0);
    gPad->SetLeftMargin(0.1);
    //gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.13);
    //gPad->SetTopMargin(0.08);
    gPad->RedrawAxis();
    
    if(hierarchy ==0 ) {
        
        c_dm2->Draw();
        c_dm2->cd();
        h_dm32_ncor->SetLineColor(kGreen+3);
        h_dm32_ncor->SetLineWidth(2);
        h_dm32_corr->SetLineColor(kBlue);
        h_dm32_corr->SetLineWidth(2);
        h_dm32_acor->SetLineColor(kRed);
        h_dm32_acor->SetLineWidth(2);
        
        TPad *p_NH = new TPad("p_NH","p_NH",0.5,0.0,1.0,1.0);
        TPad *p_IH = new TPad("p_IH","p_IH",0.0,0.0,0.5,1.0);
        p_NH->SetLeftMargin(0.01);
        p_NH->SetRightMargin(0.17);
        p_IH->SetRightMargin(0.01);
        p_IH->SetLeftMargin(0.17);
        p_NH->Draw(); p_IH->Draw();
        
        TPad *p_NH_logy = (TPad*)p_NH->Clone();
        p_NH_logy->SetLogy();
        TPad *p_IH_logy = (TPad*)p_IH->Clone();
        p_IH_logy->SetLogy();
        
        TH1D* dummyNH = (TH1D*)h_dm32_ncor->Clone("dummyNH");
        TH1D* dummyIH = (TH1D*)h_dm32_ncor->Clone("dummyIH");
        dummyNH->Reset(); dummyIH->Reset();
        dummyNH->GetXaxis()->SetRangeUser(0.00225,0.00275);
        dummyIH->GetXaxis()->SetRangeUser(-0.00275,-0.00225);
        //dummyNH->GetYaxis()->SetRangeUser(0.,h_dm32_ncor->GetMaximum()*1.5);
        //dummyIH->GetYaxis()->SetRangeUser(0.,h_dm32_ncor->GetMaximum()*1.5);
        
        dummyNH->GetXaxis()->SetLabelSize(0.07);
        dummyIH->GetXaxis()->SetLabelSize(0.07);
        dummyNH->GetXaxis()->SetTitleSize(0.08);
        dummyIH->GetXaxis()->SetTitleSize(0.08);
        dummyIH->GetYaxis()->SetTitleSize(0.08);
        dummyNH->GetXaxis()->SetTitleOffset(0.6);
        dummyIH->GetXaxis()->SetTitleOffset(1.5);
        dummyIH->GetYaxis()->SetTitleOffset(0.8);
        dummyNH->GetXaxis()->SetLabelOffset(-0.01);
        dummyIH->GetXaxis()->SetLabelOffset(-0.01);
        dummyIH->GetXaxis()->SetNdivisions(512);
        dummyNH->GetXaxis()->SetNdivisions(512);
        TGaxis::SetExponentOffset(0.01, 0.0, "x");
        
        p_NH->cd();
        dummyNH->Draw();
        h_dm32_ncor->Draw("hist same");
        h_dm32_corr->Draw("hist same");
        h_dm32_acor->Draw("hist same");
        for(unsigned int a=0;a<dm32_cred1sig_regions_corr.size();a++) {
            drawCredInterval(p_NH,dm32_cred1sig_regions_corr[a][0],dm32_cred1sig_regions_corr[a][1],h_dm32_ncor->GetMaximum()*0.25,h_dm32_ncor->GetMaximum()*0.50,h_dm32_corr->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_corr.size();a++) {
            drawCredInterval(p_NH,dm32_cred2sig_regions_corr[a][0],dm32_cred2sig_regions_corr[a][1],h_dm32_ncor->GetMaximum()*0.25,h_dm32_ncor->GetMaximum()*0.50,h_dm32_corr->GetLineColor(),0.25,2);
        }
        for(unsigned int a=0;a<dm32_cred1sig_regions_ncor.size();a++) {
            drawCredInterval(p_NH,dm32_cred1sig_regions_ncor[a][0],dm32_cred1sig_regions_ncor[a][1],h_dm32_ncor->GetMaximum()*0.50,h_dm32_ncor->GetMaximum()*0.75,h_dm32_ncor->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_ncor.size();a++) {
            drawCredInterval(p_NH,dm32_cred2sig_regions_ncor[a][0],dm32_cred2sig_regions_ncor[a][1],h_dm32_ncor->GetMaximum()*0.50,h_dm32_ncor->GetMaximum()*0.75,h_dm32_ncor->GetLineColor(),0.25,2);
        }
        for(unsigned int a=0;a<dm32_cred1sig_regions_acor.size();a++) {
            drawCredInterval(p_NH,dm32_cred1sig_regions_acor[a][0],dm32_cred1sig_regions_acor[a][1],h_dm32_ncor->GetMaximum()*0.00,h_dm32_ncor->GetMaximum()*0.25,h_dm32_acor->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_acor.size();a++) {
            drawCredInterval(p_NH,dm32_cred2sig_regions_acor[a][0],dm32_cred2sig_regions_acor[a][1],h_dm32_ncor->GetMaximum()*0.00,h_dm32_ncor->GetMaximum()*0.25,h_dm32_acor->GetLineColor(),0.25,2);
        }
        TLine* l_asim_dm2 = new TLine(asim_dm2,0.,asim_dm2,0.77*h_dm32_ncor->GetMaximum());
        l_asim_dm2->SetLineWidth(3);
        if(drawAsimovPoint) l_asim_dm2->Draw();
        c_dm2->cd();
        leg1->Draw();
        
        p_NH->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        gStyle->SetLegendBorderSize(0);
        
        p_IH->cd();
        dummyIH->Draw();
        h_dm32_ncor->Draw("hist same");
        h_dm32_corr->Draw("hist same");
        h_dm32_acor->Draw("hist same");
        for(unsigned int a=0;a<dm32_cred1sig_regions_corr.size();a++) {
            drawCredInterval(p_IH,dm32_cred1sig_regions_corr[a][0],dm32_cred1sig_regions_corr[a][1],h_dm32_ncor->GetMaximum()*0.25,h_dm32_ncor->GetMaximum()*0.50,h_dm32_corr->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_corr.size();a++) {
            drawCredInterval(p_IH,dm32_cred2sig_regions_corr[a][0],dm32_cred2sig_regions_corr[a][1],h_dm32_ncor->GetMaximum()*0.25,h_dm32_ncor->GetMaximum()*0.50,h_dm32_corr->GetLineColor(),0.25,2);
        }
        for(unsigned int a=0;a<dm32_cred1sig_regions_ncor.size();a++) {
            drawCredInterval(p_IH,dm32_cred1sig_regions_ncor[a][0],dm32_cred1sig_regions_ncor[a][1],h_dm32_ncor->GetMaximum()*0.50,h_dm32_ncor->GetMaximum()*0.75,h_dm32_ncor->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_ncor.size();a++) {
            drawCredInterval(p_IH,dm32_cred2sig_regions_ncor[a][0],dm32_cred2sig_regions_ncor[a][1],h_dm32_ncor->GetMaximum()*0.50,h_dm32_ncor->GetMaximum()*0.75,h_dm32_ncor->GetLineColor(),0.25,2);
        }
        for(unsigned int a=0;a<dm32_cred1sig_regions_acor.size();a++) {
            drawCredInterval(p_IH,dm32_cred1sig_regions_acor[a][0],dm32_cred1sig_regions_acor[a][1],h_dm32_ncor->GetMaximum()*0.00,h_dm32_ncor->GetMaximum()*0.25,h_dm32_acor->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_acor.size();a++) {
            drawCredInterval(p_IH,dm32_cred2sig_regions_acor[a][0],dm32_cred2sig_regions_acor[a][1],h_dm32_ncor->GetMaximum()*0.00,h_dm32_ncor->GetMaximum()*0.25,h_dm32_acor->GetLineColor(),0.25,2);
        }
        if(drawAsimovPoint) l_asim_dm2->Draw();
        gPad->SetTickx();
        gPad->SetTicky();
        gStyle->SetLegendBorderSize(0);
        
        /*
         // log y-axis version
         c_dm2_logy->Draw(); c_dm2_logy->cd();
         p_NH_logy->Draw();
         p_IH_logy->Draw();
         
         TH1D* dummyNH_logy = (TH1D*)dummyNH->Clone("dummNH_logy");
         TH1D* dummyIH_logy = (TH1D*)dummyIH->Clone("dummIH_logy");
         
         dummyNH_logy->GetYaxis()->SetRangeUser(0.1,h_dm32_ncor->GetMaximum()*1000.);
         dummyIH_logy->GetYaxis()->SetRangeUser(0.1,h_dm32_ncor->GetMaximum()*1000.);
         
         p_NH_logy->cd();
         dummyNH_logy->Draw();
         h_dm32_ncor->Draw("hist same");
         h_dm32_corr->Draw("hist same");
         h_dm32_acor->Draw("hist same");
         leg2->Draw();
         p_NH_logy->SetLogy();
         gPad->SetTickx();
         gPad->SetTicky();
         gPad->RedrawAxis();
         
         p_IH_logy->cd();
         dummyIH_logy->Draw();
         h_dm32_ncor->Draw("hist same");
         h_dm32_corr->Draw("hist same");
         h_dm32_acor->Draw("hist same");
         p_IH_logy->SetLogy();
         gPad->SetTickx();
         gPad->SetTicky();
         gPad->RedrawAxis();
         */
        
    }
    else {
        if(hierarchy == 1) h_dm32_ncor->GetXaxis()->SetRangeUser(0.00225,0.00275);
        else if(hierarchy == -1) h_dm32_ncor->GetXaxis()->SetRangeUser(-0.00275,-0.00225);
        c_dm2->Draw();
        c_dm2->cd();
        h_dm32_ncor->SetLineColor(kGreen+3);
        h_dm32_ncor->SetLineWidth(2);
        h_dm32_corr->SetLineColor(kBlue);
        h_dm32_corr->SetLineWidth(2);
        h_dm32_acor->SetLineColor(kRed);
        h_dm32_acor->SetLineWidth(2);
        h_dm32_ncor->Draw("hist");
        h_dm32_corr->Draw("hist same");
        h_dm32_acor->Draw("hist same");
        for(unsigned int a=0;a<dm32_cred1sig_regions_corr.size();a++) {
            drawCredInterval(c_dm2,dm32_cred1sig_regions_corr[a][0],dm32_cred1sig_regions_corr[a][1],h_dm32_ncor->GetMaximum()*0.25,h_dm32_ncor->GetMaximum()*0.50,h_dm32_corr->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_corr.size();a++) {
            drawCredInterval(c_dm2,dm32_cred2sig_regions_corr[a][0],dm32_cred2sig_regions_corr[a][1],h_dm32_ncor->GetMaximum()*0.25,h_dm32_ncor->GetMaximum()*0.50,h_dm32_corr->GetLineColor(),0.25,2);
        }
        for(unsigned int a=0;a<dm32_cred1sig_regions_ncor.size();a++) {
            drawCredInterval(c_dm2,dm32_cred1sig_regions_ncor[a][0],dm32_cred1sig_regions_ncor[a][1],h_dm32_ncor->GetMaximum()*0.50,h_dm32_ncor->GetMaximum()*0.75,h_dm32_ncor->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_ncor.size();a++) {
            drawCredInterval(c_dm2,dm32_cred2sig_regions_ncor[a][0],dm32_cred2sig_regions_ncor[a][1],h_dm32_ncor->GetMaximum()*0.50,h_dm32_ncor->GetMaximum()*0.75,h_dm32_ncor->GetLineColor(),0.25,2);
        }
        for(unsigned int a=0;a<dm32_cred1sig_regions_acor.size();a++) {
            drawCredInterval(c_dm2,dm32_cred1sig_regions_acor[a][0],dm32_cred1sig_regions_acor[a][1],h_dm32_ncor->GetMaximum()*0.00,h_dm32_ncor->GetMaximum()*0.25,h_dm32_acor->GetLineColor(),0.25,9);
        }
        for(unsigned int a=0;a<dm32_cred2sig_regions_acor.size();a++) {
            drawCredInterval(c_dm2,dm32_cred2sig_regions_acor[a][0],dm32_cred2sig_regions_acor[a][1],h_dm32_ncor->GetMaximum()*0.00,h_dm32_ncor->GetMaximum()*0.25,h_dm32_acor->GetLineColor(),0.25,2);
        }
        TLine* l_asim_dm2 = new TLine(asim_dm2,0.,asim_dm2,0.77*h_dm32_ncor->GetMaximum());
        l_asim_dm2->SetLineWidth(3);
        if(drawAsimovPoint && asim_dm2*hierarchy>=0.) l_asim_dm2->Draw();
        leg1->Draw();
    }
    std::cout << "5" << std::endl;
    
    std::ostringstream outBaseName;
    if(outTag=="sens") outBaseName << "sensitivityComp/t2knova"<<asimov<<"_v03Sep2021_sensContours_";
    else if(outTag=="stat") outBaseName << "statComp/t2knova"<<asimov<<"_v03Sep2021_statContours_";
    else if(outTag=="knob") outBaseName << "data/v03Sep2021_nuenumucorr/nova2v3knob/" << "t2knova"<<asimov<<"_v03Sep2021_nova2v3knob__";
    else outBaseName << outTag << "SystCorrComp/t2knova"<< asimov << "_v03Sep2021_" << outTag << "CorrAnticorrComp_";
    string tempName;
    if(hierarchy == 1) {
        tempName = outBaseName.str() + "dcp_NH.pdf";
        c_dcp->Print(tempName.c_str());
        tempName = outBaseName.str() + "dm2_NH.pdf";
        c_dm2->Print(tempName.c_str());
        //tempName = outBaseName.str() + "dm2LogY_NH.pdf";
        //c_dm2_logy->Print(tempName.c_str());
        tempName = outBaseName.str() + "th13_NH.pdf";
        c_th13->Print(tempName.c_str());
        tempName = outBaseName.str() + "th23_NH.pdf";
        c_th23->Print(tempName.c_str());
    }
    else if(hierarchy == 0) {
        tempName = outBaseName.str() + "dcp_both.pdf";
        c_dcp->Print(tempName.c_str());
        tempName = outBaseName.str() + "dm2_both.pdf";
        c_dm2->Print(tempName.c_str());
        //tempName = outBaseName.str() + "dm2LogY_both.pdf";
        //c_dm2_logy->Print(tempName.c_str());
        tempName = outBaseName.str() + "th13_both.pdf";
        c_th13->Print(tempName.c_str());
        tempName = outBaseName.str() + "th23_both.pdf";
        c_th23->Print(tempName.c_str());
    }
    else if(hierarchy == -1) {
        tempName = outBaseName.str() + "dcp_IH.pdf";
        c_dcp->Print(tempName.c_str());
        tempName = outBaseName.str() + "dm2_IH.pdf";
        c_dm2->Print(tempName.c_str());
        //tempName = outBaseName.str() + "dm2LogY_IH.pdf";
        //c_dm2_logy->Print(tempName.c_str());
        tempName = outBaseName.str() + "th13_IH.pdf";
        c_th13->Print(tempName.c_str());
        tempName = outBaseName.str() + "th23_IH.pdf";
        c_th23->Print(tempName.c_str());
    }
    
}

void draw1Dcreds_3comp() {
    
    /*
     // t2k vs nova vs t2k+nova Asimov0
     makePlot(0,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov0-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov0_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov0-03Sep2021_job1_red.root","sens","Asimov0");
     makePlot(1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov0-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov0_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov0-03Sep2021_job1_red.root","sens","Asimov0");
     makePlot(-1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov0_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov0-03Sep2021_job1_red.root","sens","Asimov0");
     
     // t2k vs nova vs t2k+nova Asimov1
     makePlot(0,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov1-03Sep2021_job1_red.root","sens","Asimov1");
     makePlot(1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov1-03Sep2021_job1_red.root","sens","Asimov1");
     makePlot(-1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov1-03Sep2021_job1_red.root","sens","Asimov1");
     
     // t2k vs nova vs t2k+nova Asimov4
     makePlot(0,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov4_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov4-03Sep2021_job1_red.root","sens","Asimov4");
     makePlot(1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov4_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov4-03Sep2021_job1_red.root","sens","Asimov4");
     makePlot(-1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021/MaCh3_t2konly_Asimov4_03Sep2021_job3_ch0-999_red.root","data/v03Sep2021/MaCh3-novaAsimov4-03Sep2021_job1_red.root","sens","Asimov4");
     */
    
    // stat vs. stat + det syst. vs. stat + all syst
    //!! not really!!! hacking here...
    //  makePlot(0,"data/v03Sep2021/MaCh3-t2knovaAsimov0-03Sep2021_narval_job0-2_red.root","data/v03Sep2021/t2k_nova_Asimov0_sensitivity_v4_all_ND_data_wRC_STAT_ONLY_sum_red.root","data/v03Sep2021/t2k_nova_Asimov0_sensitivity_v4_container_all_ND_data_wRC_STAT_AND_SKDET_sum_red.root","stat","Asimov0");
    //  makePlot(0,"data/v03Sep2021/MaCh3-t2knovaAsimov1-03Sep2021-grahamJob0-2-seawulfJob1-6_red.root","data/v03Sep2021/t2k_nova_Asimov1_sensitivity_v4_all_ND_data_wRC_STAT_ONLY_sum_red.root","data/v03Sep2021/t2k_nova_Asimov1_sensitivity_v4_container_all_ND_data_wRC_STAT_AND_SKDET_sum_red.root","stat","Asimov1");
    //  makePlot(0,"data/v03Sep2021/MaCh3-t2knovaAsimov4-03Sep2021_narval_job0-2_red.root","data/v03Sep2021/t2k_nova_Asimov4_sensitivity_v4_all_ND_data_wRC_STAT_ONLY_sum_red.root","data/v03Sep2021/t2k_nova_Asimov4_sensitivity_v4_container_all_ND_data_wRC_STAT_AND_SKDET_sum_red.root","stat","Asimov4");
    
    // nova knob comp
    //makePlot(0,"data/v03Sep2021/MaCh3-t2knovaAsimov1-03Sep2021-grahamJob0-2-seawulfJob1-6_red.root","data/v03Sep2021_nuenumucorr/nova2v3knob/MaCh3-t2knovaAsimov1-03Sep2021_novaOnly_2knob_job1_red_reweighted.root","data/v03Sep2021_nuenumucorr/nova2v3knob/MaCh3-t2knovaAsimov1-03Sep2021_novaOnly_3knob_job1_red_reweighted.root","knob","Asimov1");
    
    // dm2 correlations:
    
    makePlot(0,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2nocorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2corr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2anticorr_woRC_jobs4-5_red_reweighted.root","dm2","Asimov4");
    makePlot(1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2nocorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2corr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2anticorr_woRC_jobs4-5_red_reweighted.root","dm2","Asimov4");
    makePlot(-1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2nocorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2corr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_dm2anticorr_woRC_jobs4-5_red_reweighted.root","dm2","Asimov4");
    
    //makePlot(0,"data/MaCh3-T2KNOvA-08Mar2021update_reduced.root","data/MaCh3-t2knovaAsimov1-24Mar2021_dm32Corr_job1-6_ch1-20_reduced.root","data/MaCh3-t2knovaAsimov1-24Mar2021_dm2AntiCorr_reduced.root","dm2");
    
    // dcp correlations:
    ///*
    makePlot(0,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_NOnuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumuanticorr_woRC_jobs4-5_red_reweighted.root","dcp","Asimov1");
    makePlot(1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_NOnuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumuanticorr_woRC_jobs4-5_red_reweighted.root","dcp","Asimov1");
    makePlot(-1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_NOnuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumuanticorr_woRC_jobs4-5_red_reweighted.root","dcp","Asimov1");
    //*/
    
    //makePlot(0,"data/MaCh3-T2KNOvA-08Mar2021update_reduced.root","data/MaCh3-t2knovaAsimov1-nue_numu-2ndclasscurr_corr_job1-6_ch1-20_reduced.root","data/MaCh3-t2knovaAsimov1-nue_numu-2ndclasscurr_antiCorr_job1-6_ch1-20_reduced.root","dcp");
    //makePlot(0,"data/v03Sep2021/MaCh3-t2knovaAsimov1-03Sep2021-grahamJob0-2-seawulfJob1-6_red.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumuanticorr_woRC_jobs4-5_red_reweighted.root","dcp","Asimov1");
    //makePlot(1,"data/v03Sep2021/MaCh3-t2knovaAsimov1-03Sep2021-grahamJob0-2-seawulfJob1-6_red.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumuanticorr_woRC_jobs4-5_red_reweighted.root","dcp","Asimov1");
    //makePlot(-1,"data/v03Sep2021/MaCh3-t2knovaAsimov1-03Sep2021-grahamJob0-2-seawulfJob1-6_red.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumuanticorr_woRC_jobs4-5_red_reweighted.root","dcp","Asimov1");
    
    // th23 correlations:
    ///*
    makePlot(0,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_th23corr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_th23anticorr_woRC_jobs4-5_red_reweighted.root","th23","Asimov4");
    makePlot(1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_th23corr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_th23anticorr_woRC_jobs4-5_red_reweighted.root","th23","Asimov4");
    makePlot(-1,"data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_th23corr_woRC_jobs4-5_red_reweighted.root","data/v03Sep2021_nuenumucorr/MaCh3-t2knovaAsimov4-v03Sep2021_nuenumucorr_th23anticorr_woRC_jobs4-5_red_reweighted.root","th23","Asimov4");
    //*/
    //makePlot(0,"data/MaCh3-T2KNOvA-08Mar2021update_reduced.root","data/MaCh3-t2knovaAsimov1-24Mar2021_th23Corr_ICLGPUandSeaWulf_reduced.root","data/MaCh3-t2knovaAsimov1-24Mar2021_th23AntiCorr_job1-4_ch1-20_plusICLGPU_reduced.root","th23");
    
}
