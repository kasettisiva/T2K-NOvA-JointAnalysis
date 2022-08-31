//To compare MaCh3 and P-theta sensitivities for t2k-nova joint analysis
#include <iostream>
#include <string.h>

void Compare_MaCh3_Ptheta_t2knova_28June2022() {
	
	//TFile *fm=new TFile("MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
	TFile *fm=new TFile("MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
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
	
	TH1D* hm=new TH1D("hm","hm",51,-3.2044245,3.2044245); //-TMath::Pi(),TMath::Pi());
	TH1D* hmdchi2 =new TH1D("hmdchi2","hmdchi2",51,-3.2044245,3.2044245);
	for(int i=0;i<nentries;++i) {
		trm->GetEntry(i);
		hm->Fill(dcp,RCreweight);
		if(dm2>0) hmdchi2->Fill(dcp,RCreweight); //Normal hierarchy
	}
	hm->Scale(1/hm->Integral()); //normalize
	hmdchi2->Scale(1/hmdchi2->Integral()); //normalize
	
	TCanvas *c=new TCanvas("c","c",700,700);
	hm->Draw("HIST");
	
	//TFile *fp = TFile::Open("credible_T2KNOvA_both.root");
	TFile *fp = TFile::Open("credible_T2K_both.root");
	TH1D* hp=(TH1D*)fp->Get("LL_both");
	//hp->SetLineColor(kOrange);
	//hp->Draw("same");
	hp->Rebin(51);
	
	TH1D* htemp=new TH1D("htemp", "hemp",51,-3.2044245,3.2044245);
	for(int i=0;i<51;++i) {
		htemp->SetBinContent(i+1, hp->GetBinContent(i+1));
	}
	
	//TH1D* htemp=(TH1D*)hp->Clone();
	htemp->SetLineColor(kRed);
	htemp->Draw("HIST same");
	
	//Plot Dchi2
	TCanvas *c1=new TCanvas("c1","c1",700,700);
	Double_t bincontmax=hmdchi2->GetMaximum();
	
	for(int i=0;i<hm->GetNbinsX();++i) {
		hmdchi2->SetBinContent(i+1, -2*TMath::Log(hm->GetBinContent(i+1)/bincontmax));
	}
	
	TFile *fpdchi2 = new TFile("contour_dCP_100k_T2K.root");
	TH1D *hpdchi2=(TH1D*)fpdchi2->Get("SA_dcp");

	hpdchi2->SetLineColor(kRed);
	hmdchi2->Draw("HIST");
	hpdchi2->Draw("Csame");

	

}
