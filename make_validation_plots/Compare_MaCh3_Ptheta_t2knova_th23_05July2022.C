//To compare MaCh3 and P-theta sensitivities for t2k-nova joint analysis
#include <iostream>
#include <string.h>

struct Binning {
	Int_t nbins;
	Double_t bin_min;
	Double_t bin_max;
};

const char* oscparname;

Binning dcp_binning = {51,-3.2044245,3.2044245};
Binning dm2_binning = {81,2.25e-3,2.75e-3};
Binning dm2_io_binning = {81,-2.75e-3,-2.25e-3};
Binning th23_binning = {81,0.38,0.62};
Binning binning;
Binning binning_io;

void Compare_MaCh3_Ptheta_t2knova_th23_05July2022(const char* expt, const char* osc_par) {

	gStyle->SetOptStat(0);
	
	if(strstr(osc_par,"dCP") != NULL) {
		std::cout << " dcp " << std::endl;
		binning = dcp_binning;
		binning_io = dcp_binning;
		oscparname = "dcp";
	}
	else if(strstr(osc_par,"dm2") != NULL) {
		std::cout << " dm2 " << std::endl;
		binning = dm2_binning;
		binning_io = dm2_io_binning;
		oscparname = "dm2";
	}
	else if(strstr(osc_par,"th23") != NULL) {
		std::cout << " th23 " << std::endl;
		binning = th23_binning;
		binning_io = th23_binning;
		oscparname = "sin223";
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
		fm=new TFile("MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
	}
	else {
		std::cout << " MaCh3 T2K only input file: " << std::endl;
		fm=new TFile("MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
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
	
	TH1D* hm=new TH1D("hm","hm",           binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* hm_no =new TH1D("hm_no","hm_no", binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* hm_io =new TH1D("hm_io","hm_io", binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
	
	for(int i=0;i<nentries;++i) {
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
	hm->Scale(1/hm->Integral()); //normalize
	hm_no->Scale(1/hm_no->Integral()); //normalize
	hm_io->Scale(1/hm_io->Integral()); //normalize

	TLegend* leg;
	//leg=new TLegend(0.6,0.67,0.75,0.85);
	leg=new TLegend(0.2,0.7,0.35,0.81);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	
	TCanvas *c_both=new TCanvas("c_both","c_both",700,700);
	
	//TFile *fp = TFile::Open("credible_T2KNOvA_both.root");
	TFile *fp_both = TFile::Open(Form("credible_%s_%s_both.root",osc_par,expt));
	TH1D* hp=(TH1D*)fp_both->Get("LL_both");
	//hp->SetLineColor(kOrange);
	//hp->Draw("same");
	if(strstr(osc_par,"dCP") != NULL) hp->Rebin(51);
	else hp->Rebin(81);
	
	TH1D* htemp=new TH1D("htemp", "hemp",binning.nbins,binning.bin_min,binning.bin_max);
	for(int i=0;i<binning.nbins;++i) {
		htemp->SetBinContent(i+1, hp->GetBinContent(i+1));
	}
	
	leg->AddEntry(hm,"MaCh3","l");
	leg->AddEntry(htemp,"Ptheta","l");
	
	hm->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
	hm->GetYaxis()->SetTitle("Posterior probability");
	hm->SetTitle("Both hierarchies");
	htemp->SetLineColor(kRed);
	
	hm->Draw("HIST");
	htemp->Draw("HIST same");
	leg->Draw("same");
	c_both->SaveAs(Form("compare_credible_mach3_ptheta_%s_both_%s_05July2022.pdf",osc_par,expt));
	
	//Normal Ordering
	TCanvas *c_no=new TCanvas("c_no","c_no",700,700);
	TFile *fp_no = TFile::Open(Form("credible_%s_%s.root",osc_par,expt));
	TH1D* hp_no=(TH1D*)fp_no->Get("LL");
	
	if(strstr(osc_par,"dCP") != NULL) hp_no->Rebin(51);
	else hp_no->Rebin(81);
	
	TH1D* htemp_no=new TH1D("htemp", "hemp",binning.nbins,binning.bin_min,binning.bin_max);
	for(int i=0;i<binning.nbins;++i) {
		htemp_no->SetBinContent(i+1, hp_no->GetBinContent(i+1));
	}
	
	hm_no->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
	hm_no->GetYaxis()->SetTitle("Posterior probability");
	hm_no->SetTitle("NO");
	hm_no->Draw("HIST");
	htemp_no->SetLineColor(kRed);
	htemp_no->Draw("HIST same");
	leg->Draw("same");
	c_no->SaveAs(Form("compare_credible_mach3_ptheta_%s_NO_%s_05July2022.pdf",osc_par,expt));
	
	//Inverted Ordering
	TCanvas *c_io=new TCanvas("c_io","c_io",700,700);
	TFile *fp_io = TFile::Open(Form("credible_%s_%s_IH.root",osc_par,expt));
	TH1D* hp_io=(TH1D*)fp_io->Get("LL");
	
	if(strstr(osc_par,"dCP") != NULL) hp_io->Rebin(51);
	else hp_io->Rebin(81);
	
	TH1D* htemp_io=new TH1D("htemp_io", "hemp_io",binning.nbins,binning.bin_min,binning.bin_max);
	for(int i=0;i<binning.nbins;++i) {
		htemp_io->SetBinContent(i+1, hp_io->GetBinContent(i+1));
	}
	
	hm_io->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
	hm_io->GetYaxis()->SetTitle("Posterior probability");
	hm_io->SetTitle("IO");
	hm_io->Draw("HIST");
	htemp_io->SetLineColor(kRed);
	htemp_io->Draw("HIST same");
	leg->Draw("same");
	c_io->SaveAs(Form("compare_credible_mach3_ptheta_%s_IO_%s_05July2022.pdf",osc_par,expt));
	//std::cout << " here " << std::endl;
	
	//Plot Dchi2
	Double_t bincontmax_no=hm_no->GetMaximum();
	Double_t bincontmax_io=hm_io->GetMaximum();
	
	std::cout << bincontmax_no << " " << bincontmax_io << " " << hm->GetNbinsX() << std::endl;
	
	TH1D* hmdchi2_no=(TH1D*)hm_no->Clone();
	TH1D* hmdchi2_io=(TH1D*)hm_io->Clone();
	
	for(int i=0;i<hm->GetNbinsX();++i) {
		//std::cout << i+1 << " " << hm_no->GetBinContent(i+1) << " " << -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no) << " " << -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io) << std::endl;
		if(hm_no->GetBinContent(i+1)!=0) hmdchi2_no->SetBinContent(i+1, -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no));
		if(hm_io->GetBinContent(i+1)!=0) hmdchi2_io->SetBinContent(i+1, -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io));
	}
	
	// NH
	TCanvas *c1=new TCanvas("c1","c1",700,700);
	
	TFile *fpdchi2 = new TFile(Form("contour_%s_100k_%s.root",osc_par,expt));
	TH1D *hpdchi2_no=(TH1D*)fpdchi2->Get(Form("SA_%s",oscparname));
	TH1D *hpdchi2_io=(TH1D*)fpdchi2->Get(Form("SA_%s_IH",oscparname));
	TH1D *hpdchi2_temp=(TH1D*)hpdchi2_io->Clone();
	
	TH1D *htemp1=(TH1D*)hmdchi2_io->Clone();

	gStyle->SetOptStat(0);
	hmdchi2_no->SetTitle("NO");
	hmdchi2_no->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
	hmdchi2_no->GetYaxis()->SetTitle("#Delta#chi^{2}");
	//hmdchi2_no->GetXaxis()->SetMaxDigits(2);
	hmdchi2_no->GetXaxis()->SetLabelSize(0.035);
	
	hpdchi2_no->SetLineColor(kRed);
	hmdchi2_no->Draw("HIST");
	hpdchi2_no->Draw("HIST same");
	
	leg->Draw("same");
	c1->SaveAs(Form("compare_dchi2_mach3_ptheta_%s_NO_%s_05July2022.pdf",osc_par,expt));
	
	//Plot Dchi2
	//IO
	TCanvas *c2=new TCanvas("c2","c2",700,700);

	int ix_min = hpdchi2_io->GetMinimumBin()-1;
	double locmin  = hpdchi2_io->GetBinContent(ix_min+1);

	for(int i=0;i<binning.nbins;++i) hpdchi2_io->SetBinContent(i+1, hpdchi2_temp->GetBinContent(i+1)-locmin);
	
	hmdchi2_io->SetTitle("IO");
	hmdchi2_io->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
	hmdchi2_io->GetYaxis()->SetTitle("#Delta#chi^{2}");
	//hmdchi2_io->GetXaxis()->SetMaxDigits(2);
	hmdchi2_io->GetXaxis()->SetLabelSize(0.035);
	
	hpdchi2_io->SetLineColor(kRed);
	hmdchi2_io->Draw("HIST");
	hpdchi2_io->Draw("HIST same");
	
	leg->Draw("same");
	c2->SaveAs(Form("compare_dchi2_mach3_ptheta_%s_IO_%s_05July2022.pdf",osc_par,expt));

}
