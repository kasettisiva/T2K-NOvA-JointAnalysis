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

void Compare_MaCh3_Ptheta_Aria_MaCh3postBANFF_t2knova_th23_17July2022(const char* expt, const char* osc_par) {

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
	//************************//
	//== MaCh3 ==//
	
	//TFile *fm=new TFile("MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
	//TFile *fm=new TFile("MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
	TFile* fm;
	if(strstr(expt,"T2KNOvA") != NULL) {
		std::cout << " MaCh3 T2KNOvA input file: " << std::endl;
		fm=new TFile("MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
	}
	else {
		std::cout << " MaCh3 T2K only input file: " << std::endl;
		fm=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
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
	
	hm->SetLineWidth(2);
	hm_no->SetLineWidth(2);
	hm_io->SetLineWidth(2);
	
	//************************//
	//== Aria ==//
	
	TFile* fa;
	if(strstr(expt,"T2KNOvA") != NULL) {
		std::cout << " Aria T2KNOvA input file: " << std::endl;
		fa=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root");
	}
	else {
		std::cout << " Aria T2K only input file: " << std::endl;
		fa=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root");
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
	
	TH1D* ha=new TH1D("ha","ha",           binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* ha_no =new TH1D("ha_no","ha_no", binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* ha_io =new TH1D("ha_io","ha_io", binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
	
	for(int i=0;i<anentries;++i) {
		tra->GetEntry(i);
		
		Double_t ssth13 = TMath::Power(TMath::Sin(ath13),2);
		Double_t ssth23 = TMath::Power(TMath::Sin(ath23),2);
		double ssth23_weight = 2.*TMath::Sqrt((1-ssth23)*ssth23); // To go from flat prior in th23 -> flat prior in sin^2(th23)
		double ssth13_weight = 2.*TMath::Sqrt((1-ssth13)*ssth13); // To go from flat prior in th13 -> flat prior in sin^2(th13)
		Double_t rc_weight = TMath::Gaus(ssth13, 2.18e-2, 0.07e-2);
		
		if(strstr(osc_par,"dCP") != NULL) {
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
	ha->Scale(1/ha->Integral()); //normalize
	ha_no->Scale(1/ha_no->Integral()); //normalize
	ha_io->Scale(1/ha_io->Integral()); //normalize

	ha->SetLineWidth(2);
	ha_no->SetLineWidth(2);
	ha_io->SetLineWidth(2);
	
	TLegend* leg;
	//leg=new TLegend(0.6,0.67,0.75,0.85);
	leg=new TLegend(0.2,0.7,0.35,0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	
	TCanvas *c_both=new TCanvas("c_both","c_both",700,700);
	
	//************************//
	//== P-theta ==//
	
	//TFile *fp = TFile::Open("credible_T2KNOvA_both.root");
	TFile *fp_both = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_30June2022/credible_%s_both.root",osc_par,expt));
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
	leg->AddEntry(ha,"Aria","l");
	
	hm->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
	hm->GetYaxis()->SetTitle("Posterior probability");
	hm->SetTitle("Both hierarchies");

	htemp->SetLineWidth(2);
	htemp->SetLineColor(kRed);
	ha->SetLineColor(kGreen);
	
	hm->Draw("HIST");
	ha->Draw("HIST same");
	htemp->Draw("HIST same");
	leg->Draw("same");
	c_both->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_both_%s_07July2022.pdf",osc_par,expt));

	//Normal Ordering
	TCanvas *c_no=new TCanvas("c_no","c_no",700,700);
	TFile *fp_no = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_30June2022/credible_%s.root",osc_par,expt));
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
	
	htemp_no->SetLineWidth(2);
	
	htemp_no->SetLineColor(kRed);
	ha_no->SetLineColor(kGreen);
	hm_no->SetLineColor(kBlue);
	
	hm_no->Draw("HIST");
	ha_no->Draw("HIST same");
	htemp_no->Draw("HIST same");
	leg->Draw("same");
	c_no->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_NO_%s_07July2022.pdf",osc_par,expt));
	
	//Inverted Ordering
	TCanvas *c_io=new TCanvas("c_io","c_io",700,700);
	TFile *fp_io = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_30June2022/credible_%s_IH.root",osc_par,expt));
	TH1D* hp_io=(TH1D*)fp_io->Get("LL");
	
	if(strstr(osc_par,"dCP") != NULL) hp_io->Rebin(51);
	else hp_io->Rebin(81);
	
	TH1D* htemp_io=new TH1D("htemp_io", "hemp_io",binning.nbins,binning.bin_min,binning.bin_max);
	for(int i=0;i<binning.nbins;++i) {
		htemp_io->SetBinContent(i+1, hp_io->GetBinContent(i+1));
	}
	
	htemp_io->SetLineWidth(2);
	
	htemp_io->SetLineColor(kRed);
	ha_io->SetLineColor(kGreen);
	hm_io->SetLineColor(kBlue);
	
	hm_io->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
	hm_io->GetYaxis()->SetTitle("Posterior probability");
	hm_io->SetTitle("IO");
	
	hm_io->Draw("HIST");
	ha_io->Draw("HIST same");
	htemp_io->Draw("HIST same");
	leg->Draw("same");
	c_io->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_IO_%s_07July2022.pdf",osc_par,expt));
	//std::cout << " here " << std::endl;
	
	//Plot Dchi2
	Double_t bincontmax_no=hm_no->GetMaximum();
	Double_t bincontmax_io=hm_io->GetMaximum();
	
	Double_t abincontmax_no=ha_no->GetMaximum();
	Double_t abincontmax_io=ha_io->GetMaximum();
	
	std::cout << bincontmax_no << " " << bincontmax_io << abincontmax_no << " " << abincontmax_io << " " << hm->GetNbinsX() << std::endl;
	
	TH1D* hmdchi2_no=(TH1D*)hm_no->Clone();
	TH1D* hmdchi2_io=(TH1D*)hm_io->Clone();
	
	TH1D* hadchi2_no=(TH1D*)ha_no->Clone();
	TH1D* hadchi2_io=(TH1D*)ha_io->Clone();
	
	for(int i=0;i<hm->GetNbinsX();++i) {
		//std::cout << i+1 << " " << hm_no->GetBinContent(i+1) << " " << -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no) << " " << -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io) << std::endl;
		if(hm_no->GetBinContent(i+1)!=0) hmdchi2_no->SetBinContent(i+1, -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no));
		if(hm_io->GetBinContent(i+1)!=0) hmdchi2_io->SetBinContent(i+1, -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io));
		
		if(ha_no->GetBinContent(i+1)!=0) hadchi2_no->SetBinContent(i+1, -2*TMath::Log(ha_no->GetBinContent(i+1)/abincontmax_no));
		if(ha_io->GetBinContent(i+1)!=0) hadchi2_io->SetBinContent(i+1, -2*TMath::Log(ha_io->GetBinContent(i+1)/abincontmax_io));
	}
	
	// NH
	TCanvas *c1=new TCanvas("c1","c1",700,700);
	
	TFile *fpdchi2 = new TFile(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_30June2022/contour_%s_100k_%s.root",osc_par,osc_par,expt));
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

	hmdchi2_no->SetLineWidth(2);
	hpdchi2_no->SetLineWidth(2);
	hadchi2_no->SetLineWidth(2);
	
	hmdchi2_no->SetLineColor(kBlue);
	hpdchi2_no->SetLineColor(kRed);
	hadchi2_no->SetLineColor(kGreen);
	
	hmdchi2_no->Draw("HIST");
	hadchi2_no->Draw("HIST same");
	hpdchi2_no->Draw("HIST same");
	
	leg->Draw("same");
	c1->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_NO_%s_07July2022.pdf",osc_par,expt));
	
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
	
	hmdchi2_io->SetLineWidth(2);
	hpdchi2_io->SetLineWidth(2);
	hadchi2_io->SetLineWidth(2);
	
	hmdchi2_io->SetLineColor(kBlue);
	hpdchi2_io->SetLineColor(kRed);
	hadchi2_io->SetLineColor(kGreen);
	
	hmdchi2_io->Draw("HIST");
	hadchi2_io->Draw("HIST same");
	hpdchi2_io->Draw("HIST same");
	
	leg->Draw("same");
	c2->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_IO_%s_07July2022.pdf",osc_par,expt));

}
