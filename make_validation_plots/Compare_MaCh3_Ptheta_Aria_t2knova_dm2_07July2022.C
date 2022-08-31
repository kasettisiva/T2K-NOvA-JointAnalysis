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

void Compare_MaCh3_Ptheta_Aria_t2knova_dm2_07July2022(const char* expt, const char* osc_par, const char* outstring) {

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
		fm=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/MaCh3-t2knovaAsimov1-v03Sep2021_nuenumucorr_woRC_jobs4-5_red_reweighted.root");
	}
	else {
		std::cout << " MaCh3 T2K only input file: " << std::endl;
		//fm=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
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
	
	TH1D* hm=new TH1D("hm","hm",           binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* hm_no =new TH1D("hm_no","hm_no", binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* hm_io =new TH1D("hm_io","hm_io", binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
	
	for(int i=0;i<100000;++i) { //nentries
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
	
	hm_no->SetLineColor(kBlue);
	hm_io->SetLineColor(kBlue);
	
	TLegend* leg;
	//leg=new TLegend(0.6,0.67,0.75,0.85);
	leg=new TLegend(0.25,0.7,0.4,0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	
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
        fa=new TFile("~/Downloads/mcmcsamples_t2kcafana_joint_both_systs_noreactor_nddata_thinned1_postburnin20000.root"); //asimov0
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
	
	for(int i=0;i<100000;++i) { //anentries
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
	ha->Scale(1/ha->Integral()); //normalize
	ha_no->Scale(1/ha_no->Integral()); //normalize
	ha_io->Scale(1/ha_io->Integral()); //normalize

	ha->SetLineWidth(2);
	ha_no->SetLineWidth(2);
	ha_io->SetLineWidth(2);
	
	ha_no->SetLineColor(kGreen);
	ha_io->SetLineColor(kGreen);
	
	TCanvas *c=new TCanvas("c","c",700,400);
	//c->Divide(2,1);
	//hm->Draw("HIST");
	
	//TFile *fp = TFile::Open("credible_T2KNOvA_both.root");
	//TFile *fp_no = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_30June2022/credible_%s.root",osc_par,expt));
    TFile *fp_no = TFile::Open(Form("credible_%s_%s_asimov0.root",osc_par,expt));
	TH1D* hp_no=(TH1D*)fp_no->Get("LL");
	hp_no->SetLineColor(kRed);

	leg->AddEntry(hm,"MaCh3","l");
	leg->AddEntry(hp_no,"Ptheta","l");
	leg->AddEntry(ha_no,"Aria","l");
	
	if(strstr(osc_par,"dCP") != NULL) hp_no->Rebin(51);
	else hp_no->Rebin(81);
	
	TH1D* htemp_no=new TH1D("htemp_no", "hemp_no",binning.nbins,binning.bin_min,binning.bin_max);
	for(int i=0;i<binning.nbins;++i) {
		htemp_no->SetBinContent(i+1, hp_no->GetBinContent(i+1));
	}
	
	//TFile *fp_io = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_30June2022/credible_%s_IH.root",osc_par,expt));
    TFile *fp_io = TFile::Open(Form("credible_%s_%s_IH_asimov0.root",osc_par,expt));
	TH1D* hp_io=(TH1D*)fp_io->Get("LL");

	if(strstr(osc_par,"dCP") != NULL) hp_io->Rebin(51);
	else hp_io->Rebin(81);
	
	TSpline3 *sp_io=new TSpline3(hp_io);
	
	TH1D* htemp_io=new TH1D("htemp_io", "hemp_io",binning_io.nbins,binning_io.bin_min,binning_io.bin_max);
	for(int i=0;i<binning.nbins;++i) {
		//htemp_io->SetBinContent(i+1, hp_io->GetBinContent(i+1));
		double dm32 = hm_io->GetBinCenter(i+1);
		double dm21 = 7.53e-5;
		double dm31 = dm32 + dm21;
		double abs_dm31 = fabs(dm31);
		//std::cout << " abs_dm31 " << abs_dm31 << " sp_io->Eval(abs_dm31) " << sp_io->Eval(abs_dm31) << std::endl;
		htemp_io->SetBinContent(i+1, sp_io->Eval(abs_dm31));
	}
	
	//TH1D* htemp=(TH1D*)hp->Clone();
	//c->cd(1);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 0.5, 1.0);
	//pad1->SetBottomMargin(0.05);
	//pad1->SetRightMargin(0.05);
	pad1->SetLeftMargin(0.1919771);
	pad1->SetRightMargin(0.005730659);
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();
	
	hm_no->GetXaxis()->SetMaxDigits(2);
	hm_no->SetTitle("NO");
	hm_no->GetYaxis()->SetTitle("Posterior probability");
	//hm_no->GetXaxis()->SetTitle("#Delta m^{2}_{32}");
	hm_no->GetYaxis()->SetTitleSize(0.05);
	hm_no->GetYaxis()->SetLabelSize(0.05);
	hm_no->GetXaxis()->SetLabelSize(0.05);
	
	htemp_no->SetLineColor(kRed);
	
	htemp_no->SetLineWidth(2);
	hm_no->SetLineWidth(2);
	ha_no->SetLineWidth(2);
	
	hm_no->Draw("HIST");
	ha_no->Draw("HIST same");
	htemp_no->Draw("HIST same");
	leg->Draw("same");
	
	//c->cd(2);
	c->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0.5,0,1,1);
	//TLine* line=new TLine(hold[i][j]->GetXaxis()->GetXmin(),1,hold[i][j]->GetXaxis()->GetXmax(),1);
	//line->SetLineColor(kRed);
	//pad2->SetBottomMargin(0.05);
	pad2->SetLeftMargin(0.005730659);
	pad2->SetRightMargin(0.1919771);
	pad2->SetTickx(1);
	//pad2->SetTopMargin(0.05649718);
	//pad2->SetLeftMargin(0.1638418);
	pad2->SetTicky();
	pad2->Draw();
	pad2->cd();
		
	hm_io->SetTitle("IO");
	hm_io->GetXaxis()->SetTitle("#Delta m^{2}_{32}");
	hm_io->GetXaxis()->SetTitleSize(0.05);
	hm_io->GetXaxis()->SetTitleOffset(0.90);
	hm_io->GetYaxis()->SetLabelSize(0);
	hm_io->GetXaxis()->SetLabelSize(0.05);
	hm_io->GetXaxis()->SetMaxDigits(2);
	
	htemp_io->SetLineColor(kRed);
	
	htemp_io->SetLineWidth(2);
	hm_io->SetLineWidth(2);
	ha_io->SetLineWidth(2);
	
	hm_io->Draw("HIST");
	ha_io->Draw("HIST same");
	htemp_io->Draw("HIST same");
	
	c->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_both_%s_%s.pdf",osc_par,expt,outstring));

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
	
	// NO
	TCanvas *c1=new TCanvas("c1","c1",700,700);
	
	//TFile *fpdchi2 = new TFile(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_30June2022/contour_%s_100k_%s.root",osc_par,osc_par,expt));
    TFile *fpdchi2 = new TFile(Form("contour_%s_100k_%s_asimov0.root",osc_par,expt));
	TH1D *hpdchi2_no=(TH1D*)fpdchi2->Get(Form("SA_%s",oscparname));
	TH1D *hpdchi2_io=(TH1D*)fpdchi2->Get(Form("SA_%s_IH",oscparname));
	TH1D *hpdchi2_temp=(TH1D*)hpdchi2_io->Clone();
	
	
	TH1D *htemp1=(TH1D*)hmdchi2_io->Clone();

	hadchi2_no->SetTitle("NO");
	hadchi2_no->GetYaxis()->SetTitle("#Delta#chi^{2}");
	hadchi2_no->GetXaxis()->SetTitle("#Delta m^{2}_{32}");
	hadchi2_no->GetXaxis()->SetTitleSize(0.04);
	hadchi2_no->GetXaxis()->SetTitleOffset(1.14);
	hadchi2_no->GetYaxis()->SetTitleSize(0.04);
	hadchi2_no->GetXaxis()->SetMaxDigits(2);
	hadchi2_no->GetXaxis()->SetLabelSize(0.035);
	hadchi2_no->GetYaxis()->SetLabelSize(0.035);
	
	hadchi2_no->GetYaxis()->SetRangeUser(0,40);
	
	hpdchi2_no->SetLineColor(kRed);
	hadchi2_no->SetLineColor(kGreen);
	hmdchi2_no->SetLineColor(kBlue);
	
	hpdchi2_no->SetLineWidth(2);
	hadchi2_no->SetLineWidth(2);
	hmdchi2_no->SetLineWidth(2);
	
	hadchi2_no->Draw("HIST");
	hmdchi2_no->Draw("HIST same");
	hpdchi2_no->Draw("HIST same");
	
	leg->Draw("same");
	
	//c1->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_NO_%s_07July2022.pdf",osc_par,expt));
    c1->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_NO_%s_%s.pdf",osc_par,expt,outstring));
	
	/////////////////////////////
	//Plot Dchi2
	//IO
	TCanvas *c2=new TCanvas("c2","c2",700,700);

	int ix_min = hpdchi2_io->GetMinimumBin()-1;
	double locmin  = hpdchi2_io->GetBinContent(ix_min+1);

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
	
	hadchi2_io->SetTitle("IO");
	hadchi2_io->GetYaxis()->SetTitle("#Delta#chi^{2}");
	hadchi2_io->GetXaxis()->SetTitle("#Delta m^{2}_{32}");
	hadchi2_io->GetXaxis()->SetTitleSize(0.04);
	hadchi2_io->GetXaxis()->SetTitleOffset(1.14);
	hadchi2_io->GetXaxis()->SetMaxDigits(2);
	hadchi2_io->GetYaxis()->SetLabelSize(0.035);
	hadchi2_io->GetXaxis()->SetLabelSize(0.035);
	//hpdchi2_io->SetLineColor(kRed);
	
	hadchi2_io->GetYaxis()->SetRangeUser(0,40);
	
	htemp1->SetLineColor(kRed);
	hadchi2_io->SetLineColor(kGreen);
	hmdchi2_io->SetLineColor(kBlue);
	
	htemp1->SetLineWidth(2);
	hadchi2_io->SetLineWidth(2);
	hmdchi2_io->SetLineWidth(2);
	
	hadchi2_io->Draw("HIST");
	hmdchi2_io->Draw("HIST same");
	htemp1->Draw("HIST same");
	
	leg->Draw("same");
	
	//c2->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_IO_%s_07July2022.pdf",osc_par,expt));
    c2->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_IO_%s_%s.pdf",osc_par,expt,outstring));

}
