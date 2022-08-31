//To compare MaCh3 and P-theta sensitivities for t2k-nova joint analysis
#include <iostream>
#include <string.h>

struct Binning {
	Int_t nbins;
	Double_t bin_min;
	Double_t bin_max;
    const char* xtitle;
    const char* date;
};

const char* oscparname;

//Binning dcp_binning = {51,-3.2044245,3.2044245};
Binning dcp_binning = {51,-TMath::Pi(),TMath::Pi()};
Binning dm2_binning = {81,2.25e-3,2.75e-3};
Binning dm2_io_binning = {81,-2.75e-3,-2.25e-3};
Binning th23_binning = {81,0.38,0.62};
Binning binning;
Binning binning_io;

void Compare_MaCh3_Ptheta_Aria_MaCh3postBANFF_t2knova_02Aug2022(const char* expt, const char* osc_par, const char* outstring) {

	gStyle->SetOptStat(0);
	
	if(strstr(osc_par,"dCP") != NULL) {
		std::cout << " dcp " << std::endl;
		binning = dcp_binning;
		binning_io = dcp_binning;
		oscparname = "dcp";
        binning.xtitle = "#delta_{CP}";
        binning.date = "28June2022";
	}
	else if(strstr(osc_par,"dm2") != NULL) {
		std::cout << " dm2 " << std::endl;
		binning = dm2_binning;
		binning_io = dm2_io_binning;
		oscparname = "dm2";
        binning.xtitle = "#Delta m^{2}_{32}";
        binning.date = "30June2022";
	}
	else if(strstr(osc_par,"th23") != NULL) {
		std::cout << " th23 " << std::endl;
		binning = th23_binning;
		binning_io = th23_binning;
		oscparname = "sin223";
        binning.xtitle = "sin^{2}#theta_{23}";
        binning.date = "30June2022";
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
		//fm=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root");
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
	
    TH1D* hmdchi2_no = new TH1D("hmdchi2_no","hmdchi2_no", binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hmdchi2_io = new TH1D("hmdchi2_io","hmdchi2_io", binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    
	for(int i=0;i<nentries;++i) { //nentries
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
		//fa=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root");
        fa=new TFile("~/Downloads/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root");
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
    
    TH1D* hadchi2_no = new TH1D("hadchi2_no","hadchi2_no", binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hadchi2_io = new TH1D("hadchi2_io","hadchi2_io", binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    
	for(int i=0;i<anentries;++i) { //anentries
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
			
			ha->Fill(adcp, rc_weight*ssth23_weight);
			if(amh>0) ha_no->Fill(adcp, rc_weight*ssth23_weight); //Normal hierarchy
			else if (amh<0) ha_io->Fill(adcp, rc_weight*ssth23_weight); //ssth13_weight*ssth23_weight
		}
		else if(strstr(osc_par,"dm2") != NULL) {
			ha->Fill(adm2, rc_weight*ssth23_weight);
			if(amh>0) ha_no->Fill(adm2, rc_weight*ssth23_weight); //Normal hierarchy
			else if (amh<0) ha_io->Fill(adm2, rc_weight*ssth23_weight);
		}
		else if(strstr(osc_par,"th23") != NULL) {
			ha->Fill(ssth23, rc_weight*ssth23_weight);
			if(amh>0) ha_no->Fill(ssth23, rc_weight*ssth23_weight); //Normal hierarchy
			else if (amh<0) ha_io->Fill(ssth23, rc_weight*ssth23_weight);
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
	leg=new TLegend(0.2,0.65,0.35,0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	
	TLegend* leg1;
	//leg=new TLegend(0.6,0.67,0.75,0.85);
	leg1=new TLegend(0.5,0.65,0.75,0.85);
	leg1->SetBorderSize(0);
	leg1->SetTextSize(0.04);
    
    //************************//
    //== MaCh3 post-BANFF ==//
    
    TFile* fmbanff;
    if(strstr(expt,"T2KNOvA") != NULL) {
        std::cout << " MaCh3 T2KNOvA input file: " << std::endl;
        fmbanff=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/t2konly_banff_fit_job0-8_reweighted.root");
    }
    else {
        std::cout << " MaCh3 T2K only input file: " << std::endl;
        //fmbanff=new TFile("/cshare/vol2/users/kasetti/T2K-NOvA/t2konly_banff_fit_job0-8_reweighted.root"); //RC 2018
        //fmbanff=new TFile("~/Downloads/t2konly_banff_fit_job0-8_reweighted.root"); //RC 2018
        //fmbanff=new TFile("~/Downloads/t2konly_banff_job0-8_reweighted.root"); //RC 2019
        fmbanff=new TFile("t2konly_postBANFF_UPDATED_STEPSCALE_reweighted_BURN_IN_CUT.root");
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
    
    TH1D* hmbanff=new TH1D("hmbanff","hmbanff",           binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hmbanff_no =new TH1D("hmbanff_no","hmbanff_no", binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hmbanff_io =new TH1D("hmbanff_io","hmbanff_io", binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    
    TH1D* hmbanffdchi2_no = new TH1D("hmbanffdchi2_no","hmbanffdchi2_no", binning.nbins, binning.bin_min, binning.bin_max);
    TH1D* hmbanffdchi2_io = new TH1D("hmbanffdchi2_io","hmbanffdchi2_io", binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
    
    hmbanffdchi2_no->SetLineColor(kViolet);
    hmbanffdchi2_io->SetLineColor(kViolet);
    
    hmbanffdchi2_no->SetLineWidth(2);
    hmbanffdchi2_io->SetLineWidth(2);
    
    for(int i=0;i<nentries_mbanff;++i) { //nentries_mbanff
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
    hmbanff->Scale(1/hmbanff->Integral()); //normalize
    hmbanff_no->Scale(1/hmbanff_no->Integral()); //normalize
    hmbanff_io->Scale(1/hmbanff_io->Integral()); //normalize
    
    hmbanff->SetLineWidth(2);
    hmbanff_no->SetLineWidth(2);
    hmbanff_io->SetLineWidth(2);
    
    hmbanff->SetLineColor(kViolet);
    hmbanff_no->SetLineColor(kViolet);
    hmbanff_io->SetLineColor(kViolet);
	
	//************************//
	//== P-theta ==//
	
	//TFile *fp = TFile::Open("credible_T2KNOvA_both.root");
	//TFile *fp_both = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_%s/credible_%s_both.root",osc_par,binning.date,expt));
	//TFile *fp_both = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_Asimov0_t2knova_15July2022/credible_%s_both.root",osc_par,expt));
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
    
    TCanvas *c_both=new TCanvas("c_both","c_both",700,700);
    c_both->SetLeftMargin(0.14);
	
	leg->AddEntry(hm,"MaCh3","l");
	leg->AddEntry(htemp,"Ptheta","l");
	leg->AddEntry(ha,"Aria","l");
    leg->AddEntry(hmbanff,"MaCh3 post-BANFF","l");
	
	leg1->AddEntry(hm,"MaCh3","l");
	leg1->AddEntry(htemp,"Ptheta","l");
	leg1->AddEntry(ha,"Aria","l");
    leg1->AddEntry(hmbanff,"MaCh3 post-BANFF","l");
	
	hmbanff->GetXaxis()->SetTitle(binning.xtitle);
	hmbanff->GetYaxis()->SetTitle("Posterior probability");
	hmbanff->SetTitle("Both hierarchies");

	htemp->SetLineWidth(2);
	htemp->SetLineColor(kRed);
	ha->SetLineColor(kGreen);
	
    (strstr(osc_par,"th23") != NULL)?hmbanff->GetYaxis()->SetRangeUser(0,0.04):hmbanff->GetYaxis()->SetRangeUser(0,0.065);
    
    hmbanff->Draw("HIST");
	hm->Draw("HIST same");
	ha->Draw("HIST same");
	htemp->Draw("HIST same");
    (strstr(osc_par,"th23") != NULL)?leg->Draw("same"):leg1->Draw("same");
	//c_both->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_both_%s_07July2022.pdf",osc_par,expt));
	//c_both->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_both_%s_%s.pdf",osc_par,expt,outstring));
    c_both->SaveAs(Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_both_%s_%s.pdf",osc_par,expt,outstring));

	//Normal Ordering
	TCanvas *c_no=new TCanvas("c_no","c_no",700,700);
    c_no->SetLeftMargin(0.14);
	//TFile *fp_no = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_%s/credible_%s.root",osc_par,binning.date,expt));
    //TFile *fp_no = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_Asimov0_t2knova_15July2022/credible_%s.root",osc_par,expt));
    TFile *fp_no = TFile::Open(Form("credible_%s_%s.root",osc_par,expt));
	TH1D* hp_no=(TH1D*)fp_no->Get("LL");
	
	if(strstr(osc_par,"dCP") != NULL) hp_no->Rebin(51);
	else hp_no->Rebin(81);
	
	TH1D* htemp_no=new TH1D("htemp", "hemp",binning.nbins,binning.bin_min,binning.bin_max);
	for(int i=0;i<binning.nbins;++i) {
		htemp_no->SetBinContent(i+1, hp_no->GetBinContent(i+1));
	}
	
	hmbanff_no->GetXaxis()->SetTitle(binning.xtitle);
	hmbanff_no->GetYaxis()->SetTitle("Posterior probability");
	hmbanff_no->SetTitle("NO");
    
    (strstr(osc_par,"th23") != NULL)?hmbanff_no->GetYaxis()->SetRangeUser(0,0.04):hmbanff_no->GetYaxis()->SetRangeUser(0,0.06);
	
	htemp_no->SetLineWidth(2);
	
	htemp_no->SetLineColor(kRed);
	ha_no->SetLineColor(kGreen);
	hm_no->SetLineColor(kBlue);
	
    hmbanff_no->Draw("HIST");
	hm_no->Draw("HIST same");
	ha_no->Draw("HIST same");
	htemp_no->Draw("HIST same");
	//leg1->Draw("same");
    (strstr(osc_par,"th23") != NULL)?leg->Draw("same"):leg1->Draw("same");
	//c_no->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_NO_%s_07July2022.pdf",osc_par,expt));
	//c_no->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_NO_%s_%s.pdf",osc_par,expt,outstring));
    c_no->SaveAs(Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_NO_%s_%s.pdf",osc_par,expt,outstring));
	
	//Inverted Ordering
	TCanvas *c_io=new TCanvas("c_io","c_io",700,700);
    c_io->SetLeftMargin(0.14);
	//TFile *fp_io = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_%s/credible_%s_IH.root",osc_par,binning.date,expt));
    //TFile *fp_io = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_Asimov0_t2knova_15July2022/credible_%s_IH.root",osc_par,expt));
    TFile *fp_io = TFile::Open(Form("credible_%s_%s_IH.root",osc_par,expt));
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
	
	hmbanff_io->GetXaxis()->SetTitle(binning.xtitle);
	hmbanff_io->GetYaxis()->SetTitle("Posterior probability");
	hmbanff_io->SetTitle("IO");
    (strstr(osc_par,"th23") != NULL)?hmbanff_io->GetYaxis()->SetRangeUser(0,0.04):hmbanff_io->GetYaxis()->SetRangeUser(0,0.1);
	
    hmbanff_io->Draw("HIST");
	hm_io->Draw("HIST same");
	ha_io->Draw("HIST same");
	htemp_io->Draw("HIST same");
	//leg1->Draw("same");
    (strstr(osc_par,"th23") != NULL)?leg->Draw("same"):leg1->Draw("same");
	//c_io->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_IO_%s_07July2022.pdf",osc_par,expt));
	//c_io->SaveAs(Form("compare_credible_mach3_ptheta_aria_%s_IO_%s_%s.pdf",osc_par,expt,outstring));
    c_io->SaveAs(Form("compare_credible_mach3_ptheta_aria_mach3postbanff_%s_IO_%s_%s.pdf",osc_par,expt,outstring));
	//std::cout << " here " << std::endl;
	
	//Plot Dchi2
    //	Double_t bincontmax_no=hm_no->GetMaximum();
    //	Double_t bincontmax_io=hm_io->GetMaximum();
    //
    //  Double_t mbanffbincontmax_no=hmbanff_no->GetMaximum();
    //  Double_t mbanffbincontmax_io=hmbanff_io->GetMaximum();
    
    //	Double_t abincontmax_no=ha_no->GetMaximum();
    //	Double_t abincontmax_io=ha_io->GetMaximum();
    
    Double_t bincontmax_no=hm_no->GetBinContent(hm_no->GetMaximumBin());
    Double_t bincontmax_io=hm_io->GetBinContent(hm_io->GetMaximumBin());
    
    Double_t abincontmax_no=ha_no->GetBinContent(ha_no->GetMaximumBin());
    Double_t abincontmax_io=ha_io->GetBinContent(ha_io->GetMaximumBin());
    
    Double_t mbanffbincontmax_no=hmbanff_no->GetBinContent(hmbanff_no->GetMaximumBin());
    Double_t mbanffbincontmax_io=hmbanff_io->GetBinContent(hmbanff_no->GetMaximumBin());
	
	std::cout << bincontmax_no << " " << bincontmax_io << abincontmax_no << " " << abincontmax_io << " " << hm->GetNbinsX() << std::endl;
    
	
	for(int i=0;i<binning.nbins;++i) {
		//std::cout << i+1 << " " << hm_no->GetBinContent(i+1) << " " << -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no) << " " << -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io) << std::endl;
		if(hm_no->GetBinContent(i+1)!=0) hmdchi2_no->SetBinContent(i+1, -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no));
		if(hm_io->GetBinContent(i+1)!=0) hmdchi2_io->SetBinContent(i+1, -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io));
		
		if(ha_no->GetBinContent(i+1)!=0) hadchi2_no->SetBinContent(i+1, -2*TMath::Log(ha_no->GetBinContent(i+1)/abincontmax_no));
		if(ha_io->GetBinContent(i+1)!=0) hadchi2_io->SetBinContent(i+1, -2*TMath::Log(ha_io->GetBinContent(i+1)/abincontmax_io));
        
        if(hmbanff_no->GetBinContent(i+1)!=0) hmbanffdchi2_no->SetBinContent(i+1, -2*TMath::Log(hmbanff_no->GetBinContent(i+1)/mbanffbincontmax_no));
        if(hmbanff_io->GetBinContent(i+1)!=0) hmbanffdchi2_io->SetBinContent(i+1, -2*TMath::Log(hmbanff_io->GetBinContent(i+1)/mbanffbincontmax_io));
	}
	
	// NH
	TCanvas *c1=new TCanvas("c1","c1",700,700);
	
	//TFile *fpdchi2 = new TFile(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_%s/contour_%s_100k_%s.root",osc_par,binning.date,osc_par,expt));
    TFile *fpdchi2 = new TFile(Form("contour_%s_100k_%s.root",osc_par,expt));
	TH1D *hpdchi2_no=(TH1D*)fpdchi2->Get(Form("SA_%s",oscparname));
	TH1D *hpdchi2_io=(TH1D*)fpdchi2->Get(Form("SA_%s_IH",oscparname));
	TH1D *hpdchi2_temp=(TH1D*)hpdchi2_io->Clone();
	
	TH1D *htemp1=(TH1D*)hmdchi2_io->Clone();

	gStyle->SetOptStat(0);
	hadchi2_no->SetTitle("NO");
	hadchi2_no->GetXaxis()->SetTitle(binning.xtitle);
	hadchi2_no->GetYaxis()->SetTitle("#Delta#chi^{2}");
    (strstr(osc_par,"th23") != NULL)?hadchi2_no->GetYaxis()->SetRangeUser(0,30):hadchi2_no->GetYaxis()->SetRangeUser(0,15);
	//hadchi2_no->GetXaxis()->SetMaxDigits(2);
	hadchi2_no->GetXaxis()->SetLabelSize(0.035);

	hmdchi2_no->SetLineWidth(2);
	hpdchi2_no->SetLineWidth(2);
	hadchi2_no->SetLineWidth(2);
	
	hmdchi2_no->SetLineColor(kBlue);
	hpdchi2_no->SetLineColor(kRed);
	hadchi2_no->SetLineColor(kGreen);
	
	hadchi2_no->Draw("HIST");
	hmdchi2_no->Draw("HIST same");
    hmbanffdchi2_no->Draw("HIST same");
	hpdchi2_no->Draw("HIST same");
	
	leg->Draw("same");
	//c1->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_NO_%s_07July2022.pdf",osc_par,expt));
	//c1->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_NO_%s_%s.pdf",osc_par,expt,outstring));
    c1->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_mach3postbanff_%s_NO_%s_%s.pdf",osc_par,expt,outstring));
	
	//Plot Dchi2
	//IO
	TCanvas *c2=new TCanvas("c2","c2",700,700);

	int ix_min = hpdchi2_io->GetMinimumBin()-1;
	double locmin  = hpdchi2_io->GetBinContent(ix_min+1);

	for(int i=0;i<binning.nbins;++i) hpdchi2_io->SetBinContent(i+1, hpdchi2_temp->GetBinContent(i+1)-locmin);
	
	hadchi2_io->SetTitle("IO");
	hadchi2_io->GetXaxis()->SetTitle(binning.xtitle);
	hadchi2_io->GetYaxis()->SetTitle("#Delta#chi^{2}");
    (strstr(osc_par,"th23") != NULL)?hadchi2_io->GetYaxis()->SetRangeUser(0,30):hadchi2_io->GetYaxis()->SetRangeUser(0,30);
	//hadchi2_io->GetXaxis()->SetMaxDigits(2);
	hadchi2_io->GetXaxis()->SetLabelSize(0.035);
	
	hmdchi2_io->SetLineWidth(2);
	hpdchi2_io->SetLineWidth(2);
	hadchi2_io->SetLineWidth(2);
	
	hmdchi2_io->SetLineColor(kBlue);
	hpdchi2_io->SetLineColor(kRed);
	hadchi2_io->SetLineColor(kGreen);
	
	hadchi2_io->Draw("HIST");
	hmdchi2_io->Draw("HIST same");
    hmbanffdchi2_io->Draw("HIST same");
	hpdchi2_io->Draw("HIST same");
	
	leg->Draw("same");
	//c2->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_IO_%s_07July2022.pdf",osc_par,expt));
	//c2->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_%s_IO_%s_%s.pdf",osc_par,expt,outstring));
    c2->SaveAs(Form("compare_dchi2_mach3_ptheta_aria_mach3postbanff_%s_IO_%s_%s.pdf",osc_par,expt,outstring));

}
