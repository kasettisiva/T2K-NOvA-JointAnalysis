//To compare MaCh3 and P-theta sensitivities for t2k-nova joint analysis
#include <iostream>
#include <string.h>

struct Binning {
	Int_t nbins;
	Double_t bin_min;
	Double_t bin_max;
    const char* xtitle;
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
    h->GetYaxis()->SetNdivisions(309);
    std::cout << 0.8*h->GetBinContent(h->GetMinimumBin()) << " " << 1.2*h->GetBinContent(h->GetMaximumBin()) << std::endl;
    h->GetYaxis()->SetRangeUser(0.7*h->GetBinContent(h->GetMinimumBin()),std::min(3.0, 1.3*h->GetBinContent(h->GetMaximumBin())));
}

void SetTopPadStyle(TPad* pad) {
    //pad->SetBottomMargin(0.01);
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

void GetRatioHist(TH1* hp, TH1* hm, TH1* ha, TLegend *leg, Binning binning, bool isDchi2, const char* output_name) {
    TCanvas *c=new TCanvas("","",700,700);
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.346, 1, 1.0);
    SetTopPadStyle(pad1);
    
    Double_t bm = hm->GetBinContent(hm->GetMaximumBin());
    Double_t bp = hp->GetBinContent(hp->GetMaximumBin());
    Double_t ba = ha->GetBinContent(ha->GetMaximumBin());
    
    std::vector<Double_t> cont_vec={bm, bp, ba};
    Double_t bcmax = *std::max_element(cont_vec.begin(), cont_vec.end());
    std::cout << bm << " " << bp << " " << ba << " " << bcmax << std::endl;
    
    if(isDchi2)  hm->GetYaxis()->SetRangeUser(0, std::min(10.0, 1.2*hm->GetBinContent(hm->GetMaximumBin())));
    else hm->GetYaxis()->SetRangeUser(std::max(0.0, 0.8*hm->GetBinContent(hm->GetMinimumBin())), 1.2*bcmax);
    
    SetHistStyle(hm);
    
    hm->Draw("HIST");
    ha->Draw("HIST same");
    hp->Draw("HIST same");
    
    leg->Draw("same");
    
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.358);
    SetBottomPadStyle(pad2);
    
    TH1D* hratio_ma = (TH1D*)hm->Clone();
    TH1D* hratio_pa = (TH1D*)hp->Clone();
    
    hratio_ma->SetMarkerStyle(20);
    hratio_ma->SetMarkerColor(kBlue);
    hratio_pa->SetMarkerStyle(20);
    hratio_pa->SetMarkerColor(kRed);
    
    hratio_ma->Divide(ha);
    hratio_pa->Divide(ha);
    
    SetRatioStyle(hratio_ma, binning.xtitle);
    
    hratio_ma->Draw("phist");
    hratio_pa->Draw("phistsame");
    
    TLine *line = new TLine(binning.bin_min,1,binning.bin_max,1);
    line->Draw("same");
    
    c->SaveAs(output_name);
}

const char* oscparname;

Binning dcp_binning = {51,-TMath::Pi(),TMath::Pi()}; //-3.2044245,3.2044245
Binning dm2_binning = {81,2.25e-3,2.75e-3};
Binning dm2_io_binning = {81,-2.75e-3,-2.25e-3};
Binning th23_binning = {81,0.38,0.62};
Binning binning;
Binning binning_io;

Int_t kNOvAColor=TColor::GetColor("#044c94");
Int_t kNOvAColor1Sigma=TColor::GetColor("#3670a9");
Int_t kNOvAColor2Sigma=TColor::GetColor("#82a6ca");
Int_t kNOvAColor3Sigma=TColor::GetColor("#cddbea");

Int_t kT2KColor=TColor::GetColor("#840404");
Int_t kT2KColor1Sigma=TColor::GetColor("#9d3636");
Int_t kT2KColor2Sigma=TColor::GetColor("#b56868");
Int_t kT2KColor3Sigma=TColor::GetColor("#dab4b4");

Int_t kJointColor=TColor::GetColor("#44284C");

void Compare_MaCh3_Ptheta_Aria_t2knova_11Aug2022(const char* expt, const char* osc_par, const char* outstring) {

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
		oscparname = "dm2";
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
        //fm=new TFile("MaCh3_t2konly_Asimov1_03Sep2021_job3_ch0-999_red_reweighted.root"); //Asimov1
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
	
	TH1D* hm =new TH1D("hm",      Form("Both Hierarchies;%s;Posterior probability", binning.xtitle),
                      binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* hm_no =new TH1D("hm_no",Form("Normal Hierarchy;%s;Posterior probability", binning.xtitle),
                          binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* hm_io =new TH1D("hm_io",Form("Inverted Hierarchy;%s;Posterior probability", binning.xtitle),
                          binning_io.nbins, binning_io.bin_min, binning_io.bin_max);

    TH1D* hmdchi2_no=(TH1D*)hm_no->Clone();
    TH1D* hmdchi2_io=(TH1D*)hm_io->Clone();
    
    
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
    hmdchi2_no->SetLineWidth(2);
    hmdchi2_io->SetLineWidth(2);
    
    hm->SetLineColor(kBlue);
    hm_no->SetLineColor(kBlue);
    hm_io->SetLineColor(kBlue);
    hmdchi2_no->SetLineColor(kBlue);
    hmdchi2_io->SetLineColor(kBlue);
    
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
        //fa=new TFile("~/Downloads/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root"); //asimov1
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
	
	TH1D* ha=new TH1D("ha",       Form("Both Hierarchies;%s;Posterior probability", binning.xtitle),
                      binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* ha_no =new TH1D("ha_no",Form("Normal Hierarchy;%s;Posterior probability", binning.xtitle),
                          binning.nbins, binning.bin_min, binning.bin_max);
	TH1D* ha_io =new TH1D("ha_io",Form("Inverted Hierarchy;%s;Posterior probability", binning.xtitle),
                          binning_io.nbins, binning_io.bin_min, binning_io.bin_max);
	
    TH1D* hadchi2_no=(TH1D*)ha_no->Clone();
    TH1D* hadchi2_io=(TH1D*)ha_io->Clone();
    
    ha->SetLineWidth(2);
    ha_no->SetLineWidth(2);
    ha_io->SetLineWidth(2);
    hadchi2_no->SetLineWidth(2);
    hadchi2_io->SetLineWidth(2);
    
    ha->SetLineColor(kGreen);
    ha_no->SetLineColor(kGreen);
    ha_io->SetLineColor(kGreen);
    hadchi2_no->SetLineColor(kGreen);
    hadchi2_io->SetLineColor(kGreen);
    
    
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
    
	//************************//
	//== P-theta ==//
	
	//TFile *fp = TFile::Open("credible_T2KNOvA_both.root");
	//TFile *fp_both = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_28June2022/credible_%s_both.root",osc_par,expt));
	//TFile *fp_both = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_Asimov0_t2knova_15July2022/credible_%s_both.root",osc_par,expt));
    TFile *fp_both = TFile::Open(Form("credible_%s_%s_both_asimov0.root",osc_par,expt));
	TH1D* hp=(TH1D*)fp_both->Get("LL_both");
	//hp->SetLineColor(kOrange);
	//hp->Draw("same");
	if(strstr(osc_par,"dCP") != NULL) hp->Rebin(51);
	else hp->Rebin(81);
	
	TH1D* htemp=new TH1D("htemp", Form("Both Hierarchies;%s;Posterior probability", binning.xtitle),
                         binning.nbins,binning.bin_min,binning.bin_max);
    
	for(int i=0;i<binning.nbins;++i) {
		htemp->SetBinContent(i+1, hp->GetBinContent(i+1));
	}
    
    htemp->SetLineWidth(2);
    htemp->SetLineColor(kRed);
    
    TLegend* leg;
    //leg=new TLegend(0.6,0.67,0.75,0.85);
    leg=new TLegend(0.2,0.67,0.35,0.85);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.045);
    
    leg->AddEntry(hm,"MaCh3","l");
    leg->AddEntry(htemp,"Ptheta","l");
    leg->AddEntry(ha,"Aria","l");
    
    TLegend* leg1;
    //leg=new TLegend(0.6,0.67,0.75,0.85);
    leg1=new TLegend(0.65,0.67,0.82,0.85);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.045);
	
	leg1->AddEntry(hm,"MaCh3","l");
	leg1->AddEntry(htemp,"Ptheta","l");
	leg1->AddEntry(ha,"Aria","l");
    
	//hm->GetXaxis()->SetTitle(binning.xtitle);
	//hm->GetYaxis()->SetTitle("Posterior probability");
	//hm->SetTitle("Both hierarchies");
    //hm->GetXaxis()->SetLabelSize(0);
    //(strstr(osc_par,"th23") != NULL)?hm->GetYaxis()->SetRangeUser(0.008,0.04):hm->GetYaxis()->SetRangeUser(0.008,0.04);
    
    GetRatioHist(htemp, hm, ha, leg, binning, false, Form("compare_credible_mach3_ptheta_aria_%s_both_%s_%s.pdf",osc_par,expt,outstring));
    
    //Get credible intervals for p-theta
	//Normal Ordering
	//TFile *fp_no = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_28June2022/credible_%s.root",osc_par,expt));
    //TFile *fp_no = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_Asimov0_t2knova_15July2022/credible_%s.root",osc_par,expt));
    TFile *fp_no = TFile::Open(Form("credible_%s_%s_asimov0.root",osc_par,expt));
	TH1D* hp_no=(TH1D*)fp_no->Get("LL");
	
	if(strstr(osc_par,"dCP") != NULL) hp_no->Rebin(51);
	else hp_no->Rebin(81);
	
	TH1D* htemp_no=new TH1D("htemp_no", Form("Normal Hierarchy;%s;Posterior probability", binning.xtitle),
                            binning.nbins,binning.bin_min,binning.bin_max);
    
	for(int i=0;i<binning.nbins;++i) {
		htemp_no->SetBinContent(i+1, hp_no->GetBinContent(i+1));
	}
	
	//hm_no->GetXaxis()->SetTitle(binning.xtitle);
	//hm_no->GetYaxis()->SetTitle("Posterior probability");
	//hm_no->SetTitle("NO");
    //(strstr(osc_par,"th23") != NULL)?hm_no->GetYaxis()->SetRangeUser(0.002,0.035):hm_no->GetYaxis()->SetRangeUser(0.002,0.04);
	
    htemp_no->SetLineColor(kRed);
	htemp_no->SetLineWidth(2);
	
	ha_no->SetLineColor(kGreen);
	hm_no->SetLineColor(kBlue);
 
    ha_no->SetLineWidth(2);
    hm_no->SetLineWidth(2);
    
    GetRatioHist(htemp_no, hm_no, ha_no, leg, binning, false, Form("compare_credible_mach3_ptheta_aria_%s_NO_%s_%s.pdf",osc_par,expt,outstring));
 
	//Inverted Ordering
	//TFile *fp_io = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_28June2022/credible_%s_IH.root",osc_par,expt));
    //TFile *fp_io = TFile::Open(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_Asimov0_t2knova_15July2022/credible_%s_IH.root",osc_par,expt));
    TFile *fp_io = TFile::Open(Form("credible_%s_%s_IH_asimov0.root",osc_par,expt));
	TH1D* hp_io=(TH1D*)fp_io->Get("LL");
	
	if(strstr(osc_par,"dCP") != NULL) hp_io->Rebin(51);
	else hp_io->Rebin(81);
	
	TH1D* htemp_io=new TH1D("htemp_io", Form("Inverted Hierarchy;%s;Posterior probability", binning.xtitle),
                            binning.nbins,binning.bin_min,binning.bin_max);
    
	for(int i=0;i<binning.nbins;++i) {
		htemp_io->SetBinContent(i+1, hp_io->GetBinContent(i+1));
	}

    htemp_io->SetLineColor(kRed);
	ha_io->SetLineColor(kGreen);
	hm_io->SetLineColor(kBlue);
    
    htemp_io->SetLineWidth(2);
    ha_io->SetLineWidth(2);
    hm_io->SetLineWidth(2);
	
	//hm_io->GetXaxis()->SetTitle(binning.xtitle);
	//hm_io->GetYaxis()->SetTitle("Posterior probability");
	//hm_io->SetTitle("IO");
    //(strstr(osc_par,"th23") != NULL)?hm_io->GetYaxis()->SetRangeUser(0.002,0.035):hm_io->GetYaxis()->SetRangeUser(0.002,0.04);
 
    GetRatioHist(htemp_io, hm_io, ha_io, (strstr(osc_par,"th23") != NULL)?leg:leg1, binning, false, Form("compare_credible_mach3_ptheta_aria_%s_IO_%s_%s.pdf",osc_par,expt,outstring));
    

	//Plot Dchi2
    
    Double_t bincontmax_no=hm_no->GetBinContent(hm_no->GetMaximumBin());
    Double_t bincontmax_io=hm_io->GetBinContent(hm_io->GetMaximumBin());
    
    Double_t abincontmax_no=ha_no->GetBinContent(ha_no->GetMaximumBin());
    Double_t abincontmax_io=ha_io->GetBinContent(ha_io->GetMaximumBin());
	
	std::cout << bincontmax_no << " " << bincontmax_io << " " << abincontmax_no << " " << abincontmax_io << " " << hm->GetNbinsX() << std::endl;
	
	for(int i=0;i<hm->GetNbinsX();++i) {
		//std::cout << i+1 << " " << hm_no->GetBinCenter(i+1) << " MaCh3 " << hm_no->GetBinContent(i+1) << " " << -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no) << " " << hm_io->GetBinContent(i+1) << " " << -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io)
        //<< " Aria " << ha_no->GetBinContent(i+1) << " " << -2*TMath::Log(ha_no->GetBinContent(i+1)/abincontmax_no) << " " << ha_io->GetBinContent(i+1) << " " << -2*TMath::Log(ha_io->GetBinContent(i+1)/abincontmax_io) << std::endl;
        
		if(hm_no->GetBinContent(i+1)!=0) hmdchi2_no->SetBinContent(i+1, -2*TMath::Log(hm_no->GetBinContent(i+1)/bincontmax_no));
		if(hm_io->GetBinContent(i+1)!=0) hmdchi2_io->SetBinContent(i+1, -2*TMath::Log(hm_io->GetBinContent(i+1)/bincontmax_io));
		
		if(ha_no->GetBinContent(i+1)!=0) hadchi2_no->SetBinContent(i+1, -2*TMath::Log(ha_no->GetBinContent(i+1)/abincontmax_no));
		if(ha_io->GetBinContent(i+1)!=0) hadchi2_io->SetBinContent(i+1, -2*TMath::Log(ha_io->GetBinContent(i+1)/abincontmax_io));
	}
	
	// NH
	
	//TFile *fpdchi2 = new TFile(Form("/home/kasetti/T2K_on_grid/job_outputs/margTemp_%s_wRC_t2knova_28June2022/contour_%s_100k_%s.root",osc_par,osc_par,expt));
    TFile *fpdchi2 = new TFile(Form("contour_%s_100k_%s_asimov0.root",osc_par,expt));
	TH1D *hpdchi2_no=(TH1D*)fpdchi2->Get(Form("SA_%s",oscparname));
	TH1D *hpdchi2_io=(TH1D*)fpdchi2->Get(Form("SA_%s_IH",oscparname));
	TH1D *hpdchi2_temp=(TH1D*)hpdchi2_io->Clone();
	
	TH1D *htemp1=(TH1D*)hmdchi2_io->Clone();

	gStyle->SetOptStat(0);
	//hmdchi2_no->SetTitle("NO");
	//hmdchi2_no->GetXaxis()->SetTitle(binning.xtitle);
	hmdchi2_no->GetYaxis()->SetTitle("#Delta#chi^{2}");
    (strstr(osc_par,"th23") != NULL)?hmdchi2_no->GetYaxis()->SetRangeUser(0,25):hmdchi2_no->GetYaxis()->SetRangeUser(0,7);
	//hadchi2_no->GetXaxis()->SetMaxDigits(2);
	hmdchi2_no->GetXaxis()->SetLabelSize(0.035);

	hmdchi2_no->SetLineWidth(2);
	hpdchi2_no->SetLineWidth(2);
	hadchi2_no->SetLineWidth(2);
	
	hmdchi2_no->SetLineColor(kBlue);
	hpdchi2_no->SetLineColor(kRed);
	hadchi2_no->SetLineColor(kGreen);
	
    GetRatioHist(hpdchi2_no, hmdchi2_no, hadchi2_no, (strstr(osc_par,"th23") != NULL)?leg:leg1, binning, true, Form("compare_dchi2_mach3_ptheta_aria_%s_NO_%s_%s.pdf",osc_par,expt,outstring));

	//Plot Dchi2
	//IO

	int ix_min = hpdchi2_io->GetMinimumBin()-1;
	double locmin  = hpdchi2_io->GetBinContent(ix_min+1);

	for(int i=0;i<binning.nbins;++i) hpdchi2_io->SetBinContent(i+1, hpdchi2_temp->GetBinContent(i+1)-locmin);
    
	//hmdchi2_io->SetTitle("IO");
	//hmdchi2_io->GetXaxis()->SetTitle(binning.xtitle);
	hmdchi2_io->GetYaxis()->SetTitle("#Delta#chi^{2}");
    (strstr(osc_par,"th23") != NULL)?hmdchi2_io->GetYaxis()->SetRangeUser(0,25):hmdchi2_io->GetYaxis()->SetRangeUser(0,7);
	//hadchi2_io->GetXaxis()->SetMaxDigits(2);
	hadchi2_io->GetXaxis()->SetLabelSize(0.035);
	
	hmdchi2_io->SetLineWidth(2);
	hpdchi2_io->SetLineWidth(2);
	hadchi2_io->SetLineWidth(2);
	
	hmdchi2_io->SetLineColor(kBlue);
	hpdchi2_io->SetLineColor(kRed);
	hadchi2_io->SetLineColor(kGreen);
    
    GetRatioHist(hpdchi2_io, hmdchi2_io, hadchi2_io, leg, binning, true, Form("compare_dchi2_mach3_ptheta_aria_%s_IO_%s_%s.pdf",osc_par,expt,outstring));

     
}

