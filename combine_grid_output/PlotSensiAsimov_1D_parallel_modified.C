#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TMath.h"

#include <iostream>
#include <string>

using namespace std;

enum OscParam {
    not_defined = -1,
    sin2213 = 0,
    deltaCP,
    dm32,
    sin223,
    OscParamCount
};

const char *OscParamName[int(OscParam::OscParamCount)] = {
    "sin2213",
    "delta",
    "dm2",
    "sin223"
};

const char *OscParamTitle[int(OscParam::OscParamCount)] = {
    "sin^{2}#theta_{13}",
    "#delta_{CP}",
    "#Delta m^{2}_{32} [eV^{2}]",
    "sin^{2}#theta_{23}"
};

typedef int MassHierarchy;
const int kNormalHierarchy   = 0;
const int kInvertedHierarchy = 1;

void plot(const char *mcdir, const char *expt, MassHierarchy plotHierarchy, const char *mhSuffix, bool margAsimov);

void PlotSensiAsimov_1D_parallel_modified(const char *mcdir, const char *expt, bool margAsimov = false) {
// const char *mcdir = "/home/t2k/lukasb/ptheta2/Minimal/inputs/outputs/MargTemplates_2019_100k_noEb_20191007T2K2_th13_dCP";
  plot(mcdir, expt, kNormalHierarchy  , ""   , margAsimov);
  plot(mcdir, expt, kInvertedHierarchy, "_IH", margAsimov);
}

void plot(const char *mcdir, const char *expt, MassHierarchy plotHierarchy, const char *mhSuffix, bool margAsimov) {

TChain *SensiTree = new TChain(margAsimov?"MargAsimov":"MargTemplate");
if (gSystem->AccessPathName(Form("%s/allfits.root", mcdir)) == 0) {
    // result of BestFit (same format)
    SensiTree->Add(Form("%s/allfits.root", mcdir));
}
else if (gSystem->AccessPathName(Form("%s/bins*.root", mcdir)) == 0) {
    // result of MargTemplates
    SensiTree->Add(Form("%s/bins*.root", mcdir));
}
else {
    // backwards compatibility
    //SensiTree->Add(Form(margAsimov?"%s/margasimov_rc_100throws*.root":"%s/margtemplates_rc_100throws*.root", mcdir));
    SensiTree->Add(Form("%s/bins*.root", mcdir));
}
char outfile[1024];
sprintf(outfile, "%s/hist_%s%s.root", mcdir, expt, mhSuffix);


OscParam gridParam = OscParam::not_defined;
for (int io = 0; io < int(OscParam::OscParamCount); io++) {
    if (SensiTree->GetBranch(OscParamName[io]) != NULL) {
        gridParam = (OscParam)io;
        break;
    }
}

if (gridParam == OscParam::not_defined) {
    cerr << "Could not determine grid param" << endl;
    exit(39);
}

Double_t AvNLLtot;
Double_t mh = 0.;
Double_t gridParamValue;
Double_t plotParamValue; // in case of sin2213 we need to transform to sin213
if(strstr(expt,"NOvA") != NULL)
  SensiTree->SetBranchAddress(Form("AvNLLtot_%s",expt), &AvNLLtot);
else
  SensiTree->SetBranchAddress("AvNLLtot", &AvNLLtot);
SensiTree->SetBranchAddress(OscParamName[int(gridParam)], &gridParamValue);

if (SensiTree->GetBranch("mh") != NULL) {
    SensiTree->SetBranchAddress("mh", &mh); // normal if < 0.5, inverted if > 0.5
}
else {
    // backwards compatibility
    if (plotHierarchy == kInvertedHierarchy) {
        cout << "Skipping to plot inverted hierarchy since mh branch is missing (backwards compatibility)" << endl;
        return;
    }
}

int Nentries_max = SensiTree->GetEntries();
double *u_gridParamValues = new double[Nentries_max]; // u_ for unsorted
double *u_plotParamValues = new double[Nentries_max];

int Nentries = 0;
for (int ie = 0; ie < Nentries_max; ie++) {
  SensiTree->GetEntry(ie);
  MassHierarchy thisHiearchy = (mh < 0.5 ? kNormalHierarchy : kInvertedHierarchy);
  if (thisHiearchy != plotHierarchy) { continue; }

  plotParamValue = gridParamValue;
  if (gridParam == OscParam::sin2213) {
      // sqrt[sin2(2th)] = sin(2th) = sin
      // sin2(x)
      // cos(2x) = cos2(x) - sin2(x) = 1 - 2sin2(x)
      // sin2(x) = 0.5 (1 - cos(2x))
      // cos(2x) = sqrt[1 - sin2(2x)]
      double sin2213 = gridParamValue;
      double sin213 = 0.5 * (1. - sqrt(1. - sin2213));
      plotParamValue = sin213;
  }

  u_gridParamValues[Nentries] = gridParamValue;
  u_plotParamValues[Nentries] = plotParamValue;
  Nentries++;
}

if (gridParam == OscParam::deltaCP && u_gridParamValues[Nentries-1] == -TMath::Pi()) {
    // fix a bug where I was accidentally clamping the grid delta as well
    u_gridParamValues[Nentries-1] = TMath::Pi();
    u_plotParamValues[Nentries-1] = TMath::Pi();
}

double gridParamValues_max = TMath::MaxElement(Nentries, u_gridParamValues);
double gridParamValues_min = TMath::MinElement(Nentries, u_gridParamValues);

int *ord_gridParamValues = new int[Nentries];
TMath::Sort(Nentries, u_gridParamValues, ord_gridParamValues, false);

//std::cout << " Nentries " << Nentries << " Nentries_max " << Nentries_max << std::endl;

double *gridParamValues = new double[Nentries];
double *plotParamValues = new double[Nentries];
for (int ie = 0; ie < Nentries; ie++) {
  gridParamValues[ie] = u_gridParamValues[ord_gridParamValues[ie]];
  plotParamValues[ie] = u_plotParamValues[ord_gridParamValues[ie]];

  //printf("%3d: %5f -> %5f\n", ie, gridParamValues[ie], plotParamValues[ie]);
}

double *plotParamBins    = new double[Nentries+1];
plotParamBins[0] = plotParamValues[0] - (plotParamValues[1]-plotParamValues[0])/2.;
for (int i = 1; i < Nentries; i++) {
  plotParamBins[i] = (plotParamValues[i]+plotParamValues[i-1])/2.;
}
plotParamBins[Nentries] = plotParamValues[Nentries-1] + (plotParamValues[Nentries-1]-plotParamValues[Nentries-2])/2.;
for (int i = 0; i <= Nentries; i++) {
    printf("plotParamBins[%d] = %g\n", i, plotParamBins[i]);
}


TH1D *cont = new TH1D("cont1", "", Nentries, plotParamBins);
cont->GetXaxis()->SetTitle(OscParamTitle[int(gridParam)]);
cont->GetYaxis()->SetTitle("-2lnL");

TH1D *cont1;
TGraph *gr = new TGraph(Nentries);
gr->GetXaxis()->SetTitle(OscParamTitle[int(gridParam)]);
gr->SetName("graph");

if (gridParam == OscParam::deltaCP) {
  cont1 = new TH1D("cont", "", Nentries, -TMath::Pi(), TMath::Pi());
  cont1->GetXaxis()->SetTitle(OscParamTitle[int(gridParam)]);
  cont1->GetYaxis()->SetTitle("-2lnL");
}
else cont->SetName("cont");

double minnll = 1e6;
for (int ie = 0; ie < Nentries_max; ie++) {
  SensiTree->GetEntry(ie);
  MassHierarchy thisHiearchy = (mh < 0.5 ? kNormalHierarchy : kInvertedHierarchy);
  if (thisHiearchy != plotHierarchy) { continue; }

  minnll = TMath::Min(minnll, AvNLLtot);
}

int iter = 0;
for (int ie = 0; ie < Nentries_max; ie++) {
  SensiTree->GetEntry(ie);
  MassHierarchy thisHiearchy = (mh < 0.5 ? kNormalHierarchy : kInvertedHierarchy);
  if (thisHiearchy != plotHierarchy) { continue; }

  double v = (u_gridParamValues[iter] - gridParamValues_min)/(gridParamValues_max - gridParamValues_min)*double(Nentries-1);
  
  // +0.5 for proper rounding
  int ix = v+0.5;
  // printf("u_gridParamValues[%d] = %g, v = %g, ix = %d, AvNLLtot = %g\n", iter, u_gridParamValues[iter], v, ix, AvNLLtot);
  cont->SetBinContent(ix+1, 2.*AvNLLtot);
  //std::cout << " plotParamValues[" << iter << "] "<< plotParamValues[iter] << std::endl;
  //if (gridParam == OscParam::deltaCP) gr->SetPoint(iter, plotParamValues[iter], 2.*AvNLLtot);
  gr->SetPoint(iter, plotParamValues[iter], 2.*AvNLLtot);
  iter++;
}

if (gridParam == OscParam::deltaCP) {
  for (int ie = 0; ie < cont1->GetNbinsX(); ie++) {
    std::cout << ie << " bin center " << cont->GetBinCenter(ie+1) << " cont1 bin center " << cont1->GetBinCenter(ie+1) << std::endl;
    cont1->SetBinContent(ie+1, gr->Eval(cont1->GetBinCenter(ie+1)));
  }
}
/*
if (gridParam == OscParam::dm2 && plotHierarchy == kInvertedHierarchy) {
    TH1D *htemp=new TH1D("htemp","htemp",81,-2.75e-3,-2.25e-3);
    TSpline3 *sp=new TSpline3(cont);
    for(int i=0;i<cont->GetNbinsX();++i) {
		double dm32 = cont->GetBinCenter(i+1);
		double dm21 = 7.53e-5;
		double dm31 = dm32 + dm21;
		double abs_dm31 = fabs(dm31);
		htemp->SetBinContent(i+1, sp->Eval(abs_dm31));
     }
}
*/

TFile *fout = new TFile(outfile, "RECREATE");
cont->Write();
if (gridParam == OscParam::deltaCP) cont1->Write();
if (gridParam == OscParam::deltaCP) gr->Write();
/*
if (gridParam == OscParam::dm32 && plotHierarchy == kInvertedHierarchy) {
    TH1D *htemp=new TH1D("htemp","htemp",Nentries,-2.75e-3,-2.25e-3);
    TSpline3 *sp=new TSpline3(cont);
    for(int i=0;i<cont->GetNbinsX();++i) {
                double dm32 = cont->GetBinCenter(i+1);
                double dm21 = 7.53e-5;
                double dm31 = dm32 + dm21;
                double abs_dm31 = fabs(dm32);
                htemp->SetBinContent(i+1, sp->Eval(abs_dm31));
     }
     htemp->SetName("cont");
     htemp->Write();
}
else cont->Write();
*/
TParameter<int> *pGridParam = new TParameter<int>("gridParam", int(gridParam));
pGridParam->Write();
TObjString *gridParamName = new TObjString(OscParamName[int(gridParam)]);
gridParamName->Write("gridParamName");
fout->Close();
cout << "Wrote to " << fout->GetName() << endl;
delete cont;
}
