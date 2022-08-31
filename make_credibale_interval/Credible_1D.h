/*-------------------------------
  Code to make credible intervals from a TH2D containing distribution of -log(L)
  C. Bronner - 2014
  -------------------------------- */

#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TColor.h>
#include <TChain.h>
#include <TMath.h>
#include <TROOT.h>
#include <TGaxis.h>
#include <TList.h>
#include <TObjArray.h>
#include <TLegend.h>

//Members
TH1D *Smoothed;

//Output


//Functions
void MakeCredible(char *iFile, char *oFile, char *type);
bool CompareBins(int b1, int b2);
void PrintHisto(TCanvas *canvas, TString name, TString filename);
void PrintCredibleRange(TH1D *Cred, TString name, TString filename, int params);
void PreparePlot(TH1D* plot, int choice);
void PrepareBestPoint(TGraph* gr);
void SetStyleVariables(TStyle *t2kStyle);

