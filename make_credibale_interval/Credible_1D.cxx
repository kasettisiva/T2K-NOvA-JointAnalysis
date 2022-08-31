#include "Credible_1D.h"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm> 
#include <cmath>

using namespace std;

//####################################################
int main(int argc,char** argv) 
{
  char *Fout=NULL, *Fin=NULL;
  char *Type=NULL;  
  
  for(int i=1; i<argc; i++){
    if((strcmp(argv[i],"-f")==0) && i+1<argc){
      i++;
      if(Fin!=NULL) delete[] Fin;
      Fin = new char[strlen(argv[i])+1];
      strcpy(Fin,argv[i]);
      continue;
    }

    if((strcmp(argv[i],"-o")==0) && i+1<argc){
      i++;
      if(Fout!=NULL) delete[] Fout;
      Fout = new char[strlen(argv[i])+1];
      strcpy(Fout,argv[i]);
      continue;
    }

    if((strcmp(argv[i],"-t")==0) && i+1<argc){
      i++;
      if(Type!=NULL) delete[] Type;
      Type = new char[strlen(argv[i])+1];
      strcpy(Type,argv[i]);
      continue;
    }
  }
  
  if(Fout==NULL || Fin==NULL || Type==NULL){
    std::cout <<"Usage:"<<std::endl;
    std::cout << "Credible -f <input file> -o <Output file> -t <surface type = chi2 or nll>"<<std::endl;
    exit(-1);
  }

  MakeCredible(Fin, Fout, Type);
}
//###################################################
void MakeCredible(char *iFile, char *oFile, char* type)
{

  //**** Set Style for Plots ****                                                                                   
  TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
  SetStyleVariables(t2kstyle);
  gROOT->SetStyle("T2K");
  gROOT->ForceStyle();

  TColor *color95 = new TColor(9999, 1.0, 0.9, 0.9);
  


  //  TString gr_names[4] = {"nue", "numu", "nuebar", "numubar"};
  TString gr_numbers[10] = {"0","1","2","3","4","5","6","7","8","9"};

  TString filename = TString(oFile);

  //Read initial histo into arrays
  TFile *f=new TFile(iFile,"read");
  if(f->GetListOfKeys()->Contains("Smoothed")){
    Smoothed=(TH1D*)f->Get("Smoothed");
  }
  else if(f->GetListOfKeys()->Contains("cont")){
    Smoothed=(TH1D*)f->Get("cont");
  }
  else if(f->GetListOfKeys()->Contains("SA_dcp")){
    Smoothed=(TH1D*)f->Get("SA_dcp");
  }
  else if(f->GetListOfKeys()->Contains("SA_t13")){
    Smoothed=(TH1D*)f->Get("SA_t13");
  }
  else if(f->GetListOfKeys()->Contains("SA_t23")){
    Smoothed=(TH1D*)f->Get("SA_t23");
  }
  else if(f->GetListOfKeys()->Contains("SA_dm32")){
    Smoothed=(TH1D*)f->Get("SA_dm32");
  }
  else{ Printf("ERROR - input file is not suitable.\n"); return; }

  TString stype = TString(type);
  double factor;
  if(stype.Contains("chi2", TString::kIgnoreCase)){ factor = 0.5; }
  else if(stype.Contains("nll", TString::kIgnoreCase)){ factor = 1.0; }
  else{ std::cerr << "set surface type to chi-squared (chi2) or negative log likelihood (nll)" << std::endl; return; }

  TString *inFile = new TString(iFile);
  int params;
  int theta;
  if(inFile->Contains("_dcp", TString::kIgnoreCase)){
    params = 0;
    printf("dcp - oscillation parameter.\n");
  }
  else if(inFile->Contains("_th23", TString::kIgnoreCase)){
    params = 1;
    printf("t23 - oscillation parameter.\n");
  }
  else if(inFile->Contains("_th13", TString::kIgnoreCase)){
    params = 4;
    printf("t13 - oscillation parameter.\n");
  }
  else if(inFile->Contains("_dm2", TString::kIgnoreCase) && inFile->Contains("_IH")){
    params = 5;
    printf("dm32 - oscillation parameter.\n");
  }
  else if(inFile->Contains("_dm2", TString::kIgnoreCase)){
    params = 3;
    printf("dm32 - oscillation parameter.\n");
  }
  else{ printf("WARNING - oscillation parameter unclear from filename.\n"); exit(119); }


  //Find best fit value  
  int Xmax, Ymax, Zmax;
  Smoothed->GetMinimumBin(Xmax,Ymax,Zmax);
  double MaxX=Smoothed->GetXaxis()->GetBinCenter(Xmax);
  double MaxY=Smoothed->GetBinContent(Xmax);
  TGraph *gBF=new TGraph(1,&MaxX, &MaxY);
  
  TH1D *LL = new TH1D();
  Smoothed->Copy(*LL);
  LL->SetName("Likelihood");
  

  std::vector<int> BinList;
  double Total=0.;
  for (int i=0;i<Smoothed->GetNbinsX();i++)
    {
      LL->SetBinContent(i+1,TMath::Exp(-factor*Smoothed->GetBinContent(i+1)));
      BinList.push_back(i+1);
      Total+=TMath::Exp(-factor*Smoothed->GetBinContent(i+1));
    }
  
  std::cout << "Number of bins: "<< Smoothed->GetNbinsX()<<std::endl;
  std::cout << "Total Value: "<<Total<<std::endl;

  std::stable_sort (BinList.begin(), BinList.end(),CompareBins);


  Smoothed->SetMinimum(0);
  Smoothed->SetLineColor(kBlack);
  Smoothed->SetLineWidth(2);
  PreparePlot(Smoothed, params);


  //LL->Scale(1./LL->Integral());
  LL->SetMinimum(0);
  LL->SetLineColor(kBlack);
  LL->SetLineWidth(2);
  PreparePlot(LL, params);
  LL->GetYaxis()->SetTitle("Likelihood Density");

  bool is_dCP = inFile->Contains("_dcp", TString::kIgnoreCase);
  if (is_dCP) {
      LL->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
  }

  double Added=0.;
  int iter=0;
  TH1D *Cred68=new TH1D();
  Smoothed->Copy(*Cred68);
  Cred68->Reset();
  Cred68->SetName("Cred68");
  Cred68->SetLineColor(kRed-9);
  Cred68->SetLineWidth(1);
  //  Cred68->SetLineStyle(7);
  Cred68->SetFillColor(kRed-9);
  //  Cred68->SetFillStyle(3844);

  TH1D *Cred90=new TH1D();
  Smoothed->Copy(*Cred90);
  Cred90->Reset();
  Cred90->SetName("Cred90");
  Cred90->SetLineColor(kRed-10);
  Cred90->SetLineWidth(1);
  //  Cred90->SetLineStyle(7);
  Cred90->SetFillColor(kRed-10);
  //  Cred90->SetFillStyle(3444);

  TH1D *Cred95=new TH1D();
  Smoothed->Copy(*Cred95);
  Cred95->Reset();
  Cred95->SetName("Cred95");
  Cred95->SetLineColor(9999);
  Cred95->SetLineWidth(1);
  //  Cred95->SetLineStyle(7);
  Cred95->SetFillColor(9999);
  //  Cred95->SetFillStyle(3244);


  TH1D *Cred95_LL = new TH1D();
  Cred95->Copy(*Cred95_LL);
  Cred95_LL->SetName("Cred95_LL");

  TH1D *Cred90_LL = new TH1D();
  Cred90->Copy(*Cred90_LL);
  Cred90_LL->SetName("Cred90_LL");

  TH1D *Cred68_LL = new TH1D();
  Cred68->Copy(*Cred68_LL);
  Cred68_LL->SetName("Cred68_LL");


  for (int i=0;i<Smoothed->GetNbinsX();i++){ 
    Cred68->SetBinContent(i+1, 0);
    Cred90->SetBinContent(i+1, 0);
    Cred95->SetBinContent(i+1, 0);
    Cred68_LL->SetBinContent(i+1, 0);
    Cred90_LL->SetBinContent(i+1, 0);
    Cred95_LL->SetBinContent(i+1, 0);
  }

  while (iter<fabs(BinList.size()) && Added/Total*100.<68.27)
    {
      if(!is_dCP || TMath::Abs(Smoothed->GetBinCenter(BinList.at(iter))) < TMath::Pi()){
      Added+=TMath::Exp(-factor*Smoothed->GetBinContent(BinList.at(iter)));
      Cred68->SetBinContent(BinList.at(iter), Smoothed->GetBinContent(BinList.at(iter)));
      Cred68_LL->SetBinContent(BinList.at(iter), LL->GetBinContent(BinList.at(iter)));
      }
      iter++;
    }

  std::cout << (Added/Total)*100. << std::endl;

  Added=0;
  iter=0;

  while (iter<fabs(BinList.size()) && Added/Total*100.<90.)
    {
      if(!is_dCP || TMath::Abs(Smoothed->GetBinCenter(BinList.at(iter))) < TMath::Pi()){
      Added+=TMath::Exp(-factor*Smoothed->GetBinContent(BinList.at(iter)));
      Cred90->SetBinContent(BinList.at(iter), Smoothed->GetBinContent(BinList.at(iter)));
      Cred90_LL->SetBinContent(BinList.at(iter), LL->GetBinContent(BinList.at(iter)));
      }
      iter++;
    }

  std::cout << (Added/Total)*100. << std::endl;


  Added=0;
  iter=0;

  while (iter<fabs(BinList.size()) && Added/Total*100.<95.)
    {
      if(!is_dCP || TMath::Abs(Smoothed->GetBinCenter(BinList.at(iter))) < TMath::Pi()){
      Added+=TMath::Exp(-factor*Smoothed->GetBinContent(BinList.at(iter)));
      Cred95->SetBinContent(BinList.at(iter), Smoothed->GetBinContent(BinList.at(iter)));
      Cred95_LL->SetBinContent(BinList.at(iter), LL->GetBinContent(BinList.at(iter)));
      }
      iter++;
    }

	//      printf("%f\n", TMath::Abs(Smoothed->GetBinCenter(BinList.at(iter))));

  std::cout << (Added/Total)*100. << std::endl;


  TFile *out1=new TFile(filename+".root","recreate");
  out1->cd();

  TCanvas *c1 = new TCanvas();
  c1->cd();
  Smoothed->Draw();
  Cred95->Draw("same");
  Cred90->Draw("same");
  Cred68->Draw("same");
  Smoothed->Draw("same");
  c1->RedrawAxis();
  c1->RedrawAxis("G");

  PrintHisto(c1, "SensiAsimov_all", filename);

  LL->SetMaximum(1.3*LL->GetBinContent(LL->GetMaximumBin()));

  TLegend *leg2 = new TLegend(0.7,0.8,0.85,0.95);
  leg2->AddEntry(Cred68_LL, "68.3%", "F");
  leg2->AddEntry(Cred90_LL, "90%", "F");
  leg2->AddEntry(Cred95_LL, "95%", "F");


  TCanvas *c2 = new TCanvas();
  c2->cd();
  LL->Draw("C");
  Cred95_LL->Draw("same");
  Cred90_LL->Draw("same");
  Cred68_LL->Draw("same");
  LL->Draw("HIST C same");
  c2->RedrawAxis();
  c2->RedrawAxis("G");
  leg2->Draw("same");

  PrintHisto(c2, "SensiAsimov_ll_all", filename);
  PrintCredibleRange(Cred95_LL, "CredRange_95", filename, params);
  PrintCredibleRange(Cred90_LL, "CredRange_90", filename, params);
  PrintCredibleRange(Cred68_LL, "CredRange_68", filename, params);


  //NLL2D->Write();
  //  Smoothed->Write();
  Smoothed->Write("NLL");
  LL->Write("LL");
  gBF->Write("gBF");
  Cred68->Write();
  Cred90->Write();
  Cred95->Write();
  Cred68_LL->Write();
  Cred90_LL->Write();
  Cred95_LL->Write();
  out1->Close();
}
//###################################################
bool CompareBins(int b1, int b2)
{
  return (TMath::Exp(-0.5*Smoothed->GetBinContent(b1)) > TMath::Exp(-0.5*Smoothed->GetBinContent(b2)));
}
//#############################################################################
void PrintHisto(TCanvas *canvas, TString name, TString filename){
  canvas->Print(filename+"/"+name+".eps");
  canvas->Print(filename+"/"+name+".pdf");
  canvas->Print(filename+"/"+name+".png");
  canvas->Print(filename+"/"+name+".C");
}
//#############################################################################
void PrintCredibleRange(TH1D *Cred, TString name, TString filename, int params){
  std::vector<std::pair<double,double> > ranges;
  bool inRange = false;
  double rangeStart;
  double rangeEnd;
  for (int ix = 0; ix < Cred->GetNbinsX(); ix++) {
      double val = Cred->GetBinContent(ix+1);
      if (val != 0.) {
          if (!inRange) {
              rangeStart = Cred->GetXaxis()->GetBinLowEdge(ix+1);
              inRange = true;
          }
          // if already in range, just keep going
          // but for convenience (and especially for the last bin
          // we keep extending the edge of the range
          // note that we intentionally don't do an "else" here,
          // since we might have a range consisting of a single bin only
          rangeEnd = Cred->GetXaxis()->GetBinUpEdge(ix+1);
      }
      else {
          // not in range anymore
          if (inRange) {
              ranges.push_back(std::make_pair(rangeStart, rangeEnd));
              inRange = false;
          }
      }
  }

  if (inRange) {
      ranges.push_back(std::make_pair(rangeStart, rangeEnd));
      inRange = false;
  }

  int ixmax = Cred->GetMaximumBin()-1;
  double xmax = Cred->GetXaxis()->GetBinCenter(ixmax+1);

  std::string units = "";
  if (params == 3 || params == 5) {
      // dm2, convert to [10^{-3} eV^2]
      xmax *= 1000.;
      for (unsigned i = 0; i < ranges.size(); i++) {
          ranges.at(i).first  *= 1000.;
          ranges.at(i).second *= 1000.;
      }
      units = "\\times 10^{-3}\\,\\mathrm{eV^2}";
  }
  if (params == 4) { // sin2213
      xmax *= 1000.;
      for (unsigned i = 0; i < ranges.size(); i++) {
          ranges.at(i).first  *= 1000.;
          ranges.at(i).second *= 1000.;
      }
      units = "\\times 10^{-3}";
  }

  // now we have all the ranges, so let's print them
  string outname = string(filename+"/"+name+".txt");
  ofstream of(outname);
  of << std::fixed;
  if (params == 0) { // dcp
      of << std::setprecision(2);
  }
  else if (params == 4) { // sin2213
      of << std::setprecision(1);
  }
  else { // sin223, dm1
      of << std::setprecision(3);
  }
  of << "$" << xmax << units << "$ & ";
  for (unsigned i = 0; i < ranges.size(); i++) {
      if (i > 0) {
          of << " and ";
      }
      of << "$[" << ranges.at(i).first << ", " << ranges.at(i).second << "]" << units << "$";
  }
  of << " \\\\" << endl;
  of.close();
  cout << "Wrote to " << outname << endl;
}
//#############################################################################
void PreparePlot(TH1D* plot, int choice){
  plot->GetYaxis()->SetTitle("#Delta #chi^{2}");
  if     (choice == 0){ plot->GetXaxis()->SetTitle("#delta_{CP}"); }
  else if(choice == 1){ plot->GetXaxis()->SetTitle("sin^{2}#theta_{23}"); }
  else if(choice == 4){ plot->GetXaxis()->SetTitle("sin^{2}#theta_{13}"); }
  else if(choice == 2){ plot->GetXaxis()->SetTitle("sin^{2}2#theta_{13}"); }
  else if(choice == 3){ plot->GetXaxis()->SetTitle("#Deltam^{2}_{32} [eV^{2}]"); }
  else if(choice == 5){ plot->GetXaxis()->SetTitle("|#Deltam^{2}_{31}| [eV^{2}]"); }
  TGaxis* gx = (TGaxis*)plot->GetXaxis();
  gx->SetMaxDigits(3);
  TGaxis* gy = (TGaxis*)plot->GetYaxis();
  gy->SetMaxDigits(3);
}
//#############################################################################
void PrepareBestPoint(TGraph* gr){
  gr->SetMarkerStyle(34);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerSize(0.5);
}
//#############################################################################
void SetStyleVariables(TStyle *t2kStyle){

  // use plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  // set the paper & margin sizes
  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.05);
  t2kStyle->SetPadRightMargin(0.15); //0.05 
  t2kStyle->SetPadBottomMargin(0.16);
  t2kStyle->SetPadLeftMargin(0.13);

  // use large Times-Roman fonts
  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.05,"x");
  t2kStyle->SetTitleSize(0.06,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.06,"y");
  t2kStyle->SetLabelSize(0.05,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.25);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  t2kStyle->SetPadGridX(true);
  t2kStyle->SetPadGridY(true);

  // use bold lines and markers
  //  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //  t2kStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                   NCont);
  t2kStyle->SetNumberContours(NCont);

  // End of definition of t2kStyle
}
