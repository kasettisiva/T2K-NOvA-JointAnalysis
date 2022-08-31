#include "Bicubic_1D.h"
#include <vector>
#include <iostream>
#include <algorithm> 
#include <cmath>

//####################################################
int main(int argc,char** argv) 
{
  char *Fout=NULL, *Fin=NULL, *Hist=NULL;
  int nReBin=-1;
  
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

    if((strcmp(argv[i],"-n")==0) && i+1<argc){
      i++;
      nReBin=atoi(argv[i]);
      continue;
    }

    if((strcmp(argv[i],"-h")==0) && i+1<argc){
      i++;
      if(Hist!=NULL) delete[] Hist;
      Hist = new char[strlen(argv[i])+1];
      strcpy(Hist,argv[i]);
      continue;
    }


  }
  
  if(Fout==NULL || Fin==NULL || Hist==NULL || nReBin<0){
    std::cout <<"Usage:"<<std::endl;
    std::cout << "Bicubic -f <input file> -o <Output file> -n <Rebin factor> -h <input histogram name>"<<std::endl;
    exit(-1);
  }

  Resize(Fin, Fout, nReBin, Hist);
}
//###################################################
TH1D *Resize(int nReBin, TH1D *NLL1D)
{
  std::vector<double> in_vec(NLL1D->GetNbinsX());
  for (int i=0;i<NLL1D->GetNbinsX();i++)
    {
      in_vec[i]=NLL1D->GetBinContent(i+1);
    }      
  std::vector<double> out_vec=bicubicresize(in_vec, NLL1D->GetNbinsX(), nReBin);
 
  TH1D * Smoothed=new TH1D("Smoothed","Test interpolation",nReBin*NLL1D->GetNbinsX(),NLL1D->GetXaxis()->GetXmin(),NLL1D->GetXaxis()->GetXmax());
  
  for (int i=0;i<nReBin*NLL1D->GetNbinsX();i++)
    {
      Smoothed->SetBinContent(i+1, out_vec.at(i));
    }

  return Smoothed;
}
//###################################################
void Resize(char *iFile, char *oFile, int nReBin, char *nHist)
{
  //Read initial histo into vector
  TFile *f=new TFile(iFile,"read");
  //  TH2D *NLL2D=(TH2D*)f->Get("NLL2D_IH");
  TH1D *NLL2D=(TH1D*)f->Get(nHist);
  TH1D *Smoothed=Resize(nReBin, NLL2D);

  TFile *out1=new TFile(oFile,"recreate");
  NLL2D->Write();
  Smoothed->Write();
  out1->Close();
}
//###################################################
inline double getpixel(const std::vector<double>& in_vec, int src_width, int x)
{
  if (x>=0 && x < src_width) return in_vec[x];
  else if (x<0) return in_vec[0];
  else if (x>=src_width) return in_vec[src_width-1];

  return 0;
}
//###################################################
std::vector<double> bicubicresize(const std::vector<double>& in_vec, int src_width, int Rebin)
{
  int dest_width=src_width*Rebin;

  std::cout<<"Original number of bins: x=" <<src_width<<std::endl;
  std::cout<<"Target number of bins: x=" <<dest_width<<std::endl;
 
  std::vector<double> out(dest_width);

  const double tx=1./((double)Rebin);

  std::cout << tx<<std::endl;

  for (int i=0;i<dest_width;i++)
    {
      int p=int((i-(double)Rebin/2.)*tx);
      double a=(i-(double)Rebin/2.)*tx-p;
      double Val=0;
      int OK=1;
      for (int m=-1;m<3;m++)
	    {
	      if (OK==0)break;
	      double ashift=-((double)m-a);
	      // Val+=getpixel(in_vec,src_width,src_height,p+m,q+n)*Rc(ashift)*Rc(bshift);
	      Val+=getpixel(in_vec,src_width,p+m)*CatMullRom(ashift);
	      if (Val>1e6)
		{
		  std::cout<<getpixel(in_vec,src_width,p+m)<<" "<<p+m<<std::endl;
		  std::cout<<ashift<<" "<<Rc(ashift)<<std::endl;
		  OK=0;
		}
	    }
      out[i]=Val;
      //  if (Val>0)std::cout<<i<<" "<<j<<" "<<Val<<std::endl;
    }
  return out;
}
//###################################################
double Rc(double x)
{
  double f = x;
  if( f < 0.0 )
    {
      f = -f;
    }
  
  if( f >= 0.0 && f <= 1.0 )
    {
      return ( 2.0 / 3.0 ) + ( 0.5 ) * ( f* f * f ) - (f*f);
    }
  else if( f > 1.0 && f <= 2.0 )
    {
      return 1.0 / 6.0 *TMath::Power( ( 2.0 - f  ), 3.0 );
    }
  return 1.0; 
}
//###################################################
double CatMullRom( double x )
{
  const double B = 0.0;
  const double C = 0.5;
  double f = x;
  if( f < 0.0 )
    {
      f = -f;
    }
  if( f < 1.0 )
    {
      return ( ( 12 - 9 * B - 6 * C ) * ( f * f * f ) +
	       ( -18 + 12 * B + 6 *C ) * ( f * f ) +
	       ( 6 - 2 * B ) ) / 6.0;
    }
  else if( f >= 1.0 && f < 2.0 )
    {
      return ( ( -B - 6 * C ) * ( f * f * f )
	       + ( 6 * B + 30 * C ) * ( f *f ) +
	       ( - ( 12 * B ) - 48 * C  ) * f +
	       8 * B + 24 * C)/ 6.0;
    }
  else
    {
      return 0.0;
    }
} 
//###################################################
