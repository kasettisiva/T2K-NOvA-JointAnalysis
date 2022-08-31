/*-------------------------------
  Code using bicubic interpolation to increase the number of bins in a TH2D
  C. Bronner - 2014
  -------------------------------- */

#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"

//Members


//Output


//Functions
void Resize(char *iFile, char *oFile, int nReBin, char *nHist);
inline double getpixel(const std::vector<double>& in_vec, int src_width, int x);

std::vector<double> bicubicresize(const std::vector<double>& in_vec, int src_width, int Rebin);
double Rc(double x);
double CatMullRom( double x );
