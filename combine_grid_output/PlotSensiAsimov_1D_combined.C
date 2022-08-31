#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TParameter.h"

using namespace std;

enum OscParam {
    not_defined = -1,
    sin2213 = 0,
    deltaCP,
    dm32,
    sin223,
    OscParamCount
};

const char *OscParamNameShort[OscParamCount] = {
    "sin213",
    "dcp",
    "dm2",
    "sin223"
};

const char *OscParamNameForOutput[OscParamCount] = {
    "th13",
    "dCP",
    "dm2",
    "th23"
};

static int kUnusedColor = 11310;

void plotForOneRC(const char *mcdir, const char *expt, const char *runversion, bool wRC, const char *plotdir, const char *CriticalChi2File=NULL, int CriticalChi2LevelIndex=1, bool drawErrorBand=false);

void PlotSensiAsimov_1D_combined(const char *mcdir, const char *expt, const char *runversion, const char *plotdir, const char *CriticalChi2File=NULL, int CriticalChi2LevelIndex=1) {
    if (CriticalChi2File) { // we only compute CriticalChi2 for wRC
        if (CriticalChi2LevelIndex >= 0) {
            plotForOneRC(mcdir, expt, runversion, true,  plotdir, CriticalChi2File, CriticalChi2LevelIndex,0);
        }
        plotForOneRC(mcdir, expt, runversion, true,  plotdir, CriticalChi2File, CriticalChi2LevelIndex,1);
    }
    else {
        plotForOneRC(mcdir, expt, runversion, true,  plotdir);
        plotForOneRC(mcdir, expt, runversion, false, plotdir);
    }
}

void GetConfidenceInterval(std::vector<double> &xmin, std::vector<double> &xmax, std::vector<double> &xedges, TGraph *gCriticalChi2, TH1D *cont) {
    double xstart = gCriticalChi2->GetX()[0];
    double xend   = gCriticalChi2->GetX()[gCriticalChi2->GetN()-1];
    double xinc   = (xend-xstart)/700.;
    bool inBand = false;
    for (double x = xstart; x < xend; x += xinc) {
        double criticalChi2 = gCriticalChi2->Eval(x);
        double chi2 = cont->Interpolate(x);
        if (!inBand && chi2 < criticalChi2) {
            // start a band
            printf("Start band at %g\n", x);
            xmin.push_back(x);
            if (x != xstart) {
                xedges.push_back(x);
            }
            inBand = true;
        }
        else if (inBand && chi2 > criticalChi2) {
            printf("End band at %g\n", x);
            xmax.push_back(x);
            xedges.push_back(x);
            inBand = false;
        }
    }
    if (inBand) {
        xmax.push_back(xend);
    }

}

std::string GetLatexInterval(std::vector<double> xmin, std::vector<double> xmax, std::vector<double> xedges, OscParam gridParam) {
    std::stringstream latex;
    latex << "$";
    for (unsigned ib = 0; ib < xmin.size(); ib++) {
        if (ib > 0) {
            latex << " \\cup ";
        }
        switch (gridParam) {
            case deltaCP:
            latex << std::fixed << std::setprecision(2);
            break;

            case sin223:
            default:
            latex << std::fixed << std::setprecision(3);
            break;
        }
        latex << "[";
        if (xmin.at(ib) == -TMath::Pi()) {
            latex << "-\\pi";
        }
        else {
            latex << xmin.at(ib);
        }
        latex << ",";
        if (xmax.at(ib) == +TMath::Pi()) {
            latex << "\\pi";
        }
        else {
            latex << xmax.at(ib);
        }
        latex << "]";
    }
    latex << "$";
    return latex.str();
}

int *MakeLevelColors(int nlevels, int baseColorCode) {
    int *levelColorCodes = new int[nlevels];
    for (int il = 0; il < nlevels; il++) {
        levelColorCodes[il] = kUnusedColor++;
        TColor *baseColor = gROOT->GetColor(baseColorCode);
        Float_t h,l,s;
        baseColor->GetHLS(h,l,s);
        l = l + (1.-l)*double(il)/double(nlevels); // make the higher levels lighter
        Float_t r,g,b;
        TColor::HLS2RGB(h,l,s,r,g,b);
        TColor *newColor = new TColor(levelColorCodes[il],r,g,b);
    }
    return levelColorCodes;
}

void plotForOneRC(const char *mcdir, const char *expt, const char *runversion, bool wRC, const char *plotdir, const char *CriticalChi2File, int CriticalChi2LevelIndex, bool drawErrorBand) {
    const char *RCsuffix = "";
    if (!wRC) {
        RCsuffix = "_noRC";
    }

    char prefix[1024];
    bool backwardsCompatibility = false;
    if (strstr(mcdir, "/") != NULL) {
        strcpy(prefix, mcdir);
    }
    else {
        // backwards compatibility
        // before we had suffix instead of mcdir
        sprintf(prefix, "/home/t2k/lukasb/ptheta2/Minimal/inputs/outputs/MargTemplates_OA2019_100k_%s_%s_dCP", mcdir, runversion);
        //sprintf(prefix, "/home/t2k/lukasb/ptheta2/Minimal/inputs/outputs/MargTemplates_OA2019_100k_%s_%s_dCP", mcdir, runversion);
        backwardsCompatibility = true;
    }

    const int nc = 2;
    const char *confSuffix[nc] = { "", "_IH" };
    const char *MHname[nc] = { "NH", "IH" };
    char fname[nc][1024];
    for (int ic = 0; ic < nc; ic++) {
        //sprintf(fname[ic], "%s/hist_%s%s.root", prefix, RCsuffix, confSuffix[ic]);
        sprintf(fname[ic], "%s/hist_%s%s.root", prefix, expt, confSuffix[ic]);
    }

    const char *title[nc] = { "Normal ordering", "Inverted ordering" };
    int colors[nc] = { kAzure+1, kOrange+1 };

    TFile *fin[nc];
    TH1D *cont[nc];
    OscParam gridParam = deltaCP; // backwards compatibility
    double minchi2[nc]; // here separate minima for NO and IO, later we take the global minimum (different than disappearance)
    for (int ic = 0; ic < nc; ic++) {
      minchi2[ic] = 1e9;
      fin[ic] = new TFile(fname[ic]);
      cont[ic] = (TH1D *)fin[ic]->Get("cont");
      cont[ic]->SetName(Form("cont_%d", ic));
      for (int ix = 1; ix <= cont[ic]->GetNbinsX(); ix++) {
        double val = cont[ic]->GetBinContent(ix);
        if (val < minchi2[ic]) {
          minchi2[ic] = val;
        }
      }

      TParameter<int> *pGridParam = (TParameter<int> *)fin[ic]->Get("gridParam");
      if (pGridParam != NULL) {
          gridParam = static_cast<OscParam>(int(pGridParam->GetVal()));
      }
    }

    TFile *fCriticalChi2 = NULL;
    TGraph *gCriticalChi2[nc];
    TObjArray *CriticalChi2LevelLabels = NULL;

    int nlevels = 0;
    TGraph **gCriticalChi2All[nc]; // all levels [ic][il]
    int *levelStyles; // [il]
    int *levelColors[nc]; // [ic][il]

    Int_t pastelWhite = kUnusedColor++; // actual white, not transparent (this is just an arbitrary large number)
    TColor *pastelWhiteColor = new TColor(pastelWhite, 1., 1., 1.);

    if (CriticalChi2File) {
        fCriticalChi2 = new TFile(CriticalChi2File);
        if (fCriticalChi2->IsZombie()) {
            cerr << "Could not open " << CriticalChi2File << endl;
            return;
        }

        CriticalChi2LevelLabels = (TObjArray *)fCriticalChi2->Get("levelLabels");
        if (CriticalChi2LevelLabels == NULL) {
            cerr << "Could not get level labels" << endl;
            fCriticalChi2->ls();
            return;
        }
        nlevels = CriticalChi2LevelLabels->GetEntries();
        levelStyles = new int[nlevels];
        for (int il = 0; il < nlevels; il++) {
            if (il == 0) {
                levelStyles[il] = 3454; // solid
            }
            else if (il == 1) {
                levelStyles[il] = 3245; // 45 deg diagonal in both directions
            }
            else if (il == 2) {
                levelStyles[il] = 3244; // 45 deg diagonal
            }
            else if (il == 3) {
                levelStyles[il] = 3190; // dots
            }
            else {
                levelStyles[il] = 3003+il;
            }
            // levelStyles[il] = 1001;
        }

        for (int ic = 0; ic < nc; ic++) {
            // we use Upper (which has larger Dchi2 than Central) as a more conservative limit
            if (CriticalChi2LevelIndex >= 0) {
                gCriticalChi2[ic] = (TGraph *)fCriticalChi2->Get(Form("gCritical_%s_%d_Upper", MHname[ic], CriticalChi2LevelIndex));
                if (gCriticalChi2[ic] == NULL) {
                    cerr << "Could not get the CriticalChi2 graph" << endl;
                    fCriticalChi2->ls();
                    return;
                }
            }
            else {
                gCriticalChi2[ic] = NULL;
            }

            gCriticalChi2All[ic] = new TGraph*[nlevels];
            levelColors[ic] = MakeLevelColors(nlevels, colors[ic]);
            for (int il = 0; il < nlevels; il++) {
                gCriticalChi2All[ic][il] = (TGraph *)fCriticalChi2->Get(Form("gCritical_%s_%d_Upper", MHname[ic], il));
            }
        }
    }

    int *levelLegendColors = MakeLevelColors(nlevels, kGray+2);

    double global_minchi2 = TMath::MinElement(nc, minchi2);
    for (int ic = 0; ic < nc; ic++) {
      for (int ix = 1; ix <= cont[ic]->GetNbinsX(); ix++) {
        double val = cont[ic]->GetBinContent(ix);
        //cont[ic]->SetBinContent(ix,val - global_minchi2);
        cont[ic]->SetBinContent(ix,val - minchi2[ic]);
      }
    }

    const int ncont = 3;
    double contLevels[ncont] = { 1., 4., 9. };
    const char *contTitle[ncont] = { "1#sigma C.L.", "2#sigma C.L.", "3#sigma C.L." };
    const char *contName [ncont] = { "1s", "2s", "3s" };
    int styles[ncont] = { 2, 2, 2 };

    TH2D *contPerLine[nc][ncont];
    TH2D *dummyForLegend[ncont];
    double legendXoffset = 0.;

    TCanvas *c1 = new TCanvas;

    // get axes with nicer limits
    TH1D *contAxes = NULL;
    switch (gridParam) {
        case deltaCP:
            contAxes = new TH1D("cont_axes", "", 10, -TMath::Pi(), TMath::Pi());
            contAxes->SetMaximum(wRC ? 27.5 : 16.);
            break;

        case sin2213:
            if (wRC) {
                contAxes = new TH1D("cont_axes", "", 10, 0.018, 0.026);
                legendXoffset = 0.2;
            }
            else {
                double xmin = 0.;
                double xmax = 0.05;
                if (cont[0]->GetXaxis()->GetXmax() > 0.08) {
                    xmax = 0.10;
                }
                contAxes = new TH1D("cont_axes", "", 10, xmin, xmax);
                legendXoffset = 0.3;
            }
            contAxes->SetMaximum(25.);
            break;

        case dm32:
            {
                double xmin = 2.3e-3;
                double xmax = 2.7e-3;
                if (cont[0]->GetXaxis()->GetXmin() < 2.2e-3) {
                    xmin = 1.5e-3;
                    xmax = 3.5e-3;
                }
                contAxes = new TH1D("cont_axes", "", 10, xmin, xmax);
                contAxes->SetMaximum(10.);
                legendXoffset = 0.25;
            }
            break;

        case sin223:
            contAxes = new TH1D("cont_axes", "", 10, 0.3, 0.7);
            contAxes->SetMaximum(20.);
            legendXoffset = 0.25;
            break;

        default:
            std::cerr << "Not implemented." << std::endl;
            exit(189);
            return;
    }
    contAxes->GetXaxis()->SetTitle(cont[0]->GetXaxis()->GetTitle());
    if (gridParam == dm32) {
        contAxes->GetXaxis()->SetTitle("#Deltam^{2}_{32} (NO) / |#Deltam^{2}_{31}| (IO) [eV^{2}]");
        TGaxis::SetMaxDigits(3);
        c1->SetRightMargin(0.10);
    }
    contAxes->GetYaxis()->SetTitle("#Delta#chi^{2}");
    double largeNum = 1e6;
    contAxes->SetLineWidth(0);
    contAxes->SetMinimum(0.);
    contAxes->GetYaxis()->SetRangeUser(0,25);
    contAxes->Draw("");

    TLegend *l1 = new TLegend(legendXoffset + 0.17, 0.69, legendXoffset + 0.45, 0.85);
    l1->SetBorderSize(0);

    TLegend *l2 = new TLegend(legendXoffset + 0.17, 0.52, legendXoffset + 0.45, 0.68);
    l2->SetBorderSize(0);
    TLegend *l2o = new TLegend(legendXoffset + 0.17, 0.52, legendXoffset + 0.45, 0.68);
    l2o->SetBorderSize(0);

    TF1 *f1[ncont];
    std::vector<std::string> CIlatex(nc);
    if (CriticalChi2File) {
        if (drawErrorBand) {
            if (CriticalChi2LevelIndex >= 0) {
                // draw one CI
                for (int ic = 0; ic < nc; ic++) {
                    vector<double> xmin;
                    vector<double> xmax;
                    vector<double> xedges;
                    GetConfidenceInterval(xmin, xmax, xedges, gCriticalChi2[ic], cont[ic]);
                    CIlatex.at(ic) = GetLatexInterval(xmin, xmax, xedges, gridParam);
                    for (unsigned ib = 0; ib < xmin.size(); ib++) {
                        double  x = (xmax.at(ib) + xmin.at(ib))/2.;
                        double ex = (xmax.at(ib) - xmin.at(ib))/2.;
                        double  y = 0.;
                        double ey = 2.*contAxes->GetMaximum();
                        TGraphAsymmErrors *band = new TGraphAsymmErrors(1,&x,&y,&ex,&ex,&ey,&ey);
                        band->SetFillColor(colors[ic]);
                        band->SetFillStyle(3245+9*ic);
                        band->Draw("2");
                    }
                    for (unsigned ib = 0; ib < xedges.size(); ib++) {
                        double xx[2];
                        xx[0] = xedges.at(ib);
                        xx[1] = xedges.at(ib);
                        double yy[2];
                        yy[0] = -1.;
                        yy[1] = 2.*contAxes->GetMaximum();
                        TGraph *line = new TGraph(2,xx,yy);
                        line->SetLineColor(colors[ic]);
                        line->SetLineWidth(2);
                        line->Draw("L");
                    }
                }
            }
            else {
                vector<double> xmin;
                vector<double> xmax;
                vector<double> xedges;
                for (int ic = nc-1; ic >= 0; ic--) {
                    // draw all CI
                    for (int il = nlevels-1-(CriticalChi2LevelIndex==-2 ? 1 : 0); il >= 0; il--) {
                        xmin.clear();
                        xmax.clear();
                        xedges.clear();
                        GetConfidenceInterval(xmin, xmax, xedges, gCriticalChi2All[ic][il], cont[ic]);
                        TH1D *ymaxCont = cont[ic];
                        TH1D *yminCont = NULL;
                        if (ic > 0) {
                            // for IH we don't want the band to reach into the NH one
                            yminCont = cont[ic-1];
                        }
                        for (unsigned ib = 0; ib < xmin.size(); ib++) {
                            double dx = cont[ic]->GetXaxis()->GetBinWidth(1)/4.;
                            int nx = TMath::CeilNint((xmax.at(ib)-xmin.at(ib))/dx);
                            double *xval = new double[nx];
                            double *yval = new double[nx];
                            double *yerr = new double[nx];
                            for (int ix = 0; ix < nx; ix++) {
                                xval[ix] = xmin.at(ib) + double(ix)/double(nx-1) * (xmax.at(ib) - xmin.at(ib));
                                double ymax = ymaxCont->Interpolate(xval[ix]);
                                double ymin = 0.;
                                if (yminCont) {
                                    ymin = yminCont->Interpolate(xval[ix]);
                                }
                                yval[ix] = (ymax+ymin)/2.;
                                yerr[ix] = (ymax-ymin)/2.;
                            }
                            TGraphErrors *hband = new TGraphErrors(nx, xval, yval, 0, yerr);
                            hband->SetLineWidth(0);

                            TGraphErrors *hbandClear = (TGraphErrors *)hband->Clone();
                            hbandClear->SetFillStyle(1001); // solid
                            hbandClear->SetFillColor(levelColors[ic][il]);
                            hbandClear->Draw("3"); // filled

                            hband->SetFillColor(pastelWhite);
                            hband->SetFillStyle(levelStyles[il]);
                            hband->Draw("3"); // filled
                        }
                        for (unsigned ib = 0; ib < xedges.size(); ib++) {
                            double xx[2];
                            xx[0] = xedges.at(ib);
                            xx[1] = xedges.at(ib);
                            double yy[2];
                            yy[0] = -1.;
                            yy[1] = ymaxCont->Interpolate(xedges.at(ib));
                            TGraph *line = new TGraph(2,xx,yy);
                            line->SetLineColor(levelColors[ic][TMath::Max(0,il-2)]);
                            // line->SetLineColor(colors[ic]);
                            line->SetLineWidth(1);
                            line->Draw("L");
                        }
                    }
                }
                for (int il = 0; il < nlevels-(CriticalChi2LevelIndex==-2?1:0); il++) {
                    TGraph *gForLeg = new TGraph();
                    gForLeg->SetFillStyle(levelStyles[il]);
                    // gForLeg->SetFillColor(kGray);
                    // gForLeg->SetFillColor(kGray+nlevels-il-1);
                    gForLeg->SetFillColor(pastelWhite);
                    gForLeg->SetLineColor(kGray+3);
                    // gForLeg->SetLineColor(kGray+nlevels-TMath::Max(0,il-2)-1);
                    TObjString *lab = (TObjString *)CriticalChi2LevelLabels->At(il);
                    l2->AddEntry(gForLeg, lab->GetName(), "F");

                    TGraph *gForLegOverlay = (TGraph *)gForLeg->Clone();
                    gForLegOverlay->SetFillStyle(1001);
                    gForLegOverlay->SetFillColor(levelLegendColors[il]);
                    gForLegOverlay->SetLineColor(kWhite); // hide
                    l2o->AddEntry(gForLegOverlay, Form("#color[0]{%s}", lab->GetName()), "F"); // 0=white=hide
                }
            }
        }
        else {
            for (int ic = 0; ic < nc; ic++) {
                if (CriticalChi2LevelIndex >= 0) {
                    gCriticalChi2[ic]->SetLineColor(colors[ic]);
                    gCriticalChi2[ic]->SetLineWidth(2);
                    gCriticalChi2[ic]->SetLineStyle(2); // dashed
                    gCriticalChi2[ic]->Draw("L");
                    if (ic == 0) {
                        TGraph *gForLegend = (TGraph *)gCriticalChi2[ic]->Clone();
                        gForLegend->SetLineColor(kBlack);
                        TObjString *levelLabel = (TObjString *)CriticalChi2LevelLabels->At(CriticalChi2LevelIndex);
                        l2->AddEntry(gForLegend, Form("%s critical #Delta#chi^{2}", levelLabel->GetName()));
                    }
                }
            }
        }
    }
    else {
        for (int il = 0; il < ncont; il++) {
          f1[il] = new TF1(Form("f1_%d", il), Form("%g", contLevels[il]), contAxes->GetXaxis()->GetXmin(), contAxes->GetXaxis()->GetXmax());
          f1[il]->SetLineColor(kBlack);
          f1[il]->SetLineStyle(styles[il]);
          f1[il]->SetLineWidth(1);
          f1[il]->Draw("same");
          // l2->AddEntry(f1[il], contTitle[il], "L"); // don't show 1sigma etc. unless we actually do FC
        }
    }

    char plotname[1024];
    if (backwardsCompatibility) {
        sprintf(plotname, "%s/asimovA_dCP%s_%s", plotdir, RCsuffix, runversion);
    }
    else {
        sprintf(plotname, "%s/contour_%s%s_%s_%s", plotdir, OscParamNameForOutput[int(gridParam)], RCsuffix, runversion, expt);
    }

    if (CriticalChi2File != NULL) {
        if (CriticalChi2LevelIndex >= 0) {
            sprintf(plotname+strlen(plotname), "_wCI_%d", CriticalChi2LevelIndex);
        }
        else {
            sprintf(plotname+strlen(plotname), "_wCI_all");
        }
        if (drawErrorBand) {
            sprintf(plotname+strlen(plotname), "_errorBand");
        }
    }

    for (int ic = 0; ic < nc; ic++) {
      cont[ic]->SetLineColor(colors[ic]);
      cont[ic]->Draw("C SAME");
      l1->AddEntry(cont[ic], title[ic], "L");
    }
    if (drawErrorBand && CriticalChi2LevelIndex >= 0) {
        for (int ic = 0; ic < nc; ic++) {
            string filename = Form("%s%s.tex", plotname, confSuffix[ic]);
            ofstream of(filename.c_str());
            of << CIlatex.at(ic);
            of.close();
            cout << "Wrote to " << filename << endl;
        }
    }
    if (!drawErrorBand) {
        l1->Draw();
        l2->Draw();
    }
    else if (CriticalChi2LevelIndex < 0) {
        l1->Draw();
        l2o->Draw();
        l2->Draw();
    }
    gPad->RedrawAxis(); // redraw axes since they can get hidden by the error bands
    c1->SaveAs(Form("%s.pdf", plotname));

    TFile *fout = new TFile(Form("%s.root", plotname), "RECREATE");
    for (int ic = 0; ic < nc; ic++) {
      cont[ic]->Write(Form("SA_%s%s", OscParamNameShort[int(gridParam)], confSuffix[ic]));
    }
    fout->Write();
    cout << "Wrote to " << fout->GetName() << endl;
}
