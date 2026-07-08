#ifndef ANALYSISLIBRARY_H
#define ANALYSISLIBRARY_H

// ============================================================
// ROOT Includes
// ============================================================

//#if !defined(**CLING**) || defined(**ROOTCLING**)
#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TGraph.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TProfile.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#endif

#include <fmt/format.h>

// ============================================================
// GLOBAL NAMESPACE
// ============================================================

namespace pthistos
{
extern TH2D* PtBin[100];
extern std::vector<std::string> ptbinsK0sh;
extern std::vector<std::string> ptbinslambda;
extern std::vector<double> centBinsK0sh;
extern std::vector<double> centBinsLambda;
}

// ============================================================
// AO2D CONFIGURATION
// ============================================================

struct AO2DConfig
{
std::string datarootpath;
std::string mcrootpath;
std::string mcgenrootpath;
std::string datarootsubfolder;
std::string mcrootsubfolder;
std::string mcgenrootsubfolder;
bool hyperloopFiles;
};

// ============================================================
// AO2D INITIALIZATION FILES
// ============================================================
void InitAO2DFiles();

// ============================================================
// GLOBAL OBJECTS
// ============================================================

extern AO2DConfig AO2DFiles;
extern TFile* F1;

// ============================================================
// UTILITY FUNCTIONS
// ============================================================

void NormalizeByBinWidth(TH1* hist);
TH1* ReadAndSaveHistogram(const std::string& readHistName, const std::string& savedHistName);
TH1* ReadAndSaveHistogramProjection(const std::string& readHistName, const std::string& savedHistName, double cent_low, double cent_high);
TH2D* RebinTH2D(TH2D* hInput,const char* newName, std::vector<double>* bins, bool rebinX = true);

// ============================================================
// HISTOGRAM LOADERS
// ============================================================

void LoadCentralityHistograms(TFile* file, const std::string& particle, const std::string& histName, const std::vector<double>& centBins, std::vector<TH1F*>& outList);
void LoadCentralityFunctions(TFile* file, const std::string& particle, const std::string& funcName, const std::vector<double>& centBins, std::vector<TF1*>& outList);

// ============================================================
// RATIO COMPUTATION
// ============================================================

void ComputeCentralityRatios(TFile* file, const std::string& histName, const std::vector<double>& centBins, std::vector<TH1D*>& ratioList, TH1D*& hRatioSum);

// ============================================================
// DRAWING FUNCTIONS
// ============================================================

TCanvas* DrawCentralityHistoSet( const std::vector<TH1F*>& hList, const std::vector<TString>& roman, const std::vector<Color_t>& colors, const std::vector<double>& centBins, const char* canvasName, bool logy = false, bool fixYAxis = false, bool markerline = false, bool drawLevyTsallis = false);
TCanvas* DrawCentralityFunctionSet(const std::vector<TF1*>& fList, const std::vector<TString>& roman, const std::vector<Color_t>& colors, const std::vector<double>& centBins, const char* canvasName, bool logy = false, bool fixYAxis = false);
TCanvas* DrawLambdaAntiLambdaComparison( const std::vector<TH1F*>& hList1, const std::vector<TH1F*>& hList2, const std::vector<TString>& roman, const std::vector<Color_t>& colors, const std::vector<double>& centBins, const char* canvasName, bool logy = false, bool fixYAxis = false);

// ============================================================
// STORAGE HELPERS
// ============================================================

std::vector<TH1F*> CloneStorage( const std::vector<TH1F*>& input,const TString& tag);

#endif
