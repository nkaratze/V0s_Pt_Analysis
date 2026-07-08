#include "AnalysisLibrary.h"

// ============================================================
// GLOBAL VARIABLES
// ============================================================

AO2DConfig AO2DFiles;
TFile* F1 = nullptr;

// ============================================================
// INITIALIZATION OF FILE PATHS
// ============================================================

void InitAO2DFiles()
{
    AO2DFiles.datarootpath   = "PtSpectrumResults/Data/AnalysisResultsLHC22opass7Medium683739.root";
    AO2DFiles.mcrootpath     = "PtSpectrumResults/MC/AnalysisResultsLHC24f3d683254.root";
    AO2DFiles.mcgenrootpath  = "PtSpectrumResults/MC/AnalysisResultsLHC24f3d683254.root";

    AO2DFiles.hyperloopFiles = true;

    if (AO2DFiles.hyperloopFiles) {
        AO2DFiles.datarootsubfolder    = "v0-pt-inv-mass-plots";
        AO2DFiles.mcrootsubfolder      = "v0-pt-inv-mass-plots_id17569";
        AO2DFiles.mcgenrootsubfolder   = "v0-pt-inv-mass-plots_id18385";
    } else {
        AO2DFiles.datarootsubfolder    = "v0-pt-inv-mass-plots";
        AO2DFiles.mcrootsubfolder    = "v0-pt-inv-mass-plots";
        AO2DFiles.mcgenrootsubfolder = "v0-pt-inv-mass-plots";
    }
}

// ============================================================
// NAMESPACE DEFINITIONS
// ============================================================

namespace pthistos
{
    TH2D* PtBin[100];

    std::vector<std::string> ptbinsK0sh =
    {
        "0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9",
        "1.0","1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9",
        "2.0","2.2","2.4","2.6","2.8","3.0","3.3","3.6","3.9","4.2",
        "4.6","5.0","5.4","5.9","6.5","7.0","7.5","8.0","8.5","9.2",
        "10.0","11.0","12.0","13.5","15.0"
    };

    std::vector<std::string> ptbinslambda =
    {
        "0.5","0.6","0.7","0.8","0.9","1.0","1.1","1.2","1.3","1.4",
        "1.5","1.6","1.7","1.8","1.9","2.0","2.2","2.4","2.6","2.8",
        "3.0","3.2","3.4","3.6","3.8","4.0","4.4","4.9","5.5","6.0",
        "7.0","8.0","10.0"
    };

    std::vector<double> centBinsK0sh = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    std::vector<double> centBinsLambda = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
}

// ============================================================
// NORMALIZE HISTOGRAM BY BIN WIDTH
// ============================================================

void NormalizeByBinWidth(TH1* hist)
{
    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
    double binContent = hist->GetBinContent(i);
    double binError   = hist->GetBinError(i);
    double binWidth   = hist->GetBinWidth(i);

        if (binWidth > 0)
        {
            hist->SetBinContent(i, binContent / binWidth);
            hist->SetBinError(i, binError / binWidth);
        }
    }

}

// ============================================================
// READ AND SAVE HISTOGRAM
// ============================================================

TH1* ReadAndSaveHistogram(const std::string& readHistName, const std::string& savedHistName)
{
    if (TH1* oldHist = dynamic_cast<TH1*>(gDirectory->Get(savedHistName.c_str())))
    {delete oldHist;}

    TH1* hInput = dynamic_cast<TH1*>(gDirectory->Get(readHistName.c_str()));

    if (!hInput)
    {
        std::cerr << "ERROR: Histogram '" << readHistName << "' not found.\n";
        return nullptr;
    }

    TH1* hOutput =
        dynamic_cast<TH1*>(hInput->Clone(savedHistName.c_str()));

    if (!hOutput)
    {
        std::cerr << "ERROR: Failed to clone histogram '" << readHistName << "'\n";
        return nullptr;
    }
    hOutput->SetDirectory(nullptr);
    return hOutput;
}

// ============================================================
// READ HISTOGRAM PROJECTION
// ============================================================

TH1* ReadAndSaveHistogramProjection(const std::string& readHistName, const std::string& savedHistName, double cent_low, double cent_high)
{
    if (TH1* oldHist =
    dynamic_cast<TH1*>(gDirectory->Get(savedHistName.c_str())))
    {
    delete oldHist;
    }

    TObject* obj = gDirectory->Get(readHistName.c_str());

    if (!obj)
    {
        std::cerr << "Error: Histogram '"<< readHistName << "' not found.\n";
        return nullptr;
    }

    // ========================================================
    // TH1 CASE
    // ========================================================

    if (TH1* h1 = dynamic_cast<TH1*>(obj))
    {
        if (h1->GetDimension() == 1)
        {
            int nBins = h1->GetNbinsX();
            const TAxis* xAxis = h1->GetXaxis();
            std::vector<double> binEdges(nBins + 1);

            for (int i = 1; i <= nBins; ++i)
            {
                binEdges[i - 1] = xAxis->GetBinLowEdge(i);
            }

            binEdges[nBins] = xAxis->GetBinUpEdge(nBins);
            TH1D* hOut = new TH1D(savedHistName.c_str(), h1->GetTitle(), nBins, binEdges.data());

            for (int bin = 1; bin <= nBins; ++bin)
            {
                double c = h1->GetBinContent(bin);
                hOut->SetBinContent(bin, c);
                hOut->SetBinError(bin,std::sqrt(std::max(c, 0.0)));
                const char* label = xAxis->GetBinLabel(bin);

                if (label && strlen(label) > 0)
                {
                    hOut->GetXaxis()->SetBinLabel(bin, label);
                }
            }

            return hOut;
        }
    }

    // ========================================================
    // TH2 CASE
    // ========================================================

    if (TH2* h2 = dynamic_cast<TH2*>(obj))
    {
        TAxis* yAxis = h2->GetYaxis();
        int yBinLow = yAxis->FindFixBin(cent_low + 1e-8);
        int yBinHigh = yAxis->FindFixBin(cent_high - 1e-8);

        yBinLow  = std::max(1, yBinLow);
        yBinHigh = std::min(yAxis->GetNbins(), yBinHigh);
        TH1D* hProj = h2->ProjectionX(savedHistName.c_str(), yBinLow, yBinHigh);

        hProj->SetTitle(Form("%s (%.1f #leq cent < %.1f)",h2->GetTitle(), cent_low, cent_high));
        return hProj;
    }

    // ========================================================
    // TH3 CASE
    // ========================================================

    if (TH3* h3 = dynamic_cast<TH3*>(obj))
    {
        TAxis* zAxis = h3->GetZaxis();

        int zBinLow  = zAxis->FindFixBin(cent_low  + 1e-8);
        int zBinHigh = zAxis->FindFixBin(cent_high - 1e-8);

        zBinLow  = std::max(1, zBinLow);
        zBinHigh = std::min(zAxis->GetNbins(), zBinHigh);

        h3->GetZaxis()->SetRange(zBinLow, zBinHigh);

        TH2D* hProj = dynamic_cast<TH2D*>(
            h3->Project3D("yx")
        );

        if (!hProj)
        {
            std::cerr << "Error: TH3 projection failed.\n";
            return nullptr;
        }

        hProj->SetName(savedHistName.c_str());
        hProj->SetTitle(
            Form("%s (%.1f #leq cent < %.1f)",
                h3->GetTitle(),
                cent_low,
                cent_high));

        h3->GetZaxis()->SetRange(0, 0); // reset axis range

        return hProj;
    }

    std::cerr << "Error: unsupported histogram type.\n";
    return nullptr;
}

// ============================================================
// REBIN TH2D
// ============================================================

TH2D* RebinTH2D(TH2D* hInput, const char* newName, std::vector<double>* bins, bool rebinX)
{
    int nBinsX, nBinsY;
    const double* xEdges = nullptr;
    const double* yEdges = nullptr;

    if (rebinX)
    {
        nBinsX = bins->size() - 1;
        nBinsY = hInput->GetYaxis()->GetNbins();
        xEdges = bins->data();
    }
    else
    {
        nBinsX = hInput->GetXaxis()->GetNbins();
        nBinsY = bins->size() - 1;
        yEdges = bins->data();
    }

    TH2D* hOut = nullptr;

    if (rebinX)
    {
        hOut = new TH2D(newName, hInput->GetTitle(), nBinsX, xEdges, nBinsY, hInput->GetYaxis()->GetXmin(), hInput->GetYaxis()->GetXmax());
    }
    else
    {
        hOut = new TH2D( newName, hInput->GetTitle(), nBinsX, hInput->GetXaxis()->GetXmin(), hInput->GetXaxis()->GetXmax(), nBinsY, yEdges);
    }

    for (int i = 1; i <= hInput->GetNbinsX(); ++i)
    {
        double x = hInput->GetXaxis()->GetBinCenter(i);

        for (int j = 1; j <= hInput->GetNbinsY(); ++j)
        {
            double y   = hInput->GetYaxis()->GetBinCenter(j);
            double val = hInput->GetBinContent(i, j);
            hOut->Fill(x, y, val);
        }
    }
    return hOut;
}
// ============================================================
// LOAD HISTOGRAMS FROM CENTRALITY DIRECTORIES
// ============================================================
void LoadCentralityHistograms(TFile* file, const std::string& particle, const std::string& histName, const std::vector<double>& centBins, std::vector<TH1F*>& outList)
{
    outList.clear();

    for (size_t i = 0; i + 1 < centBins.size(); ++i) {
        double cent_low  = centBins[i];
        double cent_high = centBins[i + 1];

        TString centDirName = Form("%s/Cent_%.0f_%.0f", particle.c_str(), cent_low, cent_high);
        TDirectory* centDir = (TDirectory*)file->Get(centDirName);
        if (!centDir) {
            std::cout << "Directory " << centDirName << " not found!\n";
            continue;
        }

        TString fullHistName = Form("%s %s", particle.c_str(), histName.c_str());
        TH1F* h = (TH1F*)centDir->Get(fullHistName);
        if (!h) {
            std::cout << "Histogram " << histName << " not found in " << centDirName << "\n";
            continue;
        }

        TH1F* hClone = (TH1F*)h->Clone(Form("draw_%s_%.0f_%.0f", histName.c_str(), cent_low, cent_high));
        hClone->SetDirectory(0);
        outList.push_back(hClone);
    }
}
// ============================================================
// LOAD FUNCTIONS FROM CENTRALITY DIRECTORIES
// ============================================================
void LoadCentralityFunctions(TFile* file, const std::string& particle,const std::string& funcName, const std::vector<double>& centBins, std::vector<TF1*>& outList) 
{
    outList.clear();

    for (size_t i = 0; i + 1 < centBins.size(); ++i) {
        double cent_low  = centBins[i];
        double cent_high = centBins[i + 1];

        // Build the directory name for this centrality bin
        TString dirName = Form("%s/Cent_%.0f_%.0f", particle.c_str(), cent_low, cent_high);
        TDirectory* centDir = dynamic_cast<TDirectory*>(file->Get(dirName));
        if (!centDir) {
            std::cout << "Directory " << dirName << " not found!\n";
            continue;
        }

        // Get the object safely
        // Function Object Name
        TString objFitName = Form("%s %s", particle.c_str(), funcName.c_str());

      // Get the object safely
        TObject* obj = centDir->Get(objFitName);
        if (!obj) {
            std::cout << "Object " << objFitName << " not found in " << dirName << "\n";
            continue;
        }

        // Ensure it is a TF1
        TF1* f = dynamic_cast<TF1*>(obj);
        if (!f) {
            std::cout << "Object " << objFitName << " in " << dirName << " is not a TF1 (found: " << obj->ClassName() << ")\n";
            continue;
        }

        // Clone the TF1 safely
        TF1* fClone = dynamic_cast<TF1*>(f->Clone(Form("draw_%s_%.0f_%.0f", funcName.c_str(), cent_low, cent_high)));
        //fClone->SetDirectory(nullptr); // Prevent ROOT from auto-deleting

        // Store in output vector
        outList.push_back(fClone);
    }

    std::cout << "Loaded " << outList.size() << " TF1 objects for function " << funcName << "\n";
}
// ============================================================
// Compute Baryon to Meson Ratio L_AL / 2K0sh
// ============================================================
void ComputeCentralityRatios(TFile* file, const std::string& histName, const std::vector<double>& centBins, std::vector<TH1D*>& ratioList, TH1D*& hRatioSum)
{
    ratioList.clear();
    hRatioSum = nullptr;

    // ============================================================
    // Load histograms
    // ============================================================

    std::vector<TH1F*> hK_list;
    std::vector<TH1F*> hL_list;
    std::vector<TH1F*> hAL_list;

    LoadCentralityHistograms(file, "K0sh",       histName, centBins, hK_list);
    LoadCentralityHistograms(file, "Lambda",     histName, centBins, hL_list);
    LoadCentralityHistograms(file, "AntiLambda", histName, centBins, hAL_list);

    if (hK_list.size() != hL_list.size() || hK_list.size() != hAL_list.size()) {
        std::cout << "Mismatch in centrality bins!" << std::endl;
        return;
    }

    // ============================================================
    // Prepare summed histograms
    // ============================================================

    TH1D* hK_sum  = nullptr;
    TH1D* hL_sum  = nullptr;
    TH1D* hAL_sum = nullptr;

    // ============================================================
    // Loop over centrality classes
    // ============================================================

    for (size_t c = 0; c < hK_list.size(); ++c) {

        TH1D* hK  = (TH1D*)hK_list[c];
        TH1D* hL  = (TH1D*)hL_list[c];
        TH1D* hAL = (TH1D*)hAL_list[c];

        if (!hK || !hL || !hAL)
            continue;

        // --------------------------------------------------------
        // Create ratio histogram
        // --------------------------------------------------------

        TH1D* hRatio = (TH1D*)hK->Clone(Form("hRatio_cent_%zu", c));
        hRatio->Reset();

        // ========================================================
        // Compute:
        // (Lambda + AntiLambda) / (2*K0s)
        // ========================================================

        for (int i = 1; i <= hK->GetNbinsX(); ++i) {

            double k_low  = hK->GetBinLowEdge(i);
            double k_high = hK->GetBinLowEdge(i + 1);

            if (k_low < 0.5 || k_high > 3.0 + 1e-6)
                continue;

            // ----------------------------------------------------
            // Find matching Lambda bin
            // ----------------------------------------------------

            int jL = -1;

            for (int j = 1; j <= hL->GetNbinsX(); ++j) {

                double l_low  = hL->GetBinLowEdge(j);
                double l_high = hL->GetBinLowEdge(j + 1);

                if (std::abs(k_low  - l_low)  < 1e-6 && std::abs(k_high - l_high) < 1e-6) {

                    jL = j;
                    break;
                }
            }

            // ----------------------------------------------------
            // Find matching AntiLambda bin
            // ----------------------------------------------------

            int jAL = -1;

            for (int j = 1; j <= hAL->GetNbinsX(); ++j) {

                double al_low  = hAL->GetBinLowEdge(j);
                double al_high = hAL->GetBinLowEdge(j + 1);

                if (std::abs(k_low  - al_low)  < 1e-6 && std::abs(k_high - al_high) < 1e-6) {
                    jAL = j;
                    break;
                }
            }

            if (jL < 0 || jAL < 0)
                continue;

            // ----------------------------------------------------
            // Values
            // ----------------------------------------------------

            double k_val  = hK->GetBinContent(i);
            double l_val  = hL->GetBinContent(jL);
            double al_val = hAL->GetBinContent(jAL);

            double numerator   = l_val + al_val;
            double denominator = 2.0 * k_val;

            if (numerator <= 0 || denominator <= 0)
                continue;

            // ----------------------------------------------------
            // Ratio
            // ----------------------------------------------------

            double ratio = numerator / denominator;

            // ----------------------------------------------------
            // Errors
            // ----------------------------------------------------

            double k_err  = hK->GetBinError(i);
            double l_err  = hL->GetBinError(jL);
            double al_err = hAL->GetBinError(jAL);
            double numerator_err = std::sqrt(l_err*l_err + al_err*al_err);
            double denominator_err = 2.0 * k_err;
            double err = ratio * std::sqrt((numerator_err*numerator_err) / (numerator*numerator) + (denominator_err*denominator_err) / (denominator*denominator));

            hRatio->SetBinContent(i, ratio);
            hRatio->SetBinError(i, err);
        }

        ratioList.push_back(hRatio);

        // ========================================================
        // Initialize summed histograms
        // ========================================================

        if (!hK_sum) {

            hK_sum = (TH1D*)hK->Clone("hK_sum");
            hK_sum->Reset();

            hL_sum = (TH1D*)hL->Clone("hL_sum");
            hL_sum->Reset();

            hAL_sum = (TH1D*)hAL->Clone("hAL_sum");
            hAL_sum->Reset();
        }

        // ========================================================
        // Sum matching bins
        // ========================================================

        for (int i = 1; i <= hK->GetNbinsX(); ++i) {

            double k_low  = hK->GetBinLowEdge(i);
            double k_high = hK->GetBinLowEdge(i + 1);

            if (k_low < 0.5 || k_high > 3.0 + 1e-6)
                continue;

            int jL  = -1;
            int jAL = -1;

            // Match Lambda
            for (int j = 1; j <= hL->GetNbinsX(); ++j) {

                double l_low  = hL->GetBinLowEdge(j);
                double l_high = hL->GetBinLowEdge(j + 1);

                if (std::abs(k_low  - l_low)  < 1e-6 && std::abs(k_high - l_high) < 1e-6) {
                    jL = j;
                    break;
                }
            }

            // Match AntiLambda
            for (int j = 1; j <= hAL->GetNbinsX(); ++j) {

                double al_low  = hAL->GetBinLowEdge(j);
                double al_high = hAL->GetBinLowEdge(j + 1);

                if (std::abs(k_low  - al_low)  < 1e-6 &&
                    std::abs(k_high - al_high) < 1e-6) {

                    jAL = j;
                    break;
                }
            }

            if (jL < 0 || jAL < 0)
                continue;

            // ----------------------------------------------------
            // Sum contents
            // ----------------------------------------------------

            hK_sum->SetBinContent(i, hK_sum->GetBinContent(i) + hK->GetBinContent(i));
            hL_sum->SetBinContent(jL, hL_sum->GetBinContent(jL) + hL->GetBinContent(jL));
            hAL_sum->SetBinContent(jAL, hAL_sum->GetBinContent(jAL) + hAL->GetBinContent(jAL));

            // ----------------------------------------------------
            // Sum errors in quadrature
            // ----------------------------------------------------

            hK_sum->SetBinError(i, std::sqrt(std::pow(hK_sum->GetBinError(i), 2) + std::pow(hK->GetBinError(i), 2)));
            hL_sum->SetBinError(jL, std::sqrt(std::pow(hL_sum->GetBinError(jL), 2) + std::pow(hL->GetBinError(jL), 2)));
            hAL_sum->SetBinError(jAL, std::sqrt(std::pow(hAL_sum->GetBinError(jAL), 2) + std::pow(hAL->GetBinError(jAL), 2)));
        }
    }

    // ============================================================
    // Compute summed ratio
    // ============================================================

    if (!hK_sum || !hL_sum || !hAL_sum)
        return;

    hRatioSum = (TH1D*)hK_sum->Clone("hRatio_sum");
    hRatioSum->Reset();

    for (int i = 1; i <= hK_sum->GetNbinsX(); ++i) {

        double k_low  = hK_sum->GetBinLowEdge(i);
        double k_high = hK_sum->GetBinLowEdge(i + 1);

        if (k_low < 0.5 || k_high > 3.0 + 1e-6)
            continue;

        int jL  = -1;
        int jAL = -1;

        // Match Lambda
        for (int j = 1; j <= hL_sum->GetNbinsX(); ++j) {

            double l_low  = hL_sum->GetBinLowEdge(j);
            double l_high = hL_sum->GetBinLowEdge(j + 1);

            if (std::abs(k_low  - l_low)  < 1e-6 && std::abs(k_high - l_high) < 1e-6) {
                jL = j;
                break;
            }
        }

        // Match AntiLambda
        for (int j = 1; j <= hAL_sum->GetNbinsX(); ++j) {

            double al_low  = hAL_sum->GetBinLowEdge(j);
            double al_high = hAL_sum->GetBinLowEdge(j + 1);

            if (std::abs(k_low  - al_low)  < 1e-6 && std::abs(k_high - al_high) < 1e-6) {
                jAL = j;
                break;
            }
        }

        if (jL < 0 || jAL < 0)
            continue;

        double k_val  = hK_sum->GetBinContent(i);
        double l_val  = hL_sum->GetBinContent(jL);
        double al_val = hAL_sum->GetBinContent(jAL);
        double numerator   = l_val + al_val;
        double denominator = 2.0 * k_val;

        if (numerator <= 0 || denominator <= 0)
            continue;

        // --------------------------------------------------------
        // Ratio
        // --------------------------------------------------------

        double ratio = numerator / denominator;

        // --------------------------------------------------------
        // Errors
        // --------------------------------------------------------

        double k_err  = hK_sum->GetBinError(i);
        double l_err  = hL_sum->GetBinError(jL);
        double al_err = hAL_sum->GetBinError(jAL);
        double numerator_err = std::sqrt(l_err*l_err + al_err*al_err);
        double denominator_err = 2.0 * k_err;
        double err = ratio * std::sqrt((numerator_err*numerator_err) / (numerator*numerator) + (denominator_err*denominator_err) / (denominator*denominator));
        hRatioSum->SetBinContent(i, ratio);
        hRatioSum->SetBinError(i, err);
    }
}

// ============================================================
// DRAWING FUNCTIONS
// ============================================================

TCanvas* DrawCentralityHistoSet( const std::vector<TH1F*>& hList, const std::vector<TString>& roman, const std::vector<Color_t>& colors, const std::vector<double>& centBins, const char* canvasName, bool logy, bool fixYAxis, bool markerline, bool drawLevyTsallis)
{
    if(hList.empty()) {
        std::cout << "No histograms to draw for " << canvasName << "\n";
        return nullptr;
    }

    if(hList.size() != centBins.size() - 1) {
        std::cout << "Mismatch: hList has " << hList.size()
        << " entries but centBins defines "
        << (centBins.size() - 1) << " intervals.\n";
        return nullptr;
    }

    TCanvas* c = new TCanvas(canvasName, canvasName, 800, 600);
    if(logy) c->SetLogy();

    TLegend* leg = new TLegend(0.65, 0.5, 0.88, 0.88);

    // ----- Sum histogram (first set) -----
    TH1F* hSum = (TH1F*)hList[0]->Clone(Form("Sum_%s", canvasName));
    hSum->Reset();
    hSum->SetDirectory(0);

    // Find global maximum
    double globalMax = 0.0;

    for(size_t i = 0; i < hList.size(); ++i){

        double cent_low  = centBins[i];
        double cent_high = centBins[i + 1];

        hList[i]->SetLineColor(colors[i % colors.size()]);
        hList[i]->SetMarkerColor(colors[i % colors.size()]);
        hList[i]->SetMarkerSize(1.2);

        double hmax = hList[i]->GetMaximum();
        if(hmax > globalMax)
        globalMax = hmax;

        if(markerline){
            hList[i]->Draw(i == 0 ? "PEL" : "PEL SAME");
        } 
        else {
            hList[i]->Draw(i == 0 ? "PEL" : "PEL SAME");
        }

        TH1F* hClone = (TH1F*)hList[i]->Clone(Form("%s_clone_%zu", hList[i]->GetName(), i));
        hClone->SetDirectory(nullptr);
        hClone->Scale((cent_high-cent_low)/100.);
        hSum->Add(hClone);
        delete hClone;

        TString label = Form("%s (%.0f%% - %.0f%%)",roman[i % roman.size()].Data(),cent_low,cent_high);
        leg->AddEntry(hList[i], label, "lep");
    }

    // Set y-axis range
    if(fixYAxis) {
        hList[0]->GetYaxis()->SetRangeUser(1e-6, 1.0);
    } 
    else {
        hList[0]->GetYaxis()->SetRangeUser(1e-6, globalMax * 1.2);
    }

    hSum->SetMarkerStyle(kFullDiamond);
    hSum->SetMarkerSize(1.5);
    hSum->SetMarkerColor(kBlack);
    hSum->SetLineColor(kBlack);
    hSum->Draw("PEL SAME");

    leg->AddEntry(hSum, "Sum (0%% - 100%%)", "lep");
    leg->Draw();

    c->Modified();
    c->Update();

    return c;
}

TCanvas* DrawCentralityFunctionSet(const std::vector<TF1*>& fList,const std::vector<TString>& roman,const std::vector<Color_t>& colors,const std::vector<double>& centBins,const char* canvasName, bool logy,bool fixYAxis)
{
    if(fList.empty()) {
        std::cout << "No functions to draw for " << canvasName << "\n";
        return nullptr;
    }

    if(fList.size() != centBins.size() - 1) {
        std::cout << "Mismatch: fList has " << fList.size()
        << " entries but centBins defines "
        << (centBins.size() - 1) << " intervals.\n";
        return nullptr;
    }

    TCanvas* c = new TCanvas(canvasName, canvasName, 800, 600);
    if(logy) c->SetLogy();

    TLegend* leg = new TLegend(0.65, 0.5, 0.88, 0.88);

    for(size_t i = 0; i < fList.size(); ++i){
        TF1* f = fList[i];
        if(!f) continue;

        double cent_low  = centBins[i];
        double cent_high = centBins[i + 1];

        f->SetLineColor(colors[i % colors.size()]);
        f->SetLineWidth(2);

        // Draw
        f->Draw(i == 0 ? "L" : "L SAME");

        // Optional Y-axis control (only meaningful for first drawn object)
        if(i == 0 && fixYAxis) {
            f->SetMinimum(1e-10);  // tweak as needed
            f->SetMaximum(1.0);
        }

        TString label = Form("%s (%.0f%% - %.0f%%)",
        roman[i % roman.size()].Data(),
        cent_low, cent_high);

        leg->AddEntry(f, label, "l");
    }

    leg->Draw();

    return c;
}

// ============================================================
// DRAWING FUNCTIONS FOR LAMBDA VS ANTI-LAMBDA COMPARISON
// ============================================================


TCanvas* DrawLambdaAntiLambdaComparison(const std::vector<TH1F*>& hList1 /*e.g. Lambda*/, const std::vector<TH1F*>& hList2 /*e.g. AntiLambda*/, const std::vector<TString>& roman, const std::vector<Color_t>& colors, const std::vector<double>& centBins, const char* canvasName, bool logy, bool fixYAxis)
{
    if (hList1.empty() || hList2.empty()) {
        std::cout << "Empty histogram list in " << canvasName << "\n";
        return nullptr;
    }

    if (hList1.size() != hList2.size()) {
        std::cout << "Mismatch between the two histogram sets!\n";
        return nullptr;
    }

    // =========================================================
    // Canvas + pads
    // =========================================================
    TCanvas* c = new TCanvas(canvasName, canvasName, 1200, 800);

    // LEFT COLUMN
    TPad* pad1 = new TPad("pad1", "main left", 0.0, 0.3, 0.5, 1.0);
    TPad* pad2 = new TPad("pad2", "ratio left", 0.0, 0.0, 0.5, 0.3);

    // RIGHT COLUMN
    TPad* pad3 = new TPad("pad3", "main right", 0.5, 0.3, 1.0, 1.0);
    TPad* pad4 = new TPad("pad4", "ratio right", 0.5, 0.0, 1.0, 0.3);

    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();

    // =========================================================
    // LEFT TOP PAD : Individual histograms
    // =========================================================
    pad1->cd();

    if (logy)
    pad1->SetLogy();

    pad1->SetBottomMargin(0.01);

    gStyle->SetOptStat(0);

    TLegend* leg = new TLegend(0.60, 0.50, 0.88, 0.88);
    leg->SetBorderSize(0);

    // =========================================================
    // SUM HISTOGRAMS
    // =========================================================
    TH1F* hSum1 = (TH1F*)hList1[0]->Clone(Form("Sum1_%s", canvasName));
    TH1F* hSum2 = (TH1F*)hList2[0]->Clone(Form("Sum2_%s", canvasName));

    hSum1->Reset();
    hSum2->Reset();

    hSum1->SetDirectory(0);
    hSum2->SetDirectory(0);

    // =========================================================
    // DRAW INDIVIDUAL HISTOGRAMS
    // =========================================================
    for (size_t i = 0; i < hList1.size(); ++i) {

        if (!hList1[i] || !hList2[i]) {
            std::cout << "Null histogram at index " << i << "\n";
            continue;
        }

        if (hList1[i]->GetNbinsX() != hList2[i]->GetNbinsX()) {
            std::cout << "Binning mismatch at index " << i << "\n";
            continue;
        }

        double cent_low  = centBins[i];
        double cent_high = centBins[i + 1];

        Color_t col = colors[i % colors.size()];

        // =========================
        // Styling
        // =========================
        hList1[i]->SetMarkerStyle(kFullCircle);
        hList1[i]->SetMarkerColor(col);
        hList1[i]->SetLineColor(col);

        hList2[i]->SetMarkerStyle(kOpenCircle);
        hList2[i]->SetMarkerColor(col);
        hList2[i]->SetLineColor(col);

        hList1[i]->Draw(i == 0 ? "PE" : "PE SAME");
        hList2[i]->Draw("PE SAME");

        if (i == 0 && fixYAxis)
        hList1[i]->GetYaxis()->SetRangeUser(0.0, 1.0);

        // =========================
        // Safe clones for sums
        // =========================
        TH1F* h1 = (TH1F*)hList1[i]->Clone(Form("tmp_h1_%zu_%s", i, canvasName));
        TH1F* h2 = (TH1F*)hList2[i]->Clone(Form("tmp_h2_%zu_%s", i, canvasName));

        h1->Scale((cent_high - cent_low) / 100.);
        h2->Scale((cent_high - cent_low) / 100.);

        hSum1->Add(h1);
        hSum2->Add(h2);

        delete h1;
        delete h2;

        TString label = Form("%s (%.0f-%.0f%%)",
        roman[i].Data(),
        cent_low,
        cent_high);

        leg->AddEntry(hList1[i], label + " #Lambda", "lep");
        leg->AddEntry(hList2[i], label + " #bar{#Lambda}", "lep");
    }

    leg->Draw();

    // =========================================================
    // RIGHT TOP PAD : Summed histograms only
    // =========================================================
    pad3->cd();

    if (logy)
    pad3->SetLogy();

    pad3->SetBottomMargin(0.01);

    hSum1->SetMarkerStyle(kFullDiamond);
    hSum1->SetMarkerColor(kBlack);
    hSum1->SetLineColor(kBlack);

    hSum2->SetMarkerStyle(kOpenDiamond);
    hSum2->SetMarkerColor(kBlack);
    hSum2->SetLineColor(kBlack);

    hSum1->Draw("PE");
    hSum2->Draw("PE SAME");

    if (fixYAxis)
    hSum1->GetYaxis()->SetRangeUser(0.0, 1.0);

    TLegend* legSum = new TLegend(0.50, 0.75, 0.88, 0.88);
    legSum->SetBorderSize(0);

    legSum->AddEntry(hSum1, "Sum #Lambda", "lep");
    legSum->AddEntry(hSum2, "Sum #bar{#Lambda}", "lep");

    legSum->Draw();

    // =========================================================
    // LEFT BOTTOM PAD : Individual ratios
    // =========================================================
    pad2->cd();

    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.3);
    gPad->SetGridy();

    bool firstDraw = true;

    for (size_t i = 0; i < hList1.size(); ++i) {

        TH1F* hRatio =
        (TH1F*)hList1[i]->Clone(Form("hRatio_cent%zu", i));

        hRatio->Divide(hList2[i]);

        hRatio->SetMarkerStyle(kFullCircle);
        hRatio->SetMarkerColor(colors[i % colors.size()]);
        hRatio->SetLineColor(colors[i % colors.size()]);

        hRatio->GetYaxis()->SetTitle("#Lambda / #bar{#Lambda}");
        hRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");

        hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);

        hRatio->Draw(firstDraw ? "PE" : "PE SAME");

        firstDraw = false;
    }

    // Unity line for left ratio pad
    TLine* line1 = new TLine(
    hList1[0]->GetXaxis()->GetXmin(), 1.0,
    hList1[0]->GetXaxis()->GetXmax(), 1.0);

    line1->SetLineStyle(2);
    line1->Draw("SAME");

    // =========================================================
    // RIGHT BOTTOM PAD : Summed ratio only
    // =========================================================
    pad4->cd();

    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.3);
    gPad->SetGridy();

    TH1F* hSumRatio = (TH1F*)hSum1->Clone("hSumRatio");

    hSumRatio->Divide(hSum2);

    hSumRatio->SetMarkerStyle(kFullDiamond);
    hSumRatio->SetMarkerColor(kBlack);
    hSumRatio->SetLineColor(kBlack);

    hSumRatio->GetYaxis()->SetTitle("#Lambda / #bar{#Lambda}");
    hSumRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    hSumRatio->GetYaxis()->SetRangeUser(0.5, 1.5);

    hSumRatio->Draw("PE");

    // Unity line for summed ratio
    TLine* line2 = new TLine(
    hSumRatio->GetXaxis()->GetXmin(), 1.0,
    hSumRatio->GetXaxis()->GetXmax(), 1.0);

    line2->SetLineStyle(2);
    line2->Draw("SAME");

    c->Update();

    return c;
}

// ============================================================
// STORAGE HELPERS
// ============================================================

std::vector<TH1F*> CloneStorage(const std::vector<TH1F*>& input, const TString& tag)
{
    std::vector<TH1F*> output;

    for (size_t i = 0; i < input.size(); ++i) {
        TH1F* hist = input[i];
        if (!hist) continue;

        TH1F* hClone = (TH1F*)hist->Clone(Form("%s_%zu", tag.Data(), i));
        hClone->SetDirectory(nullptr);
        output.push_back(hClone);
    }

    return output;
}