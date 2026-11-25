#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TObject.h>
#include <TIterator.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

struct SystematicVariation {
    std::string file;
    std::string label;
};

struct ParticleConfig {
    std::string name;
    std::string dirName;
    std::string histNameCent;
    std::string overallDir;
    std::string canvasNameAll;
    std::string histNameOverall;
    std::string canvasNameMCChi2;
    std::string histNameMCChi2;
    std::string refFitLabel;
    std::vector<SystematicVariation> fitVariations;
    std::vector<SystematicVariation> cutVariations;
    std::vector<SystematicVariation> feeddownVariations;
};

// Helper function to find histograms inside a TCanvas/TPad
TH1D* FindHistogramInCanvas(TCanvas* c, const std::string& histName) {
    if (!c) return nullptr;

    TIter next(c->GetListOfPrimitives());
    TObject* obj = nullptr;
    while ((obj = next())) {
        if (obj->InheritsFrom(TH1::Class())) {
            TH1* h = (TH1*)obj;
            if (std::string(h->GetName()) == histName) {
                return dynamic_cast<TH1D*>(h);
            }
        }
        if (obj->InheritsFrom("TPad")) {
            TPad* pad = (TPad*)obj;
            TIter padNext(pad->GetListOfPrimitives());
            TObject* padObj = nullptr;
            while ((padObj = padNext())) {
                if (padObj->InheritsFrom(TH1::Class())) {
                    TH1* h = (TH1*)padObj;
                    if (std::string(h->GetName()) == histName) {
                        return dynamic_cast<TH1D*>(h);
                    }
                }
            }
        }
    }
    return nullptr;
}

// Enhanced Barlow Criterion Helper Function returning the numerical criterion value
double CalculateBarlowCriterionValue(const std::string& particle, const std::string& context, const std::string& varLabel,
                                     double valRef, double errRef, double valVar, double errVar, double sigmaBarlow = 2.0, bool* passesOut = nullptr) {
    double delta = std::abs(valRef - valVar);
    double sigmaSqDiff = std::abs((errRef * errRef) - (errVar * errVar));

    double criterion = 0.0;
    bool passes = false;

    if (sigmaSqDiff == 0.0) {
        passes = (delta != 0.0);
        criterion = passes ? 999.0 : 0.0;
        if (passesOut) *passesOut = passes;
        return criterion;
    }

    double sigmaDelta = std::sqrt(sigmaSqDiff);
    criterion = delta / sigmaDelta;
    passes = (criterion >= sigmaBarlow);

    if (passesOut) *passesOut = passes;
    return criterion;
}

void CalculateSystematics()
{
    // ============================================================
    // USER INPUT
    // ============================================================

    //const std::string finalFile = "V0s_Pt_Spectrum_Cent_5.36_TeVLHC24ppref.root";
    //const std::string finalFile = "V0s_Pt_Spectrum_Cent_5.36_TeV.root";
    const std::string finalFile = "PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2.root";
    std::vector<std::string> centralities = {"Cent_0_1", "Cent_1_5", "Cent_5_10","Cent_10_15", "Cent_15_20","Cent_20_30", "Cent_30_40", "Cent_40_50","Cent_50_70", "Cent_70_100"};

    // Material budget systematic value (1% quoted)
    const double materialBudgetFrac = 0.02;

    // ============================================================
    // PARTICLE-SPECIFIC VARIATION LISTS
    // ============================================================

    std::vector<SystematicVariation> k0shFitVariations = {
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol1.root", "Double-Gaussian + pol1"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2.root", "Double-Gaussian + pol2"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_GaussianPol1.root", "Gaussian + pol1"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_GaussianPol2.root", "Gaussian + pol2"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_VoigtPol1.root", "Voigt + pol1"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_VoigtPol2.root", "Voigt + pol2"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_BreitWignerPol1.root", "Breit-Wigner + pol1"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_BreitWignerPol2.root", "Breit-Wigner + pol2"}
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2Range09.root", "Double-Gaussian + pol2 + 0.9*Range"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2Range11.root", "Double-Gaussian + pol2 + 1.1*Range"}
    };

    std::vector<SystematicVariation> lambdaFitVariations = {
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol1.root", "Double-Gaussian + pol1"},
       // {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2.root", "Double-Gaussian + pol2"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_GaussianPol1.root", "Gaussian + pol1"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_GaussianPol2.root", "Gaussian + pol2"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_VoigtPol1.root", "Voigt + pol1"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_VoigtPol2.root", "Voigt + pol2"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_BreitWignerPol1.root", "Breit-Wigner + pol1"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_BreitWignerPol2.root", "Breit-Wigner + pol2"}
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2Range09.root", "Double-Gaussian + pol2 + 0.9*Range"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2Range11.root", "Double-Gaussian + pol2 + 1.1*Range"}
    };

    std::vector<SystematicVariation> antiLambdaFitVariations = {
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol1.root", "Double-Gaussian + pol1"},        
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2.root", "Double-Gaussian + pol2"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_GaussianPol1.root", "Gaussian + pol1"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_GaussianPol2.root", "Gaussian + pol2"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_VoigtPol1.root", "Voigt + pol1"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_VoigtPol2.root", "Voigt + pol2"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_BreitWignerPol1.root", "Breit-Wigner + pol1"},
        //{"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_BreitWignerPol2.root", "Breit-Wigner + pol2"}
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2Range09.root", "Double-Gaussian + pol2 + 0.9*Range"},
        {"PtSpectrumResults/Data/Systematics/V0s_Pt_Spesctrum_5.36_TeV_DoubleGaussianPol2Range11.root", "Double-Gaussian + pol2 + 1.1*Range"}
    };

    std::vector<ParticleConfig> particles = {
        {
            "K0sh", "K0sh", "K0sh Corrected Pt Spectrum",
            "K0sh/All_Centralities",
            "Corrected_PtSpectra_AllCentralities_K0sh",
            "Sum_Corrected_PtSpectra_AllCentralities_K0sh",
            "Reduced_Chi2_MC_AllCentralities_K0sh",
            "Sum_Reduced_Chi2_MC_AllCentralities_K0sh",
            "Double-Gaussian + pol1",
            k0shFitVariations, {}, {}
        },
        {
            "Lambda", "Lambda", "Lambda Corrected Pt Spectrum",
            "Lambda/All_Centralities",
            "Corrected_PtSpectra_AllCentralities_Lambda",
            "Sum_Corrected_PtSpectra_AllCentralities_Lambda",
            "Reduced_Chi2_MC_AllCentralities_Lambda",
            "Sum_Reduced_Chi2_MC_AllCentralities_Lambda",
            "Double-Gaussian + pol1",
            lambdaFitVariations, {}, {}
        },
        {
            "AntiLambda", "AntiLambda", "AntiLambda Corrected Pt Spectrum",
            "AntiLambda/All_Centralities",
            "Corrected_PtSpectra_AllCentralities_AntiLambda",
            "Sum_Corrected_PtSpectra_AllCentralities_AntiLambda",
            "Reduced_Chi2_MC_AllCentralities_AntiLambda",
            "Sum_Reduced_Chi2_MC_AllCentralities_AntiLambda",
            "Double-Gaussian + pol1",
            antiLambdaFitVariations, {}, {}
        }
    };

    // ============================================================
    // OPEN FILES & PREPARE OUTPUT ROOT FILE
    // ============================================================

    TFile* fFinal = TFile::Open(finalFile.c_str());
    if (!fFinal || fFinal->IsZombie()) return;

    TFile* fOut = TFile::Open("SystematicVariations.root", "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        fFinal->Close();
        return;
    }

    std::vector<Color_t> colors = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1,kCyan+2, kViolet+1, kAzure+2, kPink+1, kTeal+3, kYellow+2};

    // ============================================================
    // LOOP OVER PARTICLE SPECIES
    // ============================================================

    for (const auto& part : particles) {
        std::cout << "\n\n############################################################\n";
        std::cout << "PROCESSING PARTICLE: " << part.name << "\n";
        std::cout << "############################################################\n";

        fOut->mkdir(part.dirName.c_str());
        fOut->cd(part.dirName.c_str());
        TDirectory* dirPart = gDirectory;

        // ============================================================
        // 1. PROCESS INDIVIDUAL CENTRALITIES
        // ============================================================

        for (const auto& cent : centralities) {
            std::string centPath = part.dirName + "/" + cent + "/" + part.histNameCent;
            TH1* hFinalCent = dynamic_cast<TH1*>(fFinal->Get(centPath.c_str()));
            if (!hFinalCent) continue;

            dirPart->mkdir(cent.c_str());
            dirPart->cd(cent.c_str());
            TDirectory* dirCent = gDirectory;

            TH1* hRefCentClone = (TH1*)hFinalCent->Clone(Form("hRef_%s_%s", part.name.c_str(), cent.c_str()));
            hRefCentClone->SetDirectory(dirCent);
            hRefCentClone->Write();

            TH1* hSystCent = (TH1*)hFinalCent->Clone(Form("hSystematic_%s_%s", part.name.c_str(), cent.c_str()));
            hSystCent->Reset();
            hSystCent->SetDirectory(dirCent);

            auto processVariations = [&](const std::vector<SystematicVariation>& varList, const std::string& subDirPrefix, std::vector<TH1*>& outVarHists, std::vector<TH1*>& outRatios) {
                if (varList.empty()) return;
                dirCent->mkdir(subDirPrefix.c_str());
                dirCent->cd(subDirPrefix.c_str());
                TDirectory* dirSub = gDirectory;

                for (size_t i = 0; i < varList.size(); ++i) {
                    TFile* fVar = TFile::Open(varList[i].file.c_str());
                    if (!fVar || fVar->IsZombie()) continue;

                    TH1* hVarCent = dynamic_cast<TH1*>(fVar->Get(centPath.c_str()));
                    if (!hVarCent) {
                        fVar->Close();
                        continue;
                    }

                    dirSub->cd();
                    TH1* hVarClone = (TH1*)hVarCent->Clone(Form("hVar_%s_%zu_%s_%s", subDirPrefix.c_str(), i, part.name.c_str(), cent.c_str()));
                    hVarClone->Write();
                    outVarHists.push_back(hVarClone);

                    TH1* hRatioCent = (TH1*)hVarCent->Clone(Form("hRatio_%s_%zu_%s_%s", subDirPrefix.c_str(), i, part.name.c_str(), cent.c_str()));
                    hRatioCent->Divide(hFinalCent);
                    hRatioCent->Write();
                    outRatios.push_back(hRatioCent);

                    fVar->Close();
                }
                dirCent->cd();
            };

            std::vector<TH1*> centFitVarHists, centFitRatios;
            std::vector<TH1*> centCutVarHists, centCutRatios;
            std::vector<TH1*> centFeedVarHists, centFeedRatios;

            processVariations(part.fitVariations, "FitVars", centFitVarHists, centFitRatios);
            processVariations(part.cutVariations, "CutVars", centCutVarHists, centCutRatios);
            processVariations(part.feeddownVariations, "FeedVars", centFeedVarHists, centFeedRatios);

            auto saveBarlows = [&](const std::vector<SystematicVariation>& varList, const std::vector<TH1*>& varHists, const std::string& prefix) {
                for (size_t i = 0; i < varHists.size(); ++i) {
                    dirCent->cd();
                    TH1* hBarlowCent = dynamic_cast<TH1*>(hFinalCent->Clone(Form("hBarlowCriterion_%s_%zu_%s_%s", prefix.c_str(), i, part.name.c_str(), cent.c_str())));
                    hBarlowCent->Reset();

                    for (int bin = 1; bin <= hBarlowCent->GetNbinsX(); ++bin) {
                        double valRef = hFinalCent->GetBinContent(bin);
                        double errRef = hFinalCent->GetBinError(bin);
                        double valVar = varHists[i]->GetBinContent(bin);
                        double errVar = varHists[i]->GetBinError(bin);

                        bool passes = false;
                        double critVal = CalculateBarlowCriterionValue(part.name, cent, varList[i].label, valRef, errRef, valVar, errVar, 2.0, &passes);
                        hBarlowCent->SetBinContent(bin, critVal);
                    }

                    hBarlowCent->Write();
                    delete hBarlowCent;
                }
            };

            saveBarlows(part.fitVariations, centFitVarHists, "Fit");
            saveBarlows(part.cutVariations, centCutVarHists, "Cut");
            saveBarlows(part.feeddownVariations, centFeedVarHists, "Feed");

            for (int bin = 1; bin <= hSystCent->GetNbinsX(); ++bin) {
                double sumSquared = materialBudgetFrac * materialBudgetFrac;
                double valRef = hFinalCent->GetBinContent(bin);
                double errRef = hFinalCent->GetBinError(bin);

                auto accumulateVariations = [&](const std::vector<SystematicVariation>& varList, const std::vector<TH1*>& varHists, const std::vector<TH1*>& ratioHists) {
                    for (size_t i = 0; i < ratioHists.size(); ++i) {
                        double valVar = varHists[i]->GetBinContent(bin);
                        double errVar = varHists[i]->GetBinError(bin);

                        bool passes = false;
                        CalculateBarlowCriterionValue(part.name, cent, varList[i].label, valRef, errRef, valVar, errVar, 2.0, &passes);

                        if (passes) {
                            double ratio = ratioHists[i]->GetBinContent(bin);
                            double systematic = std::abs(ratio - 1.0);
                            sumSquared += systematic * systematic;
                        }
                    }
                };

                accumulateVariations(part.fitVariations, centFitVarHists, centFitRatios);
                accumulateVariations(part.cutVariations, centCutVarHists, centCutRatios);
                accumulateVariations(part.feeddownVariations, centFeedVarHists, centFeedRatios);

                hSystCent->SetBinContent(bin, std::sqrt(sumSquared) * 100.0);
            }

            dirCent->cd();
            hSystCent->Write();

            for (auto* h : centFitVarHists) delete h;
            for (auto* h : centFitRatios) delete h;
            for (auto* h : centCutVarHists) delete h;
            for (auto* h : centCutRatios) delete h;
            for (auto* h : centFeedVarHists) delete h;
            for (auto* h : centFeedRatios) delete h;
        }

        // ============================================================
        // 2. PROCESS OVERALL SPECTRUM & INDIVIDUAL SYSTEMATICS
        // ============================================================

        dirPart->cd();
        TCanvas* cRefCanvas = dynamic_cast<TCanvas*>(fFinal->Get((part.overallDir + "/" + part.canvasNameAll).c_str()));
        TH1* hFinalOverall = nullptr;
        if (cRefCanvas) {
            TH1* hTemp = FindHistogramInCanvas(cRefCanvas, part.histNameOverall);
            if (hTemp) {
                hFinalOverall = dynamic_cast<TH1*>(hTemp->Clone(Form("hFinalOverall_%s", part.name.c_str())));
                hFinalOverall->SetDirectory(nullptr);
            }
        }

        if (!hFinalOverall) continue;

        auto loadOverallVariations = [&](const std::vector<SystematicVariation>& varList, std::vector<TH1*>& outVar, std::vector<TH1*>& outRatio, std::vector<std::string>& outLabels, bool appendRefLabel) {
            for (size_t i = 0; i < varList.size(); ++i) {
                TFile* fVar = TFile::Open(varList[i].file.c_str());
                if (!fVar || fVar->IsZombie()) continue;

                TCanvas* cVarCanvas = dynamic_cast<TCanvas*>(fVar->Get((part.overallDir + "/" + part.canvasNameAll).c_str()));
                TH1* hVar = nullptr;
                if (cVarCanvas) {
                    TH1* hTempVar = FindHistogramInCanvas(cVarCanvas, part.histNameOverall);
                    if (hTempVar) {
                        hVar = dynamic_cast<TH1*>(hTempVar->Clone(Form("hVariation_%s_%zu_%s", varList[i].label.c_str(), i, part.name.c_str())));
                        hVar->SetDirectory(nullptr);
                    }
                }

                if (!hVar) {
                    fVar->Close();
                    continue;
                }

                TH1* hRatio = dynamic_cast<TH1*>(hVar->Clone(Form("hRatio_%s_%zu_%s", varList[i].label.c_str(), i, part.name.c_str())));
                hRatio->SetDirectory(nullptr);
                hRatio->Divide(hFinalOverall);

                outVar.push_back(hVar);
                outRatio.push_back(hRatio);

                std::string currentLabel = varList[i].label;
                if (appendRefLabel && currentLabel == part.refFitLabel) {
                    currentLabel += " (reference)";
                }
                outLabels.push_back(currentLabel);

                fVar->Close();
            }
        };

        std::vector<TH1*> hFitVars, hFitRatios;
        std::vector<std::string> fitLabels;
        loadOverallVariations(part.fitVariations, hFitVars, hFitRatios, fitLabels, true);

        std::vector<TH1*> hCutVars, hCutRatios;
        std::vector<std::string> cutLabels;
        loadOverallVariations(part.cutVariations, hCutVars, hCutRatios, cutLabels, false);

        std::vector<TH1*> hFeedVars, hFeedRatios;
        std::vector<std::string> feedLabels;
        loadOverallVariations(part.feeddownVariations, hFeedVars, hFeedRatios, feedLabels, false);

        dirPart->mkdir("IndividualSystematics");
        dirPart->cd("IndividualSystematics");
        TDirectory* dirIndiv = gDirectory;

        // A. Save each Fit Systematic on its own Canvas
        for (size_t i = 0; i < hFitRatios.size(); ++i) {
            TCanvas* cSingle = new TCanvas(Form("cFitSys_%s_%zu", part.name.c_str(), i), Form("Fit Systematic: %s", part.fitVariations[i].label.c_str()), 800, 600);
            cSingle->SetLeftMargin(0.12);
            cSingle->SetBottomMargin(0.12);

            hFitRatios[i]->SetTitle("");
            hFitRatios[i]->GetYaxis()->SetTitle("Variation / Reference");
            hFitRatios[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            hFitRatios[i]->GetYaxis()->SetRangeUser(0.85, 1.15);
            hFitRatios[i]->SetLineWidth(2);
            hFitRatios[i]->SetLineColor(kBlue+1);
            hFitRatios[i]->Draw("HIST");

            TLine* lineOne = new TLine(hFitRatios[i]->GetXaxis()->GetXmin(), 1.0, hFitRatios[i]->GetXaxis()->GetXmax(), 1.0);
            lineOne->SetLineStyle(2);
            lineOne->SetLineColor(kGray+2);
            lineOne->Draw();

            dirIndiv->cd();
            cSingle->Write(Form("cFitSys_%s", part.fitVariations[i].label.c_str()));
            delete cSingle;
        }

        // B. Save Material Budget on its own Canvas
        {
            TCanvas* cMat = new TCanvas(Form("cMatBudget_%s", part.name.c_str()), "Material Budget Systematic", 800, 600);
            cMat->SetLeftMargin(0.12);
            cMat->SetBottomMargin(0.12);

            TH1* hMat = (TH1*)hFinalOverall->Clone(Form("hMatBudgetPlot_%s", part.name.c_str()));
            hMat->Reset();
            for (int bin = 1; bin <= hMat->GetNbinsX(); ++bin) {
                hMat->SetBinContent(bin, materialBudgetFrac);
            }
            hMat->SetTitle("");
            hMat->GetYaxis()->SetTitle("Relative Systematic Uncertainty");
            hMat->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            hMat->GetYaxis()->SetRangeUser(0.0, materialBudgetFrac * 2.0);
            hMat->SetLineWidth(2);
            hMat->SetLineColor(kRed+1);
            hMat->Draw("HIST");

            dirIndiv->cd();
            cMat->Write("cMaterialBudget");
            delete cMat;
            delete hMat;
        }

        dirPart->cd();

        // Build Final Systematic Uncertainty Histogram for Overall Spectrum
        TH1* hSystematic = dynamic_cast<TH1*>(hFinalOverall->Clone(Form("hSystematic_%s", part.name.c_str())));
        hSystematic->Reset();

        for (int bin = 1; bin <= hSystematic->GetNbinsX(); ++bin) {
            double sumSquared = materialBudgetFrac * materialBudgetFrac;
            double valRef = hFinalOverall->GetBinContent(bin);
            double errRef = hFinalOverall->GetBinError(bin);

            auto accumulateOverall = [&](const std::vector<SystematicVariation>& varList, const std::vector<TH1*>& varHists, const std::vector<TH1*>& ratioHists) {
                for (size_t i = 0; i < ratioHists.size(); ++i) {
                    double valVar = varHists[i]->GetBinContent(bin);
                    double errVar = varHists[i]->GetBinError(bin);

                    bool passes = false;
                    CalculateBarlowCriterionValue(part.name, "Overall", varList[i].label, valRef, errRef, valVar, errVar, 2.0, &passes);

                    if (passes) {
                        double ratio = ratioHists[i]->GetBinContent(bin);
                        double systematic = std::abs(ratio - 1.0);
                        sumSquared += systematic * systematic;
                    }
                }
            };

            accumulateOverall(part.fitVariations, hFitVars, hFitRatios);
            accumulateOverall(part.cutVariations, hCutVars, hCutRatios);
            accumulateOverall(part.feeddownVariations, hFeedVars, hFeedRatios);

            hSystematic->SetBinContent(bin, std::sqrt(sumSquared) * 100.0);
        }

        // C. Combined Systematics Summary Canvas (Top + Ratios)
        if (!hFitRatios.empty() || true) {
            TCanvas* c = new TCanvas(Form("cSystematics_%s", part.name.c_str()), Form("Systematic variations - %s", part.name.c_str()), 1000, 900);
            TPad* padTop = new TPad("padTop", "Spectra", 0, 0.30, 1, 1);
            TPad* padBottom = new TPad("padBottom", "Ratios", 0, 0, 1, 0.30);

            padTop->SetBottomMargin(0.02);
            padBottom->SetTopMargin(0.02);
            padBottom->SetBottomMargin(0.30);

            padTop->Draw();
            padBottom->Draw();

            padTop->cd();
            hFinalOverall->SetTitle("");
            hFinalOverall->GetYaxis()->SetTitle(Form("%s p_{T} spectrum (Overall)", part.name.c_str()));
            hFinalOverall->GetXaxis()->SetLabelSize(0);
            hFinalOverall->SetLineWidth(3);
            hFinalOverall->SetLineColor(kBlack);
            hFinalOverall->Draw("HIST");

            TLegend* leg = new TLegend(0.52, 0.50, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);
            leg->AddEntry(hFinalOverall, Form("Final / %s (reference)", part.refFitLabel.c_str()), "l");

            int colorIdx = 0;
            for (size_t i = 0; i < hFitVars.size(); ++i) {
                Color_t col = colors[colorIdx++ % colors.size()];
                hFitVars[i]->SetLineWidth(2);
                hFitVars[i]->SetLineColor(col);
                hFitVars[i]->Draw("HIST SAME");

                hFitRatios[i]->SetLineWidth(2);
                hFitRatios[i]->SetLineColor(col);
                leg->AddEntry(hFitVars[i], fitLabels[i].c_str(), "l");
            }
            leg->Draw();

            padBottom->cd();
            if (!hFitRatios.empty()) {
                hFitRatios[0]->SetTitle("");
                hFitRatios[0]->GetYaxis()->SetTitle("Variation / Final");
                hFitRatios[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hFitRatios[0]->GetYaxis()->SetTitleSize(0.08);
                hFitRatios[0]->GetYaxis()->SetLabelSize(0.07);
                hFitRatios[0]->GetXaxis()->SetTitleSize(0.10);
                hFitRatios[0]->GetXaxis()->SetLabelSize(0.08);
                hFitRatios[0]->GetYaxis()->SetRangeUser(0.9, 1.1);
                hFitRatios[0]->Draw("HIST");

                for (size_t i = 1; i < hFitRatios.size(); ++i) {
                    hFitRatios[i]->Draw("HIST SAME");
                }
            }

            TH1* hRefRatioRef = (!hFitRatios.empty()) ? hFitRatios[0] : hFinalOverall;
            TLine* lineOne = new TLine(hRefRatioRef->GetXaxis()->GetXmin(), 1.0, hRefRatioRef->GetXaxis()->GetXmax(), 1.0);
            lineOne->SetLineStyle(2);
            lineOne->SetLineColor(kGray+2);
            lineOne->Draw();

            dirPart->cd();
            c->Write("cSystematics");
            delete c;
        }

        // D. SAVE THE FINAL SYSTEMATIC UNCERTAINTY HISTOGRAM AND CANVAS
        TCanvas* cSyst = new TCanvas(Form("cFinalSystematic_%s", part.name.c_str()), Form("Final systematic uncertainty - %s", part.name.c_str()), 1000, 700);
        hSystematic->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hSystematic->GetYaxis()->SetTitle("Systematic uncertainty (%)");
        hSystematic->SetLineWidth(3);
        hSystematic->Draw("HIST");

        // Write ALL outputs back to output ROOT file
        dirPart->cd();
        cSyst->Write("cFinalSystematic");
        hFinalOverall->Write("hFinalOverall");
        hSystematic->Write("hSystematicUncertainty");


    // ============================================================
    // E. ALL SYSTEMATICS COMPARISON CANVAS
    // ============================================================

    // Create individual subsystematic histograms
    TH1* hFitSystematic = dynamic_cast<TH1*>(
        hFinalOverall->Clone(Form("hFitSystematic_%s", part.name.c_str())));
    hFitSystematic->Reset();

    TH1* hMaterialSystematic = dynamic_cast<TH1*>(
        hFinalOverall->Clone(Form("hMaterialSystematic_%s", part.name.c_str())));
    hMaterialSystematic->Reset();

    // Calculate individual subsystematics
    for (int bin = 1; bin <= hSystematic->GetNbinsX(); ++bin) {

        double fitSumSquared = 0.0;

        double valRef = hFinalOverall->GetBinContent(bin);
        double errRef = hFinalOverall->GetBinError(bin);

        // Fit systematic
        for (size_t i = 0; i < hFitRatios.size(); ++i) {

            double valVar = hFitVars[i]->GetBinContent(bin);
            double errVar = hFitVars[i]->GetBinError(bin);

            bool passes = false;

            CalculateBarlowCriterionValue(
                part.name, "Overall", fitLabels[i],
                valRef, errRef, valVar, errVar,
                2.0, &passes);

            if (passes) {
                double ratio = hFitRatios[i]->GetBinContent(bin);
                double systematic = std::abs(ratio - 1.0);

                fitSumSquared += systematic * systematic;
            }
        }

        hFitSystematic->SetBinContent(bin, std::sqrt(fitSumSquared) * 100.0);

        // Material budget systematic
        hMaterialSystematic->SetBinContent(bin, materialBudgetFrac * 100.0);
        }

        // Create canvas
        TCanvas* cAll = new TCanvas(Form("cAllSystematics_%s", part.name.c_str()),Form("All systematic uncertainties - %s", part.name.c_str()),1000, 700);

        cAll->SetLeftMargin(0.12);
        cAll->SetBottomMargin(0.12);

        double maxY = std::max({hSystematic->GetMaximum(),hFitSystematic->GetMaximum(),hMaterialSystematic->GetMaximum()});
        hSystematic->GetYaxis()->SetRangeUser(0, maxY * 1.1);

        // Total systematic
        hSystematic->SetTitle("");
        hSystematic->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hSystematic->GetYaxis()->SetTitle("Systematic uncertainty (%)");
        hSystematic->SetLineWidth(3);
        hSystematic->SetLineColor(kBlack);
        hSystematic->Draw("HIST");

        // Fit systematic
        hFitSystematic->SetLineWidth(3);
        hFitSystematic->SetLineColor(kGreen+1);
        hFitSystematic->Draw("HIST SAME");

        // Material budget systematic
        hMaterialSystematic->SetLineWidth(3);
        hMaterialSystematic->SetLineColor(kOrange+1);
        hMaterialSystematic->Draw("HIST SAME");

        // Legend
        TLegend* legAll = new TLegend(0.55, 0.70, 0.88, 0.88);
        legAll->SetBorderSize(0);
        legAll->SetTextSize(0.035);

        legAll->AddEntry(hSystematic, "Total systematic uncertainty", "l");
        legAll->AddEntry(hFitSystematic, "Fit systematic", "l");
        legAll->AddEntry(hMaterialSystematic, "Material budget systematic", "l");

        legAll->Draw();

        // Save canvas and subsystematic histograms
        dirPart->cd();

        cAll->Write("cAllSystematics");

        hFitSystematic->Write("hFitSystematic");
        hMaterialSystematic->Write("hMaterialSystematic");

        // Cleanup
        delete cAll;
        delete legAll;

        delete hFitSystematic;
        delete hMaterialSystematic;

        delete cSyst;
        delete hSystematic;
        delete hFinalOverall;

        for (auto h : hFitVars) delete h;
        for (auto h : hFitRatios) delete h;
        for (auto h : hCutVars) delete h;
        for (auto h : hCutRatios) delete h;
        for (auto h : hFeedVars) delete h;
        for (auto h : hFeedRatios) delete h;

        // ============================================================
        // 3. PROCESS MCChi2 CANVASES & hSum HISTOGRAMS
        // ============================================================

        dirPart->cd();
        TCanvas* cRefMCChi2 = dynamic_cast<TCanvas*>(fFinal->Get((part.overallDir + "/" + part.canvasNameMCChi2).c_str()));
        TH1* hFinalMCChi2 = nullptr;
        if (cRefMCChi2) {
            TH1* hTemp = FindHistogramInCanvas(cRefMCChi2, part.histNameMCChi2);
            if (hTemp) {
                hFinalMCChi2 = dynamic_cast<TH1*>(hTemp->Clone(Form("hFinalMCChi2_%s", part.name.c_str())));
                hFinalMCChi2->SetDirectory(nullptr);
            }
        }

        if (hFinalMCChi2) {
            std::vector<TH1*> hChi2Variations;
            std::vector<std::string> chi2Labels;

            for (size_t i = 0; i < part.fitVariations.size(); ++i) {
                TFile* fVar = TFile::Open(part.fitVariations[i].file.c_str());
                if (!fVar || fVar->IsZombie()) continue;

                TCanvas* cVarMCChi2 = dynamic_cast<TCanvas*>(fVar->Get((part.overallDir + "/" + part.canvasNameMCChi2).c_str()));
                TH1* hVar = nullptr;
                if (cVarMCChi2) {
                    TH1* hTempVar = FindHistogramInCanvas(cVarMCChi2, part.histNameMCChi2);
                    if (hTempVar) {
                        hVar = dynamic_cast<TH1*>(hTempVar->Clone(Form("hChi2Variation_%zu_%s", i, part.name.c_str())));
                        hVar->SetDirectory(nullptr);
                    }
                }

                if (!hVar) {
                    fVar->Close();
                    continue;
                }

                hChi2Variations.push_back(hVar);
                std::string currentLabel = part.fitVariations[i].label;
                if (currentLabel == part.refFitLabel) currentLabel += " (reference)";
                chi2Labels.push_back(currentLabel);
                fVar->Close();
            }

            double yMax = hFinalMCChi2->GetMaximum();
            for (const auto* hVar : hChi2Variations) {
                if (hVar->GetMaximum() > yMax) yMax = hVar->GetMaximum();
            }

            TCanvas* cChi2 = new TCanvas(Form("cMCChi2_Variations_%s", part.name.c_str()), Form("MCChi2 hSum variations - %s", part.name.c_str()), 1000, 700);
            cChi2->SetLeftMargin(0.12);
            cChi2->SetBottomMargin(0.12);

            hFinalMCChi2->SetTitle("");
            hFinalMCChi2->GetYaxis()->SetTitle(Form("%s hSum (MCChi2)", part.name.c_str()));
            hFinalMCChi2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            hFinalMCChi2->GetYaxis()->SetRangeUser(0.0, yMax * 1.15);
            hFinalMCChi2->SetLineWidth(2);
            hFinalMCChi2->SetLineColor(kBlack);
            hFinalMCChi2->Draw("HIST");

            TLegend* legChi2 = new TLegend(0.52, 0.50, 0.88, 0.88);
            legChi2->SetBorderSize(0);
            legChi2->SetTextSize(0.03);
            legChi2->AddEntry(hFinalMCChi2, Form("Final / %s (reference)", part.refFitLabel.c_str()), "l");

            for (size_t i = 0; i < hChi2Variations.size(); ++i) {
                Color_t col = colors[i % colors.size()];
                hChi2Variations[i]->SetLineWidth(2);
                hChi2Variations[i]->SetLineColor(col);
                hChi2Variations[i]->Draw("HIST SAME");
                legChi2->AddEntry(hChi2Variations[i], chi2Labels[i].c_str(), "l");
            }
            legChi2->Draw();

            dirPart->cd();
            cChi2->Write("cMCChi2_Comparison");
            hFinalMCChi2->Write("hFinal_MCChi2_hSum");

            delete cChi2;
            delete legChi2;
            for (auto* h : hChi2Variations) delete h;
            delete hFinalMCChi2;
        }
    }

    fOut->Close();
    fFinal->Close();
    std::cout << "\nSuccess!!\n";
}