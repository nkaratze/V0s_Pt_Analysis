#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TText.h>
#include <TMultiGraph.h>
#include <vector>
#include <iostream>


void PlotTestFit(TGraph* g, TF1* fit, const char* canvasName)
{
    if (!g || !fit) {
        std::cout << "Cannot plot test: graph or fit is null." << std::endl;
        return;
    }

    TCanvas* c = new TCanvas(canvasName,Form("Fit test: %s", g->GetTitle()),900, 700);

    c->cd();

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.0);
    g->Draw("AP");

    // Make the fit visible on top of the data.
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->Draw("SAME");

    c->Modified();
    c->Update();
}

TGraphErrors* MakeGraphErrors(TGraph* original, const char* name)
{
    if (!original) return nullptr;

    auto* sourceErrors = dynamic_cast<TGraphErrors*>(original);
    auto* result = new TGraphErrors(original->GetN());

    result->SetName(name);
    result->SetTitle(original->GetTitle());

    for (int i = 0; i < original->GetN(); ++i) {
        double x, y;
        original->GetPoint(i, x, y);

        double ex = 0.0;
        double ey = 0.0;

        if (sourceErrors) {
            ex = sourceErrors->GetErrorX(i);
            ey = sourceErrors->GetErrorY(i);
        }

        result->SetPoint(i, x, y);
        result->SetPointError(i, ex, ey);
    }

    return result;
} 

void OverlayFinalCanvases()
{
    // Open input files
    TFile *file136 = TFile::Open("V0s_Nch_EnergySpectrumLHC24pp.root", "READ");
    TFile *file536 = TFile::Open("V0s_Nch_EnergySpectrumLHC24ppref.root", "READ");
    TFile *file502 = TFile::Open("HEPData502Paper.root","READ");
    TFile *file130 = TFile::Open("HEPData130Paper.root","READ");

    // Official Nch at 5.36 TeV
    double nchOfficial[10] = {15.4103, 13.2991, 11.6304, 10.4861, 9.62378, 8.57903, 7.44574, 6.48026, 5.30030, 3.22057};
    double nchOfficialErr[10] = {0.0403394, 0.0127358, 0.00892575, 0.00787672, 0.0103695, 0.0042811, 0.00662173, 0.00353035, 0.00440069, 0.00124722};

    double minNch = 0.0;
    double maxNch = 30.0;

    if (!file136 || file136->IsZombie()) {
        std::cout << "Error opening 13.6 TeV file." << std::endl;
        return;
    }

    if (!file536 || file536->IsZombie()) {
        std::cout << "Error opening 5.36 TeV file." << std::endl;
        return;
    }

    if (!file502 || file502->IsZombie()) {
        std::cout << "Error opening 5.02 TeV file." << std::endl;
        return;
    }

    if (!file130 || file130->IsZombie()) {
        std::cout << "Error opening 13.0 TeV file." << std::endl;
        return;
    }

    // Retrieve canvases
    TCanvas *c136Yield = (TCanvas*)file136->Get("cYieldVsNchAllParticles");
    TCanvas *c536Yield = (TCanvas*)file536->Get("cYieldVsNchAllParticles");
    TCanvas *c136MeanPt = (TCanvas*)file136->Get("cMeanPtVsNchAllParticles");
    TCanvas *c536MeanPt = (TCanvas*)file536->Get("cMeanPtVsNchAllParticles");
    TCanvas *c536YieldOfficial = (TCanvas*)file536->Get("cYieldVsNchAllParticles");
    TCanvas *c536MeanPtOfficial = (TCanvas*)file536->Get("cMeanPtVsNchAllParticles");

    if (!c136Yield || !c536Yield || !c136MeanPt || !c536MeanPt) {
        std::cout << "Could not retrieve one or more input canvases." << std::endl;
        return;
    }

    std::vector<TGraphErrors*> graphs136Yield;
    std::vector<TGraphErrors*> graphs536Yield;
    std::vector<TGraphErrors*> graphs536YieldOfficial;
    std::vector<TGraphErrors*> graphs136MeanPt;
    std::vector<TGraphErrors*> graphs536MeanPt;
    std::vector<TGraphErrors*> graphs536MeanPtOfficial;

    std::vector<TF1*> fits136Yield;
    std::vector<TF1*> fits536Yield;
    std::vector<TF1*> fits536YieldOfficial;
    std::vector<TF1*> fits136MeanPt;
    std::vector<TF1*> fits536MeanPt;
    std::vector<TF1*> fits536MeanPtOfficial;

    // ---------------------------------------
    // Remove original legends and text labels
    // ---------------------------------------
    TList *list136Yield = c136Yield->GetListOfPrimitives();
    TList *list536Yield = c536Yield->GetListOfPrimitives();
    TList *list536YieldOfficial = c536Yield->GetListOfPrimitives();
    TList *list136MeanPt = c136MeanPt->GetListOfPrimitives();
    TList *list536MeanPt = c536MeanPt->GetListOfPrimitives();
    TList *list536MeanPtOfficial = c536MeanPt->GetListOfPrimitives();

    auto cleanCanvasList = [](TList *list) {
        if (!list) return;
        for (int i = list->GetSize() - 1; i >= 0; i--) {
            TObject *obj = list->At(i);
            if (obj && (obj->InheritsFrom(TLegend::Class()) ||
                        obj->InheritsFrom(TLatex::Class()) ||
                        obj->InheritsFrom(TPaveText::Class()) ||
                        obj->InheritsFrom(TText::Class()))) {
                list->Remove(obj);
            }
        }
    };

    cleanCanvasList(list136Yield);
    cleanCanvasList(list136MeanPt);
    cleanCanvasList(list536Yield);
    cleanCanvasList(list536MeanPt);
    cleanCanvasList(list536YieldOfficial);
    cleanCanvasList(list536MeanPtOfficial);

    // Get 5.02 TeV Results
    TDirectory *dirK0sYield = (TDirectory*)file502->Get("Figure 4a");
    TDirectory *dirLambdaYield = (TDirectory*)file502->Get("Figure 4b");
    TDirectory *dirK0sMeanPt = (TDirectory*)file502->Get("Figure 2a");
    TDirectory *dirLambdaMeanPt = (TDirectory*)file502->Get("Figure 2b");

    if (!dirK0sYield || !dirLambdaYield || !dirK0sMeanPt || !dirLambdaMeanPt) {
        std::cout << "Could not find 5.02 TeV directories." << std::endl;
        return;
    }

    TGraph *g502K0sYield = (TGraph*)dirK0sYield->Get("Graph1D_y1");
    TGraph *g502LambdaYield = (TGraph*)dirLambdaYield->Get("Graph1D_y1");
    TGraph *g502K0sMeanPt = (TGraph*)dirK0sMeanPt->Get("Graph1D_y1");
    TGraph *g502LambdaMeanPt = (TGraph*)dirLambdaMeanPt->Get("Graph1D_y1");

    // Get 13.0 TeV Results
    TDirectory *dirK0sYield130 = (TDirectory*)file130->Get("Table 8c");
    TDirectory *dirLambdaYield130 = (TDirectory*)file130->Get("Table 8c");
    TDirectory *dirK0sMeanPt130 = (TDirectory*)file130->Get("Table 5c");
    TDirectory *dirLambdaMeanPt130 = (TDirectory*)file130->Get("Table 5c");

    if (!dirK0sYield130 || !dirLambdaYield130 || !dirK0sMeanPt130 || !dirLambdaMeanPt130) {
        std::cout << "Could not find 13.0 TeV directories." << std::endl;
        return;
    }

    TGraph *g130K0sYield = (TGraph*)dirK0sYield130->Get("Graph1D_y1");
    TGraph *g130LambdaYield = (TGraph*)dirLambdaYield130->Get("Graph1D_y2");
    TGraph *g130K0sMeanPt = (TGraph*)dirK0sMeanPt130->Get("Graph1D_y1");
    TGraph *g130LambdaMeanPt = (TGraph*)dirLambdaMeanPt130->Get("Graph1D_y2");

    // ---------------------------------------
    // Style 13.6 TeV (GREEN)
    // ---------------------------------------
    TIter next136Yield(c136Yield->GetListOfPrimitives());
    TObject *obj136Yield = nullptr;
    int counter136Yield = 0;

    while ((obj136Yield = next136Yield())) {
        if (obj136Yield->InheritsFrom(TGraph::Class())) {
            TGraph *gOriginal = (TGraph*)obj136Yield;
            TGraphErrors *g = MakeGraphErrors(gOriginal, Form("g136Yield_%d", counter136Yield));
            if (!g) continue;

            graphs136Yield.push_back(g);
            g->SetLineColor(kGreen+2);
            g->SetMarkerColor(kGreen+2);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);
            g->SetMarkerStyle(counter136Yield == 0 ? 21 : 20);

            std::cout << "Fitting graph 13.6 TeV Yield: " << g->GetName()<< " | " << g->GetTitle() << std::endl;

            TF1* fit136Yield = new TF1(Form("fit136_%d", counter136Yield),"pol1",minNch,maxNch);
            g->Fit(fit136Yield, "L");

            fit136Yield->SetLineColor(kGreen+2);
            fit136Yield->SetLineWidth(1);
            fits136Yield.push_back(fit136Yield);

            counter136Yield++;
        }
    }

    TIter next136MeanPt(c136MeanPt->GetListOfPrimitives());
    TObject *obj136MeanPt = nullptr;
    int counter136MeanPt = 0;

    while ((obj136MeanPt = next136MeanPt())) {
        if (obj136MeanPt->InheritsFrom(TGraph::Class())) {
            TGraph *gOriginal = (TGraph*)obj136MeanPt;
            TGraphErrors *g = MakeGraphErrors(gOriginal, Form("g136MeanPt_%d", counter136MeanPt));
            if (!g) continue;

            graphs136MeanPt.push_back(g);
            g->SetLineColor(kGreen+2);
            g->SetMarkerColor(kGreen+2);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);
            g->SetMarkerStyle(counter136MeanPt == 0 ? 21 : 20);

            //TF1 *fit136MeanPt = new TF1(Form("fit136_meanpt_%d", counter136MeanPt), "pol1");
            //g->Fit(fit136MeanPt, "L");
            //fit136MeanPt->SetLineColor(kGreen+2);
            //fit136MeanPt->SetLineWidth(1);
            //fits136MeanPt.push_back(fit136MeanPt);

            counter136MeanPt++;
        }
    }

    // ---------------------------------------
    // Style 5.36 TeV (ORANGE)
    // ---------------------------------------
    TIter next536Yield(c536Yield->GetListOfPrimitives());
    TObject *obj536Yield = nullptr;
    int counter536Yield = 0;

    while ((obj536Yield = next536Yield())) {
        if (obj536Yield->InheritsFrom(TGraph::Class())) {
            TGraph *gOriginal = (TGraph*)obj536Yield;
            TGraphErrors *g = MakeGraphErrors(gOriginal, Form("g536Yield_%d", counter536Yield));
            if (!g) continue;

            graphs536Yield.push_back(g);
            g->SetLineColor(kOrange+1);
            g->SetMarkerColor(kOrange+1);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);
            g->SetMarkerStyle(counter536Yield == 0 ? 29 : 33);

            std::cout << "Fitting graph 5.36 TeV Yield: " << g->GetName()<< " | " << g->GetTitle() << std::endl;
            TF1* fit536Yield = new TF1(Form("fit536_%d", counter536Yield),"pol1",minNch,maxNch);
            g->Fit(fit536Yield, "L");
            fit536Yield->SetLineColor(kOrange+1);
            fit536Yield->SetLineWidth(1);
            fits536Yield.push_back(fit536Yield);

            counter536Yield++;
        }
    }

    TIter next536MeanPt(c536MeanPt->GetListOfPrimitives());
    TObject *obj536MeanPt = nullptr;
    int counter536MeanPt = 0;

    while ((obj536MeanPt = next536MeanPt())) {
        if (obj536MeanPt->InheritsFrom(TGraph::Class())) {
            TGraph *gOriginal = (TGraph*)obj536MeanPt;
            TGraphErrors *g = MakeGraphErrors(gOriginal, Form("g536MeanPt_%d", counter536MeanPt));
            if (!g) continue;

            graphs536MeanPt.push_back(g);
            g->SetLineColor(kOrange+1);
            g->SetMarkerColor(kOrange+1);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);
            g->SetMarkerStyle(counter536MeanPt == 0 ? 29 : 33);

            //TF1 *fit536MeanPt = new TF1(Form("fit536_meanpt_%d", counter536MeanPt), "pol1");
            //g->Fit(fit536MeanPt, "L");
            //fit536MeanPt->SetLineColor(kOrange+1);
            //fit536MeanPt->SetLineWidth(1);
            //fits536MeanPt.push_back(fit536MeanPt);

            counter536MeanPt++;
        }
    }

    // ------------------------------------------------
    // Style 5.36 TeV (BLACK) with official Nch values
    // ------------------------------------------------
    TIter next536YieldOfficial(c536YieldOfficial->GetListOfPrimitives());
    TObject *obj536YieldOfficial = nullptr;
    int counter536YieldOfficial = 0;

    while ((obj536YieldOfficial = next536YieldOfficial())) {
        if (obj536YieldOfficial->InheritsFrom(TGraph::Class())) {

            TGraph *gOriginal = (TGraph*)obj536YieldOfficial;

            // Use the helper to copy the original graph and its errors
            TGraphErrors *g = MakeGraphErrors(gOriginal,Form("g536YieldOfficial_%d", counter536YieldOfficial));

            if (!g) continue;

            // Style the official-Nch graph
            g->SetLineColor(kBlack);
            g->SetMarkerColor(kBlack);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);
            g->SetMarkerStyle(counter536YieldOfficial == 0 ? 29 : 33);

            // Replace x-values and horizontal errors with official Nch values.
            // Preserve the y-values and vertical errors from MakeGraphErrors().
            for (int i = 0; i < g->GetN() && i < 10; i++) {

                double x, y;
                g->GetPoint(i, x, y);

                double ey = g->GetErrorY(i);

                g->SetPoint(i, nchOfficial[i], y);
                //g->SetPointError(i, nchOfficialErr[i], ey);
            }

            graphs536YieldOfficial.push_back(g);

            std::cout << "Fitting graph 5.36 TeV Yield Official: "<< g->GetName() << " | " << g->GetTitle() << std::endl;

            // Print the points and uncertainties for verification
            std::cout << "\n--- " << g->GetName() << " ---\n";

            for (int i = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                std::cout << "Point " << i<< ": x = " << x<< ", y = " << y << ", ex = " << g->GetErrorX(i) << ", ey = " << g->GetErrorY(i) << '\n';
            }

            // Create and fit the linear function
            TF1* fit536YieldOfficial = new TF1(Form("fit536Official_%d", counter536YieldOfficial),"pol1",minNch,maxNch);
            g->Fit(fit536YieldOfficial, "L");
            fit536YieldOfficial->SetLineColor(kBlack);
            fit536YieldOfficial->SetLineWidth(1);
            fits536YieldOfficial.push_back(fit536YieldOfficial);

            counter536YieldOfficial++;
        }
    }

    TIter next536MeanPtOfficial(c536MeanPtOfficial->GetListOfPrimitives());
    TObject *obj536MeanPtOfficial = nullptr;
    int counter536MeanPtOfficial = 0;

    while ((obj536MeanPtOfficial = next536MeanPtOfficial())) {
        if (obj536MeanPtOfficial->InheritsFrom(TGraph::Class())) {
            TGraph *gOriginal = (TGraph*)obj536MeanPtOfficial;
            TGraphErrors *g = new TGraphErrors(gOriginal->GetN());

            g->SetName(Form("g536MeanPtOfficial_%d", counter536MeanPtOfficial));
            g->SetTitle(gOriginal->GetTitle());
            g->SetLineColor(kBlack);
            g->SetMarkerColor(kBlack);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);
            g->SetMarkerStyle(counter536MeanPtOfficial == 0 ? 29 : 33);

            for (int i = 0; i < gOriginal->GetN() && i < 10; i++) {
                double x, y;
                gOriginal->GetPoint(i, x, y);
                double ey = 0.0;
                if (auto *ge = dynamic_cast<TGraphErrors*>(gOriginal)) {
                    ey = ge->GetErrorY(i);
                }
                g->SetPoint(i, nchOfficial[i], y);
                g->SetPointError(i, nchOfficialErr[i], ey);
            }

            graphs536MeanPtOfficial.push_back(g);

            //TF1 *fit = new TF1(Form("fit536MeanPtOfficial_%d", counter536MeanPtOfficial), "pol1");
            //g->Fit(fit, "L");
            //fit->SetLineColor(kBlack);
            //fit->SetLineWidth(1);
            //fits536MeanPtOfficial.push_back(fit);

            counter536MeanPtOfficial++;
        }
    }

    // ---------------------------------------
    // Style 5.02 TeV & 13.0 TeV
    // ---------------------------------------
    g502K0sYield->SetMarkerStyle(22); g502K0sYield->SetMarkerColor(kBlue+1); g502K0sYield->SetLineColor(kBlue+1); g502K0sYield->SetMarkerSize(1.2);
    g502LambdaYield->SetMarkerStyle(23); g502LambdaYield->SetMarkerColor(kBlue+1); g502LambdaYield->SetLineColor(kBlue+1); g502LambdaYield->SetMarkerSize(1.2);
    g502K0sMeanPt->SetMarkerStyle(22); g502K0sMeanPt->SetMarkerColor(kBlue+1); g502K0sMeanPt->SetLineColor(kBlue+1); g502K0sMeanPt->SetMarkerSize(1.2);
    g502LambdaMeanPt->SetMarkerStyle(23); g502LambdaMeanPt->SetMarkerColor(kBlue+1); g502LambdaMeanPt->SetLineColor(kBlue+1); g502LambdaMeanPt->SetMarkerSize(1.2);

    g130K0sYield->SetMarkerStyle(34); g130K0sYield->SetMarkerColor(kRed+1); g130K0sYield->SetLineColor(kRed+1); g130K0sYield->SetMarkerSize(1.2);
    g130LambdaYield->SetMarkerStyle(47); g130LambdaYield->SetMarkerColor(kRed+1); g130LambdaYield->SetLineColor(kRed+1); g130LambdaYield->SetMarkerSize(1.2);
    g130K0sMeanPt->SetMarkerStyle(34); g130K0sMeanPt->SetMarkerColor(kRed+1); g130K0sMeanPt->SetLineColor(kRed+1); g130K0sMeanPt->SetMarkerSize(1.2);
    g130LambdaMeanPt->SetMarkerStyle(47); g130LambdaMeanPt->SetMarkerColor(kRed+1); g130LambdaMeanPt->SetLineColor(kRed+1); g130LambdaMeanPt->SetMarkerSize(1.2);

    // ---------------------------------------
    // Setup Legend
    // ---------------------------------------
    TLegend *leg = new TLegend(0.1, 0.65, 0.35, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);

    TGraph *gRedCross = new TGraph(); gRedCross->SetMarkerStyle(34); gRedCross->SetMarkerColor(kRed+1);
    TGraph *gRedPlus = new TGraph(); gRedPlus->SetMarkerStyle(47); gRedPlus->SetMarkerColor(kRed+1);
    TGraph *gGreenSquare = new TGraph(); gGreenSquare->SetMarkerStyle(21); gGreenSquare->SetMarkerColor(kGreen+2);
    TGraph *gGreenCircle = new TGraph(); gGreenCircle->SetMarkerStyle(20); gGreenCircle->SetMarkerColor(kGreen+2);
    TGraph *gOrangeStar = new TGraph(); gOrangeStar->SetMarkerStyle(29); gOrangeStar->SetMarkerColor(kOrange+1);
    TGraph *gBlackStar = new TGraph(); gBlackStar->SetMarkerStyle(29); gBlackStar->SetMarkerColor(kBlack);
    TGraph *gOrangeDiamond = new TGraph(); gOrangeDiamond->SetMarkerStyle(33); gOrangeDiamond->SetMarkerColor(kOrange+1);
    TGraph *gBlackDiamond = new TGraph(); gBlackDiamond->SetMarkerStyle(33); gBlackDiamond->SetMarkerColor(kBlack);
    TGraph *gBlueTriangle = new TGraph(); gBlueTriangle->SetMarkerStyle(22); gBlueTriangle->SetMarkerColor(kBlue+1);
    TGraph *gBlueTriangleDown = new TGraph(); gBlueTriangleDown->SetMarkerStyle(23); gBlueTriangleDown->SetMarkerColor(kBlue+1);

    leg->AddEntry(gOrangeStar,"5.36 TeV: K^{0}_{S} This Analysis","p");
    leg->AddEntry(gOrangeDiamond,"5.36 TeV: #Lambda+#bar{#Lambda} This Analysis","p");
    leg->AddEntry(gBlackStar,"5.36 TeV: K^{0}_{S} Official Nch","p");
    leg->AddEntry(gBlackDiamond,"5.36 TeV: #Lambda+#bar{#Lambda} Official Nch","p");
    leg->AddEntry(gGreenSquare,"13.6 TeV: K^{0}_{S} Test","p");
    leg->AddEntry(gGreenCircle,"13.6 TeV: #Lambda+#bar{#Lambda} Test","p");
    leg->AddEntry(gRedCross,"13.0 TeV: K^{0}_{S} Run 2 ALICE","p");
    leg->AddEntry(gRedPlus,"13.0 TeV: #Lambda+#bar{#Lambda} Run 2 ALICE","p");
    leg->AddEntry(gBlueTriangle,"5.02 TeV: K^{0}_{S} Run 2 ALICE","p");
    leg->AddEntry(gBlueTriangleDown,"5.02 TeV: #Lambda+#bar{#Lambda} Run 2 ALICE","p");

    // ---------------------------------------
    // Create output canvases
    // ---------------------------------------
    TCanvas *cYieldAllResults = new TCanvas("cYield_FinalAllResults_Overlay","cYield_FinalAllResults_Overlay",800,600);
    TCanvas *cYieldAllResultsFits = new TCanvas("cYield_FinalAllResults_OverlayWithFits","cYield_FinalAllResults_OverlayWithFits",800,600);
    TCanvas *cMeanPtAllResults = new TCanvas("cMeanPt_FinalAllResults_Overlay","cMeanPt_FinalAllResults_Overlay",800,600);
    //TCanvas *cMeanPtAllResultsFits = new TCanvas("cMeanPt_FinalAllResults_OverlayWithFits","cMeanPt_FinalAllResults_OverlayWithFits",800,600);

    // Right before drawing your first plain yield/meanpt multigraphs, clear their internal functions:
    for (auto g : graphs136Yield) g->GetListOfFunctions()->Clear();
    for (auto g : graphs536Yield) g->GetListOfFunctions()->Clear();
    for (auto g : graphs536YieldOfficial) g->GetListOfFunctions()->Clear();
    for (auto g : graphs136MeanPt) g->GetListOfFunctions()->Clear();
    for (auto g : graphs536MeanPt) g->GetListOfFunctions()->Clear();
    //for (auto g : graphs536MeanPtOfficial) g->GetListOfFunctions()->Clear();


    // ======================================================
    // 1. Yields Multigraph Canvas (Points Only)
    // ======================================================
    cYieldAllResults->cd();
    TMultiGraph *mgYield = new TMultiGraph();
    mgYield->SetTitle("; <dN_{ch}/dy>; <dN/dy>");

    // Use "p" here so fits are NEVER drawn on this canvas
    for (auto g : graphs136Yield) mgYield->Add(g, "p");
    for (auto g : graphs536Yield) mgYield->Add(g, "p");
    for (auto g : graphs536YieldOfficial) mgYield->Add(g, "p");
    mgYield->Add(g502K0sYield, "p");
    mgYield->Add(g502LambdaYield, "p");
    mgYield->Add(g130K0sYield, "p");
    mgYield->Add(g130LambdaYield, "p");

    mgYield->Draw("AP");
    mgYield->GetXaxis()->SetTitleSize(0.05);
    mgYield->GetYaxis()->SetTitleSize(0.05);
    mgYield->GetXaxis()->SetLabelSize(0.045);
    mgYield->GetYaxis()->SetLabelSize(0.045);
    mgYield->GetYaxis()->SetTitleOffset(1.1);
    cYieldAllResults->Modified();

    leg->Draw();

    // ======================================================
    // 2. Mean Pt Multigraph Canvas (Points Only)
    // ======================================================
    cMeanPtAllResults->cd();
    TMultiGraph *mgMeanPt = new TMultiGraph();
    mgMeanPt->SetTitle("; <dN_{ch}/dy>; <p_{T}> (GeV/c)");

    for (auto g : graphs136MeanPt) mgMeanPt->Add(g, "p");
    for (auto g : graphs536MeanPt) mgMeanPt->Add(g, "p");
    for (auto g : graphs536MeanPtOfficial) mgMeanPt->Add(g, "p");
    mgMeanPt->Add(g502K0sMeanPt, "p");
    mgMeanPt->Add(g502LambdaMeanPt, "p");
    mgMeanPt->Add(g130K0sMeanPt, "p");
    mgMeanPt->Add(g130LambdaMeanPt, "p");

    mgMeanPt->Draw("AP");
    mgMeanPt->GetXaxis()->SetTitleSize(0.05);
    mgMeanPt->GetYaxis()->SetTitleSize(0.05);
    mgMeanPt->GetXaxis()->SetLabelSize(0.045);
    mgMeanPt->GetYaxis()->SetLabelSize(0.045);
    mgMeanPt->GetYaxis()->SetTitleOffset(1.1);
    cMeanPtAllResults->Modified();

    leg->Draw();

    // ======================================================
    // 3. Yields Plus Fits Canvas
    // ======================================================
    cYieldAllResultsFits->cd();

    TMultiGraph *mgYieldFits = new TMultiGraph();
    mgYieldFits->SetTitle("; <dN_{ch}/dy>; <dN/dy>");

    for (auto g : graphs136Yield) mgYieldFits->Add(g, "p");
    for (auto g : graphs536Yield) mgYieldFits->Add(g, "p");
    for (auto g : graphs536YieldOfficial) mgYieldFits->Add(g, "p");
    mgYieldFits->Add(g502K0sYield, "p");
    mgYieldFits->Add(g502LambdaYield, "p");
    mgYieldFits->Add(g130K0sYield, "p");
    mgYieldFits->Add(g130LambdaYield, "p");

    mgYieldFits->Draw("AP");

    mgYieldFits->GetXaxis()->SetTitleSize(0.05);
    mgYieldFits->GetYaxis()->SetTitleSize(0.05);
    mgYieldFits->GetXaxis()->SetLabelSize(0.045);
    mgYieldFits->GetYaxis()->SetLabelSize(0.045);
    mgYieldFits->GetYaxis()->SetTitleOffset(1.1);

    cYieldAllResultsFits->Modified();
    cYieldAllResultsFits->Update();

    // Draw the stored 13.6 TeV fit functions
    for (auto f : fits136Yield) {
        if (!f) continue;
        f->SetLineColor(kGreen+2);
        f->SetLineWidth(2);
        f->Draw("SAME");
    }

    // Draw the stored 5.36 TeV fit functions
    for (auto f : fits536Yield) {
        if (!f) continue;
        f->SetLineColor(kOrange+1);
        f->SetLineWidth(2);
        f->Draw("SAME");
    }

    for (auto f : fits136Yield) {
        if (f) f->Draw("SAME");
    }

    for (auto f : fits536Yield) {
        if (f) f->Draw("SAME");
    }

    for (auto f : fits536YieldOfficial) {
        if (f) f->Draw("SAME");
    }

    TLegend *leg2 = (TLegend*)leg->Clone();
    leg2->Draw();

    cYieldAllResultsFits->Modified();
    cYieldAllResultsFits->Update();

    // ---------------------------------------
    // Save output
    // ---------------------------------------
    TFile *outFile = new TFile("FinalComparisonOverlay.root","RECREATE");

    cYieldAllResults->Write();
    cYieldAllResultsFits->Write();
    cMeanPtAllResults->Write();

    outFile->Close();

    file136->Close();
    file536->Close();
    file502->Close();
    file130->Close();
}