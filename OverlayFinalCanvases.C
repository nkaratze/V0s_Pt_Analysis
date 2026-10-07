#include <vector>

void OverlayFinalCanvases()
{
    // Open input files
    TFile *file136 = TFile::Open("V0s_Nch_EnergySpectrumLHC24pp.root", "READ");
    TFile *file536 = TFile::Open("V0s_Nch_EnergySpectrumLHC24ppref.root", "READ");
    TFile *file502 = TFile::Open("HEPData502Paper.root","READ");
    TFile *file130 = TFile::Open("HEPData130Paper.root","READ");
    //7 TeV

    //Official Nch at 5.36 TeV
    double nchOfficial[10] = {15.4103, 13.2991, 11.6304, 10.4861, 9.62378, 8.57903, 7.44574, 6.48026,5.30030,3.22057};
    double nchOfficialErr[10] = {0.0403394,0.0127358,0.00892575,0.00787672,0.0103695,0.0042811,0.00662173,0.00353035,0.00440069,0.00124722};

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

    if (!c136Yield || !c536Yield) {
        std::cout << "Could not retrieve one or both Yield canvases." << std::endl;
        return;
    }

    if (!c136MeanPt || !c536MeanPt) {
        std::cout << "Could not retrieve one or both Mean Pt canvases." << std::endl;
        return;
    }

    std::vector<TGraph*> graphs136Yield;
    std::vector<TGraph*> graphs536Yield;
    std::vector<TGraph*> graphs536YieldOfficial;
    std::vector<TGraph*> graphs136MeanPt;
    std::vector<TGraph*> graphs536MeanPt;
    std::vector<TGraph*> graphs536MeanPtOfficial;

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

    //Get 136 Results
    for (int i = list136Yield->GetSize()-1; i >= 0; i--) {

        TObject *obj = list136Yield->At(i);

        if (obj->InheritsFrom(TLegend::Class()) ||
            obj->InheritsFrom(TLatex::Class()) ||
            obj->InheritsFrom(TPaveText::Class()) ||
            obj->InheritsFrom(TText::Class())) {

            list136Yield->Remove(obj);
        }
    }

    for (int i = list136MeanPt->GetSize()-1; i >= 0; i--) {

        TObject *obj = list136MeanPt->At(i);

        if (obj->InheritsFrom(TLegend::Class()) ||
            obj->InheritsFrom(TLatex::Class()) ||
            obj->InheritsFrom(TPaveText::Class()) ||
            obj->InheritsFrom(TText::Class())) {

            list136MeanPt->Remove(obj);
        }
    }

    //Get 536 Results
    for (int i = list536Yield->GetSize()-1; i >= 0; i--) {

        TObject *obj = list536Yield->At(i);

        if (obj->InheritsFrom(TLegend::Class()) ||
            obj->InheritsFrom(TLatex::Class()) ||
            obj->InheritsFrom(TPaveText::Class()) ||
            obj->InheritsFrom(TText::Class())) {

            list536Yield->Remove(obj);
        }
    }

    for (int i = list536MeanPt->GetSize()-1; i >= 0; i--) {

        TObject *obj = list536MeanPt->At(i);

        if (obj->InheritsFrom(TLegend::Class()) ||
            obj->InheritsFrom(TLatex::Class()) ||
            obj->InheritsFrom(TPaveText::Class()) ||
            obj->InheritsFrom(TText::Class())) {

            list536MeanPt->Remove(obj);
        }
    }

    //Get 536 Results Official
    for (int i = list536YieldOfficial->GetSize()-1; i >= 0; i--) {

        TObject *obj = list536YieldOfficial->At(i);

        if (obj->InheritsFrom(TLegend::Class()) ||
            obj->InheritsFrom(TLatex::Class()) ||
            obj->InheritsFrom(TPaveText::Class()) ||
            obj->InheritsFrom(TText::Class())) {

            list536YieldOfficial->Remove(obj);
        }
    }

    for (int i = list536MeanPtOfficial->GetSize()-1; i >= 0; i--) {

        TObject *obj = list536MeanPtOfficial->At(i);

        if (obj->InheritsFrom(TLegend::Class()) ||
            obj->InheritsFrom(TLatex::Class()) ||
            obj->InheritsFrom(TPaveText::Class()) ||
            obj->InheritsFrom(TText::Class())) {

            list536MeanPtOfficial->Remove(obj);
        }
    }

    //Get 502 Results
    TDirectory *dirK0sYield = (TDirectory*)file502->Get("Figure 4a");
    TDirectory *dirLambdaYield = (TDirectory*)file502->Get("Figure 4b");
    TDirectory *dirK0sMeanPt = (TDirectory*)file502->Get("Figure 2a");
    TDirectory *dirLambdaMeanPt = (TDirectory*)file502->Get("Figure 2b");

    if (!dirK0sYield || !dirLambdaYield) {
        std::cout << "Could not find Yield directories." << std::endl;
        return;
    }
    if (!dirK0sMeanPt || !dirLambdaMeanPt) {
        std::cout << "Could not find Mean Pt directories." << std::endl;
        return;
    }

    TGraph *g502K0sYield = (TGraph*)dirK0sYield->Get("Graph1D_y1");
    TGraph *g502LambdaYield = (TGraph*)dirLambdaYield->Get("Graph1D_y1");
    TGraph *g502K0sMeanPt = (TGraph*)dirK0sMeanPt->Get("Graph1D_y1");
    TGraph *g502LambdaMeanPt = (TGraph*)dirLambdaMeanPt->Get("Graph1D_y1");

    if (!g502K0sYield || !g502LambdaYield) {
        std::cout << "Could not retrieve one or both Yield graphs." << std::endl;
        return;
    }
    if (!g502K0sMeanPt || !g502LambdaMeanPt) {
        std::cout << "Could not retrieve one or both Mean Pt graphs." << std::endl;
        return;
    }

    //
    //Get 130 Results
    //
    TDirectory *dirK0sYield130 = (TDirectory*)file130->Get("Table 8c");
    TDirectory *dirLambdaYield130 = (TDirectory*)file130->Get("Table 8c");
    TDirectory *dirK0sMeanPt130 = (TDirectory*)file130->Get("Table 5c");
    TDirectory *dirLambdaMeanPt130 = (TDirectory*)file130->Get("Table 5c");

    if (!dirK0sYield130 || !dirLambdaYield130) {
        std::cout << "Could not find Yield directories." << std::endl;
        return;
    }
    if (!dirK0sMeanPt130 || !dirLambdaMeanPt130) {
        std::cout << "Could not find Mean Pt directories." << std::endl;
        return;
    }

    TGraph *g130K0sYield = (TGraph*)dirK0sYield130->Get("Graph1D_y1");
    TGraph *g130LambdaYield = (TGraph*)dirLambdaYield130->Get("Graph1D_y2");
    TGraph *g130K0sMeanPt = (TGraph*)dirK0sMeanPt130->Get("Graph1D_y1");
    TGraph *g130LambdaMeanPt = (TGraph*)dirLambdaMeanPt130->Get("Graph1D_y2");

    if (!g130K0sYield || !g130LambdaYield) {
        std::cout << "Could not retrieve one or both Yield graphs." << std::endl;
        return;
    }
    if (!g130K0sMeanPt || !g130LambdaMeanPt) {
        std::cout << "Could not retrieve one or both Mean Pt graphs." << std::endl;
        return;
    }


    // ---------------------------------------
    // Style 13.6 TeV (GREEN)
    // ---------------------------------------

    TIter next136Yield(c136Yield->GetListOfPrimitives());
    TObject *obj136Yield = nullptr;

    int counter136Yield = 0;

    while ((obj136Yield = next136Yield())) {

        if (obj136Yield->InheritsFrom(TGraph::Class())) {

            TGraph *g = (TGraph*)obj136Yield;
            graphs136Yield.push_back(g);
            g->SetLineColor(kGreen+2);
            g->SetMarkerColor(kGreen+2);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);

            if (counter136Yield == 0)
                g->SetMarkerStyle(21);   // K0s square
            else
                g->SetMarkerStyle(20);   // Lambda circle

            TF1 *fit136Yield = new TF1(Form("fit136_%d", counter136Yield), "pol1");
            g->Fit(fit136Yield, "Q0");

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

            TGraph *g = (TGraph*)obj136MeanPt;
            graphs136MeanPt.push_back(g);
            g->SetLineColor(kGreen+2);
            g->SetMarkerColor(kGreen+2);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);

            if (counter136MeanPt == 0)
                g->SetMarkerStyle(21);   // K0s square
            else
                g->SetMarkerStyle(20);   // Lambda circle

            TF1 *fit136MeanPt = new TF1(Form("fit136_%d", counter136MeanPt), "pol1");
            g->Fit(fit136MeanPt, "Q0");

            fit136MeanPt->SetLineColor(kGreen+2);
            fit136MeanPt->SetLineWidth(1);

            fits136MeanPt.push_back(fit136MeanPt);

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

            TGraph *g = (TGraph*)obj536Yield;
            graphs536Yield.push_back(g);
            g->SetLineColor(kOrange+1);
            g->SetMarkerColor(kOrange+1);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);

            if (counter536Yield == 0)
                g->SetMarkerStyle(29);   // K0s star
            else
                g->SetMarkerStyle(33);   // Lambda diamond

            TF1 *fit536Yield = new TF1(Form("fit536_%d", counter536Yield), "pol1");
            g->Fit(fit536Yield, "Q0");

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

            TGraph *g = (TGraph*)obj536MeanPt;
            graphs536MeanPt.push_back(g);
            g->SetLineColor(kOrange+1);
            g->SetMarkerColor(kOrange+1);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);

            if (counter536MeanPt == 0)
                g->SetMarkerStyle(29);   // K0s star
            else
                g->SetMarkerStyle(33);   // Lambda diamond

            TF1 *fit536MeanPt = new TF1(Form("fit536_%d", counter536MeanPt), "pol1");
            g->Fit(fit536MeanPt, "Q0");

            fit536MeanPt->SetLineColor(kOrange+1);
            fit536MeanPt->SetLineWidth(1);

            fits536MeanPt.push_back(fit536MeanPt);

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

            TGraph *g = (TGraph*)obj536YieldOfficial;

            graphs536Yield.push_back(g);

            g->SetLineColor(kBlack);
            g->SetMarkerColor(kBlack);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);

            if (counter536YieldOfficial == 0)
                g->SetMarkerStyle(29);   // K0s star
            else
                g->SetMarkerStyle(33);   // Lambda diamond

            // -----------------------------------------
            // Replace X coordinates with official <Nch>
            // -----------------------------------------
            int nPoints = g->GetN();

            for (int i = 0; i < nPoints && i < 10; i++) {

                double x, y;
                g->GetPoint(i, x, y);

                g->SetPoint(i, nchOfficial[i], y);
            }

            TF1 *fit536Yield =
                new TF1(Form("fit536Yield_%d", counter536YieldOfficial), "pol1");

            g->Fit(fit536Yield, "Q0");

            fit536Yield->SetLineColor(kBlack);
            fit536Yield->SetLineWidth(1);

            fits536YieldOfficial.push_back(fit536Yield);

            counter536YieldOfficial++;
        }
    }

    TIter next536MeanPtOfficial(c536MeanPtOfficial->GetListOfPrimitives());
    TObject *obj536MeanPtOfficial = nullptr;

    int counter536MeanPtOfficial = 0;

    while ((obj536MeanPtOfficial = next536MeanPtOfficial())) {

        if (obj536MeanPtOfficial->InheritsFrom(TGraph::Class())) {

            TGraph *g = (TGraph*)obj536MeanPtOfficial;

            graphs536MeanPt.push_back(g);

            g->SetLineColor(kBlack);
            g->SetMarkerColor(kBlack);
            g->SetLineWidth(0);
            g->SetMarkerSize(1.2);

            if (counter536MeanPtOfficial == 0)
                g->SetMarkerStyle(29);   // K0s star
            else
                g->SetMarkerStyle(33);   // Lambda diamond

            // -----------------------------------------
            // Replace X coordinates with official <Nch>
            // -----------------------------------------
            int nPoints = g->GetN();

            for (int i = 0; i < nPoints && i < 10; i++) {
                double x, y;
                g->GetPoint(i, x, y);
                g->SetPoint(i, nchOfficial[i], y);
            }

            TF1 *fit536MeanPtOfficial =new TF1(Form("fit536MeanPt_%d", counter536MeanPtOfficial), "pol1");

            g->Fit(fit536MeanPtOfficial, "Q0");

            fit536MeanPtOfficial->SetLineColor(kBlack);
            fit536MeanPtOfficial->SetLineWidth(1);

            fits536MeanPtOfficial.push_back(fit536MeanPtOfficial);

            counter536MeanPtOfficial++;
        }
    }
















    // ---------------------------------------
    // Style 5.02 TeV (BLUE)
    // ---------------------------------------

    g502K0sYield->SetMarkerStyle(22);
    g502K0sYield->SetMarkerColor(kBlue+1);
    g502K0sYield->SetLineColor(kBlue+1);
    g502K0sYield->SetMarkerSize(1.2);

    g502LambdaYield->SetMarkerStyle(23);
    g502LambdaYield->SetMarkerColor(kBlue+1);
    g502LambdaYield->SetLineColor(kBlue+1);
    g502LambdaYield->SetMarkerSize(1.2);

    g502K0sMeanPt->SetMarkerStyle(22);
    g502K0sMeanPt->SetMarkerColor(kBlue+1);
    g502K0sMeanPt->SetLineColor(kBlue+1);
    g502K0sMeanPt->SetMarkerSize(1.2);

    g502LambdaMeanPt->SetMarkerStyle(23);
    g502LambdaMeanPt->SetMarkerColor(kBlue+1);
    g502LambdaMeanPt->SetLineColor(kBlue+1);
    g502LambdaMeanPt->SetMarkerSize(1.2);

    // ---------------------------------------
    // Style 13.0 TeV (RED)
    // ---------------------------------------

    g130K0sYield->SetMarkerStyle(34);
    g130K0sYield->SetMarkerColor(kRed+1);
    g130K0sYield->SetLineColor(kRed+1);
    g130K0sYield->SetMarkerSize(1.2);

    g130LambdaYield->SetMarkerStyle(47);
    g130LambdaYield->SetMarkerColor(kRed+1);
    g130LambdaYield->SetLineColor(kRed+1);
    g130LambdaYield->SetMarkerSize(1.2);

    g130K0sMeanPt->SetMarkerStyle(34);
    g130K0sMeanPt->SetMarkerColor(kRed+1);
    g130K0sMeanPt->SetLineColor(kRed+1);
    g130K0sMeanPt->SetMarkerSize(1.2);

    g130LambdaMeanPt->SetMarkerStyle(47);
    g130LambdaMeanPt->SetMarkerColor(kRed+1);
    g130LambdaMeanPt->SetLineColor(kRed+1);
    g130LambdaMeanPt->SetMarkerSize(1.2);

    // ---------------------------------------
    // Create output canvases
    // ---------------------------------------

    TCanvas *cYieldAllResults = new TCanvas("cYield_FinalAllResults_Overlay","cYield_FinalAllResults_Overlay",800,600);
    TCanvas *cYieldAllResultsFits = new TCanvas("cYield_FinalAllResults_OverlayWithFits","cYield_FinalAllResults_OverlayWithFits",800,600);
    TCanvas *cMeanPtAllResults = new TCanvas("cMeanPt_FinalAllResults_Overlay","cMeanPt_FinalAllResults_Overlay",800,600);
    TCanvas *cMeanPtAllResultsFits = new TCanvas("cMeanPt_FinalAllResults_OverlayWithFits","cMeanPt_FinalAllResults_OverlayWithFits",800,600);

    // Draw all objects from the 13.6 TeV canvas first
    TLegend *leg = new TLegend(0.15,0.70,0.55,0.88);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);

    //-----------------------------------------
    // Set styles for each dataset
    //-----------------------------------------

    TGraph *gRedCross = new TGraph();
    gRedCross->SetMarkerStyle(34);
    gRedCross->SetMarkerColor(kRed+1);

    TGraph *gRedPlus = new TGraph();
    gRedPlus->SetMarkerStyle(47);
    gRedPlus->SetMarkerColor(kRed+1);

    TGraph *gGreenSquare = new TGraph();
    gGreenSquare->SetMarkerStyle(21);
    gGreenSquare->SetMarkerColor(kGreen+2);

    TGraph *gGreenCircle = new TGraph();
    gGreenCircle->SetMarkerStyle(20);
    gGreenCircle->SetMarkerColor(kGreen+2);

    TGraph *gOrangeStar = new TGraph();
    gOrangeStar->SetMarkerStyle(29);
    gOrangeStar->SetMarkerColor(kOrange+1);

    TGraph *gOrangeDiamond = new TGraph();
    gOrangeDiamond->SetMarkerStyle(33);
    gOrangeDiamond->SetMarkerColor(kOrange+1);

    TGraph *gBlueTriangle = new TGraph();
    gBlueTriangle->SetMarkerStyle(22);
    gBlueTriangle->SetMarkerColor(kBlue+1);

    TGraph *gBlueTriangleDown = new TGraph();
    gBlueTriangleDown->SetMarkerStyle(23);
    gBlueTriangleDown->SetMarkerColor(kBlue+1);

    leg->AddEntry(gOrangeStar,"5.36 TeV: K^{0}_{S} This Analysis","p");
    leg->AddEntry(gOrangeDiamond,"5.36 TeV: #Lambda+#bar{#Lambda} This Analysis","p");
    leg->AddEntry(gGreenSquare,"13.6 TeV: K^{0}_{S} Test","p");
    leg->AddEntry(gGreenCircle,"13.6 TeV: #Lambda+#bar{#Lambda} Test","p");
    leg->AddEntry(gRedCross,"13.0 TeV: K^{0}_{S} Run 2 ALICE","p");
    leg->AddEntry(gRedPlus,"13.0 TeV: #Lambda+#bar{#Lambda} Run 2 ALICE","p");
    leg->AddEntry(gBlueTriangle,"5.02 TeV: K^{0}_{S} Run 2 ALICE","p");
    leg->AddEntry(gBlueTriangleDown,"5.02 TeV: #Lambda+#bar{#Lambda} Run 2 ALICE","p");

    // ======================================================
    // First canvas Yields
    // ======================================================

    cYieldAllResults->cd();

    graphs136Yield[0]->SetTitle("");
    graphs136Yield[0]->GetXaxis()->SetTitle("<dN_{ch}/dy>");
    graphs136Yield[0]->GetYaxis()->SetTitle("<dN/dy>");

    graphs136Yield[0]->GetXaxis()->SetTitleSize(0.05);
    graphs136Yield[0]->GetYaxis()->SetTitleSize(0.05);
    graphs136Yield[0]->GetXaxis()->SetLabelSize(0.045);
    graphs136Yield[0]->GetYaxis()->SetLabelSize(0.045);
    graphs136Yield[0]->GetYaxis()->SetTitleOffset(1.1);
    graphs136Yield[0]->Draw("AP");

    for (size_t i=1;i<graphs136Yield.size();++i)
        graphs136Yield[i]->Draw("P SAME");

    for (auto g : graphs536Yield)
        g->Draw("P SAME");

    g502K0sYield->Draw("P SAME");
    g502LambdaYield->Draw("P SAME");
    g130K0sYield->Draw("P SAME");
    g130LambdaYield->Draw("P SAME");
    leg->Draw();

    // ======================================================
    // Second canvas Mean Pt
    // ======================================================

    cMeanPtAllResults->cd();

    graphs136MeanPt[0]->SetTitle("");
    graphs136MeanPt[0]->GetXaxis()->SetTitle("<dN_{ch}/dy>");
    graphs136MeanPt[0]->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");

    graphs136MeanPt[0]->GetXaxis()->SetTitleSize(0.05);
    graphs136MeanPt[0]->GetYaxis()->SetTitleSize(0.05);
    graphs136MeanPt[0]->GetXaxis()->SetLabelSize(0.045);
    graphs136MeanPt[0]->GetYaxis()->SetLabelSize(0.045);
    graphs136MeanPt[0]->GetYaxis()->SetTitleOffset(1.1);

    graphs136MeanPt[0]->Draw("AP");

    for (size_t i=1;i<graphs136MeanPt.size();++i)
        graphs136MeanPt[i]->Draw("P SAME");

    for (auto g : graphs536MeanPt)
        g->Draw("P SAME");

    g502K0sMeanPt->Draw("P SAME");
    g502LambdaMeanPt->Draw("P SAME");
    g130K0sMeanPt->Draw("P SAME");
    g130LambdaMeanPt->Draw("P SAME");
    leg->Draw();

    // ======================================================
    // Third canvas Yields plus fits
    // ======================================================

    cYieldAllResultsFits->cd();

    graphs136Yield[0]->Draw("AP");

    for (size_t i=1;i<graphs136Yield.size();++i)
        graphs136Yield[i]->Draw("P SAME");

    for (auto g : graphs536Yield)
        g->Draw("P SAME");

    for (auto f : fits136Yield)
        f->Draw("SAME");

    for (auto f : fits536Yield)
        f->Draw("SAME");

    TLegend *leg2 = (TLegend*)leg->Clone();
    leg2->Draw();


    g502K0sYield->Draw("P SAME");
    g502LambdaYield->Draw("P SAME");
    g130K0sYield->Draw("P SAME");
    //fit502K0s->Draw("SAME");
    //fit502Lambda->Draw("SAME");

    // ======================================================
    // Fourth canvas Mean Pt plus fits
    // ======================================================
/*
    cMeanPtAllResultsFits->cd();

    graphs136MeanPt[0]->Draw("AP");

    for (size_t i=1;i<graphs136MeanPt.size();++i)
        graphs136MeanPt[i]->Draw("P SAME");

    for (auto g : graphs536MeanPt)
        g->Draw("P SAME");

    for (auto f : fits136MeanPt)
        f->Draw("SAME");

    for (auto f : fits536MeanPt)
        f->Draw("SAME");

    TLegend *leg2 = (TLegend*)leg->Clone();
    leg2->Draw();


    g502K0sMeanPt->Draw("P SAME");
    g502LambdaMeanPt->Draw("P SAME");
    */

    //fit502K0s->Draw("SAME");
    //fit502Lambda->Draw("SAME");
    
    // ---------------------------------------
    // Save output
    // ---------------------------------------

    TFile *outFile = new TFile("FinalComparisonOverlay.root","RECREATE");

    cYieldAllResults->Write();
    cYieldAllResultsFits->Write();
    cMeanPtAllResults->Write();
    cMeanPtAllResultsFits->Write();

    outFile->Close();

    file136->Close();
    file536->Close();
    file502->Close();
}
