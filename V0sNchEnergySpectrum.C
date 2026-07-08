#include "AnalysisLibrary.cxx"

void ExtractCentralityValues(TCanvas* c, std::vector<double>& values, std::vector<std::string>& names)
{
    values.clear();
    names.clear();

    if (!c) {
        std::cerr << "ExtractCentralityValues: null canvas\n";
        return;
    }

    TIter next(c->GetListOfPrimitives());
    TObject* obj = nullptr;

    while ((obj = next())) {

        TH1* h = dynamic_cast<TH1*>(obj);
        if (!h) continue;

        if (TString(h->GetName()).BeginsWith("Sum_"))
            continue;

        names.push_back(h->GetName());
        values.push_back(h->GetBinContent(1));
    }
}

double GetExtrapolatedYield(TF1* f, double xmin, double xmax)
{
    if (!f) return 0.0;
    return f->Integral(xmin, xmax);
}

void V0sNchEnergySpectrum(){

    //Particles to be processed
   //********************************************************************************/
   //********************************************************************************/
    bool K0sh = true;
    bool LambdaAntiLambda = true;
    // ============================================================
    // DEFINING THE DIFFERENT BINNING
    // ============================================================

    std::vector<double>* nchCentBinsKaons   = &pthistos::centBinsK0sh;
    std::vector<double>* nchCentBinsLambdas = &pthistos::centBinsLambda;

    // ============================================================
    // FILES
    // ============================================================
    InitAO2DFiles();

    //Getting Nch information
    TFile* dataFile  = new TFile(AO2DFiles.datarootpath.c_str());
    TFile* mcFile    = new TFile(AO2DFiles.mcrootpath.c_str());
    TFile* mcGenFile = new TFile(AO2DFiles.mcgenrootpath.c_str());

    //Getting already calculated histos    
    TFile* processedFile  = new TFile("V0s_Pt_Spectrum_Cent_13.6_TeV.root", "READ");

    // ============================================================
    // LOAD REAL DATA
    // ============================================================

    dataFile->cd(Form("%s/NchAnalysis",AO2DFiles.datarootsubfolder.c_str()));
    TH2D* hNchData = (TH2D*)ReadAndSaveHistogram("hNchCentrality","hNchCentralityData");

    dataFile->cd(Form("%s/PtAnalysis",AO2DFiles.datarootsubfolder.c_str()));
    TH2D* hNRecEventsDataNch = (TH2D*)ReadAndSaveHistogram("hNRecEvents","hNRecEventsDataNch");

    // ============================================================
    // LOAD MC / GEN
    // ============================================================

    mcFile->cd(Form("%s/NchAnalysis",AO2DFiles.mcrootsubfolder.c_str()));
    TH2D* hNchRecoMC = (TH2D*)ReadAndSaveHistogram("hNchCentrality","hNchCentralityMC");

    mcGenFile->cd(Form("%s/NchAnalysis",AO2DFiles.mcgenrootsubfolder.c_str()));
    TH2D* hNchGenMC = (TH2D*)ReadAndSaveHistogram("hNchCentralityGenerated","hNchCentralityGeneratedNchMC");

    // ============================================================
    // Nch EFFICIENCY (AND SPLITTING (TO BE DONE)
    // ============================================================

    TH2D* hNchEfficiency = (TH2D*)hNchRecoMC->Clone("hNchEfficiency");
    hNchEfficiency->Reset();

    int nRatioAboveOne = 0;
    int nRatioBetweenZeroAndOne = 0;
    int nRatioMinus999 = 0;
    int nRatioEqualOne = 0;

    for (int ix = 1; ix <= hNchRecoMC->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hNchRecoMC->GetNbinsY(); ++iy) {

            double reco = hNchRecoMC->GetBinContent(ix, iy);
            double gen  = hNchGenMC->GetBinContent(ix, iy);

            double ratio = 0.0;

            ratio = reco / gen;

            // ONLY PRINT VALID NON-ZERO RATIOS
            if (ratio != 0.0) {
                std::cout << "Bin (" << ix << "," << iy << ")" << "  Reco = " << reco << "  Gen = " << gen << "  Ratio = " << ratio << std::endl;
            }

            if (ratio > 1.0) {
                ++nRatioAboveOne;
            }
            else if (ratio > 0.0 && ratio < 1.0) {
                ++nRatioBetweenZeroAndOne;
            }
            else if (ratio == 1.0) {
                ++nRatioEqualOne;
            }
        hNchEfficiency->SetBinContent(ix, iy, ratio);
        }
    }



    ////DEBUGGING START





TGraph* gNchMatrixMean = new TGraph();
TGraph* gNchMatrix = new TGraph();

for (int ic = 1; ic <= hNchGenMC->GetNbinsX(); ++ic) {

    // Nch distribution for this centrality bin
    TH1D* hGen  = hNchGenMC ->ProjectionY(Form("hGen_%d",ic),  ic, ic);
    TH1D* hReco = hNchRecoMC->ProjectionY(Form("hReco_%d",ic), ic, ic);

    double meanGen  = hGen->GetMean();
    double meanReco = hReco->GetMean();

    gNchMatrixMean->SetPoint(ic-1, meanGen, meanReco);

    std::cout << "Cent bin " << ic
              << "  Gen=" << meanGen
              << "  Reco=" << meanReco
              << std::endl;

    delete hGen;
    delete hReco;
}


gNchMatrixMean->SetTitle("Gen <N_{ch}>/Reco <N_{ch}>");
gNchMatrixMean->SetMarkerStyle(20);
gNchMatrixMean->SetMarkerSize(1.2);







auto& cbins = *nchCentBinsKaons;

TH1D* hRecoGenRatio = new TH1D("hRecoGenRatio",";Centrality (%);Reco/Gen",cbins.size()-1,cbins.data());

for (size_t i = 0; i < cbins.size()-1; ++i) {

    double cmin = cbins[i];
    double cmax = cbins[i+1];

    int binLo = hNchGenMC->GetXaxis()->FindBin(cmin + 1e-6);
    int binHi = hNchGenMC->GetXaxis()->FindBin(cmax - 1e-6);

    TH1D* hGen  = hNchGenMC ->ProjectionY(Form("hGen_%zu",i),  binLo, binHi);
    TH1D* hReco = hNchRecoMC->ProjectionY(Form("hReco_%zu",i), binLo, binHi);

    double meanGen  = hGen->GetMean();
    double meanReco = hReco->GetMean();

    double ratio = (meanGen > 0.) ? meanReco/meanGen : 0.;

    hRecoGenRatio->SetBinContent(i+1, ratio);

    delete hGen;
    delete hReco;
}

hRecoGenRatio->Draw("HIST");






    ////DEBUGGING END





    TCanvas* cEff = new TCanvas("cEff", "Nch Efficiency", 800, 600);

    hNchEfficiency->SetStats(0);
    hNchEfficiency->SetTitle("Nch Efficiency;X axis;Y axis");
    hNchEfficiency->Draw("COLZ");
    std::cout << "\nSummary:\n" << "Ratio > 1      : " << nRatioAboveOne << '\n' << "0 < Ratio < 1  : " << nRatioBetweenZeroAndOne << '\n' << "Ratio == 1     : " << nRatioEqualOne << '\n'<< "Ratio = -999   : " << nRatioMinus999 << '\n';

    // ============================================================
    // EVENT EFFICIENCY AND SPLITTING (ALREADY CALCULATED)
    // ============================================================

    std::vector<double> eventEfficiencyKaonBins, eventEfficiencyLambdaBins;
    std::vector<std::string> eventEfficiencyKaonBinsNames, eventEfficiencyLambdaBinsNames;
    TCanvas* cEventEffKaonBins = nullptr; 
    TCanvas* cEventSplitKaonBins = nullptr;
    TCanvas* cEventEffLambdaBins = nullptr; 
    TCanvas* cEventSplitLambdaBins = nullptr;
    std::vector<double> eventSplittingKaonBins, eventSplittingLambdaBins;
    std::vector<std::string> eventSplittingKaonBinsNames, eventSplittingLambdaBinsNames;

    if(K0sh){
        TDirectory* dirResultsK0sh = processedFile->GetDirectory("K0sh/All_Centralities");

        cEventEffKaonBins = (TCanvas*)dirResultsK0sh->Get("Event_Efficiency_AllCentralities_K0sh");
        cEventSplitKaonBins =(TCanvas*)dirResultsK0sh->Get("Event_Splitting_AllCentralities_K0sh");

        ExtractCentralityValues(cEventEffKaonBins, eventEfficiencyKaonBins,eventEfficiencyKaonBinsNames);
        ExtractCentralityValues(cEventSplitKaonBins, eventSplittingKaonBins,eventSplittingKaonBinsNames);
    }

    //SAME FOR LAMBDAS
    if(LambdaAntiLambda){
        TDirectory* dirResultsLambda = processedFile->GetDirectory("Lambda/All_Centralities");

        cEventEffLambdaBins = (TCanvas*)dirResultsLambda->Get("Event_Efficiency_AllCentralities_Lambda");
        cEventSplitLambdaBins = (TCanvas*)dirResultsLambda->Get("Event_Splitting_AllCentralities_Lambda");

        ExtractCentralityValues(cEventEffLambdaBins, eventEfficiencyLambdaBins,eventEfficiencyLambdaBinsNames);
        ExtractCentralityValues(cEventSplitLambdaBins, eventSplittingLambdaBins,eventSplittingLambdaBinsNames);
    }

    // ============================================================
    // ===================== KAON BLOCK ============================
    // ============================================================
    TProfile* pNchK0sh = nullptr;
    TH2D* hNchRecoMcK0sh = nullptr;
    TH2D* hNchGenMcK0sh = nullptr;
    TH2D* hNchEfficiencyKaon = nullptr;
    TH2D* hNchDataKaon = nullptr;
    TH2D* hNchCorrectedKaonBins = nullptr;
    TH1D* hNchK0shBins = nullptr;

    if(K0sh){

        hNchEfficiencyKaon = (TH2D*)hNchEfficiency->Clone("hNchEfficiencyKaon");

        hNchDataKaon = RebinTH2D(hNchData,"hNchDataKaon",nchCentBinsKaons);
        hNchCorrectedKaonBins = (TH2D*)hNchDataKaon->Clone("hNchCorrectedKaonBins");

        TH1D* hNRecEventsKaon = hNRecEventsDataNch->ProjectionY("hRecEvKaon");
        hNRecEventsKaon = (TH1D*)hNRecEventsKaon->Rebin(nchCentBinsKaons->size() - 1,"hRecEvKaonReb",nchCentBinsKaons->data());

        for (int i = 1; i <= hNchCorrectedKaonBins->GetNbinsX(); ++i){
            double nev        = hNRecEventsKaon->GetBinContent(i);
            double eventEff   = eventEfficiencyKaonBins[i-1];
            double eventSplit = eventSplittingKaonBins[i-1];
            double eventWeight = (eventEff * eventSplit) / nev;

            std::cout << "Bin " << i << "  nev = " << nev<< "  eventEff = " << eventEff << "  eventSplit = " << eventSplit << "  eventWeight = " << eventWeight << std::endl;

            for (int j = 1; j <= hNchCorrectedKaonBins->GetNbinsY(); ++j){
                double effReco = hNchEfficiencyKaon->GetBinContent(i, j);

                double v = hNchCorrectedKaonBins->GetBinContent(i, j);
                double e = hNchCorrectedKaonBins->GetBinError(i, j);

                // 1D event-level correction last
                v *= eventWeight;
                e *= eventWeight;

                hNchCorrectedKaonBins->SetBinContent(i, j, v);
                hNchCorrectedKaonBins->SetBinError(i, j, e);
                std::cout << "Corrected Value" << "Bin [" << i << "," << j<< "]" << "  v = " << v<< "  e = " << e << std::endl;
            }
        }
        pNchK0sh = hNchCorrectedKaonBins->ProfileX("pNchK0sh");

        hNchK0shBins = pNchK0sh->ProjectionX("hNchK0shBins");

        for (int i = 1; i <= hNchK0shBins->GetNbinsX(); ++i) {
            double nchBefore = hNchK0shBins->GetBinContent(i);
            double errBefore = hNchK0shBins->GetBinError(i);
            double corr      = hRecoGenRatio->GetBinContent(i);
            double nchAfter = nchBefore / corr;
            hNchK0shBins->SetBinContent(i, nchAfter);
            std::cout<< "Bin " << i << "  <Nch> before = " << nchBefore<< "  err = " << errBefore<< "  corr = " << corr << "  <Nch> after = " << nchAfter<< std::endl;
        }
    }

    // ============================================================
    // ===================== LAMBDA BLOCK =========================
    // ============================================================
    TProfile* pNchLambdaAntiLambda = nullptr;
    TH2D* hNchRecoMCLambda = nullptr;
    TH2D* hNchGenMCLambda = nullptr;
    TH2D* hNchEfficiencyLambda = nullptr;
    TH2D* hNchDataLambda = nullptr;
    TH2D* hNchCorrectedLambdaBins = nullptr;
    TH1D* hNchLambdaAntiLambdaBins = nullptr;

    if(LambdaAntiLambda){
        
        hNchEfficiencyLambda = (TH2D*)hNchEfficiency->Clone("hNchEfficiencyLambda");

        hNchDataLambda = RebinTH2D(hNchData,"hNchDataLambda",nchCentBinsLambdas);
        hNchCorrectedLambdaBins = (TH2D*)hNchDataLambda->Clone("hNchCorrectedLambdaBins");

        TH1D* hNRecEventsLambda = hNRecEventsDataNch->ProjectionY("hRecEvLambda");
        hNRecEventsLambda = (TH1D*)hNRecEventsLambda->Rebin(nchCentBinsLambdas->size() - 1,"hRecEvLambdaReb",nchCentBinsLambdas->data());

        for (int i = 1; i <= hNchCorrectedLambdaBins->GetNbinsX(); ++i){
            double nev        = hNRecEventsLambda->GetBinContent(i);
            double eventEff   = eventEfficiencyLambdaBins[i-1];
            double eventSplit = eventSplittingLambdaBins[i-1];
            double eventWeight = (eventEff * eventSplit) / nev;

            std::cout << "Bin " << i << "  nev = " << nev<< "  eventEff = " << eventEff << "  eventSplit = " << eventSplit << "  eventWeight = " << eventWeight << std::endl;

            for (int j = 1; j <= hNchCorrectedLambdaBins->GetNbinsY(); ++j){
                double effReco = hNchEfficiencyLambda->GetBinContent(i, j);

                double v = hNchCorrectedLambdaBins->GetBinContent(i, j);
                double e = hNchCorrectedLambdaBins->GetBinError(i, j);

                // 1D event-level correction last
                v *= eventWeight;
                e *= eventWeight;

                hNchCorrectedLambdaBins->SetBinContent(i, j, v);
                hNchCorrectedLambdaBins->SetBinError(i, j, e);
                std::cout << "Corrected Value" << "Bin [" << i << "," << j<< "]" << "  v = " << v<< "  e = " << e << std::endl;
            }
        }
        pNchLambdaAntiLambda = hNchCorrectedLambdaBins->ProfileX("pNchLambdaAntiLambda");

        hNchLambdaAntiLambdaBins = pNchLambdaAntiLambda->ProjectionX("hNchLambdaAntiLambdaBins");

        for (int i = 1; i <= hNchLambdaAntiLambdaBins->GetNbinsX(); ++i) {
            double nchBefore = hNchLambdaAntiLambdaBins->GetBinContent(i);
            double errBefore = hNchLambdaAntiLambdaBins->GetBinError(i);
            double corr      = hRecoGenRatio->GetBinContent(i);
            double nchAfter = nchBefore / corr;
            hNchLambdaAntiLambdaBins->SetBinContent(i, nchAfter);
            std::cout<< "Bin " << i << "  <Nch> before = " << nchBefore<< "  err = " << errBefore<< "  corr = " << corr << "  <Nch> after = " << nchAfter<< std::endl;
        }
    }



    // ============================================================
    // OUTPUT OF NCH VS CENTRALITY RESULTS
    // ============================================================

    TFile* NchOutputFile = new TFile("V0s_Nch_EnergySpectrum.root", "RECREATE");

    if (!NchOutputFile || NchOutputFile->IsZombie()) {
        std::cerr << "ERROR: output file not created!" << std::endl;
        return;
    }

    //If we want directories
    //TDirectory* dirNch = NchOutputFile->mkdir("NchAnalysis");
    //dirNch->cd();

    if(K0sh){
        hNchK0shBins->Write();
        cEff->Write();
        cEventEffKaonBins->Write();
        cEventSplitKaonBins->Write();
        hNchEfficiency->Write();
        hNchData->Write();
        hNchEfficiencyKaon->Write();
        hNchDataKaon->Write();
        hNchCorrectedKaonBins->Write();
    }

    if(LambdaAntiLambda){
        hNchLambdaAntiLambdaBins->Write();
        cEff->Write();
        cEventEffLambdaBins->Write();
        cEventSplitLambdaBins->Write();
        hNchEfficiency->Write();
        hNchEfficiencyLambda->Write();
        hNchDataLambda->Write();
        hNchCorrectedLambdaBins->Write();
    }
    gNchMatrix->Write();
    gNchMatrixMean->Write();
    hRecoGenRatio->Write();

    // ============================================================
    // KAON & LAMBDA+ANTILAMBDA VS CENTRALITY RESULTS
    // ============================================================

    TCanvas* cK0shPtSpectra = new TCanvas("cK0shPtSpectra", "K0sh Pt Spectra", 800, 600);
    TCanvas* cK0shCentrality = new TCanvas("cK0shCentrality", "K0sh Centrality", 800, 600);
    TCanvas* cLambdaAntiLambdaPtSpectra = new TCanvas("cLambdaAntiLambdaPtSpectra", "Lambda+Antilambda Pt Spectra", 800, 600);
    TCanvas* cLambdaAntiLambdaCentrality = new TCanvas("cLambdaAntiLambdaCentrality", "Lambda+Antilambda Centrality", 800, 600);

    TH1D* hK0shYieldCent = nullptr;
    TH1D* hK0shExtrapFrac = new TH1D("hK0shExtrapFrac","K^{0}_{S} Extrapolation Fraction vs Centrality;Centrality (%);Fraction",nchCentBinsKaons->size() - 1,nchCentBinsKaons->data());
    TH1D* hLambdaAntiLambdaYieldCent = nullptr;
    TH1D* hLambdaAntiLambdaExtrapFrac = new TH1D("hLambdaAntiLambdaExtrapFrac","Lambda+Antilambda Extrapolation Fraction vs Centrality;Centrality (%);Fraction",nchCentBinsLambdas->size() - 1,nchCentBinsLambdas->data());

    TGraph* gK0shNch = nullptr;
    TGraph* gLambdaAntiLambdaNch = nullptr;

    std::map<std::string, TF1*> fitMapK0sh;
    std::map<std::string, TF1*> fitMapLambdaAntiLambda;

    if(K0sh){
        TDirectory* dirResultsK0sh = processedFile->GetDirectory("Comparisons");

        cK0shPtSpectra = (TCanvas*)dirResultsK0sh->Get("cPt_Spectrum_Fit_Comparison_K0sh");      //This will change to cPt_Spectrum_Fit_Comparison_K0sh       
        TPad* pad1 = (TPad*)cK0shPtSpectra->FindObject("pad1");
        if (!pad1) {
            std::cout << "Pad1 not found!" << std::endl;
            return;
        }

        //Getting the Histograms and Fits for yield extrapolation
        std::map<std::string, TH1F*> hMapK0sh;
        std::map<std::string, TF1*>  fMapK0sh;

        TIter next(pad1->GetListOfPrimitives());
        TObject* obj = nullptr;
        while ((obj = next())) {

            TString name = obj->GetName();

            // remove sum objects
            if (name == "hSum_combined") continue;
            if (name == "levySum") continue;

            // clean key
            TString key = name;
            key.ReplaceAll("draw_Corrected Pt Spectrum_", "");
            key.ReplaceAll("draw_Levy-Tsallis Fit_", "");

            if (obj->InheritsFrom(TH1::Class())) {
                hMapK0sh[key.Data()] = (TH1F*)obj;
                std::cout << "HIST key: " << key << std::endl;
            }
            else if (obj->InheritsFrom(TF1::Class())) {
                fMapK0sh[key.Data()] = (TF1*)obj;
                std::cout << "FIT key: " << key << std::endl;
            }
        }
        TH1D* hK0shYieldCent = new TH1D("hK0shYieldCent","K^{0}_{S} Yield vs Centrality;Centrality (%);Integral",nchCentBinsKaons->size() - 1,nchCentBinsKaons->data());
        //TIter next(cK0shPtSpectra->GetListOfPrimitives());
        //TObject* obj = nullptr;
        std::vector<TString> centKeysK0sh = {"0_1","1_5","5_10","10_15","15_20","20_30","30_40","40_50","50_70","70_100"};

        for (auto& [key, hPtSpectra] : hMapK0sh) {

            TF1* fLevyTsallis = nullptr;

            auto it = fMapK0sh.find(key);
            if (it != fMapK0sh.end())
                fLevyTsallis = it->second;

            std::cout << "Processing centrality: " << key << std::endl;

            double yield_measured = hPtSpectra->Integral("width");

            double pT_min = 0.0;
            double pT_max = 15.0;

            double yield_low = 0.0;
            double yield_high = 0.0;

            if (fLevyTsallis) {
                yield_low  = fLevyTsallis->Integral(0.0, pT_min);
                yield_high = fLevyTsallis->Integral(pT_max, 30.0); //LevyTsallis only goes up to 14 GeV/c currently
            }
            double yield_extrap = yield_low + yield_high;
            double yield_total  = yield_measured + yield_extrap;
            double frac_extrap = (yield_total > 0) ? yield_extrap / yield_total : 0.0;
            
            auto pos = std::find(centKeysK0sh.begin(), centKeysK0sh.end(), key);

            if (pos == centKeysK0sh.end())
                continue;

            int iCent = std::distance(centKeysK0sh.begin(), pos) + 1;

            hK0shYieldCent->SetBinContent(iCent, yield_total);
            hK0shExtrapFrac->SetBinContent(iCent, frac_extrap);

            if(yield_total > 0) {
                std::cout << "Centrality Bin " << iCent << "  Measured Yield = " << yield_measured << "  Extrapolated Yield = " << yield_extrap << "  Total Yield = " << yield_total << std::endl;
            }
            else {
                std::cout << "Centrality Bin " << iCent << "  Measured Yield = " << yield_measured << "  Extrapolated Yield = " << yield_extrap << "  Total Yield = " << yield_total << " (Warning: Total Yield is zero or negative!)" << std::endl;
            }
            ++iCent;
        }
        TCanvas* cK0shYield = new TCanvas("cK0shYield","Yield vs Centrality",800, 600);
        cK0shYield->Divide(1,2);
        cK0shYield->cd(1);
        hK0shYieldCent->SetMarkerStyle(kFullCircle);
        hK0shYieldCent->SetMarkerSize(1.2);
        hK0shYieldCent->Draw("P");
        cK0shYield->cd(2);
        hK0shExtrapFrac->SetMarkerStyle(kFullCircle);
        hK0shExtrapFrac->SetMarkerSize(1.2);
        hK0shExtrapFrac->Draw("P");
        cK0shYield->Write();

        if(pNchK0sh->GetNbinsX() == hK0shYieldCent->GetNbinsX()){
            int nBins = hK0shYieldCent->GetNbinsX();
            std::vector<double> x, y;

            for (int i = 1; i <= nBins; i++) {

                double xval = hNchK0shBins->GetBinContent(i);
                double yval = hK0shYieldCent->GetBinContent(i);

                double yerr = hK0shYieldCent->GetBinError(i);

                x.push_back(xval);
                y.push_back(yval);
            }

            gK0shNch =new TGraph(x.size(),x.data(), y.data());

            gK0shNch->SetTitle("K^{0}_{S} Yield vs N_{ch};N_{ch};Yield");
            gK0shNch->SetMarkerStyle(kFullCircle);

            TCanvas* c = new TCanvas("cK0shYieldVsNch","",800,600);
            gK0shNch->Draw("AP");
            gK0shNch->Write();

        }
        else {
            std::cerr << "Error: pNchK0sh and hK0shYieldCent have different number of bins!" << std::endl;
        }
    }

    // ============================================================
    // ===================== LAMBDA BLOCK =========================
    // ============================================================
    cout<< "Processing Lambda+Antilambda..." << std::endl;
    if(LambdaAntiLambda){

        TDirectory* dirResultsLambdaAntiLambda = processedFile->GetDirectory("Comparisons");
        cLambdaAntiLambdaPtSpectra = (TCanvas*)dirResultsLambdaAntiLambda->Get("cLambdaAntiLambdaSpecFit");      //This will change to cPt_Spectrum_Fit_Comparison_Lambda       

       TPad* pad1 = (TPad*)cLambdaAntiLambdaPtSpectra->FindObject("pad1");
        if (!pad1) {
            std::cout << "Pad1 not found!" << std::endl;
            return;
        }
        //Getting the Histograms and Fits for yield extrapolation
        std::map<std::string, TH1F*> hMapLambdaAntiLambda;
        std::map<std::string, TF1*>  fMapLambdaAntiLambda;

        TIter next(pad1->GetListOfPrimitives());
        TObject* obj = nullptr;
        while ((obj = next())) {

            TString name = obj->GetName();
            std::cout << "Object name = " << name << std::endl;

            // remove sum objects
            if (name == "hSum_Lambda") continue;
            if (name == "levySum") continue;

            // clean key
            TString key = name;
            key.ReplaceAll("hLambdaPlusAntiLambda_cent_", "");
            key.ReplaceAll("levy_cent_", "");

            if (obj->InheritsFrom(TH1::Class())) {
                hMapLambdaAntiLambda[key.Data()] = (TH1F*)obj;
                std::cout << "HIST key: " << key << std::endl;
            }
            else if (obj->InheritsFrom(TF1::Class())) {
                fMapLambdaAntiLambda[key.Data()] = (TF1*)obj;
                std::cout << "FIT key: " << key << std::endl;
            }
        }

        TH1D* hLambdaAntiLambdaYieldCent = new TH1D("hLambdaAntiLambdaYieldCent","#Lambda + #bar{#Lambda} Yield vs Centrality;Centrality (%);Integral",nchCentBinsLambdas->size() - 1,nchCentBinsLambdas->data());
        //TIter next(cK0shPtSpectra->GetListOfPrimitives());
        //TObject* obj = nullptr;
        std::vector<TString> centKeysLambdaAntiLambda ={"0","1","2","3","4","5","6","7","8","9"};

        for (auto& [key, hPtSpectra] : hMapLambdaAntiLambda) {
            TF1* fLevyTsallis = nullptr;

            auto it = fMapLambdaAntiLambda.find(key);
            if (it != fMapLambdaAntiLambda.end())
                fLevyTsallis = it->second;

            std::cout << "Processing centrality: " << key << std::endl;

            double yield_measured = hPtSpectra->Integral("width");

            double pT_min = 0.5;
            double pT_max = 10.0;

            double yield_low = 0.0;
            double yield_high = 0.0;

            if (fLevyTsallis) {
                yield_low  = fLevyTsallis->Integral(0.0, pT_min);
                yield_high = fLevyTsallis->Integral(pT_max, 30.0); //LevyTsallis only goes up to 14 GeV/c currently
            }
            double yield_extrap = yield_low + yield_high;
            double yield_total  = yield_measured + yield_extrap;
            double frac_extrap = (yield_total > 0) ? yield_extrap / yield_total : 0.0;
            
            auto pos = std::find(centKeysLambdaAntiLambda.begin(), centKeysLambdaAntiLambda.end(), key);

            if (pos == centKeysLambdaAntiLambda.end())
                continue;

            int iCent = std::distance(centKeysLambdaAntiLambda.begin(), pos) + 1;

            hLambdaAntiLambdaYieldCent->SetBinContent(iCent, yield_total);
            hLambdaAntiLambdaExtrapFrac->SetBinContent(iCent, frac_extrap);

            if(yield_total > 0) {
                std::cout << "Centrality Bin " << iCent << "  Measured Yield = " << yield_measured << "  Extrapolated Yield = " << yield_extrap << "  Total Yield = " << yield_total << std::endl;
            }
            else {
                std::cout << "Centrality Bin " << iCent << "  Measured Yield = " << yield_measured << "  Extrapolated Yield = " << yield_extrap << "  Total Yield = " << yield_total << " (Warning: Total Yield is zero or negative!)" << std::endl;
            }
            ++iCent;
        }
        
        TCanvas* cLambdaAntiLambdaYield = new TCanvas("cLambdaAntiLambdaYield","Yield vs Centrality",800, 600);
        cLambdaAntiLambdaYield->Divide(1,2);
        cLambdaAntiLambdaYield->cd(1);
        hLambdaAntiLambdaYieldCent->SetMarkerStyle(kFullCircle);
        hLambdaAntiLambdaYieldCent->SetMarkerSize(1.2);
        hLambdaAntiLambdaYieldCent->Draw("P");
        cLambdaAntiLambdaYield->cd(2);
        hLambdaAntiLambdaExtrapFrac->SetMarkerStyle(kFullCircle);
        hLambdaAntiLambdaExtrapFrac->SetMarkerSize(1.2);
        hLambdaAntiLambdaExtrapFrac->Draw("P");
        cLambdaAntiLambdaYield->Write();

        if(hNchLambdaAntiLambdaBins->GetNbinsX() == hLambdaAntiLambdaYieldCent->GetNbinsX()){
            int nBins = hLambdaAntiLambdaYieldCent->GetNbinsX();
            std::vector<double> x, y;

            for (int i = 1; i <= nBins; i++) {

                double xval = hNchLambdaAntiLambdaBins->GetBinContent(i);
                double yval = hLambdaAntiLambdaYieldCent->GetBinContent(i);

                double yerr = hLambdaAntiLambdaYieldCent->GetBinError(i);

                x.push_back(xval);
                y.push_back(yval);
            }

            gLambdaAntiLambdaNch =new TGraph(x.size(),x.data(), y.data());

            gLambdaAntiLambdaNch->SetTitle("#Lambda + #bar{#Lambda} Yield vs N_{ch};N_{ch};Yield");
            gLambdaAntiLambdaNch->SetMarkerStyle(kFullCircle);

            TCanvas* c = new TCanvas("cLambdaAntiLambdaYieldVsNch","",800,600);
            gLambdaAntiLambdaNch->Draw("AP");
            gLambdaAntiLambdaNch->Write();

        }
        else {
            std::cerr << "Error: hNchLambdaAntiLambdaBins and hLambdaAntiLambdaYieldCent have different number of bins!" << std::endl;
        }
    }


    // ============================================================
    // FINAL PLOTS OF KAONS AND LAMBDAS VS NCH
    // ============================================================


    if(K0sh && LambdaAntiLambda){
        TCanvas* cFinalAllParticles = new TCanvas("cFinalAllParticles","Particles vs Nch",800,600);

        // --- Kaons
        gK0shNch->SetMarkerStyle(kFullCircle);
        gK0shNch->SetMarkerColor(kBlue);
        gK0shNch->SetLineColor(kBlue);
        gK0shNch->SetMarkerSize(1.2);

        // --- Lambdas
        gLambdaAntiLambdaNch->SetMarkerStyle(kFullSquare);
        gLambdaAntiLambdaNch->SetMarkerColor(kRed);
        gLambdaAntiLambdaNch->SetLineColor(kRed);
        gLambdaAntiLambdaNch->SetMarkerSize(1.2);

        // Draw Kaons first (sets axes)
        gK0shNch->SetTitle(";N_{ch};Yield");
        gK0shNch->Draw("AP");
        gK0shNch->SetMinimum(0);

        // Overlay Lambdas
        gLambdaAntiLambdaNch->Draw("P SAME");

        TLegend* leg = new TLegend(0.15, 0.73, 0.30, 0.88);

        leg->AddEntry(gK0shNch, "K^{0}_{S}", "p");
        leg->AddEntry(gLambdaAntiLambdaNch, "#Lambda + #bar{#Lambda}", "p");
        leg->Draw();

        TLatex *latexFinal = new TLatex();
        latexFinal->SetNDC();
        latexFinal->SetTextSize(0.03);
        latexFinal->DrawLatex(0.15, 0.70, "pp");
        TString textlatexFinal = "#sqrt{#it{s}} = 13.6 TeV";
        latexFinal->DrawLatex(0.15, 0.65, textlatexFinal.Data());
        latexFinal->DrawLatex(0.15, 0.60, "|y| < 0.5");
        cFinalAllParticles->Write();
    }




    delete dataFile;
    delete mcFile;
    delete mcGenFile;
    delete processedFile;
    NchOutputFile->Close();

    //TO ADD: Kaon vs Nch Plot
    //TO ADD: Lambda vs Nch Plot
    //TO ADD: Kaon & Lambda vs Nch Plot
}


/*
cout<<"[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][][]"<<endl;
cout<<"[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][][]"<<endl;
cout<<"[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][][]"<<endl;
cout<<"[][][][][][][][][][][][][][][][][][][][]All good up to here 1[][][][][][][][][][][][][][][][][][][][]"<<endl;
cout<<"[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][][]"<<endl;
cout<<"[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][][]"<<endl;
cout<<"[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][][]"<<endl;
*/