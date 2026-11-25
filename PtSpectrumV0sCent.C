#include "AnalysisLibrary.cxx"

void CalculateSignal(const std::string& particle, const std::string& histSuffix, const std::vector<std::string>& ptbins, std::vector<double_t>& signal, TFile* file, TDirectory* centdir, int mcflag, std::vector<double_t>& mc_sig_low, std::vector<double_t>& mc_sig_high, std::vector<double_t>& peakpt, std::vector<double_t>& widthpt, int nhistos, float cent_low, float cent_high, TH2D* hStatError, TH2D* hRelStatError, std::vector<double_t>& reducedChi2, const std::string& fitFunction, const std::string& bkgOrder, float systInvMassRangeScale){
    Double_t M_low = 0.45;
    Double_t M_high = 0.55;

   // Build histname from particle base and suffix, e.g. hK0shPt_vs_Cent_from
    std::string histname = Form("h%s%s", particle.c_str(), histSuffix.c_str());

   // initialize and convert tokenized strings into vector of doubles
    std::vector<std::string> histptbinsnames(nhistos + 1);                                        
    for (int i = 0; i < nhistos+1; i++) {            //Histos won't accept "." character so converting it to "_"
        std::string histptbin = ptbins[i];          // getting the value of the bin edge
        size_t pos = histptbin.find(".");          // finding the "." character
        histptbin[pos] = '_';                      // changing the "." character of thestring-value to a "_"
        histptbinsnames[i] = histptbin;             // filling bin edges list
    }
    for(int i=0;i<nhistos;i++){

        TH2F* h2D=(TH2F*)gDirectory->Get(fmt::format("{}_{}_to_{}", histname, histptbinsnames[i],histptbinsnames[i+1]).data());
        if(h2D == nullptr) {
            std::cout << "Trying to read histogram:" << std::endl;
            std::cout << "  histname = " << histname << std::endl;
            std::cout << "  pt bin i = " << histptbinsnames[i] << std::endl;
            std::cout << "  pt bin i+1 = " << histptbinsnames[i+1] << std::endl;
            std::cout << "  Current directory = " << gDirectory->GetPath() << std::endl;
            std::cout << "h2D was not read" << std::endl; 
            exit(1);
        }
      // Convert centrality ranges to Y-axis bins of the TH2
        int y1 = h2D->GetYaxis()->FindBin(cent_low + 1e-6);
        int y2 = h2D->GetYaxis()->FindBin(cent_high - 1e-6);
   
      // Create a unique name for this slice
        std::string hname_slice =fmt::format("{}_{}_to_{}_cent{}_to{}",histname,histptbinsnames[i], histptbinsnames[i+1],(int)cent_low, (int)cent_high);
        
        TH1D* h1 = h2D->ProjectionX(hname_slice.c_str(), y1, y2);
        if (!h1 || h1->GetEntries() == 0) {
            std::cout << "Empty slice: pT "  << histptbinsnames[i] << "–" << histptbinsnames[i+1]<< ", cent " << cent_low << "–" << cent_high << "\n";
            continue;
        }

        h1->Sumw2();
        Double_t norm = h1->GetEntries();

      // Choose particle-dependent mass and label for fit
        double particleMass = 0.497611; // default K0s
        TString latexTitle = "K^{0}_{s} Invariant Mass";
        if (TString(particle).Contains("AntiLambda", TString::kIgnoreCase)) {
            particleMass = 1.115683;
            latexTitle = "\\bar{\\Lambda} Invariant Mass";
            M_low = 1.08; M_high = 1.16;
        } else if (TString(particle).Contains("Lambda", TString::kIgnoreCase)) {
            particleMass = 1.115683;
            latexTitle = "\\Lambda Invariant Mass";
            M_low = 1.08; M_high = 1.16;
        }

        // -------------------------------------------------------------
        // Construct Signal (+ Optional Background) Formula dynamically
        // -------------------------------------------------------------
        std::string sigFormula = "";
        int bkgStartIndex = 0;

        cout<< "Fit Function: " << fitFunction << " | Background Order: " << bkgOrder << " | mcflag: " << mcflag << endl;

        if (fitFunction == "Voigt") {
            sigFormula = "[0]*TMath::Voigt(x - [1], [2], [3])";
            bkgStartIndex = 4;
        } else if (fitFunction == "Gaussian") {
            sigFormula = "[0]*TMath::Gaus(x, [1], [2])";
            bkgStartIndex = 3; 
        } else if (fitFunction == "Double-Gaussian") {
            sigFormula = "[0]*( [4]*TMath::Gaus(x,[1],[2]) + (1-[4])*TMath::Gaus(x,[1],[3]) )";
            bkgStartIndex = 5;
        } else if (fitFunction == "Breit-Wigner") {
            sigFormula = "[0]*TMath::BreitWigner(x, [1], [2])";
            bkgStartIndex = 3;
        } else if (fitFunction == "Bukin") {
            sigFormula = "ROOT::Math::BukinPdf(x, [0], [1], [2], [3], [4])";
            bkgStartIndex = 5;
        } else {
            sigFormula = "[0]*TMath::Voigt(x - [1], [2], [3])";
            bkgStartIndex = 4;
        }

        // Only add background formula if mcflag == 0
        std::string fullFormula = sigFormula;
        if (mcflag == 0) {
            fullFormula += "+" + bkgOrder + "(" + std::to_string(bkgStartIndex) + ")";
        }

        TF1* ff = new TF1("fitFunction", fullFormula.c_str(), M_low, M_high);
        ff->SetNpx(1000);

        // -------------------------------------------------------------
        // Set Initial Parameters and Limits based on Fit Function
        // -------------------------------------------------------------
        bool isLambda = TString(particle).Contains("Lambda", TString::kIgnoreCase);

        if (fitFunction == "Voigt") {
            if (isLambda) {
                double init_sigma = 0.006;
                double init_gamma = 0.006;
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, init_sigma, init_gamma, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm","mean","sigma","gamma", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, init_sigma, init_gamma);
                    ff->SetParNames("norm","mean","sigma","gamma");
                }
                ff->SetParLimits(1, particleMass - 0.005, particleMass + 0.005);
                ff->SetParLimits(2, 0.001, 0.02);
                ff->SetParLimits(3, 0.002, 0.02);
            } else {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 0.003, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm","mean","sigma","gamma", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 0.003);
                    ff->SetParNames("norm","mean","sigma","gamma");
                }
                ff->SetParLimits(1, 0.4965, 0.4985);
                ff->SetParLimits(2, 0.001, 0.008);
                ff->SetParLimits(3, 0.001, 0.008);
            }
        } else if (fitFunction == "Gaussian") {
            if (isLambda) {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm","mean","sigma", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006);
                    ff->SetParNames("norm","mean","sigma");
                }
                ff->SetParLimits(1, particleMass - 0.005, particleMass + 0.005);
                ff->SetParLimits(2, 0.001, 0.02);
            } else {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm","mean","sigma", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003);
                    ff->SetParNames("norm","mean","sigma");
                }
                ff->SetParLimits(1, 0.4965, 0.4985);
                ff->SetParLimits(2, 0.001, 0.008);
            }
        } else if (fitFunction == "Double-Gaussian") {
            if (isLambda) {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.004, 0.009, 0.5, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm", "mean", "sigma1", "sigma2", "frac", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.004, 0.009, 0.5);
                    ff->SetParNames("norm", "mean", "sigma1", "sigma2", "frac");
                }
                ff->SetParLimits(1, particleMass - 0.005, particleMass + 0.005);
                ff->SetParLimits(2, 0.001, 0.015);
                ff->SetParLimits(3, 0.003, 0.03);
                ff->SetParLimits(4, 0.0, 1.0);
            } else {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.002, 0.005, 0.5, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm", "mean", "sigma1", "sigma2", "frac", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.002, 0.005, 0.5);
                    ff->SetParNames("norm", "mean", "sigma1", "sigma2", "frac");
                }
                ff->SetParLimits(1, 0.4965, 0.4985);
                ff->SetParLimits(2, 0.001, 0.008);
                ff->SetParLimits(3, 0.002, 0.015);
                ff->SetParLimits(4, 0.0, 1.0);
            }
        } else if (fitFunction == "Breit-Wigner") {
            if (isLambda) {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm", "mean", "width", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006);
                    ff->SetParNames("norm", "mean", "width");
                }
                ff->SetParLimits(1, particleMass - 0.005, particleMass + 0.005);
                ff->SetParLimits(2, 0.002, 0.05);
            } else {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm", "mean", "width", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003);
                    ff->SetParNames("norm", "mean", "width");
                }
                ff->SetParLimits(1, 0.4965, 0.4985);
                ff->SetParLimits(2, 0.001, 0.02);
            }
        } else if (fitFunction == "Bukin") {
            if (isLambda) {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006, 0.0, 0.0, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm", "mean", "sigma", "xi", "rho", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006, 0.0, 0.0);
                    ff->SetParNames("norm", "mean", "sigma", "xi", "rho");
                }
                ff->SetParLimits(1, particleMass - 0.005, particleMass + 0.005);
                ff->SetParLimits(2, 0.001, 0.02);
                ff->SetParLimits(3, -0.5, 0.5);
                ff->SetParLimits(4, -0.5, 0.5);
            } else {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 0.0, 0.0, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm", "mean", "sigma", "xi", "rho", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 0.0, 0.0);
                    ff->SetParNames("norm", "mean", "sigma", "xi", "rho");
                }
                ff->SetParLimits(1, 0.4965, 0.4985);
                ff->SetParLimits(2, 0.001, 0.008);
                ff->SetParLimits(3, -0.5, 0.5);
                ff->SetParLimits(4, -0.5, 0.5);
            }
        } else {
            if (isLambda) {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006, 0.006, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm","mean","sigma","gamma", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.006, 0.006);
                    ff->SetParNames("norm","mean","sigma","gamma");
                }
                ff->SetParLimits(1, particleMass - 0.005, particleMass + 0.005);
                ff->SetParLimits(2, 0.001, 0.02);
                ff->SetParLimits(3, 0.002, 0.02);
            } else {
                if (mcflag == 0) {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 0.003, 3e-3, 0, 1e-3);
                    ff->SetParNames("norm","mean","sigma","gamma", "p_{0}", "p_{1}", "p_{2}");
                } else {
                    ff->SetParameters(1e-2*norm, particleMass, 0.003, 0.003);
                    ff->SetParNames("norm","mean","sigma","gamma");
                }
                ff->SetParLimits(1, 0.4965, 0.4985);
                ff->SetParLimits(2, 0.001, 0.008);
                ff->SetParLimits(3, 0.001, 0.008);
            }
        }

        // Apply background parameter limits only if mcflag == 0
        if (mcflag == 0) {
            if (bkgOrder == "pol1" || bkgOrder == "pol2" || bkgOrder == "pol3") {
                ff->SetParLimits(bkgStartIndex, 0.0, 1e9);     // p0 must be >= 0
            }
        }

        cout << "Invariant Mass Plot to be fitted Pt ranges: " << histptbinsnames[i] << " to " << histptbinsnames[i+1] << " and Centrality: " << cent_low << " to " << cent_high << " with Fit Function: " << fitFunction << endl;
        TFitResultPtr fitResult = h1->Fit(ff, "S B");

        Double_t *pars = ff->GetParameters(); 
        Double_t mean = pars[1];
        Double_t sigma = std::abs(pars[2]);
        
        Double_t widthVal = 0.0;
        if (fitFunction == "Voigt") {
            Double_t gamma = std::abs(pars[3]);
            widthVal = 0.5346 * gamma + std::sqrt(0.2166 * std::pow(gamma, 2) + std::pow(sigma, 2));
        } else if (fitFunction == "Gaussian") {
            widthVal = 3.0 * sigma;
        } else if (fitFunction == "Double-Gaussian") {
            Double_t frac = pars[4];
            Double_t s1 = std::abs(pars[2]);
            Double_t s2 = std::abs(pars[3]);
            Double_t sigEff = std::sqrt(frac * std::pow(s1, 2) + (1.0 - frac) * std::pow(s2, 2));
            widthVal = 3.0 * sigEff;
        } else if (fitFunction == "Breit-Wigner") {
            widthVal = sigma; 
        } else if (fitFunction == "Bukin") {
            widthVal = 3.0 * sigma; 
        } else {
            Double_t gamma = std::abs(pars[3]);
            widthVal = 0.5346 * gamma + std::sqrt(0.2166 * std::pow(gamma, 2) + std::pow(sigma, 2));
        }

        Double_t sig_low = mean - (widthVal * systInvMassRangeScale); 
        Double_t sig_high = mean + (widthVal * systInvMassRangeScale);

        if (mcflag == 0) {
            mc_sig_low[i] = sig_low;
            mc_sig_high[i] = sig_high;
        }
        if (mcflag == 1) {
            sig_low = mc_sig_low[i];
            sig_high = mc_sig_high[i];
        }

        // Fit Diagnostics
        Double_t chi2 = fitResult->Chi2();
        Int_t ndf = fitResult->Ndf();
        Double_t fittedreducedChi2 = (ndf > 0) ? chi2 / ndf : -1.0;
        reducedChi2[i] = fittedreducedChi2;

        Int_t i_low = h1->FindBin(sig_low);
        Int_t i_high = h1->FindBin(sig_high);
        Double_t sigplusbckg = h1->Integral(i_low, i_high);

        // --- Background & Signal Draw Functions ---
        TF1* fpol = nullptr;
        if (mcflag == 0) {
            fpol = new TF1("fpol", bkgOrder.c_str(), M_low, M_high);
            if (bkgOrder == "pol1") {
                fpol->SetParameters(pars[bkgStartIndex], pars[bkgStartIndex+1]);
            } else if (bkgOrder == "pol2") {
                fpol->SetParameters(pars[bkgStartIndex], pars[bkgStartIndex+1], pars[bkgStartIndex+2]);
            } else if (bkgOrder == "pol3") {
                fpol->SetParameters(pars[bkgStartIndex], pars[bkgStartIndex+1], pars[bkgStartIndex+2], pars[bkgStartIndex+3]);
            } else if (bkgOrder == "exp") {
                fpol->SetParameters(pars[bkgStartIndex], pars[bkgStartIndex+1]);
            }
        }

        TF1* signalDraw = new TF1("signalDraw", sigFormula.c_str(), M_low, M_high);
        if (fitFunction == "Voigt") {
            signalDraw->SetParameters(pars[0], pars[1], pars[2], pars[3]);
        } else if (fitFunction == "Gaussian") {
            signalDraw->SetParameters(pars[0], pars[1], pars[2]);
        } else if (fitFunction == "Double-Gaussian") {
            signalDraw->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
        } else if (fitFunction == "Breit-Wigner") {
            signalDraw->SetParameters(pars[0], pars[1], pars[2]);
        } else if (fitFunction == "Bukin") {
            signalDraw->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
        }

        TF1* totalfit = new TF1("TotalFit", fullFormula.c_str(), M_low, M_high);
        for (int p = 0; p < ff->GetNpar(); ++p) {
            totalfit->SetParameter(p, pars[p]);
        }
      
        // --- Canvas and Drawing ---
        TCanvas *fitcheck = new TCanvas(Form("fit_check_%d_mc%d", i, mcflag), Form("fit_check_%d_mc%d", i, mcflag), 600, 600);
        h1->GetYaxis()->SetTitle("Entries");
        h1->Draw();  
        fitcheck->Update();  
        gPad->SetTicky();

        TH1F* hfill = (TH1F*)h1->Clone("hfill");
        hfill->Reset();  

        int bin_low = h1->FindBin(sig_low);
        int bin_high = h1->FindBin(sig_high);
        double bin_low_edge = h1->GetBinLowEdge(bin_low);                           
        double bin_high_edge = h1->GetBinLowEdge(bin_high) + h1->GetBinWidth(bin_high);  

        TLine *llow = new TLine(bin_low_edge, 0, bin_low_edge, fitcheck->GetUymax());
        TLine *lhigh = new TLine(bin_high_edge, 0, bin_high_edge, fitcheck->GetUymax());
        llow->SetLineColor(kOrange-9); llow->SetLineStyle(2); llow->SetLineWidth(3); llow->Draw("same");
        lhigh->SetLineColor(kOrange-9); lhigh->SetLineStyle(2); lhigh->SetLineWidth(3); lhigh->Draw("same");

        for (int b = bin_low; b <= bin_high; ++b) {
            hfill->SetBinContent(b, h1->GetBinContent(b));
        }
        hfill->SetFillColor(kOrange-9);
        hfill->SetLineColor(kOrange-9);
        hfill->SetFillStyle(3002);  
        hfill->Draw("hist same");  

        if (mcflag == 0 && fpol) {
            fpol->SetLineColor(kGreen-3);
            fpol->Draw("same");
        }

        signalDraw->SetLineColor(kMagenta+1);
        signalDraw->SetNpx(1000);
        signalDraw->Draw("same");

        totalfit->SetLineColor(kRed+1);
        totalfit->SetNpx(1000);
        totalfit->Draw("same");

        TH1F* houtline = (TH1F*)hfill->Clone("houtline");
        houtline->SetFillStyle(0);               
        houtline->SetLineColor(kOrange-9);     
        houtline->SetLineWidth(2);               
        houtline->Draw("hist same");             

        h1->SetMarkerStyle(20);   
        h1->SetMarkerSize(1.0);   
        h1->SetMarkerColor(kAzure);
        h1->Draw("same PE");        

        TLatex *latexinv = new TLatex();
        latexinv->SetNDC();  
        latexinv->SetTextSize(0.045);  
        latexinv->DrawLatex(0.15, 0.80, latexTitle);
        latexinv->SetTextSize(0.03);  
        latexinv->DrawLatex(0.15, 0.72, "pp, #sqrt{#it{s}} = " + AO2DFiles.energy);
        TString latexLabel = TString::Format("%s #leq #it{p}_{T} < %s GeV/#it{c}", ptbins[i].c_str(), ptbins[i+1].c_str());
        latexinv->DrawLatex(0.15, 0.67, latexLabel);
        latexinv->DrawLatex(0.15, 0.62, Form("Centrality %.0f%% to %.0f%%", cent_low, cent_high));

        TLegend* legend = new TLegend(0.6, 0.65, 0.88, 0.88);  
        legend->AddEntry(h1, (mcflag == 0 ? "Data" : "MC Signal"), "lp");                 
        legend->AddEntry(hfill, (mcflag == 0 ? "Signal+Background" : "MC Signal Range"), "f");     
        if (mcflag == 0 && fpol) {
            legend->AddEntry(fpol, "Background Fit", "lp");     
        }
        legend->AddEntry(signalDraw, (fitFunction + " Fit").c_str(), "lp");         
        legend->AddEntry(totalfit, (mcflag == 0 ? "Total Fit" : "Signal Fit"), "lp");      
        legend->SetBorderSize(0);  
        legend->SetFillStyle(0);   
        legend->Draw("same");

        TDirectory* prevdir = gDirectory; 
        centdir->cd(); 
        fitcheck->Write(fmt::format("{}_{}_{}", histname, ptbins[i], ptbins[i+1]).data());
        prevdir->cd(); 

        Double_t bckg = (mcflag == 0 && fpol) ? fpol->Integral(sig_low, sig_high) : 0.0;
        Double_t realbckg = bckg / h1->GetBinWidth(0);
        if (mcflag == 1) {
            realbckg = 0.;
        }
        signal[i] = sigplusbckg - realbckg;
        peakpt[i] = pars[1];
        widthpt[i] = widthVal;
        cout<< "Signal plus background: " << sigplusbckg << " | Signal: " << signal[i] << " | Background: " << realbckg << endl;

        // Error Propagation (Stat errors)
        Double_t sigplusbckgError2 = 0.0;
        for (Int_t b = i_low; b <= i_high; ++b) {
            Double_t err = h1->GetBinError(b);
            sigplusbckgError2 += err * err;
        }
        Double_t sigplusbckgError = TMath::Sqrt(sigplusbckgError2);

        Double_t signalError = sigplusbckgError; 

        if (mcflag == 0) {
            hStatError->SetBinContent(1, i + 1, signalError);
            hRelStatError->SetBinContent(1, i + 1, 100.0 * signalError / signal[i]);
            hStatError->GetXaxis()->SetBinLabel(1, "Inv. mass");
            hRelStatError->GetXaxis()->SetBinLabel(1, "Inv. mass");
        }
    }
    return;
}
//
//PlottingFunction
//
void PlotSignal(const std::string& particle, const std::vector<std::string>& ptbins, TFile* file, AO2DConfig& AO2DFiles, TDirectory* centdir, float cent_low, float cent_high, int centralityIndex, bool Systematics_Analysis) {
    //
    //
// Initializing measured and reconstructed signals to be calculated from inv mass pltos
   //
    int nhistos = ptbins.size()-1;  //warning for less than a few thousand
    std::vector<Double_t> signal_measured(nhistos);
    std::vector<Double_t> signal_reconstructed(nhistos);
    std::vector<Double_t> signal_split_reconstructed(nhistos);
    std::vector<Double_t> peak_measured(nhistos);
    std::vector<Double_t> peak_reconstructed(nhistos);
    std::vector<Double_t> width_measured(nhistos);
    std::vector<Double_t> width_reconstructed(nhistos);
    std::vector<Double_t> mc_sig_low(nhistos);
    std::vector<Double_t> mc_sig_high(nhistos);
    std::vector<Double_t> reducedChi2_data(nhistos);
    std::vector<Double_t> reducedChi2_mc(nhistos);
    std::vector<Double_t> reducedChi2_mc_split(nhistos);


    //********************************************************************************/
    //Systematics to be processed (Default is Double-Gaussian + pol2 + 1.0*Range
   //********************************************************************************/
    std::string systInvMassFitFunction = "Double-Gaussian";
    std::string systInvMassBackgroundOrder = "pol2";
    double systInvMassRangeScale = 1.0;

    // 1. Set default baseline configuration based on the particle type
    if (particle == "K0sh") {
        systInvMassFitFunction = "Double-Gaussian";
        systInvMassBackgroundOrder = "pol2";
    } else if (particle == "Lambda") {
        systInvMassFitFunction = "Double-Gaussian";
        systInvMassBackgroundOrder = "pol2";
    } else if (particle == "AntiLambda") {
        systInvMassFitFunction = "Double-Gaussian";
        systInvMassBackgroundOrder = "pol2";
    }

    // 2. Override or select specific variations if Systematics Analysis is active
    if (Systematics_Analysis == true) {
        // The same for all particles
        systInvMassFitFunction = "Double-Gaussian"; // e.g., "Voigt", "Gaussian", "Double-Gaussian", "Breit-Wigner", "Bukin Not Yet Implemented"
        systInvMassBackgroundOrder = "pol2"; // e.g., "pol1", "pol2"
        // Invariant-mass fit range scale (un-shadowed variable fix)
        systInvMassRangeScale = 1.1; // Change to 0.8, 0.9, 1.0, 1.1, 1.2 as needed
    }

    cout << "Particle: " << particle << std::endl;
    cout << "Systematics Analysis: " << (Systematics_Analysis ? "Enabled" : "Disabled") << std::endl;
    cout << "Fit Function: " << systInvMassFitFunction << std::endl;
    cout << "Background Order: " << systInvMassBackgroundOrder << std::endl;
    cout << "Range Scale: " << systInvMassRangeScale << std::endl;
    cout << "Systematics Analysis: " << (Systematics_Analysis ? "Enabled" : "Disabled") << std::endl;
    cout << "Fit Function " << systInvMassFitFunction << std::endl;
    cout << "Background Order " << systInvMassBackgroundOrder << std::endl;

    int mcflag=0;
   //make strings of pt values into double for use
    std::vector<Double_t> ptbinedges(nhistos+1);
    for (int i = 0; i < nhistos+1; i++) {
        ptbinedges[i] = std::stod(ptbins[i]);
    }
   //end initialisation
   //
    //Histograms to be filled with the error evolution
    //
    TH2D* hStatError = new TH2D("hStatError","Evolution of Statistical Uncertainty;Analysis stage;p_{T} (GeV/c)",10, 0.5, 10.5,nhistos, ptbinedges.data());
    TH2D* hRelStatError = new TH2D("hRelStatError","Evolution of Relative Statistical Uncertainty;Analysis stage;p_{T} (GeV/c)",10, 0.5, 10.5,nhistos, ptbinedges.data());
 //
//Access Reconstructed DATA to get Measured Pt Inv. Mass Plots and Ncollisions
   //
   //
    TFile *dataroot = new TFile(AO2DFiles.datarootpath.c_str());
   // Determine particle-specific directory and histogram name base
    std::string particleName;        // base name used for histograms (e.g. K0sh, Lambda, AntiLambda)
    if (particle == "K0sh") {
        particleName = "K0sh";
    } else if (particle == "Lambda") {
        particleName = "Lambda";
    } else { // lambdabar
        particleName = "AntiLambda";
    }
   // Count histogram name derived from particleName
    std::string nCountHistName = std::string("hN") + particleName;
    TString datarootSubDir = AO2DFiles.datarootsubfolder.c_str();

    TString datamassSubDir = Form("%s/%sMassPlotsPerPtBin",datarootSubDir.Data(),particleName.c_str());
    dataroot->cd(datamassSubDir);
    mcflag=0;
    std::string histname;
    CalculateSignal(particleName, "Pt_vs_Cent_from", ptbins, signal_measured, file, centdir, mcflag, mc_sig_low, mc_sig_high, peak_measured, width_measured, nhistos, cent_low, cent_high, hStatError, hRelStatError, reducedChi2_data, systInvMassFitFunction, systInvMassBackgroundOrder, systInvMassRangeScale);

    TH1* hMeasured = new TH1D(Form("%sMeasuredPtSpectrum", particleName.c_str()), Form("%sMeasuredPtSpectrum", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data());
    TH1* hPeak = new TH1D(Form("%sPeakvsPtSpectrum", particleName.c_str()), Form("%sPeakvsPtSpectrum", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data());
    TH1* hWidth = new TH1D(Form("%sWidthvsPtSpectrum", particleName.c_str()), Form("%sWidthvsPtSpectrum", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data());
    TH1* hReducedChi2 = new TH1D(Form("%sReducedChi2vsPtSpectrum", particleName.c_str()), Form("%sReducedChi2vsPtSpectrum", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data());
    for(int i=0;i<nhistos;i++){
        hMeasured->SetBinContent(i+1,signal_measured[i]);
        hPeak->SetBinContent(i+1,peak_measured[i]);
        hWidth->SetBinContent(i+1,width_measured[i]);
        hReducedChi2->SetBinContent(i+1,reducedChi2_data[i]);
    }
    TString ptAnalysisSubDir = Form("%s/PtAnalysis", datarootSubDir.Data());
    dataroot->cd(ptAnalysisSubDir);

   //Getting Number of Events and Particles after each applied cut
    TH1* hNRecEventsData = ReadAndSaveHistogramProjection("hNRecEvents", "hNRecEventsData", cent_low, cent_high);
    TH1* hNEventsData = ReadAndSaveHistogramProjection("hNEvents", "hNEventsData", cent_low, cent_high);
    TH1* hNV0s = ReadAndSaveHistogramProjection("hNV0s", "hNV0sData", cent_low, cent_high);
    TH1* hParticleCount = ReadAndSaveHistogramProjection(nCountHistName.c_str(), (nCountHistName+"Data").c_str(), cent_low, cent_high);

   //Getting Number of Collisions for Normalisation
    double nCollisionsData = hNRecEventsData->GetBinContent(1);
   //
    hNEventsData->Scale(1.0 / hNEventsData->GetMaximum());
    hNV0s->Scale(1.0 / hNV0s->GetMaximum());
    hParticleCount->Scale(1.0 / hParticleCount->GetMaximum());

    hNEventsData->SetMinimum(0);
    hNV0s->SetMinimum(0);
    hNV0s->GetXaxis()->SetRange(1, 4);
    hParticleCount->SetMinimum(0);
   //
   /*Getting Generated Pt Spectrum from "Data-input" MC for Beta Test
   //
    TFile *datarootBetaGen= new TFile("PtSpectrumResults/Data/AnalysisResultsBetaGen.root"); //DATA Root
   //TFile *dataroot= new TFile("PtSpectrumResults/Data/AnalysisResultsReal.root"); //DATA Root Local
    datarootBetaGen->cd("v0-pt-inv-mass-plots/MCCorrections");//Generated subfolder with MC Corrections  from Local Run
   // Generated Pt Spectrum for this particle
    TH1* hGeneratedAfterEventSelectionB = ReadAndSaveHistogramProjection(Form("h%sGeneratedPtSpectrum", particleName.c_str()), Form("%sGeneratedPtSpectrumB", particleName.c_str()), cent_low, cent_high);
   //
    TH1* hNEventsCorrectionsGenBeta = ReadAndSaveHistogramProjection("hNEvents_Corrections", "hNEvents_CorrectionsB", cent_low, cent_high);
    double nCollisionsGenBeta = hNEventsCorrectionsGenBeta->GetBinContent(1);
   */
   //
//Access Generated MC to get Gen Pt Spectrums
   //
   //
//Go to MC Gen Folder to get MC Corrections 
    TFile *mcrootgen = new TFile(AO2DFiles.mcgenrootpath.c_str());
    TString mcgenSubDir = AO2DFiles.mcgenrootsubfolder.c_str();
    TString mcgencorrectionsSubDir = Form("%s/MCCorrections", mcgenSubDir.Data());
    mcrootgen->cd(mcgencorrectionsSubDir);
   // K0sh Generated Pt Spectrum
   // --- Generated spectrum ---
   
    std::string genName = Form("h%sGeneratedPtSpectrum", particle.c_str());
    TH1* hGenerated = ReadAndSaveHistogramProjection(genName.c_str(),Form("%sGeneratedPtSpectrumSaved", particle.c_str()),cent_low, cent_high);

   // --- After event selection ---
    std::string afterSelName = Form("h%sGeneratedRecoPtSpectrum", particle.c_str());
    TH1* hGeneratedAfterEventSelection = ReadAndSaveHistogramProjection(afterSelName.c_str(),Form("%sGeneratedRecoPtSpectrumSaved", particle.c_str()),cent_low, cent_high);

    if (!hGenerated || !hGeneratedAfterEventSelection) {
        std::cerr << "ERROR: Missing generated or after-event-selection histogram for " << particle << std::endl;
        return;
    }

   // --- Event corrections ---
   //Getting Histo for Event Efficiency, bin 0.5 is Event Efficiency Denominator and bin 2.5 is Numerator (also Event Split Denomimator)
   //No it's 1.5 and 3.5
    TH1* hNEventsCorrections = ReadAndSaveHistogramProjection("hNEvents_Corrections","hNEvents_CorrectionsSaved",cent_low, cent_high);
    double nCollisionsGen = hNEventsCorrections->GetBinContent(2); // Because we apply cut-Z vertex to both

    TH1* hGeneratedXiMinus = nullptr;
    TH1* hGeneratedXiZero = nullptr;
    TH1* hGeneratedXiPlus = nullptr;
    TH1* hGeneratedAntiXiZero = nullptr;

   // --- Xi feed-down (only for Lambda / AntiLambda) ---
    if (particle == "Lambda") {
        hGeneratedXiMinus = ReadAndSaveHistogramProjection("hXiMinusGeneratedPtSpectrum", "hXiMinusGeneratedPtSpectrumSaved",cent_low, cent_high);
        hGeneratedXiZero = ReadAndSaveHistogramProjection("hXiZeroGeneratedPtSpectrum", "hXiZeroGeneratedPtSpectrumSaved",cent_low, cent_high);
    }
    else if (particle == "AntiLambda") {
        hGeneratedXiPlus = ReadAndSaveHistogramProjection("hXiPlusGeneratedPtSpectrum", "hXiPlusGeneratedPtSpectrumSaved",cent_low, cent_high);
        hGeneratedAntiXiZero = ReadAndSaveHistogramProjection("hAntiXiZeroGeneratedPtSpectrum", "hAntiXiZeroGeneratedPtSpectrumSaved",cent_low, cent_high);
    }
   //
   //Getting Histos for Particle Signal Loss Numerator and Denominator
    TH1D* hSignalLossNumerator = (TH1D*)hGeneratedAfterEventSelection->Clone(Form("h%sGeneratedRecoPtSpectrum", particleName.c_str()));
    TH1* hSignalLossDenominator = ReadAndSaveHistogramProjection(Form("h%sGeneratedPtSpectrum", particleName.c_str()), Form("%sSignalLossDenominatorSaved", particleName.c_str()), cent_low, cent_high);
   //
   //
//Access Reconstructed MC to get MC Reconstructed Pt for Efficiency and Feeddown Matrix of Xi
   //
   // 
   //Getting the Reconstructed MC Information: Pt Invariant Mass Plots for Efficiency and Event/Particle Corrections
    TFile *mcroot = new TFile(AO2DFiles.mcrootpath.c_str());
    TString mcSubDir = AO2DFiles.mcrootsubfolder.c_str();
    TString mcMassSubDir = Form("%s/%sMassPlotsPerPtBin",mcSubDir.Data(), particleName.c_str());
    mcroot->cd(mcMassSubDir);
    mcflag=1;
    CalculateSignal(particleName, "Pt_vs_Cent_from", ptbins, signal_reconstructed, file, centdir,mcflag, mc_sig_low, mc_sig_high, peak_reconstructed, width_reconstructed, nhistos, cent_low, cent_high, hStatError, hRelStatError, reducedChi2_mc, systInvMassFitFunction, systInvMassBackgroundOrder,systInvMassRangeScale);

    TH1* hReconstructedMC = new TH1D(Form("%sReconstructedPtSpectrum", particleName.c_str()), Form("%sReconstructedPtSpectrum", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data());
    TH1* hReducedChi2MC = new TH1D(Form("%sReducedChi2vsPtSpectrumMC", particleName.c_str()), Form("%sReducedChi2vsPtSpectrumMC", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data());
    for(int i=0;i<nhistos;i++){
        hReconstructedMC->SetBinContent(i+1,signal_reconstructed[i]); //Filling hReconstructed
        hReducedChi2MC->SetBinContent(i+1,reducedChi2_mc[i]);
    }
   //Getting Histos and Calculating Particle Signal Split
    TString mcMassSplitSubDir = Form("%s/%sSplitMassPlotsPerPtBin",mcSubDir.Data(), particleName.c_str());
    mcroot->cd(mcMassSplitSubDir);
    CalculateSignal(particleName, "SplitPt_vs_Cent_from", ptbins, signal_split_reconstructed, file, centdir, mcflag, mc_sig_low, mc_sig_high, peak_reconstructed, width_reconstructed, nhistos, cent_low, cent_high, hStatError, hRelStatError, reducedChi2_mc_split, systInvMassFitFunction, systInvMassBackgroundOrder,systInvMassRangeScale);
    histname = Form("h%sPt_vs_Cent_from", particleName.c_str());

    TH1* hParticleSignalSplitNumerator = new TH1D(Form("h%sSignalSplitNumerator", particleName.c_str()), Form("h%sSignalSplitNumerator", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data());
    for(int i=0;i<nhistos;i++){
        hParticleSignalSplitNumerator->SetBinContent(i+1,signal_reconstructed[i]); //Filling hSignalSplitNumerator
    }

    TString mcFeeddownSubDir = Form("%s/FeeddownMatrices",mcSubDir.Data());
    mcroot->cd(mcFeeddownSubDir);
   // Feeddown matrices: choose Lambda or AntiLambda matrices depending on particle
    TH2D* hFeeddownXiMinusMatrix = nullptr;
    TH2D* hFeeddownXiZeroMatrix = nullptr;
    TH2D* hFeeddownXiPlusMatrix = nullptr;
    TH2D* hFeeddownAntiXiZeroMatrix = nullptr;
    if (particle == "Lambda") {
      // For Lambda, project TH3 feeddown matrices over centrality to get TH2
        std::cout << "[DEBUG 1] Matrix Before Projection Read = " << hFeeddownXiMinusMatrix << std::endl;
        hFeeddownXiMinusMatrix = dynamic_cast<TH2D*>(ReadAndSaveHistogramProjection("hLambdaXiMinusFeeddownMatrix","hLambdaXiMinusFeeddownMatrixSaved", cent_low, cent_high));
        std::cout << "[DEBUG 2] Matrix After Projection Read = " << hFeeddownXiMinusMatrix << std::endl;
        hFeeddownXiZeroMatrix  = dynamic_cast<TH2D*>(ReadAndSaveHistogramProjection("hLambdaXiZeroFeeddownMatrix","hLambdaXiZeroFeeddownMatrixSaved", cent_low, cent_high));
    } else if (particle == "AntiLambda") {
      // For AntiLambda, project TH3 feeddown matrices over centrality to get TH2
        hFeeddownXiPlusMatrix = dynamic_cast<TH2D*>(ReadAndSaveHistogramProjection("hAntiLambdaXiPlusFeeddownMatrix","hAntiLambdaXiPlusFeeddownMatrixSaved", cent_low, cent_high));
        hFeeddownAntiXiZeroMatrix = dynamic_cast<TH2D*>(ReadAndSaveHistogramProjection("hAntiLambdaAntiXiZeroFeeddownMatrix","hAntiLambdaAntiXiZeroFeeddownMatrixSaved", cent_low, cent_high));
    }
//
//
    TString mcPtAnalysisSubDir = Form("%s/PtAnalysis",mcSubDir.Data());
    mcroot->cd(mcPtAnalysisSubDir);
    TH1* hNEventsMC = ReadAndSaveHistogramProjection("hNEvents", "hNEventsMC", cent_low, cent_high);
    TH1* hNRecEventsMC = ReadAndSaveHistogramProjection("hNRecEvents", "hNRecEventsMC",cent_low, cent_high);
    TH1* hNV0sMC = ReadAndSaveHistogramProjection("hNV0s", "hNV0sMC", cent_low, cent_high);
    TH1* hParticleCountMC = ReadAndSaveHistogramProjection(nCountHistName.c_str(), Form("%sMC", nCountHistName.c_str()), cent_low, cent_high);
    hNEventsMC->SetMinimum(0);
    hNV0sMC->SetMinimum(0);
    hParticleCountMC->SetMinimum(0);
    double nCollisionsMC = hNEventsMC->GetBinContent(1);
//
//Access XiMinus+Xiplus Corrected Pt Spectrum downloaded from HEPData for Feeddown. (Run 2 Xi: https://www.hepdata.net/record/ins1748157?version=1 )
// Only perform Run2 Xi handling for Lambda or AntiLambda (not for K0sh)


    TH1D* hXiCorrectedSpectrumInput = nullptr;
    std::vector<double> XiscaleFactors = {12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0};

    if (particle == "Lambda" || particle == "AntiLambda") {

        hXiCorrectedSpectrumInput = new TH1D("CorrectedXiPtSpectrumInput", "CorrectedXiPtSpectrumInput", ptbinedges.size()-1, ptbinedges.data());

        TString xiFileName;
        TString xiDirName;
        TString xiHistName;
        TH1F *xiHisto = nullptr;

        if (AO2DFiles.dataperiod == "pp24" || AO2DFiles.dataperiod == "pp22o") {
            xiFileName = "XiCorrectedRun3RomanRUSRenamed.root";
            xiDirName = "";
            xiHistName = Form("Yield_XiPm_m%d", centralityIndex + 1);
        }
        else if (AO2DFiles.dataperiod == "ppref24") { //for ppref
            xiFileName = "HEPData502Paper_ModifiedXi.root";  
            xiDirName = "Figure 1c";
            xiHistName = Form("Hist1D_y%d", centralityIndex + 1);
        }
        else {
            std::cout << "ERROR: Unsupported data period: "
                    << AO2DFiles.dataperiod << std::endl;
            return;
        }

        TFile *xiFile = TFile::Open(xiFileName);
        if (!xiFile || xiFile->IsZombie()) {
            std::cout << "ERROR: cannot open file " << xiFileName << std::endl;
            return;
        }

        if (xiDirName.IsNull()) {
            // Histogram is in the top level of the ROOT file
            xiHisto = dynamic_cast<TH1F*>(xiFile->Get(xiHistName));
            // Only apply the scaling if needed
            double Xiscale = std::pow(2.0, XiscaleFactors[centralityIndex]);
            xiHisto->Scale(1.0 / Xiscale);
        }
        else {
        // Histogram is inside a directory
            TDirectory *dirXiInput = dynamic_cast<TDirectory*>(xiFile->Get(xiDirName));
            if (!dirXiInput) {
                std::cout << "ERROR: Directory " << xiDirName << " not found!" << std::endl;
                return;
            }
            xiHisto = dynamic_cast<TH1F*>(dirXiInput->Get(xiHistName));
        }

        std::cout << "Using TH1F '" << xiHisto->GetName() << "' for feeddown." << std::endl;
        cout<<"Xi Histogram fitted for Centrality Index: "<<centralityIndex<<endl;
        //
        double fitMin, fitMax;
        double mXi = 1.32171; // Xi mass

        if (AO2DFiles.dataperiod == "ppref24") {
            fitMin = 0.6;
            fitMax = 6.5;
        } else {
            fitMin = 0.6;
            fitMax = 6.5;
        }

        TF1 *levyTsallisXi = new TF1("levyTsallisXi","x*[0]*((([1]-1.)*([1]-2.))/(([1]*[2])*(([1]*[2])+([3]*([1]-2.)))))*pow(1.+((sqrt(pow([3],2.)+pow(x,2.))-[3])/([1]*[2])),-[1])",fitMin,fitMax);
        levyTsallisXi->FixParameter(3, mXi);

        // Different initial parameters depending on the dataset
        if (AO2DFiles.dataperiod == "ppref24") {
            // pp reference (HEPData Figure 1c)
            levyTsallisXi->SetParameters(5e-3, 8.0, 0.40);
            levyTsallisXi->SetParLimits(0, 1e-6, 1.0);   // Normalisation
            levyTsallisXi->SetParLimits(1, 4.0, 20.0);   // n
            levyTsallisXi->SetParLimits(2, 0.10, 1.20);  // C (GeV/c)
        }
        else {
            // Run 3 Xi pp from Roman
            levyTsallisXi->SetParameters(1e-2, 10.0, 0.30);
            levyTsallisXi->SetParLimits(0, 1e-6, 1.0);
            levyTsallisXi->SetParLimits(1, 4.0, 20.0);
            levyTsallisXi->SetParLimits(2, 0.10, 1.20);
        }

        xiHisto->Fit(levyTsallisXi, "R");

        for (size_t bin = 0; bin < ptbinedges.size()-1; ++bin) {  
            double binLow = ptbinedges[bin];
            double binHigh = ptbinedges[bin+1];
            double integral = levyTsallisXi->Integral(binLow, binHigh); // Integrate Levy-Tsallis function over bin range
            double  Xierror = sqrt(integral);
            hXiCorrectedSpectrumInput->SetBinContent(bin+1, integral);  // ROOT bins start at 1
            hXiCorrectedSpectrumInput->SetBinError(bin+1, 0.0);  // Set error to 0 for now
         //hXiCorrectedSpectrumInput->SetBinError(bin+1, Xierror);
        }
    }
//
//end of accessing and saving DATA
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////                Start of Analysis                /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////
//MC Corrections
   //////////////
    TH1* hStatErrorEvolution =(TH1*)hMeasured->Clone("hStatErrorEvolution");

   //Event Efficiency
   //Reminder: Bin 0.5 is Event Efficiency Denominator and Bin 2.5 is Numenator NO it's 
   //Reminder: Bin 0.5 NO its 1.5! is Event Efficiency Denominator and Bin 2.5 is Numeneator NO its 3.5!
   //Reminder: New changes it's Bin 1.5 and Bin 5.5!!
    if(hNEventsCorrections->GetBinContent(hNEventsCorrections->FindBin(1.5))==0){
        cout<<"ERROR: Event Efficiency Denominator is 0!"<<endl;
        return;
    }
    double EventEfficiency = hNEventsCorrections->GetBinContent(hNEventsCorrections->FindBin(5.5)) / hNEventsCorrections->GetBinContent(hNEventsCorrections->FindBin(1.5));
    double EventEfficiencyError = std::sqrt(EventEfficiency * (1.0 - EventEfficiency) / hNEventsCorrections->GetBinContent(hNEventsCorrections->FindBin(1.5)));
    TH1* hParticleEventEfficiency = new  TH1D(Form("%s Event Reconstruction Efficiency", particleName.c_str()), Form("%s Event Reconstruction Efficiency", particleName.c_str()), ptbinedges.size()-1,ptbinedges.data()); //Save it for each particle to calculate spectrum but it's the same
    for (int bin = 1; bin <= hParticleEventEfficiency->GetNbinsX(); ++bin) {
        hParticleEventEfficiency->SetBinContent(bin, EventEfficiency);
        hParticleEventEfficiency->SetBinError(bin, EventEfficiencyError);
    }
    hParticleEventEfficiency->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hParticleEventEfficiency->GetYaxis()->SetTitle("Event Efficiency");
    hParticleEventEfficiency->SetMarkerStyle(20);  // full circle
    hParticleEventEfficiency->SetMarkerSize(1.2);


    hStatError->GetXaxis()->SetBinLabel(6, "Event efficiency");
    hRelStatError->GetXaxis()->SetBinLabel(6, "Event efficiency");

    //Event Spliting
    const int eventSplittingBin = hNEventsCorrections->FindBin(5.5);  //<-------- CHECK THIS
    const double nRecEventsMC = hNRecEventsMC->GetBinContent(1);
    const double nEventsCorrections = hNEventsCorrections->GetBinContent(eventSplittingBin);
    const double nEventsCorrectionsErr = hNEventsCorrections->GetBinError(eventSplittingBin);

    if (nEventsCorrections == 0.0) {
        cout << "ERROR: Event Splitting Denominator is 0!" << endl;
        return;
    }

    const double EventSplitting = nRecEventsMC / nEventsCorrections;
    const double nRecEventsMCErr = (nRecEventsMC > 0.0) ? std::sqrt(nRecEventsMC) : 0.0;
    const double EventSplittingError = (nRecEventsMC > 0.0)
        ? EventSplitting * std::sqrt(std::pow(nRecEventsMCErr / nRecEventsMC, 2) + std::pow(nEventsCorrectionsErr / nEventsCorrections, 2))
        : 0.0;

    TH1* hParticleEventSplitting = new TH1D(Form("%s Event Splitting", particleName.c_str()), Form("%s Event Splitting", particleName.c_str()), ptbinedges.size() - 1, ptbinedges.data());
    hParticleEventSplitting->Sumw2();
    for (int bin = 1; bin <= hParticleEventSplitting->GetNbinsX(); ++bin) {
        hParticleEventSplitting->SetBinContent(bin, EventSplitting);
        hParticleEventSplitting->SetBinError(bin, EventSplittingError);
    }
    hParticleEventSplitting->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hParticleEventSplitting->GetYaxis()->SetTitle("Event Splitting");
    hParticleEventSplitting->SetMarkerStyle(20);  // full circle
    hParticleEventSplitting->SetMarkerSize(1.2);

    hStatError->GetXaxis()->SetBinLabel(7, "Event splitting");
    hRelStatError->GetXaxis()->SetBinLabel(7, "Event splitting");

   //Event Reconstruction Particle Signal Loss
    TH1* hParticleSignalLoss = (TH1*)hSignalLossNumerator->Clone(Form("%sSignalLoss", particleName.c_str()));
    hParticleSignalLoss->Divide(hSignalLossNumerator, hSignalLossDenominator, 1.0, 1.0,"B");
    hParticleSignalLoss->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hParticleSignalLoss->GetYaxis()->SetTitle("Signal Loss");

    hStatErrorEvolution->Divide(hParticleSignalLoss);
    RecordStatError(hStatError,hRelStatError,4,hStatErrorEvolution);

    hStatError->GetXaxis()->SetBinLabel(4, "Signal loss");
    hRelStatError->GetXaxis()->SetBinLabel(4, "Signal loss");

   //Xi Feeddown Correction
//Warning: Our Xi Spectrum is XiMinus+XiPlus; we divide by 2 to account for XiMinus and XiZero symmetry
//We have taken into account that XiMinus(Plus) and (Anti)XiZero have different feeddown matrices
//Warning: We assume the Omega contribution is negligible
   //
    // Pointers to the appropriate feeddown matrices based on particle type
    TH1D* hFromXiSignSpectrum = nullptr;
    TH1D* hFromXiNeutralSpectrum = nullptr;
    TH1D* hTotalFeeddown = nullptr;
    TH1D* hTotalFeeddownFraction = nullptr;
    TH2D* hFeeddownXiSignMatrix = nullptr;
    TH2D* hFeeddownXiNeutralMatrix = nullptr;
    TH1* hGenXiSign = nullptr;
    TH1* hGenXiNeutral = nullptr;
    std::string xiSignName, xiNeutralName;

   //Case for K0sh: we assume any feeddown is negligible and set feeddown histograms to zero
    if (particle == "K0sh") {
        hTotalFeeddown = new TH1D(Form("%sFromXiFeeddown", particleName.c_str()), Form("%sFromXiFeeddown", particleName.c_str()), ptbinedges.size() - 1, ptbinedges.data());
        for (int bin = 1; bin <= hTotalFeeddown->GetNbinsX(); ++bin){
            hTotalFeeddown->SetBinContent(bin, 0);
            hTotalFeeddown->SetBinError(bin, 0);
        }
    }

    if (particle == "Lambda" || particle == "AntiLambda") {
        int nRecoPtBins = hXiCorrectedSpectrumInput->GetNbinsX();
      // Feeddown matrices and accumulators
      // Xi charged (XiMinus for Lambda, XiPlus for AntiLambda)
        std::vector<std::vector<double>> feeddownXiSignMatrix(nRecoPtBins, std::vector<double>(nRecoPtBins, 0.0));
        std::vector<double> totalFeeddownXiSign(nRecoPtBins, 0.0);
      // Xi neutral (XiZero for Lambda, AntiXiZero for AntiLambda)
        std::vector<std::vector<double>> feeddownXiNeutralMatrix(nRecoPtBins, std::vector<double>(nRecoPtBins, 0.0));
        std::vector<double> totalFeeddownXiNeutral(nRecoPtBins, 0.0);
        if (particle == "Lambda") {
            std::cout << "[DEBUG END] Matrix As Read = " << hFeeddownXiMinusMatrix << std::endl;
            hFeeddownXiSignMatrix = hFeeddownXiMinusMatrix;
            hFeeddownXiNeutralMatrix = hFeeddownXiZeroMatrix;
            hGenXiSign = dynamic_cast<TH1D*>(hGeneratedXiMinus);
            hGenXiNeutral = dynamic_cast<TH1D*>(hGeneratedXiZero);
            xiSignName = "XiMinus";
            xiNeutralName = "XiZero";
        } else { // lambdabar
            hFeeddownXiSignMatrix = hFeeddownXiPlusMatrix;
            hFeeddownXiNeutralMatrix = hFeeddownAntiXiZeroMatrix;
            hGenXiSign = dynamic_cast<TH1D*>(hGeneratedXiPlus);
            hGenXiNeutral = dynamic_cast<TH1D*>(hGeneratedAntiXiZero);
            xiSignName = "XiPlus";
            xiNeutralName = "AntiXiZero";
        }
      // Calculate feeddown for both Xi types
        for (int lambdaRecoPtBin = 1; lambdaRecoPtBin <= nRecoPtBins; ++lambdaRecoPtBin) {
            for (int xiGenPtBin = 1; xiGenPtBin <= nRecoPtBins; ++xiGenPtBin) {
            // Get normalized Xi generation content
                double genXiSignContent = hGenXiSign->GetBinContent(xiGenPtBin);
                double genXiNeutralContent = hGenXiNeutral->GetBinContent(xiGenPtBin);
                double xiSpectrumContent = hXiCorrectedSpectrumInput->GetBinContent(xiGenPtBin) / 2.0; // This is XiMinus+XiPlus spectrum, divide by 2 for individual Xi contribution
            
            // Avoid division by zero
                if (genXiSignContent < 1e-10) genXiSignContent = 1.0;
                if (genXiNeutralContent < 1e-10) genXiNeutralContent = 1.0;
            
            // Normalize feeddown matrices by generated Xi spectra
                if (hFeeddownXiSignMatrix) {
                    feeddownXiSignMatrix[lambdaRecoPtBin-1][xiGenPtBin-1] = hFeeddownXiSignMatrix->GetBinContent(lambdaRecoPtBin, xiGenPtBin) / genXiSignContent;
                }
                if (hFeeddownXiNeutralMatrix) {
                    feeddownXiNeutralMatrix[lambdaRecoPtBin-1][xiGenPtBin-1] = hFeeddownXiNeutralMatrix->GetBinContent(lambdaRecoPtBin, xiGenPtBin) / genXiNeutralContent;
                }
            
            // Accumulate total feeddown for this Lambda/AntiLambda pT
                totalFeeddownXiSign[lambdaRecoPtBin-1] += feeddownXiSignMatrix[lambdaRecoPtBin-1][xiGenPtBin-1] * xiSpectrumContent * nCollisionsData;
                totalFeeddownXiNeutral[lambdaRecoPtBin-1] += feeddownXiNeutralMatrix[lambdaRecoPtBin-1][xiGenPtBin-1] * xiSpectrumContent * nCollisionsData; 
            }
        }

      // Define output histograms with particle-specific names
        hFromXiSignSpectrum = new TH1D(Form("%sFrom%sFeeddown", particleName.c_str(), xiSignName.c_str()), Form("%sFrom%sFeeddown", particleName.c_str(), xiSignName.c_str()), ptbinedges.size() - 1, ptbinedges.data());
        hFromXiNeutralSpectrum = new TH1D(Form("%sFrom%sFeeddown", particleName.c_str(), xiNeutralName.c_str()), Form("%sFrom%sFeeddown", particleName.c_str(), xiNeutralName.c_str()), ptbinedges.size() - 1, ptbinedges.data());
        hTotalFeeddown = new TH1D(Form("%sTotalFeeddown", particleName.c_str()), Form("%sTotalFeeddown", particleName.c_str()), ptbinedges.size() - 1, ptbinedges.data());
        hTotalFeeddownFraction = new TH1D(Form("%sTotalFeeddownFraction", particleName.c_str()), Form("%sTotalFeeddownFraction", particleName.c_str()), ptbinedges.size() - 1, ptbinedges.data());
      
      // Fill the spectra with error
        for (int lambdaRecoPtBin = 1; lambdaRecoPtBin <= nRecoPtBins; ++lambdaRecoPtBin) {
            double feeddownSign = totalFeeddownXiSign[lambdaRecoPtBin - 1];
            double feeddownNeutral = totalFeeddownXiNeutral[lambdaRecoPtBin - 1];
            double totalFeeddown = feeddownSign + feeddownNeutral;
         
            double errSign = std::sqrt(feeddownSign);
            double errNeutral = std::sqrt(feeddownNeutral);
         
            double measuredSignal = signal_measured[lambdaRecoPtBin-1];
            double errMeasuredSignal = std::sqrt(measuredSignal);
         
         // Set individual feeddown histograms
            hFromXiSignSpectrum->SetBinContent(lambdaRecoPtBin, feeddownSign);
            hFromXiSignSpectrum->SetBinError(lambdaRecoPtBin, errSign);
         
            hFromXiNeutralSpectrum->SetBinContent(lambdaRecoPtBin, feeddownNeutral);
            hFromXiNeutralSpectrum->SetBinError(lambdaRecoPtBin, errNeutral);
         
         // Set total feeddown histogram
            hTotalFeeddown->SetBinContent(lambdaRecoPtBin, totalFeeddown);
            hTotalFeeddown->SetBinError(lambdaRecoPtBin, std::sqrt(errSign * errSign + errNeutral * errNeutral));
         
         // Compute feeddown fraction with proper error propagation
            const double EPS = 1e-10;
            if (measuredSignal > EPS && totalFeeddown > EPS) {
                double feeddownFraction = totalFeeddown / measuredSignal;
                double errFeeddownFraction = feeddownFraction * std::sqrt((errSign * errSign + errNeutral * errNeutral) / (totalFeeddown * totalFeeddown) + (errMeasuredSignal * errMeasuredSignal) / (measuredSignal * measuredSignal));
                hTotalFeeddownFraction->SetBinContent(lambdaRecoPtBin, feeddownFraction);
                hTotalFeeddownFraction->SetBinError(lambdaRecoPtBin, errFeeddownFraction);
            } else {
                hTotalFeeddownFraction->SetBinContent(lambdaRecoPtBin, 0);
                hTotalFeeddownFraction->SetBinError(lambdaRecoPtBin, 0);
            }
        }
    }

    hStatError->GetXaxis()->SetBinLabel(2, "Feeddown");
    hRelStatError->GetXaxis()->SetBinLabel(2, "Feeddown");

   //Particle Signal Splitting
    TH1* hParticleSignalSplitting = (TH1*)hParticleSignalSplitNumerator->Clone(Form("%sSignalSplitting", particleName.c_str()));
    hParticleSignalSplitting->Divide(hParticleSignalSplitting, hReconstructedMC, 1.0, 1.0,"B");
    hParticleSignalSplitting->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hParticleSignalSplitting->GetYaxis()->SetTitle("Signal Splitting");

    hStatError->GetXaxis()->SetBinLabel(5, "Signal splitting");
    hRelStatError->GetXaxis()->SetBinLabel(5, "Signal splitting");

   //V0 Reconstruction Efficiency for this particle
    TH1* hParticleEfficiency = new  TH1D(Form("%s Reconstruction Efficiency", particleName.c_str()), Form("%s Reconstruction Efficiency", particleName.c_str()), ptbinedges.size()-1,ptbinedges.data());
    hParticleEfficiency->Divide(hReconstructedMC, hGeneratedAfterEventSelection, 1.0, 1.0, "B");
    hParticleEfficiency->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hParticleEfficiency->GetYaxis()->SetTitle("Efficiency");

    hStatError->GetXaxis()->SetBinLabel(3, "V0 efficiency");
    hRelStatError->GetXaxis()->SetBinLabel(3, "V0 efficiency");

   //V0 Ultimate Efficiency
    TH1* hParticleUltimateEfficiency = new  TH1D(Form("%s Ultimate Reconstruction Efficiency", particleName.c_str()), Form("%s Ultimate Reconstruction Efficiency", particleName.c_str()), ptbinedges.size()-1,ptbinedges.data());
    hParticleUltimateEfficiency->Divide(hParticleSignalSplitNumerator, hGenerated, 1.0, 1.0, "B");

    TH1* hRawParticle = (TH1*)hMeasured->Clone(Form("%sRawPtSpectrum", particleName.c_str()));
    hRawParticle->Scale(1.0/nCollisionsData);
    NormalizeByBinWidth(hRawParticle);

//////Event Level Corrections Applied
    double EventeLvlCorrections = (EventEfficiency * EventSplitting) / nCollisionsData;
    cout<<"////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"/////////////////////////Event Level Corrections////////////////////"<<endl;
    cout<<"////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"Event Efficiency:   "<<EventEfficiency<<endl;
    cout<<"Event Splitting:    "<<EventSplitting<<endl;
    cout<<"1 / nCollisions:   "<<1.0/nCollisionsData<<endl;
    cout<<"nCollisions:   "<<nCollisionsData<<endl;
    cout<<"Generated Collisiions:  "<<nCollisionsGen<<"  Corrected Reconstructed Collisiions:  "<<1.0/EventeLvlCorrections<<endl; 
    cout<<"////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"/////////////////////////Event Level Corrections////////////////////"<<endl;
    cout<<"////////////////////////////////////////////////////////////////////"<<endl;

//////Particle Reconstruction Level Corrections Applied//////

    TH1* hCorrectedParticle = new TH1D(Form("%sCorrectedPtSpectrum", particleName.c_str()), Form("%sCorrectedPtSpectrum", particleName.c_str()), ptbinedges.size()-1, ptbinedges.data()); 
    hCorrectedParticle->Add(hMeasured, hTotalFeeddown,1.0, -1.0); 
    //hCorrectedParticle->Divide(hParticleEfficiency);     //errors for independent histos
    //hCorrectedParticle->Divide(hParticleSignalLoss);
    //hCorrectedParticle->Divide(hParticleSignalSplitting);
    hCorrectedParticle->Divide(hParticleUltimateEfficiency);

   RecordStatError(hStatError,hRelStatError,8,hCorrectedParticle);

    hStatError->GetXaxis()->SetBinLabel(8, "Final");
    hRelStatError->GetXaxis()->SetBinLabel(8, "Final");

   //DEB
    int nbins = hCorrectedParticle->GetNbinsX();
    std::cout << std::setw(10) << "Bin" 
    << std::setw(20) << "Measured " 
    << std::setw(20) << "MC Reconstructed (Split)" 
    << std::setw(20) << "MC_Reconstructed_True" 
    << std::setw(20) << "Efficiency" 
    << std::setw(20) << "Corrected"
    << std::setw(20) << "Generated"
    << std::endl;

    for (int i = 1; i <= nbins; ++i) {
        double measured = hMeasured->GetBinContent(i);
        double splitnumerator = hParticleSignalSplitNumerator->GetBinContent(i);
        double mcreconstructed = hReconstructedMC->GetBinContent(i);
        double eff      = hParticleEfficiency->GetBinContent(i);
        double corr     = hCorrectedParticle->GetBinContent(i);
        double gener     = hGeneratedAfterEventSelection->GetBinContent(i);

        std::cout << std::setw(10) << i
        << std::setw(20) << measured
        << std::setw(20) << splitnumerator
        << std::setw(20) << mcreconstructed
        << std::setw(20) << eff
        << std::setw(20) << corr
        << std::setw(20) << gener
        << std::endl;
    }

    //Final Corrected pt Spectrum
    hCorrectedParticle->Multiply(hParticleEventEfficiency);
    hCorrectedParticle->Multiply(hParticleEventSplitting);
    hCorrectedParticle->Scale(1.0 / nCollisionsData);
    NormalizeByBinWidth(hCorrectedParticle); //Normalize before fitting

   //Levy-Tsallis Fit 
    double mParticle = (particle == "K0sh") ? 0.497611 : 1.115683; // K0sh K0s mass vs Lambda mass
    TF1 *levyTsallis = new TF1("levyTsallis", "x*[0]*((([1]-1.)*([1]-2.))/(([1]*[2])*(([1]*[2])+([3]*([1]-2.))))) * pow(1.+((sqrt(pow([3], 2.)+pow(x, 2.))-[3])/([1]*[2])), -[1])", 0.0, 16.0);
    levyTsallis->SetParameter(3, mParticle);  // Set m as a parameter
    levyTsallis->FixParameter(3, mParticle);  // Fix m so it does not change during fitting as it is fixed
    levyTsallis->SetParameters(1.0, 7.0, 0.25);
    levyTsallis->SetParLimits(0, 0.01, 5.0);   // A
    levyTsallis->SetParLimits(1, 3.0, 20.0);   // n
    levyTsallis->SetParLimits(2, 0.1, 0.6);    // T
    hCorrectedParticle->Fit(levyTsallis, "RN");

   // Create a canvas for plotting the Fit
    TCanvas *FinalFit = new TCanvas("Final Fit", "Levy-Tsallis Fit", 800, 600);
    hCorrectedParticle->Draw("E");  // "E" for errors
    levyTsallis->SetLineColor(kRed);
    levyTsallis->SetLineWidth(2);
    levyTsallis->Draw("same"); 
   //legend
    TLegend *legend = new TLegend(0.6, 0.7, 0.85, 0.85);
    legend->AddEntry(hCorrectedParticle, Form("Corrected %s P_{T} Spectrum", particleName.c_str()), "lep");
    legend->AddEntry(levyTsallis, "Levy-Tsallis Fit", "l");
    legend->Draw();
    FinalFit->Update(); // Update the canvas to display everything

    hMeasured->SetOption("PE");
    hReconstructedMC->SetOption("PE");



    TH1* hStatErrorEvolutionCheck = (TH1*)hMeasured->Clone("hStatErrorEvolutionCheck");

// Stage 1: measured signal
RecordStatError(hStatError, hRelStatError, 1,hStatErrorEvolutionCheck);

// Stage 2: feeddown subtraction
hStatErrorEvolutionCheck->Add(hTotalFeeddown,-1.0);

RecordStatError(hStatError, hRelStatError, 2,hStatErrorEvolutionCheck);

// Stage 3: V0 efficiency
hStatErrorEvolutionCheck->Divide(hParticleEfficiency);
RecordStatError(hStatError, hRelStatError, 3,hStatErrorEvolutionCheck);

// Stage 4: signal loss
hStatErrorEvolutionCheck->Divide(hParticleSignalLoss);
RecordStatError(hStatError, hRelStatError, 4,hStatErrorEvolutionCheck);

// Stage 5: signal splitting
hStatErrorEvolutionCheck->Divide(hParticleSignalSplitting);
RecordStatError(hStatError, hRelStatError, 5,hStatErrorEvolutionCheck);

// Stage 6: event efficiency
hStatErrorEvolutionCheck->Divide(hParticleEventEfficiency);
RecordStatError(hStatError, hRelStatError, 6,hStatErrorEvolutionCheck);

// Stage 7: event splitting
hStatErrorEvolutionCheck->Divide(hParticleEventSplitting);
RecordStatError(hStatError, hRelStatError, 7,hStatErrorEvolutionCheck);

// Stage 8: final
RecordStatError(hStatError, hRelStatError, 8,hCorrectedParticle);

// ============================================================
// Print relative statistical uncertainty evolution
// ============================================================

std::cout << "\n============================================" << std::endl;
std::cout << "Relative statistical uncertainty (%)" << std::endl;
std::cout << "============================================" << std::endl;

for (int stage = 1; stage <= 8; ++stage) {

    std::cout << "\n=== Stage " << stage << " ===" << std::endl;

    for (int bin = 1; bin <= hRelStatError->GetNbinsY(); ++bin) {

        double ptLow  = hRelStatError->GetYaxis()->GetBinLowEdge(bin);
        double ptHigh = hRelStatError->GetYaxis()->GetBinUpEdge(bin);
        double relError =hRelStatError->GetBinContent(stage, bin);
        std::cout << "pT = [" << ptLow << ", "  << ptHigh  << "] : " << relError << " %" << std::endl;
    }
}
//
//
//Start of Alpha Test
//
//
    TCanvas *cAlphaTest = new TCanvas("cAlphaTest", "Alpha test", 800, 800);
    TPad *padAlpha1 = new TPad("padAlpha1", "Main pad", 0, 0.3, 1, 1.0);
    TPad *padAlpha2 = new TPad("padAlpha2", "Ratio pad", 0, 0.05, 1, 0.3);

    padAlpha1->Draw();
    padAlpha2->Draw();

    // -- Top pad: draw corrected and generated spectra --
    padAlpha1->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    padAlpha1->SetLogy(1);
    padAlpha1->SetBottomMargin(0.01);

    // Clone the original histogram to keep original unchanged and multiply by nCollisions so they are comparible
    TH1D *hGeneratedAlpha = (TH1D*)hGenerated->Clone("hGeneratedAlpha");
    hGeneratedAlpha->Scale(1.0/nCollisionsGen);

    TH1D *hCorrectedAlpha = (TH1D*)hCorrectedParticle->Clone("hCorrectedAlpha");
    //hCorrectedAlpha->Scale(nCollisions);
    hCorrectedAlpha->SetMarkerStyle(kFullCircle);
    hCorrectedAlpha->SetMarkerColor(kBlue+1);
    hCorrectedAlpha->SetLineColor(kBlue+1);
    hCorrectedAlpha->SetMarkerSize(1);
    hCorrectedAlpha->GetYaxis()->SetTitle("1/N_{INEL>0} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
    hCorrectedAlpha->Draw("PE");

    hGeneratedAlpha->SetMarkerStyle(22);
    hGeneratedAlpha->SetMarkerColor(kRed+1);
    hGeneratedAlpha->SetLineColor(kRed+1);
    hGeneratedAlpha->SetMarkerSize(1);
    NormalizeByBinWidth(hGeneratedAlpha);
    hGeneratedAlpha->Draw("PE SAME");

    TLegend *legAlpha = new TLegend(0.55, 0.70, 0.88, 0.88);
    legAlpha->SetTextSize(0.03);
    legAlpha->SetBorderSize(0);
    legAlpha->SetFillColor(kWhite);
    legAlpha->AddEntry(hCorrectedAlpha, "Corrected Spectrum", "lp");
    legAlpha->AddEntry(hGeneratedAlpha, "Generated Spectrum", "lp");
    legAlpha->Draw();

    TLatex latexAlpha;
    latexAlpha.SetNDC();
    latexAlpha.SetTextSize(0.03);
    latexAlpha.DrawLatex(0.6, 0.65, "MC Closure Test (Alpha)");
    latexAlpha.DrawLatex(0.6, 0.60, "pp");
    latexAlpha.DrawLatex(0.6, 0.55, "#sqrt{#it{s}} = " +AO2DFiles.energy);
    latexAlpha.DrawLatex(0.6, 0.50, "|y| < 0.5");

        // -- Bottom pad: ratio Corrected / Generated --
    padAlpha2->cd();
    padAlpha2->SetGridy();
    padAlpha2->SetTopMargin(0.01);
    padAlpha2->SetBottomMargin(0.3);

    // Checking errors
    hGeneratedAlpha->Sumw2();
    hCorrectedAlpha->Sumw2();

        // Create ratio histogram
    TH1D *hRatioAlpha = (TH1D*)hCorrectedAlpha->Clone("hRatioAlpha");
    hRatioAlpha->SetTitle("");
    hRatioAlpha->Divide(hCorrectedAlpha, hGeneratedAlpha, 1.0, 1.0);

    hRatioAlpha->GetYaxis()->SetTitle("Corrected / Generated");
    hRatioAlpha->GetYaxis()->SetTitleSize(0.08);
    hRatioAlpha->GetYaxis()->SetTitleOffset(0.5);
    hRatioAlpha->GetYaxis()->SetLabelSize(0.08);

    hRatioAlpha->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatioAlpha->GetXaxis()->SetTitleSize(0.10);
    hRatioAlpha->GetXaxis()->SetTitleOffset(1.0);
    hRatioAlpha->GetXaxis()->SetLabelSize(0.08);

    hRatioAlpha->SetMarkerStyle(33);
    hRatioAlpha->SetMarkerSize(1);
    hRatioAlpha->SetMarkerColor(kCyan+2);
    hRatioAlpha->SetLineColor(kCyan+2);
    for (int i = 1; i <= hRatioAlpha->GetNbinsX(); ++i) {
        hRatioAlpha->SetBinError(i, 0);
    }
    hRatioAlpha->Draw("PE");

        // Reference line at 1
    TLine *lineAlpha = new TLine(hRatioAlpha->GetXaxis()->GetXmin(), 1.0,
    hRatioAlpha->GetXaxis()->GetXmax(), 1.0);
    lineAlpha->SetLineStyle(2);
    lineAlpha->Draw("SAME");

    cAlphaTest->Update();
    //
    //End of Alpha Test
    //
    //
    //Save everything
    //
   //Function for safe-writting objects into intended directory and warning if they don't exist
    auto safeWriteObject = [&] (TObject* obj, const char* name, TDirectory* dir) {
        if (!obj) {
            std::cerr << "Warning: Object pointer '" << name << "' is null, skipping WriteObject." << std::endl;
            return;
        }

      // Update canvases
        if (TCanvas* canvas = dynamic_cast<TCanvas*>(obj)) {
            canvas->cd();
            canvas->Update();
        }

        if (!dir) {
            std::cerr << "Error: Directory pointer is null for object '" << name << "'!" << std::endl;
            return;
        }

      // Go to the intended directory
        dir->cd();

      // Write into the current directory (the subdirectory)
        obj->Write(name);
    };

   // Use function to safe save everything into centdir
   //safeWriteObject(cBetaTest, "K0sh Correction Beta Test", centdir);
   //safeWriteObject(FinalFit, "Final Fit of Corrected Pt Spectrum", centdir);

   //For Debugging - write particle-dependent objects using mapped names
    safeWriteObject(hPeak, Form("%s Peak vs Pt Quality Control", particle.c_str()), centdir);
    safeWriteObject(hWidth, Form("%s Width vs Pt Quality Control", particle.c_str()), centdir);
    safeWriteObject(hNEventsCorrections, "Different Event Corrections", centdir);
    safeWriteObject(hParticleCount, Form("Number of %s", particleName.c_str()), centdir);
    safeWriteObject(hNEventsMC, "Number of Events MC", centdir);
    safeWriteObject(hNV0sMC, "Number of V0s MC", centdir);
    safeWriteObject(hParticleCountMC, Form("Number of %s MC", particleName.c_str()), centdir);
    safeWriteObject(hMeasured, Form("Measured %s Pt Spectrum", particleName.c_str()), centdir);
    safeWriteObject(hGeneratedAfterEventSelection, Form("Generated %s Pt Spectrum After Event Selection", particleName.c_str()), centdir);
    safeWriteObject(hParticleUltimateEfficiency, Form("%s Ultimate Signal Efficiency", particle.c_str()), centdir);
   
   //Final Results
    safeWriteObject(hParticleEfficiency, Form("%s Signal Efficiency", particle.c_str()), centdir);
    safeWriteObject(hParticleSignalLoss, Form("%s Signal Loss due to Event Selection and Reconstruction", particle.c_str()), centdir);
    safeWriteObject(hParticleSignalSplitting, Form("%s Signal Splitting", particle.c_str()), centdir);
    safeWriteObject(hCorrectedParticle, Form("%s Corrected Pt Spectrum", particle.c_str()), centdir);
    safeWriteObject(hParticleEventEfficiency, Form("%s Event Efficiency", particle.c_str()), centdir);
    safeWriteObject(hParticleEventSplitting, Form("%s Event Splitting", particle.c_str()), centdir);
    safeWriteObject(levyTsallis, Form("%s Levy-Tsallis Fit", particle.c_str()), centdir);
    safeWriteObject(hReconstructedMC, Form("%s MC Reconstructed Pt Spectrum", particle.c_str()), centdir);
    safeWriteObject(hNEventsCorrections, Form("%s Different Event Corrections", particle.c_str()), centdir);
    safeWriteObject(hRawParticle, Form("%s Raw Pt Spectrum", particle.c_str()), centdir);
    safeWriteObject(hGenerated, Form("%s Generated Pt Spectrum", particleName.c_str()), centdir);
    safeWriteObject(cAlphaTest, Form("%s Alpha test", particle.c_str()),centdir);
    safeWriteObject(hRatioAlpha, Form("%s Ratio Alpha Test", particleName.c_str()), centdir);
    safeWriteObject(hNEventsData, Form("%s Number of Events", particle.c_str()), centdir);
    safeWriteObject(hNV0s, Form("%s Number of V0s", particle.c_str()), centdir);
    safeWriteObject(hParticleCount, Form("%s Number of Particles", particleName.c_str()), centdir);
    safeWriteObject(hPeak, Form("%s Peak vs Pt Quality Control", particle.c_str()), centdir);
    safeWriteObject(hWidth, Form("%s Width vs Pt Quality Control", particle.c_str()), centdir);
    safeWriteObject(hStatErrorEvolution, Form("%s Statistical Error Evolution", particle.c_str()), centdir);
    safeWriteObject(hStatError, Form("%s Statistical Error", particle.c_str()), centdir);
    safeWriteObject(hRelStatError, Form("%s Relative Statistical Error", particle.c_str()), centdir);
    safeWriteObject(hReducedChi2, Form("%s Reduced Chi2", particle.c_str()), centdir);
    safeWriteObject(hReducedChi2MC, Form("%s Reduced Chi2 MC", particle.c_str()), centdir);

    if (particle == "Lambda" || particle == "AntiLambda") {
        safeWriteObject(hTotalFeeddown, Form("%s Total Feeddown", particle.c_str()), centdir);
        safeWriteObject(hTotalFeeddownFraction, Form("%s Total Feeddown Fraction", particle.c_str()),  centdir);
      //Normalize before sabing for visualization purposes
        NormalizeByBinWidth(hXiCorrectedSpectrumInput);
        NormalizeByBinWidth(hGenXiSign);
        NormalizeByBinWidth(hGenXiNeutral);
      //NormalizeByBinWidth(hFeeddownXiSignMatrix);
      //NormalizeByBinWidth(hFeeddownXiNeutralMatrix);
        NormalizeByBinWidth(hFromXiSignSpectrum);
        NormalizeByBinWidth(hFromXiNeutralSpectrum);
      // Set y-axis range for better visualization
        hGenXiSign->GetYaxis()->SetRangeUser(1e-6, 400000);
        hGenXiNeutral->GetYaxis()->SetRangeUser(1e-6, 400000);
      ////////////////////////////////////////////////////
        safeWriteObject(hFromXiSignSpectrum,Form("%s Feeddown from %s Spectrum", particle.c_str(), xiSignName.c_str()),centdir);
        safeWriteObject(hFromXiNeutralSpectrum,Form("%s Feeddown from %s Spectrum", particle.c_str(), xiNeutralName.c_str()),centdir);
        safeWriteObject(hXiCorrectedSpectrumInput, Form("%s Corrected Run 3 Xi Spectrum Used for Feeddown", particle.c_str()), centdir);
        safeWriteObject(hGenXiSign, Form("%s Generated %s Spectrum Used for Feeddown", particle.c_str(), xiSignName.c_str()),centdir);
        safeWriteObject(hGenXiNeutral, Form("%s Generated %s Spectrum Used for Feeddown", particle.c_str(), xiNeutralName.c_str()),centdir);
        safeWriteObject(hFeeddownXiSignMatrix, Form("%s from %s Feeddown Matrix Used for Feeddown", particle.c_str(),  xiSignName.c_str()),centdir);
        safeWriteObject(hFeeddownXiNeutralMatrix, Form("%s from %s Feeddown Matrix Used for Feeddown", particle.c_str(), xiNeutralName.c_str()),centdir);
    }
}
//
//*********************************************************************************************************************************************************************************************************************************************************************/
//*********************************************************************************************************************************************************************************************************************************************************************/
//*********************************************************************************************************************************************************************************************************************************************************************/
/////////////////////////////////////////////////////////////////////////////////
//Start of Code
//////////////////////////////////////////////////////////////////////////////////
//*********************************************************************************************************************************************************************************************************************************************************************/
//*********************************************************************************************************************************************************************************************************************************************************************/
//*********************************************************************************************************************************************************************************************************************************************************************/
int PtSpectrumV0sCent(){  
   //AO2D Files to be accessed configuqration
    InitAO2DFiles();

   //********************************************************************************/
    //Particles to be processed
   //********************************************************************************/
   std::vector<std::string> particles = {"K0sh", "Lambda", "AntiLambda"};
   //std::vector<std::string> particles = {"K0sh"};
   //std::vector<std::string> particles = {"Lambda"};
   //std::vector<std::string> particles = {"AntiLambda"};
   //std::vector<std::string> particles = {"Lambda", "AntiLambda"};
    bool NchAnalysis = true; //Run Nch Alysis
    //********************************************************************************/
    //Systematics to be processed?
   //********************************************************************************/

   bool Systematics_Analysis = false; // Run Systematics

   //********************************************************************************/
   //********************************************************************************/
   //File Where Everything will be Saved
    TString outputRootFile = "V0s_Pt_Spectrum_Cent_"+AO2DFiles.energy+".root";
    TFile *F1 = new TFile(outputRootFile, "RECREATE");
    //********************************************************************************/

    int romanIndex = 0;
   // Color palette
    std::vector<Color_t> colors = {kRed+1, kRed-4, kOrange+7, kOrange-3, kYellow+1, kSpring+5, kGreen+2, kGreen-7, kTeal+3, kCyan+2, kAzure+2, kBlue-4, kBlue+2};
   //Different Colors for Centrality Classes
    std::vector<TString> roman = {"I", "II", "III", "IV", "V", "VI", "VII","VIII", "IX", "X", "XI", "XII", "XIII"};
    bool k0shProcessed = false;   
    bool lambdaProcessed = false;
    bool antiLambdaProcessed = false;
    TH1F* hLambdaPlusAntiLambdaForComparison = nullptr;
    TH1F* hK0shCorrectedCompare = nullptr;
    std::vector<TH1F*> hEffLambdaStorage , hEffAntiLambdaStorage;
    std::vector<TH1F*> hFeeddownFracLambdaStorage, hFeeddownFracAntiLambdaStorage;
    std::vector<TH1F*> hPtSpectrumCorrectedLambdaStorage, hPtSpectrumCorrectedAntiLambdaStorage;
    const std::vector<double>* nCentBinsKaons;
    const std::vector<double>* nCentBinsLambdas;
    nCentBinsKaons   = &pthistos::centBinsK0sh;
    nCentBinsLambdas = &pthistos::centBinsLambda;

    for (const auto& particle : particles) {

        cout<<"********************************************************************************************************"<<endl;
        cout<<"********************************************************************************************************"<<endl;
        cout<<"***********************   Processing particle: " << particle <<"   *****************************************"<<endl;
        cout<<"********************************************************************************************************"<<endl;
        cout<<"********************************************************************************************************"<<endl;

      // Create and enter particle directory
        F1->mkdir(particle.c_str());
        F1->cd(particle.c_str());
        TDirectory* particleDir = F1->GetDirectory(particle.c_str()); //Go to partice directory

      // Select pt bins depending on particle
        const std::vector<std::string>* ptBins = nullptr;
        const std::vector<double>* centBins = nullptr;
        if (particle == "K0sh"){
            ptBins = &pthistos::ptbinsK0sh;
            centBins = &pthistos::centBinsK0sh;
        } 
        else if (particle == "Lambda" || particle == "AntiLambda") { // same binning for Lambda and AntiLambda
            ptBins = &pthistos::ptbinslambda;
            centBins = &pthistos::centBinsLambda;
        }
      
        int centralityCounter = 0; 
        for (size_t i = 0; i + 1 < centBins->size(); ++i) { // use centBins pointer
            double cent_low  = (*centBins)[i];
            double cent_high = (*centBins)[i + 1];

            TString centDirName = Form("Cent_%.0f_%.0f", cent_low, cent_high);
            TDirectory* centDir = gDirectory->GetDirectory(centDirName);
            if (!centDir) {
                centDir = gDirectory->mkdir(centDirName);
            }
            centDir->cd();  // go into the newly created (or existing) directory

            int centralityIndex = -1;
            if (particle == "Lambda" || particle == "AntiLambda") { //Add the current centrality index so we know which Xi Spectrum to use for feeddown correction
                centralityIndex = centralityCounter;
                centralityCounter++;
            }
         
            PlotSignal(particle, *ptBins, F1, AO2DFiles, centDir, cent_low, cent_high, centralityIndex, Systematics_Analysis); //Plots Signal for each centrality bin within particle directory

      // Go back to particle directory
            F1->cd(particle.c_str());
        }
      //
      //
      // Start Drawing
      // Check whether any centrality directories with results were produced for this particle.
        bool anyCentFound = false;
        for (size_t j = 0; j + 1 < centBins->size(); ++j) {
            double cl = (*centBins)[j];
            double ch = (*centBins)[j + 1];
            TString checkDir = Form("%s/Cent_%.0f_%.0f", particle.c_str(), cl, ch);
            if (F1->GetDirectory(checkDir)) { anyCentFound = true; break;}
        }
        if (!anyCentFound) {
            std::cerr << "Warning: No centrality results produced for particle " << particle << ", skipping drawing and comparisons." << std::endl;
            F1->cd();
            continue; // next particle
        }

        gDirectory->mkdir("All_Centralities"); 
        gDirectory->cd("All_Centralities");

        int colorIndex = 0;
      // store histograms for drawing
        std::vector<TH1F*> hCorrectedStorage;
        std::vector<TH1F*> hGeneratedStorage;
        std::vector<TH1F*> hEfficiencyStorage;
        std::vector<TH1F*> hSignalSplittingStorage;
        std::vector<TH1F*> hSignalLossStorage;
        std::vector<TH1F*> hUltimateEfficiencyStorage;
        std::vector<TH1F*> hEventEfficiencyStorage;
        std::vector<TH1F*> hEventSplittingStorage;
        std::vector<TH1F*> hRawStorage;
        std::vector<TH1F*> hReconstructedMCStorage;
        std::vector<TF1*> fLevyTsallisStorage;
        std::vector<TH1F*> hNEventGenStorage;
        std::vector<TH1F*> hTotalFeeddownStorage;
        std::vector<TH1F*> hTotalFeeddownFractionStorage;
        std::vector<TH1F*> hFromXiSignSpectrumStorage;
        std::vector<TH1F*> hFromXiNeutralSpectrumStorage;
        std::vector<TH1F*> hXiCorrectedSpectrumInputStorage;
        std::vector<TH1F*> hGenXiSignStorage, hGenXiNeutralStorage;
        std::vector<TH1F*> hFeeddownXiSignMatrixStorage, hFeeddownXiNeutralMatrixStorage;
        std::vector<TH1F*> hRatioAlphaStorage;
        std::vector<TH1F*> hNEventsStorage;
        std::vector<TH1F*> hNV0sStorage;
        std::vector<TH1F*> hNParticlesStorage;
        std::vector<TH1F*> hPeakStorage, hWidthStorage;
        std::vector<TH1F*> hReducedChi2Storage, hReducedChi2MCStorage;
        //std::vector<TF1*> fLevyTsallisLambdaStorage, fLevyTsallisAntiLambdaStorage; NOT NEEDED?

      // Load Centrality Histograms to Superimpose
        LoadCentralityHistograms(F1, particle, "Corrected Pt Spectrum", *centBins, hCorrectedStorage);
        LoadCentralityFunctions(F1, particle, "Levy-Tsallis Fit", *centBins, fLevyTsallisStorage);
        LoadCentralityHistograms(F1, particle, "Raw Pt Spectrum", *centBins, hRawStorage);
        LoadCentralityHistograms(F1, particle, "Signal Efficiency", *centBins, hEfficiencyStorage);
        LoadCentralityHistograms(F1, particle, "Ultimate Signal Efficiency", *centBins, hUltimateEfficiencyStorage);
        LoadCentralityHistograms(F1, particle, "Signal Splitting", *centBins, hSignalSplittingStorage);
        LoadCentralityHistograms(F1, particle, "Signal Loss due to Event Selection and Reconstruction", *centBins, hSignalLossStorage);
        LoadCentralityHistograms(F1, particle, "Event Efficiency", *centBins, hEventEfficiencyStorage);
        LoadCentralityHistograms(F1, particle, "Event Splitting", *centBins, hEventSplittingStorage);
        LoadCentralityHistograms(F1, particle, "MC Reconstructed Pt Spectrum", *centBins, hReconstructedMCStorage);
        LoadCentralityHistograms(F1, particle, "Different Event Corrections", *centBins, hNEventGenStorage);
        LoadCentralityHistograms(F1, particle, "Generated Pt Spectrum", *centBins, hGeneratedStorage);
        LoadCentralityHistograms(F1, particle, "Ratio Alpha Test", *centBins, hRatioAlphaStorage);
        LoadCentralityHistograms(F1, particle, "Number of Events", *centBins, hNEventsStorage);
        LoadCentralityHistograms(F1, particle, "Number of V0s", *centBins, hNV0sStorage);
        LoadCentralityHistograms(F1, particle, "Number of Particles", *centBins, hNParticlesStorage);
        LoadCentralityHistograms(F1, particle, "Peak vs Pt Quality Control", *centBins, hPeakStorage);
        LoadCentralityHistograms(F1, particle, "Width vs Pt Quality Control", *centBins, hWidthStorage);
        LoadCentralityHistograms(F1, particle, "Reduced Chi2", *centBins, hReducedChi2Storage);
        LoadCentralityHistograms(F1, particle, "Reduced Chi2 MC", *centBins, hReducedChi2MCStorage);

        if (particle == "Lambda") {
            LoadCentralityHistograms(F1, particle, "Total Feeddown", *centBins, hTotalFeeddownStorage);
            LoadCentralityHistograms(F1, particle, "Total Feeddown Fraction", *centBins, hTotalFeeddownFractionStorage);
            LoadCentralityHistograms(F1, particle, "Feeddown from XiMinus Spectrum", *centBins, hFromXiSignSpectrumStorage);
            LoadCentralityHistograms(F1, particle, "Feeddown from XiZero Spectrum", *centBins, hFromXiNeutralSpectrumStorage);
            LoadCentralityHistograms(F1, particle, "Corrected Run 3 Xi Spectrum Used for Feeddown", *centBins, hXiCorrectedSpectrumInputStorage);
            LoadCentralityHistograms(F1, particle, "Generated XiMinus Spectrum Used for Feeddown", *centBins, hGenXiSignStorage);
            LoadCentralityHistograms(F1, particle, "Generated XiZero Spectrum Used for Feeddown", *centBins, hGenXiNeutralStorage);
            LoadCentralityHistograms(F1, particle, "from XiMinus Feeddown Matrix Used for Feeddown", *centBins, hFeeddownXiSignMatrixStorage);
            LoadCentralityHistograms(F1, particle, "from XiZero Feeddown Matrix Used for Feeddown", *centBins, hFeeddownXiNeutralMatrixStorage);
         //Store for Comparison
            hPtSpectrumCorrectedLambdaStorage = CloneStorage(hCorrectedStorage, "hPtSpectrumCorrected");
            hEffLambdaStorage  = CloneStorage(hEfficiencyStorage, "hEffLambda");
            hFeeddownFracLambdaStorage = CloneStorage(hTotalFeeddownFractionStorage, "hFeeddownFracLambda");
        }

        if (particle == "AntiLambda") { //Store for Comparison
            LoadCentralityHistograms(F1, particle, "Total Feeddown", *centBins, hTotalFeeddownStorage);
            LoadCentralityHistograms(F1, particle, "Total Feeddown Fraction", *centBins, hTotalFeeddownFractionStorage);
            LoadCentralityHistograms(F1, particle, "Feeddown from XiPlus Spectrum", *centBins, hFromXiSignSpectrumStorage);
            LoadCentralityHistograms(F1, particle, "Feeddown from AntiXiZero Spectrum", *centBins, hFromXiNeutralSpectrumStorage);
            LoadCentralityHistograms(F1, particle, "Corrected Run 3 Xi Spectrum Used for Feeddown", *centBins, hXiCorrectedSpectrumInputStorage);
            LoadCentralityHistograms(F1, particle, "Generated XiPlus Spectrum Used for Feeddown", *centBins, hGenXiSignStorage);
            LoadCentralityHistograms(F1, particle, "Generated AntiXiZero Spectrum Used for Feeddown", *centBins, hGenXiNeutralStorage);
            LoadCentralityHistograms(F1, particle, "from XiPlus Feeddown Matrix Used for Feeddown", *centBins, hFeeddownXiSignMatrixStorage);
            LoadCentralityHistograms(F1, particle, "from AntiXiZero Feeddown Matrix Used for Feeddown", *centBins, hFeeddownXiNeutralMatrixStorage);
         //Store for Comparison
            hPtSpectrumCorrectedAntiLambdaStorage = CloneStorage(hCorrectedStorage, "hPtSpectrumCorrectedAntiLambda");
            hEffAntiLambdaStorage = CloneStorage(hEfficiencyStorage, "hEffAntiLambda");
            hFeeddownFracAntiLambdaStorage = CloneStorage(hTotalFeeddownFractionStorage, "hFeeddownFracAntiLambda");
        }
      //Comparison Canvases
        TCanvas* cEffCompare = nullptr;
        TCanvas* cFracCompare = nullptr;

      // Draw canvases
        TCanvas* cPt = DrawCentralityHistoSet(hCorrectedStorage, roman, colors, *centBins, "Corrected_PtSpectra_AllCentralities_"+TString(particle), false /*logY*/, false, false , true);
        TCanvas* cLevy = DrawCentralityFunctionSet(fLevyTsallisStorage, roman, colors, *centBins, "LevyTsallis_AllCentralities_"+TString(particle), true /*logY*/, false);
        TCanvas* cRaw = DrawCentralityHistoSet(hRawStorage, roman, colors, *centBins, "RawSpectra_AllCentralities_"+TString(particle), true /*logY*/, false, false , true);
        TCanvas* cEff = DrawCentralityHistoSet(hEfficiencyStorage, roman, colors, *centBins,"Efficiency_AllCentralities_"+TString(particle), false /*logY*/, true  /*fix axis to 1.0*/, false, false);
        TCanvas* cUltEff = DrawCentralityHistoSet(hUltimateEfficiencyStorage, roman, colors, *centBins,"Ultimate_Efficiency_AllCentralities_"+TString(particle), false /*logY*/, true  /*fix axis to 1.0*/, false, false);
        TCanvas* cSigSplit = DrawCentralityHistoSet(hSignalSplittingStorage, roman, colors, *centBins,"Signal_Splitting_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cSigLoss = DrawCentralityHistoSet(hSignalLossStorage, roman, colors, *centBins,"Signal_Loss_AllCentralities_"+TString(particle), false /*logY*/, true  /*fix axis to 1.0*/, false, false);
        TCanvas* cEvEff = DrawCentralityHistoSet(hEventEfficiencyStorage, roman, colors, *centBins,"Event_Efficiency_AllCentralities_"+TString(particle), false /*logY*/, true  /*fix axis to 1.0*/, true, false);
        TCanvas* cEvSplit = DrawCentralityHistoSet(hEventSplittingStorage, roman, colors, *centBins,"Event_Splitting_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cRecoMC = DrawCentralityHistoSet(hReconstructedMCStorage, roman, colors, *centBins,"ReconstructedMC_PtSpectrum_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cEvent = DrawCentralityHistoSet(hNEventGenStorage, roman, colors, *centBins,"NEventsGen_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cGen = DrawCentralityHistoSet(hGeneratedStorage, roman, colors, *centBins,"NGenerated_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cAlphaTestRatio = DrawCentralityHistoSet(hRatioAlphaStorage, roman, colors, *centBins,"Ratio_Alpha_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cNEvents = DrawCentralityHistoSet(hNEventsStorage, roman, colors, *centBins,"Number_Events_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cNV0s = DrawCentralityHistoSet(hNV0sStorage, roman, colors, *centBins,"Number_V0s_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cNParticles = DrawCentralityHistoSet(hNParticlesStorage, roman, colors, *centBins,"Number_Particles_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cPeak = DrawCentralityHistoSet(hPeakStorage, roman, colors, *centBins,"Peak_vs_Pt_Quality_Control_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cWidth = DrawCentralityHistoSet(hWidthStorage, roman, colors, *centBins,"Width_vs_Pt_Quality_Control_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cReducedChi2 = DrawCentralityHistoSet(hReducedChi2Storage, roman, colors, *centBins,"Reduced_Chi2_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        TCanvas* cReducedChi2MC = DrawCentralityHistoSet(hReducedChi2MCStorage, roman, colors, *centBins,"Reduced_Chi2_MC_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);

      //Feeddown only for Lambda and AntiLambda
        TCanvas* cFeeddown = nullptr;
        TCanvas* cFeeddownFraction = nullptr;
        TCanvas* cFeeddownXiSign = nullptr;
        TCanvas* cFeeddownXiNeutral = nullptr;
        TCanvas* cXiInput = nullptr;
        TCanvas* cGenXiSign = nullptr;
        TCanvas* cGenXiNeutral = nullptr;
        TCanvas* cFeeddownXiSignMatrix = nullptr;
        TCanvas* cFeeddownXiNeutralMatrix = nullptr;
        if (particle == "Lambda") {
            cFeeddown = DrawCentralityHistoSet(hTotalFeeddownStorage, roman, colors, *centBins,"Total_Feeddown_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownFraction = DrawCentralityHistoSet(hTotalFeeddownFractionStorage, roman, colors, *centBins,"Total_Feeddown_Fraction_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiSign = DrawCentralityHistoSet(hFromXiSignSpectrumStorage, roman, colors, *centBins,"Feeddown_XiMinus_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiNeutral = DrawCentralityHistoSet(hFromXiNeutralSpectrumStorage, roman, colors, *centBins,"Feeddown_XiZero_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cXiInput = DrawCentralityHistoSet(hXiCorrectedSpectrumInputStorage, roman, colors, *centBins,"Xi_Input_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cGenXiSign = DrawCentralityHistoSet(hGenXiSignStorage, roman, colors, *centBins,"Generated_XiSign_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cGenXiNeutral = DrawCentralityHistoSet(hGenXiNeutralStorage, roman, colors, *centBins,"Generated_XiNeutral_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiSignMatrix = DrawCentralityHistoSet(hFeeddownXiSignMatrixStorage, roman, colors, *centBins,"Feeddown_XiSign_Matrix_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiNeutralMatrix = DrawCentralityHistoSet(hFeeddownXiNeutralMatrixStorage, roman, colors, *centBins,"Feeddown_XiNeutral_Matrix_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        }
        if (particle == "AntiLambda") {
            cFeeddown = DrawCentralityHistoSet(hTotalFeeddownStorage, roman, colors, *centBins,"Total_Feeddown_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownFraction = DrawCentralityHistoSet(hTotalFeeddownFractionStorage, roman, colors, *centBins,"Total_Feeddown_Fraction_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiSign = DrawCentralityHistoSet(hFromXiSignSpectrumStorage, roman, colors, *centBins,"Feeddown_XiPlus_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiNeutral = DrawCentralityHistoSet(hFromXiNeutralSpectrumStorage, roman, colors, *centBins,"Feeddown_AntiXiZero_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cXiInput = DrawCentralityHistoSet(hXiCorrectedSpectrumInputStorage, roman, colors, *centBins,"Xi_Input_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cGenXiSign = DrawCentralityHistoSet(hGenXiSignStorage, roman, colors, *centBins,"Generated_XiSign_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cGenXiNeutral = DrawCentralityHistoSet(hGenXiNeutralStorage, roman, colors, *centBins,"Generated_XiNeutral_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiSignMatrix = DrawCentralityHistoSet(hFeeddownXiSignMatrixStorage, roman, colors, *centBins,"Feeddown_XiSign_Matrix_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
            cFeeddownXiNeutralMatrix = DrawCentralityHistoSet(hFeeddownXiNeutralMatrixStorage, roman, colors, *centBins,"Feeddown_XiNeutral_Matrix_AllCentralities_"+TString(particle), false /*logY*/, false  /*fix axis to 1.0*/, false, false);
        }
      // 
      //
      /////////////////////////////////////////////////////////////////////////////////
      //Comparison to Run2
      //////////////////////////////////////////////////////////////////////////////////
      //
      //
      // Particle Specific Files
        TString run2File, run2Dir, run2Graph,legendName;

      // Access comparison files (lambda and lambdabar use the same)
        if (particle == "K0sh") {

            legendName = "K^{0}_{S}";

            if (AO2DFiles.dataperiod == "ppref24") {
                run2File  = "HEPData502Paper_ModifiedXi.root";
                run2Dir   = "Figure 1a";
                run2Graph = "Graph1D_y12";
            }
            else if (AO2DFiles.dataperiod == "pp24" || AO2DFiles.dataperiod == "pp22o") {
                run2File  = "HEPDataINELgt0K0sh.root";
                run2Dir   = "Table 1";
                run2Graph = "Graph1D_y11";
            }
        }
        
        else if (particle == "Lambda") {

            legendName = "#Lambda + #bar{#Lambda}";

            if (AO2DFiles.dataperiod == "ppref24") {
                run2File  = "HEPData502Paper_ModifiedXi.root";
                run2Dir   = "Figure 1b";
                run2Graph = "Graph1D_y12";
            }
            else if (AO2DFiles.dataperiod == "pp24" || AO2DFiles.dataperiod == "pp22o") {
                run2File  = "HEPDataINELgt0Lambda.root";
                run2Dir   = "Table 2";
                run2Graph = "Graph1D_y11";
            }
        }

        else if (particle == "AntiLambda") {

            legendName = "#Lambda + #bar{#Lambda}";

            if (AO2DFiles.dataperiod == "ppref24") {
                run2File  = "HEPData502Paper_ModifiedXi.root";
                run2Dir   = "Figure 1b";
                run2Graph = "Graph1D_y12";
            }
            else if (AO2DFiles.dataperiod == "pp24" || AO2DFiles.dataperiod == "pp22o") {
                run2File  = "HEPDataINELgt0AntiLambda.root";
                run2Dir   = "Table 2";
                run2Graph = "Graph1D_y11";
            }
        }

        cPt->cd();
      // Get the pad primitives
        TH1F* hCorrectedSum = (TH1F*)cPt->GetPrimitive(Form("Sum_Corrected_PtSpectra_AllCentralities_%s", particle.c_str()));
        if(!hCorrectedSum) {
            std::cerr << "hSum not found in canvas! Skipping comparison for particle " << particle << std::endl;
        } 
        else {
            TH1F* hCorrectedCompare = (TH1F*)hCorrectedSum->Clone("hCorrectedCompare");

            if(particle == "K0sh"){
                hK0shCorrectedCompare = hCorrectedCompare;
                k0shProcessed = true; 
            }

            if(particle == "Lambda"){
                if (!hLambdaPlusAntiLambdaForComparison){
                    hLambdaPlusAntiLambdaForComparison = (TH1F*)hCorrectedCompare->Clone("hLambdaPlusAntiLambdaForComparison");
                }
                else{ 
                    hLambdaPlusAntiLambdaForComparison->Add(hCorrectedCompare);
                }
                lambdaProcessed = true;
            }
            else if (particle == "AntiLambda") {
                if (!hLambdaPlusAntiLambdaForComparison) {
                    std::cerr << "Warning: Lambda histogram missing, cloning AntiLambda only!" << std::endl;
                    hLambdaPlusAntiLambdaForComparison = (TH1F*)hCorrectedCompare->Clone("hLambdaPlusAntiLambdaForComparison");
                } else {
                    hLambdaPlusAntiLambdaForComparison->Add(hCorrectedCompare);
                }
                antiLambdaProcessed = true;
            }
        }
      // Save
      //
      //
      //Corrected Centrality Pt Spectrum Ratio Comparison to hSum
      //
      //
        if(cPt) cPt->Write();
        if(cLevy)  cLevy->Write();
        if(cEff) cEff->Write();
        if(cUltEff) cUltEff->Write();
        if(cSigSplit) cSigSplit->Write();
        if(cSigLoss) cSigLoss->Write();
        if(cEvEff) cEvEff->Write();
        if(cEvSplit) cEvSplit->Write();
        if(cRaw) cRaw->Write();
        if(cRecoMC) cRecoMC->Write();
        if(cGen) cGen->Write();
        if(cEvent) cEvent->Write();
        if(cAlphaTestRatio) cAlphaTestRatio->Write();
        if(cNEvents) cNEvents->Write();
        if(cNV0s) cNV0s->Write();
        if(cNParticles) cNParticles->Write();
        if(cPeak) cPeak->Write();
        if(cWidth) cWidth->Write();
        if(cReducedChi2) cReducedChi2->Write();
        if(cReducedChi2MC) cReducedChi2MC->Write();

        if (particle == "Lambda" || particle == "AntiLambda") {
            if(cFeeddown) cFeeddown->Write();
            if(cFeeddownFraction) cFeeddownFraction->Write();
            if(cFeeddownXiSign) cFeeddownXiSign->Write();
            if(cFeeddownXiNeutral) cFeeddownXiNeutral->Write();
            if(cXiInput) cXiInput->Write();
            if(cGenXiSign) cGenXiSign->Write();
            if(cGenXiNeutral) cGenXiNeutral->Write();
            if(cFeeddownXiSignMatrix) cFeeddownXiSignMatrix->Write();
            if(cFeeddownXiNeutralMatrix) cFeeddownXiNeutralMatrix->Write();     
        } 

      // Back to F1 home directory
        F1->cd();

      // Create Comparisons directory
        if (!F1->GetDirectory("Comparisons")) {
            F1->mkdir("Comparisons");
        }
        F1->cd("Comparisons");

        if(particle == "K0sh"){
      //
      //
      /////////////////////////////////////////////////////////////////////////////////
      //Comparison Pt Spectra to Levy-Tsallis Fit
      //////////////////////////////////////////////////////////////////////////////////
      //
      //
            auto& hListK0sh = hCorrectedStorage;   // Pthistos
            auto& fListK0sh = fLevyTsallisStorage;         // levytsallis corresponding      

            TCanvas* cPt_Spectrum_Fit_Comparison_K0sh = new TCanvas("cPt_Spectrum_Fit_Comparison_K0sh", "Final Pt Spectra", 800, 800);

            TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1);
            TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.3);

            pad1->SetBottomMargin(0.02);
            pad2->SetTopMargin(0.05);
            pad2->SetBottomMargin(0.3);

            pad1->Draw();
            pad2->Draw();

      //First Pad: Draw Everything
            pad1->cd();
            pad1->SetLogy();

            TLegend* leg = new TLegend(0.65,0.5,0.88,0.88);

      // ----- Build hSum -----
            TH1F* hSum = (TH1F*)hListK0sh[0]->Clone("hSum_combined");
            hSum->Reset();
            hSum->SetDirectory(0);

            for(size_t i = 0; i < hListK0sh.size(); ++i){
                double cent_low  = (*centBins)[i];
                double cent_high = (*centBins)[i+1];

         // Draw histogram
                hListK0sh[i]->Draw(i==0 ? "PEL" : "PEL SAME");

         // Draw corresponding fit
                fListK0sh[i]->Draw("L SAME");

         //Build sum
                TH1F* tmp = (TH1F*)hListK0sh[i]->Clone();
                tmp->Scale((cent_high-cent_low)/100.);
                hSum->Add(tmp);
                delete tmp;

                TString label = Form("%s (%.0f-%.0f%%)",roman[i % roman.size()].Data(),cent_low, cent_high);

                leg->AddEntry(hListK0sh[i], label, "lep");
            }

      //Fit hSum
            double mParticle = (particle == "K0sh") ? 0.497611 : 1.115683;

            TF1 *levySum = new TF1("levySum","x*[0]*((([1]-1.)*([1]-2.))/(([1]*[2])*(([1]*[2])+([3]*([1]-2.))))) * pow(1.+((sqrt(pow([3],2.)+pow(x,2.))-[3])/([1]*[2])), -[1])",0.0, 14.0);

            levySum->FixParameter(3, mParticle);
            levySum->SetParameters(1.0, 7.0, 0.25);
            levySum->SetParLimits(0, 0.01, 5.0);
            levySum->SetParLimits(1, 3.0, 20.0);
            levySum->SetParLimits(2, 0.1, 0.6);

            hSum->Fit(levySum, "RN");

      // Draw sum
            hSum->SetMarkerStyle(kFullDiamond);
            hSum->SetMarkerColor(kBlack);
            hSum->SetLineColor(kBlack);
            hSum->Draw("PEL SAME");

            levySum->SetLineColor(kBlack);
            levySum->SetLineWidth(3);
            levySum->Draw("L SAME");

            leg->AddEntry(hSum, "Sum", "lep");
            leg->AddEntry(levySum, "Sum fit", "l");
            leg->Draw();

      // Second Pad: Ratio 
            pad2->cd();

            for(size_t i = 0; i < hListK0sh.size(); ++i){

                TH1F* ratio = (TH1F*)hListK0sh[i]->Clone(Form("ratio_%zu", i));
                ratio->SetDirectory(0);

                for(int b = 1; b <= ratio->GetNbinsX(); ++b){
                    double x = ratio->GetBinCenter(b);
                    double fval = fListK0sh[i]->Eval(x);

                    if(fval > 0){
                        ratio->SetBinContent(b, ratio->GetBinContent(b)/fval);
                        ratio->SetBinError(b, ratio->GetBinError(b)/fval);
                    } else {
                        ratio->SetBinContent(b, 0);
                        ratio->SetBinError(b, 0);
                    }
                }

                ratio->GetYaxis()->SetTitle("Ratio to Fit");
                ratio->GetYaxis()->SetRangeUser(0.5, 1.5);

                ratio->Draw(i==0 ? "PEL" : "PEL SAME");
            }

      // ----- Ratio for hSum -----
            TH1F* ratioSum = (TH1F*)hSum->Clone("ratioSum");

            for(int b = 1; b <= ratioSum->GetNbinsX(); ++b){
                double x = ratioSum->GetBinCenter(b);
                double fval = levySum->Eval(x);

                if(fval > 0){
                    ratioSum->SetBinContent(b, ratioSum->GetBinContent(b)/fval);
                    ratioSum->SetBinError(b, ratioSum->GetBinError(b)/fval);
                } else {
                    ratioSum->SetBinContent(b, 0);
                    ratioSum->SetBinError(b, 0);
                }
            }

            ratioSum->SetMarkerColor(kBlack);
            ratioSum->Draw("PEL SAME");
            cPt_Spectrum_Fit_Comparison_K0sh->Write();
      //
      //
      /////////////////////////////////////////////////////////////////////////////////
      //Comparison Centrality-Fits to hSum-Fit
      //////////////////////////////////////////////////////////////////////////////////
      //
      //

            TCanvas* cPt_Spectrum_INEL_Comparison_K0sh = new TCanvas("cPt_Spectrum_INEL_Comparison_K0sh", "Fit Ratios to Sum", 800, 800);

            TPad* pad1b = new TPad("pad1b","pad1b",0,0.3,1,1);
            TPad* pad2b = new TPad("pad2b","pad2b",0,0,1,0.3);

            pad1b->SetBottomMargin(0.02);
            pad2b->SetTopMargin(0.05);
            pad2b->SetBottomMargin(0.3);

            pad1b->Draw();
            pad2b->Draw();

            pad1b->cd();
            pad1b->SetLogy();

            TLegend* leg2 = new TLegend(0.65,0.5,0.88,0.88);

            for(size_t i = 0; i < hListK0sh.size(); ++i){

                double cent_low  = (*centBins)[i];
                double cent_high = (*centBins)[i+1];

                hListK0sh[i]->Draw(i==0 ? "PEL" : "PEL SAME");

                fListK0sh[i]->Draw("L SAME");

                TString label = Form("%s (%.0f-%.0f%%)",roman[i % roman.size()].Data(),cent_low, cent_high);

                leg2->AddEntry(hListK0sh[i], label, "lep");
            }

      // draw sum + fit again
            hSum->Draw("PEL SAME");
            levySum->Draw("L SAME");

            leg2->AddEntry(hSum, "Sum", "lep");
            leg2->AddEntry(levySum, "Sum fit", "l");
            leg2->Draw();

            pad2b->cd();

      // Use binning from first histogram as reference
            TH1F* frame = (TH1F*)hListK0sh[0]->Clone("frame_ratio_func");
            frame->Reset();
            frame->SetDirectory(0);

      //frame->SetMinimum(0.0);
      //frame->SetMaximum(5.0);
            frame->GetYaxis()->SetTitle("Levy-Tsallis to INEL>0");
      //frame->Draw("AXIS");

            for(size_t i = 0; i < fListK0sh.size(); ++i){

                TH1F* ratioFunc = (TH1F*)hListK0sh[0]->Clone(Form("ratioFunc_%zu", i));
                ratioFunc->Reset();
                ratioFunc->SetDirectory(0);

                for(size_t i = 0; i < fListK0sh.size(); ++i){

                    TH1F* ratioFunc = (TH1F*)hListK0sh  [0]->Clone(Form("ratioFunc_%zu", i));
                    ratioFunc->Reset();
                    ratioFunc->SetDirectory(0);

                    for(int b = 1; b <= ratioFunc->GetNbinsX(); ++b){
                        double x = ratioFunc->GetBinCenter(b);

                        double f_cent = fListK0sh[i]->Eval(x);
                        double f_sum  = levySum->Eval(x);

                        if(f_sum > 0){
                            ratioFunc->SetBinContent(b, f_cent / f_sum);
                        } else {
                            ratioFunc->SetBinContent(b, 0);
                        }
                    }

                    ratioFunc->SetLineColor(colors[i % colors.size()]);
                    ratioFunc->SetLineWidth(2);

            // THIS is the important part
                    ratioFunc->SetMinimum(0.0);
                    ratioFunc->SetMaximum(5.0);

                    ratioFunc->Draw(i==0 ? "HIST" : "HIST SAME");
                }

                ratioFunc->SetLineColor(colors[i % colors.size()]);
                ratioFunc->SetLineWidth(2);

                ratioFunc->Draw(i==0 ? "HIST" : "HIST SAME");
            }
            TLine* line = new TLine(frame->GetXaxis()->GetXmin(), 1.0,frame->GetXaxis()->GetXmax(),1.0);

            line->SetLineColor(kBlack);
            line->SetLineStyle(2);   // dashed
            line->SetLineWidth(2);
            line->Draw("SAME");
            cPt_Spectrum_INEL_Comparison_K0sh->Write();


      //
      //
      /////////////////////////////////////////////////////////////////////////////////
      //Comparison K0shs to Run2
      //////////////////////////////////////////////////////////////////////////////////
      //
      //

            TString legendNameK0sh = "K^{0}_{S}";

            TFile *fileGraphK0sh = TFile::Open(run2File);  
            if (!fileGraphK0sh) { std::cerr << "Error: Cannot open K0sh graph file.\n"; return 0; }

            TDirectory *dirK0sh = (TDirectory*)fileGraphK0sh->Get(run2Dir);  
            if (!dirK0sh) { std::cerr << "Error: K0sh directory not found.\n"; return 0; }

            TGraph *gExternalK0sh = (TGraph*)dirK0sh->Get(run2Graph);
            if (!gExternalK0sh) { std::cerr << "Error: K0sh TGraph not found.\n"; return 0; }

      // Prepare new canvas with ratio pad
            TCanvas *cCompareRun2K0sh = new TCanvas("cCompareRun2_K0sh", "Comparison K0shs to Run 2", 800, 800);
            TPad *pad1_Run2K0sh = new TPad("pad1_run2_K0sh", "Main pad", 0, 0.3, 1, 1.0);
            TPad *pad2_Run2K0sh = new TPad("pad2_run2_K0sh", "Ratio pad", 0, 0.05, 1, 0.3);
            pad1_Run2K0sh->Draw();
            pad2_Run2K0sh->Draw();

      // -- Draw top pad --
            pad1_Run2K0sh->cd();
            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);
            pad1_Run2K0sh->SetLogy(1);
            pad1_Run2K0sh->SetTickx(1);
            pad1_Run2K0sh->SetBottomMargin(0.01);

            if (!hK0shCorrectedCompare) { std::cerr << "Error: K0sh histogram for comparison not found.\n"; return 0; }

            hK0shCorrectedCompare->SetMarkerStyle(kFullCircle);
            hK0shCorrectedCompare->SetMarkerColor(kBlue+1);
            hK0shCorrectedCompare->SetMarkerSize(1);
            hK0shCorrectedCompare->GetYaxis()->SetTitle("1/N_{INEL>0} #times d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
            //hK0shCorrectedCompare->GetXaxis()->SetRangeUser(0.0,6.5);
            hK0shCorrectedCompare->GetXaxis()->SetLimits(0.0, 6.5);  
            hK0shCorrectedCompare->Draw("PE");

            gExternalK0sh->SetLineColor(kGreen+2);
            gExternalK0sh->SetMarkerColor(kGreen+2);
            gExternalK0sh->SetMarkerStyle(22);
            gExternalK0sh->Draw("LP SAME");

      // Legend
            TLegend* legendcompareK0sh = new TLegend(0.53, 0.70, 0.88, 0.88);
            legendcompareK0sh->SetTextSize(0.03);
            legendcompareK0sh->SetBorderSize(0);
            legendcompareK0sh->SetFillColor(kWhite);
            legendcompareK0sh->AddEntry(hK0shCorrectedCompare, legendNameK0sh + " Corrected p_{T} Spectrum", "lp");
            legendcompareK0sh->AddEntry(gExternalK0sh, "ALICE 2025 (#sqrt{#it{s}} = 5.02 TeV)", "pl");
            legendcompareK0sh->Draw();

      // Label
            TLatex *latexK0sh = new TLatex();
            latexK0sh->SetNDC();
            latexK0sh->SetTextSize(0.03);
            latexK0sh->DrawLatex(0.65, 0.65, "pp");
            TString textlatexK0sh = "#sqrt{#it{s}} = " + AO2DFiles.energy;
            latexK0sh->DrawLatex(0.65, 0.60, textlatexK0sh.Data());
            latexK0sh->DrawLatex(0.65, 0.55, "|y| < 0.5");

      // -- Bottom pad: ratio --
            pad2_Run2K0sh->cd();
            gPad->SetTopMargin(0.01);
            gPad->SetBottomMargin(0.3);
            gPad->SetGridy();
            int nK0sh = gExternalK0sh->GetN();
            double *xgK0sh = gExternalK0sh->GetX();
            double *ygK0sh = gExternalK0sh->GetY();

            TGraphErrors *gRatioK0sh = new TGraphErrors();

            for (int i = 0; i < nK0sh; ++i) {
                double x = xgK0sh[i];
                double y = ygK0sh[i];
                int bin = hK0shCorrectedCompare->FindBin(x);
                double yref = hK0shCorrectedCompare->GetBinContent(bin);

                if (yref <= 0 || y <= 0) continue;

                double sigma_y = 0;
                TGraphErrors *gExtErrors = dynamic_cast<TGraphErrors*>(gExternalK0sh);
                if (gExtErrors) {
                    sigma_y = gExtErrors->GetErrorY(i);
                }

                double sigma_ref = hK0shCorrectedCompare->GetBinError(bin);
                double ratio = yref / y;
                double ratio_err = ratio * sqrt((sigma_y > 0 ? (sigma_y / y) * (sigma_y / y) : 0) + (sigma_ref / yref) * (sigma_ref / yref));
                gRatioK0sh->SetPoint(i, x, ratio);
                gRatioK0sh->SetPointError(i, 0, ratio_err);
            }

            gRatioK0sh->SetTitle("");
            gRatioK0sh->GetYaxis()->SetTitle("Run 3 / Run 2");
            gRatioK0sh->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            gRatioK0sh->SetMarkerStyle(33);
            gRatioK0sh->SetMarkerSize(1);
            gRatioK0sh->SetMarkerColor(kCyan+2);
            gRatioK0sh->SetLineColor(kCyan+2);
            gRatioK0sh->GetYaxis()->SetTitleSize(0.08);
            gRatioK0sh->GetYaxis()->SetTitleOffset(0.5);
            gRatioK0sh->GetYaxis()->SetLabelSize(0.08);
            gRatioK0sh->GetXaxis()->SetLabelSize(0.08);
            gRatioK0sh->GetXaxis()->SetTitleSize(0.10);
            gRatioK0sh->GetXaxis()->SetTitleOffset(1.0);
            gRatioK0sh->GetXaxis()->SetRangeUser(0.0,8.0);
            gRatioK0sh->Draw("APE");

            TLine *lineK0sh = new TLine(hK0shCorrectedCompare->GetXaxis()->GetXmin(), 1.0, hK0shCorrectedCompare->GetXaxis()->GetXmax(), 1.0);
            lineK0sh->SetLineStyle(2);
            lineK0sh->Draw("SAME");

            F1->cd();
            gDirectory->cd("Comparisons");
            cCompareRun2K0sh->Update();
            cCompareRun2K0sh->Write();
        } //K0sh Comparisons
      /////////////////////////////////////////////////////////////////
      //Lambda+AntiLambda Sum Comparisons Done Together 
        if (lambdaProcessed && antiLambdaProcessed) { //First the Lambdas must be analysed and then the AntiLambdas
         

        // ============================================================
        // Lambda + AntiLambda Pt Spectra + Levy-Tsallis Fits
        // ============================================================

        // Storage: build combined spectra
        std::vector<TH1F*> hLambdaPlusAntiLambdaStorage;
        std::vector<TF1*>  fLevyTsallisLambdaPlusAntiLambdaStorage;

        // Mass
        double mParticle = 1.115683;

        // ============================================================
        // 1. BUILD (Lambda + AntiLambda) PER CENTRALITY BIN
        // ============================================================

        for (size_t i = 0; i < hPtSpectrumCorrectedLambdaStorage.size(); ++i) {

            TH1F* hLambda = hPtSpectrumCorrectedLambdaStorage[i];
            TH1F* hAntiLambda = hPtSpectrumCorrectedAntiLambdaStorage[i];

            TH1F* hSum = (TH1F*)hLambda->Clone(Form("hLambdaPlusAntiLambda_cent_%zu", i));
            hSum->SetDirectory(0);
            hSum->Reset();

            hSum->Add(hLambda);
            hSum->Add(hAntiLambda);

            hLambdaPlusAntiLambdaStorage.push_back(hSum);
        }

        // ============================================================
        // 2. CANVAS: SPECTRA + FITS + INEL SUM + RATIO TO FIT
        // ============================================================

        TCanvas* cLambdaAntiLambdaSpecFit = new TCanvas("cLambdaAntiLambdaSpecFit", "Lambda+AntiLambda Pt Spectra",800, 800);

        TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1);
        TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.3);

        pad1->SetBottomMargin(0.02);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.3);

        pad1->Draw();
        pad2->Draw();

        // ------------------------
        // Pad 1: spectra + fits
        // ------------------------
        pad1->cd();
        pad1->SetLogy();

        TLegend* leg = new TLegend(0.60,0.55,0.88,0.88);
        leg->SetBorderSize(0);

        // INEL sum
        TH1F* hSum = (TH1F*)hLambdaPlusAntiLambdaStorage[0]->Clone("hSum_Lambda");
        hSum->Reset();
        hSum->SetDirectory(0);

        // Levy function (same form as yours)
        TF1* levySum = new TF1("levySum","x*[0]*((([1]-1.)*([1]-2.))/(([1]*[2])*(([1]*[2])+([3]*([1]-2.)))))""* pow(1.+((sqrt(pow([3],2.)+pow(x,2.))-[3])/([1]*[2])), -[1])",0.0, 10.0);

        levySum->FixParameter(3, mParticle);
        levySum->SetParameters(1.0, 7.0, 0.25);
        levySum->SetParLimits(0, 0.01, 5.0);
        levySum->SetParLimits(1, 3.0, 20.0);
        levySum->SetParLimits(2, 0.1, 0.6);

        // loop centrality bins
        for (size_t i = 0; i < hLambdaPlusAntiLambdaStorage.size(); ++i) {

            double cent_low  = (*nCentBinsLambdas)[i];
            double cent_high = (*nCentBinsLambdas)[i+1];

            TH1F* h = hLambdaPlusAntiLambdaStorage[i];

            // draw spectrum
            h->SetMarkerStyle(20);
            h->SetMarkerColor(colors[i % colors.size()]);
            h->SetLineColor(colors[i % colors.size()]);

            h->Draw(i==0 ? "E" : "E SAME");

            // fit
            TF1* f = (TF1*)levySum->Clone(Form("levy_cent_%zu", i));
            f->SetLineColor(colors[i % colors.size()]);
            f->SetLineWidth(2);

            h->Fit(f, "RQ0"); // quiet fit

            fLevyTsallisLambdaPlusAntiLambdaStorage.push_back(f);

            f->Draw("L SAME");

            leg->AddEntry(h, Form("%s (%.0f-%.0f%%)",roman[i].Data(),cent_low, cent_high),"lep");

            // build INEL sum (weighted like your K0s code)
            double weight = (cent_high - cent_low) / 100.0;

            TH1F* tmp = (TH1F*)h->Clone();
            tmp->Scale(weight);
            hSum->Add(tmp);
            delete tmp;
        }

        // fit INEL sum
        hSum->SetMarkerStyle(kFullDiamond);
        hSum->SetMarkerColor(kBlack);
        hSum->SetLineColor(kBlack);
        hSum->Draw("E SAME");

        hSum->Fit(levySum, "RQ0");

        levySum->SetLineColor(kBlack);
        levySum->SetLineWidth(3);
        levySum->Draw("L SAME");

        leg->AddEntry(hSum, "INEL sum", "lep");
        leg->AddEntry(levySum, "INEL fit", "l");

        leg->Draw();

        // ------------------------
        // Pad 2: ratio to fit
        // ------------------------
        pad2->cd();

        for (size_t i = 0; i < hLambdaPlusAntiLambdaStorage.size(); ++i) {

            TH1F* ratio = (TH1F*)hLambdaPlusAntiLambdaStorage[i]->Clone(
                Form("ratio_%zu", i));

            for (int b = 1; b <= ratio->GetNbinsX(); ++b) {

                double x = ratio->GetBinCenter(b);
                double fval = fLevyTsallisLambdaPlusAntiLambdaStorage[i]->Eval(x);

                if (fval > 0) {
                    ratio->SetBinContent(b, ratio->GetBinContent(b)/fval);
                    ratio->SetBinError(b, ratio->GetBinError(b)/fval);
                } else {
                    ratio->SetBinContent(b, 0);
                    ratio->SetBinError(b, 0);
                }
            }

            ratio->SetMarkerStyle(20);
            ratio->SetMarkerColor(colors[i % colors.size()]);
            ratio->SetLineColor(colors[i % colors.size()]);

            ratio->GetYaxis()->SetTitle("Ratio to Fit");
            ratio->GetYaxis()->SetRangeUser(0.5, 1.5);

            ratio->Draw(i==0 ? "E" : "E SAME");
        }

        // INEL ratio
        TH1F* ratioSum = (TH1F*)hSum->Clone("ratioSum");

        for (int b = 1; b <= ratioSum->GetNbinsX(); ++b) {

            double x = ratioSum->GetBinCenter(b);
            double fval = levySum->Eval(x);

            if (fval > 0) {
                ratioSum->SetBinContent(b, ratioSum->GetBinContent(b)/fval);
                ratioSum->SetBinError(b, ratioSum->GetBinError(b)/fval);
            } else {
                ratioSum->SetBinContent(b, 0);
                ratioSum->SetBinError(b, 0);
            }
        }

        ratioSum->SetMarkerColor(kBlack);
        ratioSum->Draw("E SAME");

        cLambdaAntiLambdaSpecFit->Write();


        // ============================================================
        // 3. SECOND CANVAS: CENTRALITY FITS / INEL FIT
        // ============================================================

        TCanvas* cLambdaAntiLambdaFitComparison = new TCanvas("cLambdaAntiLambdaFitComparison", "Lambda+AntiLambda Fit Ratios", 800, 800);

        TPad* pad1b = new TPad("pad1b","pad1b",0,0.3,1,1);
        TPad* pad2b = new TPad("pad2b","pad2b",0,0,1,0.3);

        pad1b->SetBottomMargin(0.02);
        pad2b->SetTopMargin(0.05);
        pad2b->SetBottomMargin(0.3);

        pad1b->Draw();
        pad2b->Draw();

        // ------------------------
        // Pad 1b
        // ------------------------
        pad1b->cd();
        pad1b->SetLogy();

        TLegend* leg2 = new TLegend(0.60,0.55,0.88,0.88);

        for (size_t i = 0; i < hLambdaPlusAntiLambdaStorage.size(); ++i) {

            TH1F* h = hLambdaPlusAntiLambdaStorage[i];

            h->SetMarkerColor(colors[i % colors.size()]);
            h->SetLineColor(colors[i % colors.size()]);

            h->Draw(i==0 ? "E" : "E SAME");

            fLevyTsallisLambdaPlusAntiLambdaStorage[i]->Draw("L SAME");

            leg2->AddEntry(h, Form("%s", roman[i].Data()), "lep");
        }

        // INEL + fit
        hSum->Draw("E SAME");
        levySum->Draw("L SAME");

        leg2->AddEntry(hSum, "INEL sum", "lep");
        leg2->AddEntry(levySum, "INEL fit", "l");

        leg2->Draw();

        // ------------------------
        // Pad 2b: fit ratios
        // ------------------------
        pad2b->cd();

        for (size_t i = 0; i < hLambdaPlusAntiLambdaStorage.size(); ++i) {

            TH1F* ratioFunc = (TH1F*)hLambdaPlusAntiLambdaStorage[i]->Clone(Form("ratioFunc_%zu", i));

            ratioFunc->Reset();

            for (int b = 1; b <= ratioFunc->GetNbinsX(); ++b) {

                double x = ratioFunc->GetBinCenter(b);

                double f_cent = fLevyTsallisLambdaPlusAntiLambdaStorage[i]->Eval(x);
                double f_sum  = levySum->Eval(x);

                if (f_sum > 0)
                    ratioFunc->SetBinContent(b, f_cent / f_sum);
                else
                    ratioFunc->SetBinContent(b, 0);
            }

            ratioFunc->SetLineColor(colors[i % colors.size()]);
            ratioFunc->SetLineWidth(2);

            ratioFunc->GetYaxis()->SetRangeUser(0,5);

            ratioFunc->Draw(i==0 ? "HIST" : "HIST SAME");
        }

        // unity line
        TLine* line = new TLine(0.0,1.0,14.0,1.0);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("SAME");

        cLambdaAntiLambdaFitComparison->Write();     
         
         //
         //
         /////////////////////////////////////////////////////////////////////////////////
         //Comparison Lambda+LambdaBar to Run2
         //////////////////////////////////////////////////////////////////////////////////
         //
         // 
            TString legendNameLambda = "#Lambda + #bar{#Lambda}";

            TFile *fileGraphLambda = TFile::Open(run2File);  
            if (!fileGraphLambda) { std::cerr << "Error: Cannot open lambda graph file.\n"; return 0; }

            TDirectory *dirLambda = (TDirectory*)fileGraphLambda->Get(run2Dir);  
            if (!dirLambda) { std::cerr << "Error: Lambda directory not found.\n"; return 0; }

            TGraph *gExternalLambda = (TGraph*)dirLambda->Get(run2Graph);
            if (!gExternalLambda) { std::cerr << "Error: Lambda TGraph not found.\n"; return 0; }

         // Prepare new canvas with ratio pad
            TCanvas *cCompareRun2Lambda = new TCanvas("cCompareRun2_Lambda+AntiLambda", "Comparison Lambda+LambdaBar to Run 2", 800, 800);
            TPad *pad1_Run2Lambda = new TPad("pad1_run2_lambda", "Main pad", 0, 0.3, 1, 1.0);
            TPad *pad2_Run2Lambda = new TPad("pad2_run2_lambda", "Ratio pad", 0, 0.05, 1, 0.3);
            pad1_Run2Lambda->Draw();
            pad2_Run2Lambda->Draw();

         // -- Draw top pad --
            pad1_Run2Lambda->cd();
            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);
            pad1_Run2Lambda->SetLogy(1);
            pad1_Run2Lambda->SetTickx(1);
            pad1_Run2Lambda->SetBottomMargin(0.01);

            if (!hLambdaPlusAntiLambdaForComparison) { std::cerr << "Error: Lambda histogram for comparison not found.\n"; return 0; }

            hLambdaPlusAntiLambdaForComparison->SetMarkerStyle(kFullCircle);
            hLambdaPlusAntiLambdaForComparison->SetMarkerColor(kBlue+1);
            hLambdaPlusAntiLambdaForComparison->SetMarkerSize(1);
            hLambdaPlusAntiLambdaForComparison->GetYaxis()->SetTitle("1/N_{INEL>0} #times d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
            //hLambdaPlusAntiLambdaForComparison->GetXaxis()->SetRangeUser(0.0,8.0);
            hLambdaPlusAntiLambdaForComparison->GetXaxis()->SetLimits(0.0, 6.5); 
            hLambdaPlusAntiLambdaForComparison->Draw("PE");

            gExternalLambda->SetLineColor(kGreen+2);
            gExternalLambda->SetMarkerColor(kGreen+2);
            gExternalLambda->SetMarkerStyle(22);
            gExternalLambda->Draw("LP SAME");

         // Legend
            TLegend* legendcompareLambda = new TLegend(0.53, 0.70, 0.88, 0.88);
            legendcompareLambda->SetTextSize(0.03);
            legendcompareLambda->SetBorderSize(0);
            legendcompareLambda->SetFillColor(kWhite);
            legendcompareLambda->AddEntry(hLambdaPlusAntiLambdaForComparison, legendNameLambda + " Corrected p_{T} Spectrum", "lp");
            legendcompareLambda->AddEntry(gExternalLambda, "ALICE 2025 (#sqrt{#it{s}} = 5.02 TeV)", "pl");
            legendcompareLambda->Draw();

         // Label
            TLatex *latexLambda = new TLatex();
            latexLambda->SetNDC();
            latexLambda->SetTextSize(0.03);
            latexLambda->DrawLatex(0.65, 0.65, "pp");
            TString textLambda = "#sqrt{#it{s}} = " + AO2DFiles.energy;
            latexLambda->DrawLatex(0.65, 0.60, textLambda.Data());
            latexLambda->DrawLatex(0.65, 0.55, "|y| < 0.5");

         // -- Bottom pad: ratio --
            pad2_Run2Lambda->cd();
            gPad->SetTopMargin(0.01);
            gPad->SetBottomMargin(0.3);
            gPad->SetGridy();
            int nLambda = gExternalLambda->GetN();
            double *xgLambda = gExternalLambda->GetX();
            double *ygLambda = gExternalLambda->GetY();

            TGraphErrors *gRatioLambda = new TGraphErrors();

            for (int i = 0; i < nLambda; ++i) {
                double x = xgLambda[i];
                double y = ygLambda[i];
                int bin = hLambdaPlusAntiLambdaForComparison->FindBin(x);
                double yref = hLambdaPlusAntiLambdaForComparison->GetBinContent(bin);

                if (yref <= 0 || y <= 0) continue;

                double sigma_y = 0;
                TGraphErrors *gExtErrors = dynamic_cast<TGraphErrors*>(gExternalLambda);
                if (gExtErrors) {
                    sigma_y = gExtErrors->GetErrorY(i);
                }

                double sigma_ref = hLambdaPlusAntiLambdaForComparison->GetBinError(bin);
                double ratio = yref / y;
                double ratio_err = ratio * sqrt((sigma_y > 0 ? (sigma_y / y) * (sigma_y / y) : 0) + (sigma_ref / yref) * (sigma_ref / yref));
                gRatioLambda->SetPoint(i, x, ratio);
                gRatioLambda->SetPointError(i, 0, ratio_err);
            }

            gRatioLambda->SetTitle("");
            gRatioLambda->GetYaxis()->SetTitle("Run 3 / Run 2");
            gRatioLambda->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            gRatioLambda->SetMarkerStyle(33);
            gRatioLambda->SetMarkerSize(1);
            gRatioLambda->SetMarkerColor(kCyan+2);
            gRatioLambda->SetLineColor(kCyan+2);
            gRatioLambda->GetYaxis()->SetTitleSize(0.08);
            gRatioLambda->GetYaxis()->SetTitleOffset(0.5);
            gRatioLambda->GetYaxis()->SetLabelSize(0.08);
            gRatioLambda->GetXaxis()->SetLabelSize(0.08);
            gRatioLambda->GetXaxis()->SetTitleSize(0.10);
            gRatioLambda->GetXaxis()->SetTitleOffset(1.0);
            gRatioLambda->GetXaxis()->SetRangeUser(0.0,8.0);
            gRatioLambda->Draw("APE");

            TLine *lineLambda = new TLine(hLambdaPlusAntiLambdaForComparison->GetXaxis()->GetXmin(), 1.0, hLambdaPlusAntiLambdaForComparison->GetXaxis()->GetXmax(), 1.0);
            lineLambda->SetLineStyle(2);
            lineLambda->Draw("SAME");

            cCompareRun2Lambda->Update();
            F1->cd();
            gDirectory->cd("Comparisons");
            cCompareRun2Lambda->Write();
      //
      //
      /////////////////////////////////////////////////////////////////////////////////
      //Comparison Lambda to LambdaBar for Differences 
      //////////////////////////////////////////////////////////////////////////////////
      //
         //
         /////////////////////////////////////////////////////////////////////////////////
         //Lambda and AntiLambda Comparisons (using Function DrawLambdaAntiLambdaComparison)
         //////////////////////////////////////////////////////////////////////////////////
         //
         //
            if (hEffLambdaStorage.empty()) {
                std::cout << "ERROR: Missing Lambda efficiency histograms\n";
                return 0;
            }

            if (hEffAntiLambdaStorage.empty()) {
                std::cout << "ERROR: Missing AntiLambda efficiency histograms\n";
                return 0;
            }

            if (hEffLambdaStorage.size() != hEffAntiLambdaStorage.size()) {
                std::cout << "ERROR: Lambda/AntiLambda size mismatch: "<< hEffLambdaStorage.size() << " vs "<< hEffAntiLambdaStorage.size() << "\n";
                return 0;
            }

            if (hEffLambdaStorage.size() != centBins->size() - 1) {
                std::cout << "ERROR: Histogram size does not match centBins\n";
                return 0;
            }
               
                  //Efficiency Comparison
            TCanvas* cPtCompare = DrawLambdaAntiLambdaComparison(hPtSpectrumCorrectedLambdaStorage,hPtSpectrumCorrectedAntiLambdaStorage,roman,colors,*centBins,"PtSpesctrum_Lambda_and_AntiLambda",false/*logy*/,false/*fixYAxis*/);
         //Efficiency Comparison
            TCanvas* cEffCompare = DrawLambdaAntiLambdaComparison(hEffLambdaStorage,hEffAntiLambdaStorage,roman,colors,*centBins,"Efficiency_Lambda_and_AntiLambda",false/*logy*/,false/*fixYAxis*/);
         //Feeddown Fraction Comparison
            TCanvas* cFracCompare = DrawLambdaAntiLambdaComparison(hFeeddownFracLambdaStorage,hFeeddownFracAntiLambdaStorage,roman,colors,*centBins, "Feeddown_Lambda_and_AntiLambda",false/*logy*/,false/*fixYAxis*/);
                     
            F1->cd();
            gDirectory->cd("Comparisons");
            if (cPtCompare) cPtCompare->Write();
            if (cFracCompare) cFracCompare->Write();
            if (cEffCompare)  cEffCompare->Write();



        } //end of Lambda+AntiLambda Comparisons

        //
        //
        /////////////////////////////////////////////////////////////////////////////////
        //Baryon to Meson Ratio
        //////////////////////////////////////////////////////////////////////////////////
        //
        //

        if (k0shProcessed && lambdaProcessed && antiLambdaProcessed && (*nCentBinsKaons == *nCentBinsLambdas)) {

            std::vector<TH1D*> ratioList;
            TH1D* hRatioSum = nullptr;
            //Function that Computes Lambda+AntiLambda / 2K0s
            ComputeBaryonToMesonRatio(F1,"Corrected Pt Spectrum",*nCentBinsKaons,ratioList,hRatioSum);

            TCanvas* cBaryonMesonRatio = new TCanvas("BaryonMesonRatio","(Lambda+AntiLambda)/2K0s",900,800);
            gStyle->SetOptStat(0);

            TPad* pad1 =new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
            TPad* pad2 =new TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.30);
            pad1->SetBottomMargin(0.02);
            pad2->SetTopMargin(0.02);
            pad2->SetBottomMargin(0.30);
            pad1->Draw();
            pad2->Draw();

            //Top Pad
            pad1->cd();
            TLegend* legBaryonMesonRatio =new TLegend(0.55, 0.60, 0.88, 0.88);
            legBaryonMesonRatio->SetBorderSize(0);
            bool first = true;

            for (size_t i = 0; i < ratioList.size(); ++i) {

                TH1D* hBaryonMesonRatio = ratioList[i];

                hBaryonMesonRatio->SetLineColor(colors[i % colors.size()]);
                hBaryonMesonRatio->SetMarkerColor(colors[i % colors.size()]);
                hBaryonMesonRatio->SetMarkerStyle(20 + (i % 10));
                hBaryonMesonRatio->SetTitle("");
                hBaryonMesonRatio->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda})/2K^{0}_{s}");
                hBaryonMesonRatio->GetYaxis()->SetTitleSize(0.05);
                hBaryonMesonRatio->GetYaxis()->SetLabelSize(0.045);
                hBaryonMesonRatio->GetXaxis()->SetLabelSize(0);
                hBaryonMesonRatio->GetYaxis()->SetRangeUser(0.0, 1.0);

                if (first) {

                    hBaryonMesonRatio->Draw("E");
                    first = false;
                }
                else {
                    hBaryonMesonRatio->Draw("E SAME");
                }
                legBaryonMesonRatio->AddEntry(hBaryonMesonRatio, Form("%s (%.0f-%.0f%%)",roman[i].Data(),(*nCentBinsKaons)[i], (*nCentBinsKaons)[i+1]),"lep");
            }

            //The Sum
            if (hRatioSum) {

                hRatioSum->SetLineColor(kBlack);
                hRatioSum->SetMarkerColor(kBlack);
                hRatioSum->SetMarkerStyle(24);
                hRatioSum->SetLineWidth(3);
                hRatioSum->Draw("E SAME");
                legBaryonMesonRatio->AddEntry(hRatioSum, "Weighted Sum", "lep");
            }
            legBaryonMesonRatio->Draw();

            // Bottom Pad
            pad2->cd();
            bool firstComparison = true;

            for (size_t i = 0; i < ratioList.size(); ++i) {

                TH1D* hComparison = (TH1D*)ratioList[i]->Clone(Form("hComparison_%zu", i));
                 // Compute ratio to sum
                for (int b = 1; b <= hComparison->GetNbinsX(); ++b) {

                    double val    = ratioList[i]->GetBinContent(b);
                    double valErr = ratioList[i]->GetBinError(b);
                    double avg    = hRatioSum->GetBinContent(b);
                    double avgErr = hRatioSum->GetBinError(b);

                    if (val <= 0 || avg <= 0) {
                        hComparison->SetBinContent(b, 0);
                        hComparison->SetBinError(b, 0);
                        continue;
                    }

                    double ratio = val / avg;
                    double ratioErr = ratio * std::sqrt((valErr * valErr) / (val * val) +(avgErr * avgErr) / (avg * avg));
                    hComparison->SetBinContent(b, ratio);
                    hComparison->SetBinError(b, ratioErr);
                }

                hComparison->SetLineColor(colors[i % colors.size()]);
                hComparison->SetMarkerColor(colors[i % colors.size()]);
                hComparison->SetMarkerStyle(20 + (i % 10));
                hComparison->SetTitle("");
                hComparison->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hComparison->GetYaxis()->SetTitle("Cent./Avg.");
                hComparison->GetYaxis()->CenterTitle();
                hComparison->GetXaxis()->SetTitleSize(0.12);
                hComparison->GetXaxis()->SetLabelSize(0.10);
                hComparison->GetYaxis()->SetTitleSize(0.10);
                hComparison->GetYaxis()->SetLabelSize(0.08);
                hComparison->GetYaxis()->SetNdivisions(505);
                hComparison->GetYaxis()->SetRangeUser(0.4, 1.6);

                if (firstComparison) {
                    hComparison->Draw("E");
                    firstComparison = false;
                }
                else {
                    hComparison->Draw("E SAME");
                }
            }
            // Unity Line

            TLine* unity =new TLine(0.5, 1.0, 6.5, 1.0);
            unity->SetLineStyle(2);
            unity->SetLineWidth(2);
            unity->Draw("SAME");

            cBaryonMesonRatio->Update();
            if (cBaryonMesonRatio)  cBaryonMesonRatio->Write();

        //
        // ============================================================
        // Compare Baryon-to-Meson Ratio Sum to Run 2
        // ============================================================
        //

            TCanvas* cBaryonMesonRun2Comparison = new TCanvas("cBaryonMesonRun2", "Baryon-to-Meson Ratio: Run 3 vs Run 2",900, 800);

            TFile *fileGraphBaryontoMeson = TFile::Open(run2File);  
            if (!fileGraphBaryontoMeson) { std::cerr << "Error: Cannot open baryon-to-meson graph file.\n"; return 0; }

            TString run2DirBaryonMeson = "Figure 3";
            TDirectory *dirBaryonMeson = (TDirectory*)fileGraphBaryontoMeson->Get(run2DirBaryonMeson);  
            if (!dirBaryonMeson) { std::cerr << "Error: Baryon-to-Meson directory not found.\n"; return 0; }

            TString run2GraphBaryonMeson = "Graph1D_y10";
            TGraph *gExternalBaryonMeson = (TGraph*)dirBaryonMeson->Get(run2GraphBaryonMeson);
            if (!gExternalBaryonMeson) { std::cerr << "Error: Baryon-to-Meson TGraph not found.\n"; return 0; }

            gStyle->SetOptStat(0);

            TPad* pad1Run2 = new TPad("pad1BaryonMesonRun2","Main pad",0.0, 0.30, 1.0, 1.0);
            TPad* pad2Run2 = new TPad("pad2BaryonMesonRun2","Ratio pad",0.0, 0.00, 1.0, 0.30);

            pad1Run2->SetBottomMargin(0.02);
            pad2Run2->SetTopMargin(0.02);
            pad2Run2->SetBottomMargin(0.30);

            pad1Run2->Draw();
            pad2Run2->Draw();


            // ========================================================
            // Top pad
            // ========================================================

            pad1Run2->cd();

            hRatioSum->SetTitle("");
            hRatioSum->SetLineColor(kGreen+1);
            hRatioSum->SetMarkerColor(kGreen+1);
            hRatioSum->SetMarkerStyle(33);
            hRatioSum->SetMarkerSize(1.2);
            hRatioSum->SetLineWidth(2);

            hRatioSum->GetYaxis()->SetTitle( "(#Lambda + #bar{#Lambda})/2K^{0}_{S}");
            hRatioSum->GetYaxis()->SetTitleSize(0.05);
            hRatioSum->GetYaxis()->SetLabelSize(0.045);
            hRatioSum->GetXaxis()->SetLabelSize(0);
            hRatioSum->GetYaxis()->SetRangeUser(0.0, 0.8);
            hRatioSum->Draw("E");


            // Run 2
            gExternalBaryonMeson->SetLineColor(kBlue+1);
            gExternalBaryonMeson->SetMarkerColor(kBlue+1);
            gExternalBaryonMeson->SetMarkerStyle(23);
            gExternalBaryonMeson->SetMarkerSize(1.2);

            gExternalBaryonMeson->Draw("P SAME");


            // ========================================================
            // Legend
            // ========================================================

            TLegend* legBaryonMesonRun2 = new TLegend(0.55, 0.68, 0.88, 0.86);

            legBaryonMesonRun2->SetBorderSize(0);
            legBaryonMesonRun2->AddEntry(hRatioSum,"(Baryon-to-Meson ratio, #sqrt{#it{s}} = 5.36 TeV)","lep");
            legBaryonMesonRun2->AddEntry(gExternalBaryonMeson,"ALICE 2025 (#sqrt{#it{s}} = 5.02 TeV)","lp");
            legBaryonMesonRun2->Draw();

            // ========================================================
            // Bottom pad: Run 3 / Run 2
            // ========================================================

            pad2Run2->cd();

            gPad->SetGridy();

            int nRun2 = gExternalBaryonMeson->GetN();

            double* xRun2 = gExternalBaryonMeson->GetX();
            double* yRun2 = gExternalBaryonMeson->GetY();

            TGraphErrors* gRatioBaryonMesonRun2Comparison = new TGraphErrors();

            for (int i = 0; i < nRun2; ++i) {

                double x = xRun2[i];
                double y = yRun2[i];

                int bin = hRatioSum->FindBin(x);

                double yRun3 = hRatioSum->GetBinContent(bin);

                if (y <= 0 || yRun3 <= 0)
                    continue;

                double sigmaRun2 = 0.0;

                TGraphErrors* gRun2Errors = dynamic_cast<TGraphErrors*>(gExternalBaryonMeson);

                if (gRun2Errors)
                    sigmaRun2 = gRun2Errors->GetErrorY(i);

                double sigmaRun3 = hRatioSum->GetBinError(bin);

                double ratio = yRun3 / y;

                double ratioErr = ratio * std::sqrt((sigmaRun3 / yRun3) * (sigmaRun3 / yRun3) +(sigmaRun2 / y) * (sigmaRun2 / y));

                int point = gRatioBaryonMesonRun2Comparison->GetN();

                gRatioBaryonMesonRun2Comparison->SetPoint(point, x, ratio);
                gRatioBaryonMesonRun2Comparison->SetPointError(point, 0, ratioErr);
            }

            gRatioBaryonMesonRun2Comparison->SetTitle("");

            gRatioBaryonMesonRun2Comparison->SetMarkerStyle(33);
            gRatioBaryonMesonRun2Comparison->SetMarkerSize(1.0);
            gRatioBaryonMesonRun2Comparison->SetMarkerColor(kCyan+1);
            gRatioBaryonMesonRun2Comparison->SetLineColor(kCyan+2);

            gRatioBaryonMesonRun2Comparison->GetYaxis()->SetTitle("Run 3 / Run 2");
            gRatioBaryonMesonRun2Comparison->GetXaxis()->SetTitle("p_{T} (GeV/c)");

            gRatioBaryonMesonRun2Comparison->GetYaxis()->SetTitleSize(0.10);
            gRatioBaryonMesonRun2Comparison->GetYaxis()->SetTitleOffset(0.5);
            gRatioBaryonMesonRun2Comparison->GetYaxis()->SetLabelSize(0.08);
            gRatioBaryonMesonRun2Comparison->GetXaxis()->SetLabelSize(0.08);
            gRatioBaryonMesonRun2Comparison->GetXaxis()->SetTitleSize(0.10);
            gRatioBaryonMesonRun2Comparison->GetXaxis()->SetTitleOffset(1.0);
            gRatioBaryonMesonRun2Comparison->GetYaxis()->SetRangeUser(0.5, 1.5);
            gRatioBaryonMesonRun2Comparison->Draw("APE");

            // Unity line

            TLine* unityRun2 = new TLine(0.5, 1.0, 6.5, 1.0);

            unityRun2->SetLineStyle(2);
            unityRun2->SetLineWidth(2);
            unityRun2->Draw("SAME");

            cBaryonMesonRun2Comparison->Update();

            F1->cd();
            gDirectory->cd("Comparisons");

            cBaryonMesonRun2Comparison->Write();
        }
    }//end of particle Loop. WATCH OUT FOR THE BRACKETS!
    F1->Close();
    return 0;
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

/*
void CalculateSignalSimple(TString& histname,std::vector<TString>& ptbins, std::vector<double_t>& signal, TFile* file, centdir, int mcflag, std::vector<double_t>& mc_sig_low, std::vector<double_t>& mc_sig_high, std::vector<double_t>& peakpt, std::vector<double_t>& widthpt, int nhistos){
    Double_t M_low = 0.00;
    Double_t M_high = 100.0;

   // initialize and convert tokenized strings into vector of doubles
    std::vector<std::string> histptbinsnames(nhistos + 1);                                        
    for (int i = 0; i < nhistos+1; i++) {            //Histos won't accept "." character so converting it to "_"
        std::string histptbin = ptbins[i];          // getting the value of the bin edge
        size_t pos = histptbin.find(".");          // finding the "." character
        histptbin[pos] = '_';                      // changing the "." character of thestring-value to a "_"
        histptbinsnames[i] = histptbin;             // filling bin edges list
    }
    for(int i=0;i<nhistos;i++){
       
        TH2F* h1=(TH2F*)gDirectory->Get(fmt::format("{}_{}_to_{}", histname, histptbinsnames[i],histptbinsnames[i+1]).data());
        if(h1 == nullptr) {
            std::cout << "h was not read Simple" << std::endl; 
            exit(1);
        }
      //for MC we count the bins based on the same Pt bin integration range
        Double_t sig_low = 0.0;
        Double_t sig_high = 100.0;
       //
       // Signal+Backg from bin counts
       //
        Int_t i_low = h1->FindBin(sig_low);
        Int_t i_high = h1->FindBin(sig_high);
        cout <<"low:" << sig_low << " high:" << sig_high << std::endl;
        Double_t sigplusbckg = 0.0; //h1->Integral(i_low,i_high);
        std::cout << "sig+bckg:" << sigplusbckg << std::endl;
        double_t realbckg=0.;

        signal[i]=sigplusbckg-realbckg;
        cout<<"Signal:"<<signal[i]<<endl<<"Number of Plot:"<<i+1<<endl;
    }
    return;
}
*/

/*
/////////////////////////////////////////////////////////////////////////////////
//Beta Test
//////////////////////////////////////////////////////////////////////////////////
//
//
TCanvas *cBetaTest = new TCanvas("cBetaTest", "Beta test", 800, 800);
TPad *padBeta1 = new TPad("padBeta1", "Main pad", 0, 0.3, 1, 1.0);
TPad *padBeta2 = new TPad("padBeta2", "Ratio pad", 0, 0.05, 1, 0.3);

padBeta1->Draw();
padBeta2->Draw();

    // -- Top pad: draw corrected and generated spectra --
padBeta1->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
padBeta1->SetLogy(1);
padBeta1->SetBottomMargin(0.01);

    // Clone the original histogram to keep original unchanged and multiply by nCollisions so they are comparible
TH1D *hGeneratedK0shBeta = (TH1D*)hGeneratedK0shB->Clone("hGeneratedK0shBeta");
    //hGeneratedK0shBeta->Scale(1.0/nCollisionsGenBeta);

TH1D *hCorrectedK0shBeta = (TH1D*)hCorrectedK0sh->Clone("hCorrectedK0shBeta");
hCorrectedK0shBeta->Scale(nCollisionsGenBeta);
hCorrectedK0shBeta->SetMarkerStyle(kFullCircle);
hCorrectedK0shBeta->SetMarkerColor(kBlue+1);
hCorrectedK0shBeta->SetLineColor(kBlue+1);
hCorrectedK0shBeta->SetMarkerSize(1);
hCorrectedK0shBeta->GetYaxis()->SetTitle("1/N_{INEL>0} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
hCorrectedK0shBeta->GetXaxis()->SetRangeUser(0.0,8.0);
hCorrectedK0shBeta->Draw("PE");

hGeneratedK0shBeta->SetMarkerStyle(22);
hGeneratedK0shBeta->SetMarkerColor(kRed+1);
hGeneratedK0shBeta->SetLineColor(kRed+1);
hGeneratedK0shBeta->SetMarkerSize(1);
NormalizeByBinWidth(hGeneratedK0shBeta);
hGeneratedK0shBeta->Draw("PE SAME");

TLegend *legBeta = new TLegend(0.55, 0.70, 0.88, 0.88);
legBeta->SetTextSize(0.03);
legBeta->SetBorderSize(0);
legBeta->SetFillColor(kWhite);
legBeta->AddEntry(hCorrectedK0shBeta, "Corrected Spectrum", "lp");
legBeta->AddEntry(hGeneratedK0shBeta, "Generated Spectrum", "lp");
legBeta->Draw();

TLatex latexBeta;
latexBeta.SetNDC();
latexBeta.SetTextSize(0.03);
latexBeta.DrawLatex(0.6, 0.65, "MC Closure Test (Beta)");
latexBeta.DrawLatex(0.6, 0.60, "pp");
latexBeta.DrawLatex(0.6, 0.55, "#sqrt{#it{s}} = "+AO2DFiles.energy);
latexBeta.DrawLatex(0.6, 0.50, "|y| < 0.5");

    // -- Bottom pad: ratio Corrected / Generated --
padBeta2->cd();
padBeta2->SetGridy();
padBeta2->SetTopMargin(0.01);
padBeta2->SetBottomMargin(0.3);

   // Checking errors
hGeneratedK0sh->Sumw2();
hCorrectedK0sh->Sumw2();

    // Create ratio histogram
TH1D *hRatioBeta = (TH1D*)hCorrectedK0shBeta->Clone("hRatioBeta");
hRatioBeta->SetTitle("");
hRatioBeta->Divide(hCorrectedK0shBeta, hGeneratedK0shBeta, 1.0, 1.0);

hRatioBeta->GetYaxis()->SetTitle("Corrected / Generated");
hRatioBeta->GetYaxis()->SetTitleSize(0.08);
hRatioBeta->GetYaxis()->SetTitleOffset(0.5);
hRatioBeta->GetYaxis()->SetLabelSize(0.08);

hRatioBeta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
hRatioBeta->GetXaxis()->SetTitleSize(0.10);
hRatioBeta->GetXaxis()->SetTitleOffset(1.0);
hRatioBeta->GetXaxis()->SetLabelSize(0.08);

hRatioBeta->SetMarkerStyle(33);
hRatioBeta->SetMarkerSize(1);
hRatioBeta->SetMarkerColor(kCyan+2);
hRatioBeta->SetLineColor(kCyan+2);
hRatioBeta->Draw("PE");

    // Reference line at 1
TLine *lineBeta = new TLine(hRatioBeta->GetXaxis()->GetXmin(), 1.0, hRatioBeta->GetXaxis()->GetXmax(), 1.0);
lineBeta->SetLineStyle(2);
lineBeta->Draw("SAME");

cBetaTest->Update();
*/
//End of Beta Test
/*
TH1D* hRun2Xi = nullptr;
if (particle == "Lambda" || particle == "AntiLambda") {
    hRun2Xi = new TH1D("CorrectedRun2XiPtSpectrum", "CorrectedRun2XiPtSpectrum", ptbinedges.size()-1, ptbinedges.data());

    TFile *xi= new TFile("HEPDataXiCorrectedSpectrum.root"); //https://www.hepdata.net/record/ins1797443
    if (!xi->cd("Table 8")) {
        std::cout << "ERROR: Directory 'Table 8' not found in file!" << std::endl;
        return;
    }
    xi->cd("Table 8");
    TGraphErrors *xigraph = (TGraphErrors*)gDirectory->Get("Graph1D_y1");
    if (!xigraph) {
        std::cout << "ERROR: TGraphErrors 'Graph1D_y1' not found in gDirectory!" << std::endl;
        return;
    }
      //
    double mXi=1.321; //mass of Xi Particle
    TF1 *levyTsallisXi = new TF1("levyTsallisXi", "x*[0]*((([1]-1.)*([1]-2.))/(([1]*[2])*(([1]*[2])+([3]*([1]-2.))))) * pow(1.+((sqrt(pow([3], 2.)+pow(x, 2.))-[3])/([1]*[2])), -[1])", 0.1, 10);
    levyTsallisXi->SetParameter(3, mXi);  // Set m as a parameter
    levyTsallisXi->FixParameter(3, mXi);  // Fix m so it does not change during fitting as it is fixed
    levyTsallisXi->SetParameters(0.02, 10, 0.3);
    xigraph->Fit(levyTsallisXi, "R");

    for (size_t bin = 0; bin < ptbinedges.size()-1; ++bin) {  
        double binLow = ptbinedges[bin];
        double binHigh = ptbinedges[bin+1];
        double integral = levyTsallisXi->Integral(binLow, binHigh); // Integrate Levy-Tsallis function over bin range
        hRun2Xi->SetBinContent(bin+1, integral);  // ROOT bins start at 1
    }
}
   */


   /*Run3 Xi Feeddown
TFile *xiFile = TFile::Open("HEPDataXiCorrectedPtSpectrum.root");
if (!xiFile || xiFile->IsZombie()) {
    std::cout << "ERROR: cannot open file!" << std::endl;
    return;
}
        // Access the canvas
TCanvas *canvas = (TCanvas*)xiFile->Get("yieldSummary_XiPm");
if (!canvas) {
    std::cout << "ERROR: Canvas 'yieldSummary_XiPm' not found!" << std::endl;
    return;
}
        // Access the first pad
TPad *pad = (TPad*)canvas->GetPrimitive("padYieldUp_XiPm");
if (!pad) {
    std::cout << "ERROR: Pad 'padYieldUp_XiPm' not found!" << std::endl;
    return;
}
    // Look for the TF1 inside the pad
TH1F* xiHisto = dynamic_cast<TH1F*>(pad->GetPrimitive("Yield_XiPm"));
if (!xiHisto) {
    std::cout << "ERROR: Histogram 'Yield_XiPm' not found!" << std::endl;
    return;
}
if (!xiHisto) {
    std::cout << "ERROR: No TF1 found in pad 'padYieldUp_XiPm'!" << std::endl;
    return;
}
std::cout << "Using TH1F '" << xiHisto->GetName() << "' for feeddown." << std::endl;
        */

/* Run2 Xi Feeddown
TFile *xiFile = TFile::Open("HEPDataXiCorrectedPtSpectrum.root");
if (!xiFile || xiFile->IsZombie()) {
    std::cout << "ERROR: cannot open file!" << std::endl;
    return;
}

   // Access the directory
TDirectory *table3 = (TDirectory*)xiFile->Get("Table 3");
if (!table3) {
    std::cout << "ERROR: Directory 'Table 3' not found!" << std::endl;
    return;
}

   // Access the histogram "Hist1D_y11"
TH1F *xiHisto = dynamic_cast<TH1F*>(table3->Get("Hist1D_y11"));
if (!xiHisto) {
    std::cout << "ERROR: Histogram 'Hist1D_y11' not found in 'Table 3'!" << std::endl;
    return;
}
      */


            //
      //
      /////////////////////////////////////////////////////////////////////////////////
      //Comparison K0shs to Francesca
      //////////////////////////////////////////////////////////////////////////////////
      //
      /*
            TFile *ResultsFrancescaK0sh = TFile::Open(francescaFile);  
            if (!ResultsFrancescaK0sh || ResultsFrancescaK0sh->IsZombie()) {
                std::cerr << "Error: Cannot open Francesca K0sh file: " << francescaFile << "\n";
                return 0;
            }

            TH1D *hFrancescaK0sh = (TH1D*)ResultsFrancescaK0sh->Get(francescaHist);
            if (!hFrancescaK0sh) {
                std::cerr << "Error: Francesca Histogram '" << francescaHist << "' not found.\n";
                return 0;
            }

      // Prepare new canvas with ratio pad
            TCanvas *cCompareFranK0sh = new TCanvas("cCompareFran_K0sh", "Comparison K0shs to Francesca", 800, 800);
            TPad *pad1_FranK0sh = new TPad("pad1_Fran_K0sh", "Main pad", 0, 0.3, 1, 1.0);
            TPad *pad2_FranK0sh = new TPad("pad2_Fran_K0sh", "Ratio pad", 0, 0.05, 1, 0.3);
            pad1_FranK0sh->Draw();
            pad2_FranK0sh->Draw();

      // -- Draw top pad --
            pad1_FranK0sh->cd();
            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);
            pad1_FranK0sh->SetLogy(1);
            pad1_FranK0sh->SetTickx(1);
            pad1_FranK0sh->SetTicky(1);
            pad1_FranK0sh->SetBottomMargin(0.01);

            hK0shCorrectedCompare->Draw("PE");
            hK0shCorrectedCompare->GetXaxis()->SetRangeUser(0.0,8.0);

            hFrancescaK0sh->SetLineColor(kTeal+2);
            hFrancescaK0sh->SetMarkerColor(kTeal+2);
            hFrancescaK0sh->SetMarkerStyle(22);
            hFrancescaK0sh->Draw("SAME PE");

      // Legend
            TLegend* legendcompareFranK0sh = new TLegend(0.55, 0.70, 0.88, 0.88);
            legendcompareFranK0sh->SetTextSize(0.03);
            legendcompareFranK0sh->SetBorderSize(0);
            legendcompareFranK0sh->SetFillColor(kWhite);
            legendcompareFranK0sh->AddEntry(hK0shCorrectedCompare, legendNameK0sh + " Corrected p_{T} Spectrum", "lp");
            TString textCompareFranK0sh = "Main Analysis " + AO2DFiles.energy;
            legendcompareFranK0sh->AddEntry(hFrancescaK0sh, textCompareFranK0sh, "pl");
            legendcompareFranK0sh->Draw();

      // Labels
            TLatex *latexFranK0sh = new TLatex();
            latexFranK0sh->SetNDC();
            latexFranK0sh->SetTextSize(0.03);
            latexFranK0sh->DrawLatex(0.65, 0.65, "pp");
            latexFranK0sh->DrawLatex(0.65, 0.60, "#sqrt{#it{s}} = " + AO2DFiles.energy);
            latexFranK0sh->DrawLatex(0.65, 0.55, "|y| < 0.5");

      // -- Bottom pad: ratio --
            pad2_FranK0sh->cd();
            gPad->SetTopMargin(0.01);
            gPad->SetBottomMargin(0.3);
            gPad->SetGridy();

            TH1D *hRatioFranK0sh = (TH1D*)hK0shCorrectedCompare->Clone("hRatioFranK0sh");
            int nbinsFranK0sh = hRatioFranK0sh->GetNbinsX();
            int firstValidK0sh = 3;
            int lastValidK0sh = nbinsFranK0sh - 9;

            for (int i = 1; i <= nbinsFranK0sh; ++i) {
                double num = hRatioFranK0sh->GetBinContent(i);
                double den = hFrancescaK0sh->GetBinContent(i-2);
                if (i >= firstValidK0sh && i <= lastValidK0sh && den != 0) {
                    hRatioFranK0sh->SetBinContent(i, num / den);
                    double err = hRatioFranK0sh->GetBinError(i);
                    double denErr = hFrancescaK0sh->GetBinError(i-2);
                    double newErr = num / den * std::sqrt( std::pow(err/num,2) + std::pow(denErr/den,2) );
                    hRatioFranK0sh->SetBinError(i, newErr);
                } else {
                    hRatioFranK0sh->SetBinContent(i, 0);
                    hRatioFranK0sh->SetBinError(i, 0);
                }
            }

            hRatioFranK0sh->SetTitle("");
            hRatioFranK0sh->GetYaxis()->SetTitle("My Results / Main Analysis");
            hRatioFranK0sh->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            hRatioFranK0sh->SetMarkerStyle(33);
            hRatioFranK0sh->SetMarkerSize(1);
            hRatioFranK0sh->SetMarkerColor(kCyan+2);
            hRatioFranK0sh->SetLineColor(kCyan+2);
            hRatioFranK0sh->GetYaxis()->SetTitleSize(0.08);
            hRatioFranK0sh->GetYaxis()->SetTitleOffset(0.5);
            hRatioFranK0sh->GetYaxis()->SetLabelSize(0.08);
            hRatioFranK0sh->GetXaxis()->SetLabelSize(0.08);
            hRatioFranK0sh->GetXaxis()->SetTitleSize(0.10);
            hRatioFranK0sh->GetXaxis()->SetTitleOffset(1.0);
            hRatioFranK0sh->GetXaxis()->SetRangeUser(0.0,8.0);
            hRatioFranK0sh->Draw("PE");

            F1->cd();
            gDirectory->cd("Comparisons");
            cCompareFranK0sh->Update();
            cCompareFranK0sh->Write();
        */

                 //
         //
         /////////////////////////////////////////////////////////////////////////////////
         //Comparison Lambda+LambdaBar to Francesca
         //////////////////////////////////////////////////////////////////////////////////
         //
         /*
            TFile *fFranLambda = TFile::Open(francescaFile); //Load Francesca's Lambda histogram
            TFile *fFranAntiLambda = TFile::Open(francescaFile); //Load Francesca's AntiLambda histogram

            if (!fFranLambda || fFranLambda->IsZombie()) {
                std::cerr << "Error opening Results_Lambda.root\n";
                return 0;
            }

            if (!fFranAntiLambda || fFranAntiLambda->IsZombie()) {
                std::cerr << "Error opening Results_AntiLambda.root\n";
                return 0;
            }

            TH1D *hFranLambda = (TH1D*)fFranLambda->Get(francescaHist);
            TH1D *hFranAntiLambda = (TH1D*)fFranAntiLambda->Get(francescaHist);

            if (!hFranLambda || !hFranAntiLambda) {
                std::cerr << "Error loading Francesca histograms\n";
                return 0;
            }

            TH1D *hFrancescaLambdaPlusAntiLambda = (TH1D*)hFranLambda->Clone("hFrancescaLambdaPlusAntiLambda");
            hFrancescaLambdaPlusAntiLambda->Add(hFranAntiLambda); //Add both her results together for comparison

         // Prepare new canvas with ratio pad
            TCanvas *cCompareFranLambda = new TCanvas("cCompareFran_Lambda+AntiLambda", "Comparison Lambda+LambdaBar to Francesca", 800, 800);
            TPad *pad1_FranLambda = new TPad("pad1_Fran_lambda", "Main pad", 0, 0.3, 1, 1.0);
            TPad *pad2_FranLambda = new TPad("pad2_Fran_lambda", "Ratio pad", 0, 0.05, 1, 0.3);
            pad1_FranLambda->Draw();
            pad2_FranLambda->Draw();

         // -- Draw top pad --
            pad1_FranLambda->cd();
            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);
            pad1_FranLambda->SetLogy(1);
            pad1_FranLambda->SetTickx(1);
            pad1_FranLambda->SetTicky(1);
            pad1_FranLambda->SetBottomMargin(0.01);

            hLambdaPlusAntiLambdaForComparison->Draw("PE");
            hLambdaPlusAntiLambdaForComparison->GetXaxis()->SetRangeUser(0.0,8.0);

            hFrancescaLambdaPlusAntiLambda->SetLineColor(kTeal+2);
            hFrancescaLambdaPlusAntiLambda->SetMarkerColor(kTeal+2);
            hFrancescaLambdaPlusAntiLambda->SetMarkerStyle(22);
            hFrancescaLambdaPlusAntiLambda->Draw("SAME PE");

         // Legend
            TLegend* legendcompareFranLambda = new TLegend(0.55, 0.70, 0.88, 0.88);
            legendcompareFranLambda->SetTextSize(0.03);
            legendcompareFranLambda->SetBorderSize(0);
            legendcompareFranLambda->SetFillColor(kWhite);
            legendcompareFranLambda->AddEntry(hLambdaPlusAntiLambdaForComparison, legendNameLambda + " Corrected p_{T} Spectrum", "lp");
            TString textCompareFranLambda = "Main Analysis " + AO2DFiles.energy;
            legendcompareFranLambda->AddEntry(hFrancescaLambdaPlusAntiLambda, textCompareFranLambda, "pl");
            legendcompareFranLambda->Draw();

         // Labels
            TLatex *latexFranLambda = new TLatex();
            latexFranLambda->SetNDC();
            latexFranLambda->SetTextSize(0.03);
            latexFranLambda->DrawLatex(0.65, 0.65, "pp");
            latexFranLambda->DrawLatex(0.65, 0.60, textCompareFranLambda.Data());
            latexFranLambda->DrawLatex(0.65, 0.55, "|y| < 0.5");

         // -- Bottom pad: ratio --
            pad2_FranLambda->cd();
            gPad->SetTopMargin(0.01);
            gPad->SetBottomMargin(0.3);
            gPad->SetGridy();

            TH1D *hRatioFranLambda = (TH1D*)hLambdaPlusAntiLambdaForComparison->Clone("hRatioFranLambda");
            int nbinsFranLambda = hRatioFranLambda->GetNbinsX();
            int firstValidLambda = 3;
            int lastValidLambda = nbinsFranLambda - 9;

            for (int i = 1; i <= nbinsFranLambda; ++i) {
                double num = hRatioFranLambda->GetBinContent(i);
                double den = hFrancescaLambdaPlusAntiLambda->GetBinContent(i-2);
                if (i >= firstValidLambda && i <= lastValidLambda && den != 0) {
                    hRatioFranLambda->SetBinContent(i, num / den);
                    double err = hRatioFranLambda->GetBinError(i);
                    double denErr = hFrancescaLambdaPlusAntiLambda->GetBinError(i-2);
                    double newErr = num / den * std::sqrt( std::pow(err/num,2) + std::pow(denErr/den,2) );
                    hRatioFranLambda->SetBinError(i, newErr);
                } else {
                    hRatioFranLambda->SetBinContent(i, 0);
                    hRatioFranLambda->SetBinError(i, 0);
                }
            }

            hRatioFranLambda->SetTitle("");
            hRatioFranLambda->GetYaxis()->SetTitle("My Results / Main Analysis");
            hRatioFranLambda->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            hRatioFranLambda->SetMarkerStyle(33);
            hRatioFranLambda->SetMarkerSize(1);
            hRatioFranLambda->SetMarkerColor(kCyan+2);
            hRatioFranLambda->SetLineColor(kCyan+2);
            hRatioFranLambda->GetYaxis()->SetTitleSize(0.08);
            hRatioFranLambda->GetYaxis()->SetTitleOffset(0.5);
            hRatioFranLambda->GetYaxis()->SetLabelSize(0.08);
            hRatioFranLambda->GetXaxis()->SetLabelSize(0.08);
            hRatioFranLambda->GetXaxis()->SetTitleSize(0.10);
            hRatioFranLambda->GetXaxis()->SetTitleOffset(1.0);
            hRatioFranLambda->GetXaxis()->SetRangeUser(0.0,8.0);
            hRatioFranLambda->Draw("PE");

            cCompareFranLambda->Update();
            F1->cd();
            gDirectory->cd("Comparisons");
            if (cCompareFranLambda)cCompareFranLambda->Write();
            */