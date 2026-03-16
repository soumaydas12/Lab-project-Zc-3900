#define ZcAnalysis_cxx
#include "ZcAnalysis.h"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLine.h>
#include <TLegend.h>
#include <TFile.h>
#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

// Vector definitions
using P3 = ROOT::Math::XYZVector;       // Three-vector
using P4E = ROOT::Math::PxPyPzEVector;  // Four-vector with energy
using P4M = ROOT::Math::PxPyPzMVector;  // Four-vector with mass
using ROOT::Math::VectorUtil::CosTheta; // angle between two four-vectors


//=============================================================================
// In this Loop function:
// 1. All histograms are defined
// 2. A for-loop iterates over all events in the data
// 3. All histograms are plotted on a canvas and saved
//=============================================================================

// decide which graphs should be plotted
bool task_2_1 = false;
bool task_2_2 = false;
bool task_2_3 = false;
bool task_3_1_and_3_2 = false;
bool task_3_3 = false;
bool task_3_5 = false;
bool task_3_6 = true;

void ZcAnalysis::Loop(TString savePath)
{
    // Some settings
    gROOT->SetBatch(true); // no graphical output during execution
    gStyle->SetOptStat(0); // no statistics box on histograms

    //=============================================================================
    // All histogram definitions should be placed here:
    // 1D parameters are: ("uniqueName", "title;x-axis-label;y-axis-label", bins, min, max)
    // 2D parameters are: ("uniqueName", "title;x-axis-label;y-axis-label;z-axis-label", x-bins, x-min, x-max, y-bins, y-min, y-max)
    //=============================================================================

    // TASK 2.1
    TH1D* h1_pTracks = nullptr;
    // TASK 2.2
    TH1D* h2_eEMC = nullptr;
    TH2D* h3_eEMC_vs_MUC = nullptr;
    // TASK 2.3
    TH1D* h4_eEMC_over_pMDC_e = nullptr;
    TH1D* h4_eEMC_over_pMDC_mu = nullptr;
    TH1D* h4_eEMC_over_pMDC_pi = nullptr;
    // TASK 3.1/3.2
    TH1D* h6_Jpsi_mass_electron = nullptr;
    TH1D* h5_Jpsi_mass_muon = nullptr;
    // TASK 3.3
    TH1D* h7_mJpsi_pion_recoil = nullptr;
    // TASK 3.5
    TH1D* h8_mass_ee = nullptr;
    TH1D* h8_mass_mumu = nullptr;
    TH1D* h8_mass_pipi = nullptr;
    TH1D* h8_mass_recoil = nullptr;
    // TASK 3.6
    TH2D* h9_dalitz = nullptr;

    // PDG masses in GeV/c^2
    const double electron_mass = 0.000510998946;
    const double muon_mass = 0.1056583745;
    const double pion_mass = 0.13957039;

    // Beam parameters
    double sqrt_s = 4.25797;
    double theta = 0.011;
    double me_beam = 0.000511;

    // Beam energy and momentum
    double Ei = sqrt((sqrt_s*sqrt_s - 4*me_beam*me_beam*sin(theta)*sin(theta))) / (2*cos(theta));
    double pi_beam = sqrt(Ei*Ei - me_beam*me_beam);

    // Beam four-vectors
    P4E pe_plus(  pi_beam*sin(theta), 0,  pi_beam*cos(theta), Ei );
    P4E pe_minus( -pi_beam*sin(theta), 0, -pi_beam*cos(theta), Ei );

    // Center-of-mass four-vector
    P4E pCMS = pe_plus + pe_minus;

    // TASK 2.1: Histogram for the absolute momentum of all charged tracks
    if (task_2_1) {
        h1_pTracks = new TH1D("h1absMomentum", 
            "Absolute momentum of charged tracks;|p| [GeV];Events", 
            100, 0, 2.5
        );
    }

    // TASK 2.2
    if (task_2_2) {
        h2_eEMC = new TH1D(
        "h2_eEMC",
        "E_{EMC} of leptons;E_{EMC} [GeV];Events",
        100, 0, 2.0
        );
        h3_eEMC_vs_MUC = new TH2D(
        "h3_eEMC_vs_MUC",
        "E_{EMC} vs d_{MUC};E_{EMC} [GeV];d_{MUC} [mm]",
        100, 0, 2.0,
        100, 0, 200
        );
    }

    // TASK 2.3
    if (task_2_3) {
        h4_eEMC_over_pMDC_e = new TH1D(
        "h4_eEMC_over_pMDC_e",
        "E_{EMC}/p_{MDC} ratio;E_{EMC}/p_{MDC};Events",
        100, 0, 2.0
        );
        h4_eEMC_over_pMDC_mu = new TH1D(
        "h4_eEMC_over_pMDC_mu",
        "E_{EMC}/p_{MDC} ratio;E_{EMC}/p_{MDC};Events",
        100, 0, 2.0
        );
        h4_eEMC_over_pMDC_pi = new TH1D(
        "h4_eEMC_over_pMDC_pi",
        "E_{EMC}/p_{MDC} ratio;E_{EMC}/p_{MDC};Events",
        100, 0, 2.0
        );
    }

    // TASK 3.1/3.2
    if (task_3_1_and_3_2) {
        h6_Jpsi_mass_electron = new TH1D(
            "h6_Jpsi_mass_electron",
            "J/#psi invariant mass (e^{+}e^{-});M_{J/#psi} [GeV];Events",
            150, 2.9, 3.2
        );
        h5_Jpsi_mass_muon = new TH1D(
            "h5_Jpsi_mass_muon",
            "J/#psi invariant mass (#mu^{+}#mu^{-});M_{J/#psi} [GeV];Events",
            150, 2.9, 3.2
        );
    }

    // TASK 3.3
    if (task_3_3) {
        h7_mJpsi_pion_recoil = new TH1D(
            "h7_mJpsi_recoil",
            "J/#psi mass from #pi^{+}#pi^{-} recoil;M_{recoil}(#pi^{+}#pi^{-}) [GeV];Events",
            100, 2.9, 3.2
        );
    }

    // TASK 3.5: Combined mass spectra
    if (task_3_5) {
        h8_mass_ee = new TH1D("h8_mass_ee", "Invariant mass spectra;Mass [GeV];Events", 200,0,4);
        h8_mass_mumu = new TH1D("h8_mass_mumu", "Invariant mass spectra;Mass [GeV];Events", 200,0,4);
        h8_mass_pipi = new TH1D("h8_mass_pipi", "Invariant mass spectra;Mass [GeV];Events", 200,0,4);
        h8_mass_recoil = new TH1D("h8_mass_recoil", "Invariant mass spectra;Mass [GeV];Events", 200,0,4);
    }

    // TASK 3.6: Dalitz plot
    if (task_3_6) {
        h9_dalitz = new TH2D(
            "h_dalitz",
            "Dalitz Plot; m^{2}(#pi^{+}#pi^{-}); m^{2}(J/#psi #pi)",
            300, 0, 1.5,
            200, 9, 18
        );
    }
    


    //=============================================================================
    // For-loop over all events in the root file
    //=============================================================================

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry); nbytes += nb;

        //=============================================================================
        // Your selection for each event
        //=============================================================================

        // Declare four-vectors for the six-particles
        // pions
        P4M pion_plus, pion_minus;
        // leptons
        P4M e_plus, e_minus;
        P4M mu_plus, mu_minus;

        // Skip events that don't have exactly 2 pions and 2 leptons of the same family
        if (intNumberPions != 2) continue;
        if (!((intNumberElectrons == 2 && intNumberMuons == 0) || 
            (intNumberElectrons == 0 && intNumberMuons == 2))) continue;

        // Loop over all charged tracks in the event (there are always 4 charged tracks: l+ l- pi+ pi-)
        for (int i = 0; i < 4; i++)
        {
            // Create a momentum-vector for each charged track
            P3 pTrack(dblTracksPx[i], dblTracksPy[i], dblTracksPz[i]);
            
            double p = pTrack.R();

            // TASK 2.1
            if (task_2_1) {
                // Fill the absolute momentum of all charged tracks in the histogram
                h1_pTracks->Fill(p); // R() gives the length of the vector
            }

            // TASK 2.2
            if (task_2_2) {
                double eEMC = dblTracksECal[i];
                double dMUC = dblTracksMucDepth[i];
                // momentum cut to select leptons
                if (p > 0.8)
                {

                    h2_eEMC->Fill(eEMC);
                    h3_eEMC_vs_MUC->Fill(eEMC, dMUC);
                }
            }

            // TASK 2.3
            if (task_2_3) {
                double eEMC = dblTracksECal[i];
                double eEMC_over_p = eEMC / p;

                // differentiate between particles
                // electrons
                if (p > 1.0 && eEMC_over_p > 0.8)
                {
                    h4_eEMC_over_pMDC_e->Fill(eEMC_over_p);
                }
                // muons
                else if (p > 1.0 && eEMC_over_p < 0.8)
                {
                    h4_eEMC_over_pMDC_mu->Fill(eEMC_over_p);
                }
                else if (p < 1.0)
                {
                    h4_eEMC_over_pMDC_pi->Fill(eEMC_over_p);
                }
            }

            // TASK 2.4
            // Create four-vectors with measured momentum from MDC and PDG masses
            pion_plus = P4M(dblTracksPx[2], dblTracksPy[2], dblTracksPz[2], pion_mass);
            pion_minus = P4M(dblTracksPx[3], dblTracksPy[3], dblTracksPz[3], pion_mass);
            e_plus = P4M(dblTracksPx[0], dblTracksPy[0], dblTracksPz[0], electron_mass);
            e_minus = P4M(dblTracksPx[1], dblTracksPy[1], dblTracksPz[1], electron_mass);
            mu_plus = P4M(dblTracksPx[0], dblTracksPy[0], dblTracksPz[0], muon_mass);
            mu_minus = P4M(dblTracksPx[1], dblTracksPy[1], dblTracksPz[1], muon_mass);

        }

        // TASK 3/3.1/3.2
        if (task_3_1_and_3_2) {
            if (intNumberElectrons == 2 && intNumberMuons == 0) {
                P4M jpsi_electron = e_plus + e_minus;
                h6_Jpsi_mass_electron->Fill(jpsi_electron.M());
            }
            else if (intNumberElectrons == 0 && intNumberMuons == 2) {
                P4M jpsi_muon = mu_plus + mu_minus;
                h5_Jpsi_mass_muon->Fill(jpsi_muon.M());
            }
        }

        // TASK 3.3
        if (task_3_3) {
            P4M pipi = pion_plus + pion_minus;
            P4E recoil = pCMS - pipi;
            h7_mJpsi_pion_recoil->Fill(recoil.M());
        }

        // TASK 3.5
        if (task_3_5) {
            // Fill combined mass histograms
            if (intNumberElectrons == 2 && intNumberMuons == 0) {
                P4M ee = e_plus + e_minus;
                h8_mass_ee->Fill(ee.M());
            }
            if (intNumberElectrons == 0 && intNumberMuons == 2) {
                P4M mumu = mu_plus + mu_minus;
                h8_mass_mumu->Fill(mumu.M());
            }
            P4M pipi = pion_plus + pion_minus;
            h8_mass_pipi->Fill(pipi.M());
            // recoil mass
            P4E recoil = pCMS - pipi;
            h8_mass_recoil->Fill(recoil.M());
        }

        // TASK 3.6
        if (task_3_6) {

            P4M jpsi;

            if (intNumberElectrons == 2) {
                jpsi = e_plus + e_minus;
            }
            else {
                jpsi = mu_plus + mu_minus;
            }

            double m2_pipi = (pion_plus + pion_minus).M2();
            double m2_jpsi_pi_plus  = (jpsi + pion_plus).M2();
            double m2_jpsi_pi_minus = (jpsi + pion_minus).M2();

            h9_dalitz->Fill(m2_pipi, m2_jpsi_pi_plus);
            h9_dalitz->Fill(m2_pipi, m2_jpsi_pi_minus);
        }

        //=============================================================================
        // End of selection
        //=============================================================================
    }


    //=============================================================================
    // All canvas definitions should be placed here
    //=============================================================================

    // TASK 2.1: Absolute momentum of all charged tracks
    if (task_2_1) {
        TCanvas* canvas = new TCanvas(); // Create an empty canvas
        h1_pTracks->SetLineWidth(2); // thicker histogram line
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        h1_pTracks->GetYaxis()->SetTitleSize(0.05);
        h1_pTracks->GetXaxis()->SetTitleSize(0.05);
        double binWidth =  h1_pTracks->GetBinWidth(1);
        h1_pTracks->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));
        h1_pTracks->Draw(); // Draw the histogram on the canvas
        // Save the canvas under the save path (set in the run.sh script)
        canvas->SaveAs(savePath + "1_absoluteMomentum.png"); // You can use .png or .pdf or ...
    }

    // TASK 2.2:
    if (task_2_2) {
        {
        TCanvas* canvas = new TCanvas();
        h2_eEMC->SetLineWidth(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        h2_eEMC->GetYaxis()->SetTitleSize(0.05);
        h2_eEMC->GetXaxis()->SetTitleSize(0.05);
        double binWidth =  h1_pTracks->GetBinWidth(1);
        h2_eEMC->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));
        h2_eEMC->Draw();
        canvas->SaveAs(savePath + "2_EEMC_leptons.png");
        }
        {
        TCanvas* canvas = new TCanvas();
        h3_eEMC_vs_MUC->SetLineWidth(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.15);
        h3_eEMC_vs_MUC->GetYaxis()->SetTitleSize(0.05);
        h3_eEMC_vs_MUC->GetXaxis()->SetTitleSize(0.05);
        h3_eEMC_vs_MUC->Draw("COLZ");
        canvas->SaveAs(savePath + "3_EEMC_vs_MUC.png");
        }
    }

    // TASK 2.3
    if (task_2_3) {
        TCanvas* canvas = new TCanvas();

        h4_eEMC_over_pMDC_e->SetLineColor(kRed);
        h4_eEMC_over_pMDC_e->SetLineWidth(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        h4_eEMC_over_pMDC_e->GetYaxis()->SetTitleSize(0.05);
        h4_eEMC_over_pMDC_e->GetXaxis()->SetTitleSize(0.05);

        double binWidth =  h4_eEMC_over_pMDC_e->GetBinWidth(1);
        h4_eEMC_over_pMDC_e->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

        h4_eEMC_over_pMDC_e->Draw("COLZ");

        h4_eEMC_over_pMDC_mu->SetLineColor(kBlue);
        h4_eEMC_over_pMDC_pi->SetLineColor(kGreen+2);

        h4_eEMC_over_pMDC_e->Draw();
        h4_eEMC_over_pMDC_mu->Draw("SAME");
        h4_eEMC_over_pMDC_pi->Draw("SAME");

        canvas->SaveAs(savePath + "4_EEMC_over_PMUC.png"); 
    }

    // TASK 3.1/3.2
    if (task_3_1_and_3_2) {
       
        // Plot J/Psi invariant mass for muons
        {
            // single-Gaussian fit
            {
                TCanvas* canvas = new TCanvas();

                TF1 *f_single = new TF1("f_single", "gaus", 2.9, 3.2);

                f_single->SetParameters(25000, 3.097, 0.012);

                h5_Jpsi_mass_muon->SetLineWidth(2);
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.15);
                h5_Jpsi_mass_muon->GetYaxis()->SetTitleSize(0.05);
                h5_Jpsi_mass_muon->GetXaxis()->SetTitleSize(0.05);

                double binWidth =  h5_Jpsi_mass_muon->GetBinWidth(1);
                h5_Jpsi_mass_muon->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

                h5_Jpsi_mass_muon->Draw();

                h5_Jpsi_mass_muon->Fit(f_single, "R");

                canvas->SaveAs(savePath + "5_Jpsi_mass_muon_single_gaussian.png");
            }

            // double–Gaussian fit
            {
                TCanvas* canvas = new TCanvas();

                // parameters are
                // [0] = scale 1, [1] = mean1,  [2] = sigma1  (narrow peak)
                // [3] = scale 2, [4] = mean2,  [5] = sigma2  (broader component)
                TF1 *f_double = new TF1("f_double",
                    "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5])",
                    2.9, 3.2);

                // set reasonable starting values: mass ≃3.097, narrow width few MeV,
                f_double->SetParameters(1000, 3.097, 0.005,   200, 3.097, 0.020);

                h5_Jpsi_mass_muon->SetLineWidth(2);
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.15);
                h5_Jpsi_mass_muon->GetYaxis()->SetTitleSize(0.05);
                h5_Jpsi_mass_muon->GetXaxis()->SetTitleSize(0.05);

                double binWidth =  h5_Jpsi_mass_muon->GetBinWidth(1);
                h5_Jpsi_mass_muon->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

                h5_Jpsi_mass_muon->Draw();

                h5_Jpsi_mass_muon->Fit(f_double, "R");    // “R” keeps the fit inside [2.9,3.2]

                canvas->SaveAs(savePath + "5_Jpsi_mass_muon_double_gaussian.png");

                // read out the narrow Gaussian's parameters
                double jpsi_mass  = f_double->GetParameter(1); // mean of first (narrow) Gaussian
                double jpsi_width = f_double->GetParameter(2); // sigma of first Gaussian

                printf("J/psi mass = %.6f GeV, width = %.6f GeV\n",
                    jpsi_mass, jpsi_width);
            }
        }

        // Plot J/Psi invariant mass for electrons
        {
            TCanvas* canvas = new TCanvas();

            // crystal ball fit
            TF1 *f_fit = new TF1("f_fit",
                "crystalball", 
                2.9, 3.2);

            f_fit->SetParameters(
                20000,  // normalization
                3.097,  // mean      
                0.02,   // sigma
                1.5,    // alpha (tail transition)
                2.0     // n (tail exponent)
            );

            h6_Jpsi_mass_electron->SetLineWidth(2);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            h6_Jpsi_mass_electron->GetYaxis()->SetTitleSize(0.05);
            h6_Jpsi_mass_electron->GetXaxis()->SetTitleSize(0.05);

            double binWidth =  h6_Jpsi_mass_electron->GetBinWidth(1);
            h6_Jpsi_mass_electron->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

            h6_Jpsi_mass_electron->Draw();

            h6_Jpsi_mass_electron->Fit(f_fit, "R");

            canvas->SaveAs(savePath + "6_Jpsi_mass_electron.png");

            // read out parameters
            double jpsi_mass  = f_fit->GetParameter(3);
            double jpsi_width = f_fit->GetParameter(4);

            printf("J/psi mass = %.6f GeV, width = %.6f GeV\n",
                jpsi_mass, jpsi_width);
        }
    }

    // TASK 3.3
    if (task_3_3) {

        // single-Gaussian fit
        {
            TCanvas* canvas = new TCanvas("canvas","canvas", 1200, 900);

            TF1* f_gaus = new TF1("f_gaus",
                "gaus",
                2.95, 
                3.15
            );

            f_gaus->SetParameters(
                20000,   // amplitude
                3.097,   // mean
                0.01     // sigma
            );

            h7_mJpsi_pion_recoil->Fit(f_gaus,"R");

            h7_mJpsi_pion_recoil->SetLineWidth(2);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            h7_mJpsi_pion_recoil->GetYaxis()->SetTitleSize(0.05);
            h7_mJpsi_pion_recoil->GetXaxis()->SetTitleSize(0.05);

            double binWidth =  h7_mJpsi_pion_recoil->GetBinWidth(1);
            h7_mJpsi_pion_recoil->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

            h7_mJpsi_pion_recoil->Draw();

            canvas->SaveAs(savePath + "7_mJpsi_recoil_from_pions_gaussian.png");
        }

        // Voigt fit
        {
            TCanvas* canvas = new TCanvas("canvas","canvas", 1200, 900);

            double binWidth = h7_mJpsi_pion_recoil->GetBinWidth(1);
            double norm = h7_mJpsi_pion_recoil->GetMaximum();

            TF1* f_voigt = new TF1(
                "f_voigt",
                Form("[0]*TMath::Voigt(x-[1],[2],[3])*%f", binWidth),
                2.95,
                3.15
            );

            f_voigt->SetParameters(
                norm,     // normalization
                3.097,    // mean
                0.01,     // sigma (detector resolution)
                0.0001    // gamma (natural width)
            );

            h7_mJpsi_pion_recoil->Fit(f_voigt,"R");

            h7_mJpsi_pion_recoil->SetLineWidth(2);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            h7_mJpsi_pion_recoil->GetYaxis()->SetTitleSize(0.05);
            h7_mJpsi_pion_recoil->GetXaxis()->SetTitleSize(0.05);

            h7_mJpsi_pion_recoil->Draw();

            canvas->SaveAs(savePath + "7_mJpsi_recoil_from_pions_voigt.png");
        }
        
    }

    // TASK 3.5
    if (task_3_5) {
        TCanvas* canvas = new TCanvas("canvas_mass","canvas_mass",1200,900);
        h8_mass_ee->SetLineColor(kRed);
        h8_mass_mumu->SetLineColor(kBlue);
        h8_mass_pipi->SetLineColor(kGreen+2);
        h8_mass_recoil->SetLineColor(kMagenta);

        h8_mass_ee->SetLineWidth(2);
        h8_mass_mumu->SetLineWidth(2);
        h8_mass_pipi->SetLineWidth(2);
        h8_mass_recoil->SetLineWidth(2);

        h8_mass_ee->SetMaximum(1.3 * h8_mass_recoil->GetMaximum());
        
        double binWidth = h8_mass_ee->GetBinWidth(1);

        h8_mass_ee->SetLineWidth(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        h8_mass_ee->GetYaxis()->SetTitleSize(0.05);
        h8_mass_ee->GetXaxis()->SetTitleSize(0.05);

        h8_mass_ee->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

        h8_mass_ee->Draw();
        h8_mass_mumu->Draw("SAME");
        h8_mass_pipi->Draw("SAME");
        h8_mass_recoil->Draw("SAME");

        TLegend* legend = new TLegend(0.15,0.65,0.38,0.88);
        legend->AddEntry(h8_mass_ee,"e^{+}e^{-}","l");
        legend->AddEntry(h8_mass_mumu,"#mu^{+}#mu^{-}","l");
        legend->AddEntry(h8_mass_pipi,"#pi^{+}#pi^{-}","l");
        legend->AddEntry(h8_mass_recoil,"Recoil(#pi^{+}#pi^{-})","l");
        legend->Draw();
        
        canvas->SaveAs(savePath + "8_all_mass_spectra.png");
    }

    // TASK 3.6
    if (task_3_6) {
        TCanvas* c_dalitz = new TCanvas();
        h9_dalitz->Draw("COLZ");

        // overlay the kinematic boundary
        const double M  = sqrt_s;      // total CM energy
        const double m1 = pion_mass;
        const double m2 = pion_mass;
        const double m3 = 3.0969;

        std::vector<double> x;
        std::vector<double> y_min;
        std::vector<double> y_max;

        for (double s12 = pow(2*pion_mass,2); s12 < pow(M - m3,2); s12 += 0.005)
        {
            double sqrt_s12 = sqrt(s12);

            double E2 = (s12 - m1*m1 + m2*m2) / (2*sqrt_s12);
            double E3 = (M*M - s12 - m3*m3) / (2*sqrt_s12);

            double p2 = sqrt(E2*E2 - m2*m2);
            double p3 = sqrt(E3*E3 - m3*m3);

            double s23_max = pow(E2+E3,2) - pow(p2-p3,2);
            double s23_min = pow(E2+E3,2) - pow(p2+p3,2);

            x.push_back(s12);
            y_min.push_back(s23_min);
            y_max.push_back(s23_max);
        }

        TGraph* g_min = new TGraph(x.size(), &x[0], &y_min[0]);
        TGraph* g_max = new TGraph(x.size(), &x[0], &y_max[0]);

        g_min->SetLineColor(kRed);
        g_max->SetLineColor(kRed);

        g_min->SetLineWidth(3);
        g_max->SetLineWidth(3);

        double binWidth = h9_dalitz->GetBinWidth(1);

        h9_dalitz->SetLineWidth(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        h9_dalitz->GetYaxis()->SetTitleSize(0.05);
        h9_dalitz->GetXaxis()->SetTitleSize(0.05);

        g_min->Draw("L SAME");
        g_max->Draw("L SAME");

        c_dalitz->SaveAs(savePath + "9_dalitz_plot.png");
    }
}
