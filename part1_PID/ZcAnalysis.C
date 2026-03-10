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
bool task_2_1 = true;
bool task_2_2 = true;
bool task_2_3 = true;
bool task_3 = true;

void ZcAnalysis::Loop(TString savePath)
{
    // Some settings
    gROOT->SetBatch(true); // no graphical output during execution
    // gStyle->SetOptStat(0); // no statistics box on histograms


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
    // TASK 3
    TH1D* h5_jpsi_mass_electron = nullptr;
    TH1D* h6_jpsi_mass_muon = nullptr;

    // PDG masses in GeV/c^2
    const double electron_mass = 0.000510998946;
    const double muon_mass = 0.1056583745;
    const double pion_mass = 0.13957039;

    // TASK 2.1: Histogram for the absolute momentum of all charged tracks
    if (task_2_1) {
        h1_pTracks = new TH1D("h1absMomentum", "Absolute momentum of charged tracks;|p| [GeV];Events / 0.025 GeV", 100, 0, 2.5);
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

    if (task_3) {
        h5_jpsi_mass_electron = new TH1D(
            "h_jpsi_mass_electron",
            "J/\\psi invariant mass (e^{+}e^{-});M_{J/\\psi} [GeV];Events",
            150, 2.9, 3.2
        );
        h6_jpsi_mass_muon = new TH1D(
            "h_jpsi_mass_muon",
            "J/\\psi invariant mass (\\mu^{+}\\mu^{-});M_{J/\\psi} [GeV];Events",
            150, 2.9, 3.2
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
                // momentum cut to select leptons
                if (p > 1.0)
                {
                    double eEMC = dblTracksECal[i];
                    double dMUC = dblTracksMucDepth[i];

                    h2_eEMC->Fill(eEMC);
                    h3_eEMC_vs_MUC->Fill(eEMC, dMUC);
                }
            }

            // TASK 2.3
            if (task_2_3) {
                double eEMC = dblTracksECal[i];
                double eEMC_over_p = eEMC / p;

                // differentiate between particles
                if (dblTracksProbElectron[i] > dblTracksProbMuon[i] &&
                    dblTracksProbElectron[i] > dblTracksProbPion[i])
                {
                    h4_eEMC_over_pMDC_e->Fill(eEMC_over_p);
                }
                else if (dblTracksProbMuon[i] > dblTracksProbElectron[i] &&
                        dblTracksProbMuon[i] > dblTracksProbPion[i])
                {
                    h4_eEMC_over_pMDC_mu->Fill(eEMC_over_p);
                }
                else if (dblTracksProbPion[i] > dblTracksProbElectron[i] &&
                        dblTracksProbPion[i] > dblTracksProbMuon[i])
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

            // TASK 3
            // Create J/Psi four-vectors by adding lepton pairs
            if (intNumberElectrons == 2 && intNumberMuons == 0) {
                P4M jpsi_electron = e_plus + e_minus;
                h5_jpsi_mass_electron->Fill(jpsi_electron.M());
            }
            else if (intNumberElectrons == 0 && intNumberMuons == 2) {
                P4M jpsi_muon = mu_plus + mu_minus;
                h6_jpsi_mass_muon->Fill(jpsi_muon.M());
            }

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
        h1_pTracks->Draw(); // Draw the histogram on the canvas
        // Save the canvas under the save path (set in the run.sh script)
        canvas->SaveAs(savePath + "1_absoluteMomentum.png"); // You can use .png or .pdf or ...
    }

    // TASK 2.2
    if (task_2_2) {
        {
        TCanvas* canvas = new TCanvas();
        h2_eEMC->Draw();
        canvas->SaveAs(savePath + "2_EEMC_leptons.png");
        }
        {
        TCanvas* canvas = new TCanvas();
        h3_eEMC_vs_MUC->Draw("COLZ");
        canvas->SaveAs(savePath + "3_EEMC_vs_MUC.png");
        }
    }

    // TASK 2.3
    if (task_2_3) {
        TCanvas* canvas = new TCanvas();

        h4_eEMC_over_pMDC_e->SetLineColor(kRed);
        h4_eEMC_over_pMDC_mu->SetLineColor(kBlue);
        h4_eEMC_over_pMDC_pi->SetLineColor(kGreen+2);

        h4_eEMC_over_pMDC_e->Draw();
        h4_eEMC_over_pMDC_mu->Draw("SAME");
        h4_eEMC_over_pMDC_pi->Draw("SAME");

        canvas->SaveAs(savePath + "4_EEMC_over_PMUC.png"); 
    }

    // TASK 3
    if (task_3) {
        // Plot J/Psi invariant mass for electrons
        {
            TCanvas* canvas = new TCanvas();
            h5_jpsi_mass_electron->Draw();
            canvas->SaveAs(savePath + "5_jpsi_mass_electron.png");
        }

        // Plot J/Psi invariant mass for muons
        {
            TCanvas* canvas = new TCanvas();

            // double–Gaussian fit: parameters are
            // [0] = scale 1, [1] = mean1,  [2] = sigma1  (narrow peak)
            // [3] = scale 2, [4] = mean2,  [5] = sigma2  (broader component)
            TF1 *f_fit = new TF1("f_fit",
                "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5])",
                2.9, 3.2);

            // set reasonable starting values: mass ≃3.097, narrow width few MeV,
            // broader width maybe a few × larger
            f_fit->SetParameters(1000, 3.097, 0.005,   200, 3.097, 0.020);

            // optionally fix the second mean equal to the first if desired:
            // f_fit->FixParameter(4,3.097);

            h6_jpsi_mass_muon->Draw();

            h6_jpsi_mass_muon->Fit(f_fit, "R");    // “R” keeps the fit inside [2.9,3.2]

            canvas->SaveAs(savePath + "6_jpsi_mass_muon.png");

            // read out the narrow Gaussian's parameters
            double jpsi_mass  = f_fit->GetParameter(1); // mean of first (narrow) Gaussian
            double jpsi_width = f_fit->GetParameter(2); // sigma of first Gaussian

            printf("J/psi mass = %.6f GeV, width = %.6f GeV\n",
                jpsi_mass, jpsi_width);
        }
    }
}
