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
bool task_2_3 = true;

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

    TH1D* h1_pTracks = nullptr;
    TH1D* h2_eEMC = nullptr;
    TH2D* h3_eEMC_vs_MUC = nullptr;
    TH1D* h4_eEMC_over_pMDC_e = nullptr;
    TH1D* h4_eEMC_over_pMDC_mu = nullptr;
    TH1D* h4_eEMC_over_pMDC_pi = nullptr;

    // TASK 2.1: Histogram for the absolute momentum of all charged tracks
    if (task_2_1) {
        h1_pTracks = new TH1D("h1absMomentum", "Absolute momentum of charged tracks;|p| [GeV];Events / 0.025 GeV", 100, 0, 2.5);
    }

    // TASK 2.2
    if (task_2_2) {
        h2_eEMC = new TH1D(
        "h2_eEMC",
        "EMC energy of leptons;E_{EMC} [GeV];Events",
        100, 0, 2.0
        );
        h3_eEMC_vs_MUC = new TH2D(
        "h3_eEMC_vs_MUC",
        "E_{EMC} vs MUC depth;E_{EMC} [GeV];d_{MUC} [mm]",
        100, 0, 2.0,
        100, 0, 200
    );
    }

    // TASK 2.3
    if (task_2_3) {
        h4_eEMC_over_pMDC_e = new TH1D(
        "h4_eEMC_over_pMDC_e",
        "eEMC/pMDC ratio;E_{EMC};Events",
        100, 0, 2.0
        );
        h4_eEMC_over_pMDC_mu = new TH1D(
        "h4_eEMC_over_pMDC_mu",
        "eEMC/pMDC ratio;E_{EMC};Events",
        100, 0, 2.0
        );
        h4_eEMC_over_pMDC_pi = new TH1D(
        "h4_eEMC_over_pMDC_pi",
        "eEMC/pMDC ratio;E_{EMC};Events",
        100, 0, 2.0
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
}
