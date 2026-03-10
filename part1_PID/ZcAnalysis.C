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

    // Histogram for the absolute momentum of all charged tracks
    TH1D* h1_pTracks = new TH1D("h1absMomentum", "Absolute momentum of charged tracks;|p| [GeV];Events / 0.025 GeV", 100, 0, 2.5);
// problem 2.2
    TH1D* h1_eEMC = new TH1D(
        "h1_eEMC",
        "EMC energy of leptons;E_{EMC} [GeV];Events",
        100, 0, 2.0
        );

    TH2D* h2_eEMC_vs_MUC = new TH2D(
        "h2_eEMC_vs_MUC",
        "E_{EMC} vs MUC depth;E_{EMC} [GeV];d_{MUC} [mm]",
        100, 0, 2.0,
        100, 0, 200
       );
       // problem 2.3
TH1D* h1_EoverP_e = new TH1D(
    "h1_EoverP_e",
    "E/p for particles;E_{EMC}/|p_{MDC}|;Events",
    100, 0, 2
);

TH1D* h1_EoverP_mu = new TH1D(
    "h1_EoverP_mu",
    "E/p for particles;E_{EMC}/|p_{MDC}|;Events",
    100, 0, 2
);

TH1D* h1_EoverP_pi = new TH1D(
    "h1_EoverP_pi",
    "E/p for particles;E_{EMC}/|p_{MDC}|;Events",
    100, 0, 2
);
// problem 3.0
TH1D* h1_mJpsi_e = new TH1D(
    "h1_mJpsi_e",
    "J/#psi invariant mass (electron channel);M_{e^{+}e^{-}} [GeV];Events",
    100, 2.9, 3.2
);

TH1D* h1_mJpsi_mu = new TH1D(
    "h1_mJpsi_mu",
    "J/#psi invariant mass (muon channel);M_{#mu^{+}#mu^{-}} [GeV];Events",
    100, 2.9, 3.2
);
// problem 3.3
TH1D* h1_mJpsi_recoil = new TH1D(
    "h1_mJpsi_recoil",
    "J/#psi mass from #pi^{+}#pi^{-} recoil;M_{recoil}(#pi^{+}#pi^{-}) [GeV];Events",
    100, 2.9, 3.2
);
//Problem 2.4
double m_pi = 0.13957;   // GeV
double m_e  = 0.000511;  // GeV
double m_mu = 0.10566;   // GeV
// problem 3.3
double sqrt_s = 4.25797;
double theta  = 0.011;
double me_beam = 0.000511;

// beam energy and momentum
double Ei = sqrt((sqrt_s*sqrt_s - 4*me_beam*me_beam*sin(theta)*sin(theta))) / (2*cos(theta));
double pi = sqrt(Ei*Ei - me_beam*me_beam);

// beam four-vectors
P4E pe_plus(  pi*sin(theta), 0,  pi*cos(theta), Ei );
P4E pe_minus( -pi*sin(theta),0, -pi*cos(theta), Ei );

// center-of-mass four-vector
P4E pCMS = pe_plus + pe_minus;
int selectedEvents = 0;
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
// problem 2.4
    // containers for particles in this event
    std::vector<P4M> pi_plus;
    std::vector<P4M> pi_minus;
    std::vector<P4M> e_plus;
    std::vector<P4M> e_minus;
    std::vector<P4M> mu_plus;
    std::vector<P4M> mu_minus;
    
        // Loop over all charged tracks in the event (there are always 4 charged tracks: l+ l- pi+ pi-)
        for (int i = 0; i < 4; i++)
        {
            // Create a momentum-vector for each charged track
            P3 pTrack(dblTracksPx[i], dblTracksPy[i], dblTracksPz[i]);
          //problem 2.2  
            double p = pTrack.R();
            double eEMC = dblTracksECal[i];
            double dMUC = dblTracksMucDepth[i];
            // momentum cut to select leptons

            if (p > 1.0)
            {
                
                h1_eEMC->Fill(eEMC);
                h2_eEMC_vs_MUC->Fill(eEMC, dMUC);
            }

// problem 2.3
double Eemc = dblTracksECal[i];

if (p > 0) {

    double EoverP = Eemc / p;

    
    // PID using detector information

// electrons: large E/p
if (EoverP > 0.8)
{
    h1_EoverP_e->Fill(EoverP);
}

// muons: penetrate deeply into MUC
else if (dMUC > 40)
{
    h1_EoverP_mu->Fill(EoverP);
}

// pions: hadrons, small MUC penetration
else
{
    h1_EoverP_pi->Fill(EoverP);
}
    // problem 2.4
double px = dblTracksPx[i];
double py = dblTracksPy[i];
double pz = dblTracksPz[i];

double charge = dblTracksCharge[i];


// identify leptons vs pions
if (p < 1.0)  // pion
{
    P4M pion(px, py, pz, m_pi);

    if (charge > 0)
        pi_plus.push_back(pion);
    else
        pi_minus.push_back(pion);
}
else
{
    if (EoverP > 0.8)   // electron
    {
        P4M electron(px, py, pz, m_e);

        if (charge > 0)
            e_plus.push_back(electron);
        else
            e_minus.push_back(electron);
    }
    else if (dblTracksMucDepth[i] > 40)   // muon
    {
        P4M muon(px, py, pz, m_mu);

        if (charge > 0)
            mu_plus.push_back(muon);
        else
            mu_minus.push_back(muon);
    }
    else
    {
        P4M pion(px, py, pz, m_pi);

        if (charge > 0)
            pi_plus.push_back(pion);
        else
            pi_minus.push_back(pion);
    }
}
     // Fill the absolute momentum of all charged tracks in the histogram
            h1_pTracks->Fill(pTrack.R()); // R() gives the length of the vector
        }
bool electronEvent =
    (pi_plus.size()==1 && pi_minus.size()==1 &&
     e_plus.size()==1 && e_minus.size()==1);

bool muonEvent =
    (pi_plus.size()==1 && pi_minus.size()==1 &&
     mu_plus.size()==1 && mu_minus.size()==1);
// problem 3.3
// reconstruct J/psi from pion recoil
if (pi_plus.size()==1 && pi_minus.size()==1)
{
    P4M pions = pi_plus[0] + pi_minus[0];

    P4E pJpsi_recoil = pCMS - pions;

    h1_mJpsi_recoil->Fill(pJpsi_recoil.M());
}
if (!(electronEvent || muonEvent))
    continue;
    selectedEvents++;
    // // problem 3.0
    if (electronEvent)
{
    P4M Jpsi = e_plus[0] + e_minus[0];
    h1_mJpsi_e->Fill(Jpsi.M());
}

if (muonEvent)
{
    P4M Jpsi = mu_plus[0] + mu_minus[0];
    h1_mJpsi_mu->Fill(Jpsi.M());
} // end of muon loop
} // end of event loop
        //=============================================================================
        // End of selection
        //=============================================================================
    }
    // problem 2.4
    std::cout << "Total events: " << nentries << std::endl;
    std::cout << "Selected events: " << selectedEvents << std::endl;


    //=============================================================================
    // All canvas definitions should be placed here
    //=============================================================================

    // Absolute momentum of all charged tracks
    // problem 2.2
    {
        TCanvas* canvas = new TCanvas(); // Create an empty canvas

        double binWidth =  h1_eEMC->GetBinWidth(1);
        h1_eEMC->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

        h1_eEMC->Draw();
        canvas->SaveAs(savePath + "2_EEMC_leptons.png");

    {
        TCanvas* canvas = new TCanvas();
        h2_eEMC_vs_MUC->Draw("COLZ");
        canvas->SaveAs(savePath + "3_EEMC_vs_MUC.png");
    }
    // problem 2.3
{
    TCanvas* canvas = new TCanvas("canvas","canvas",1200,900);
    
    double binWidth =  h1_EoverP_e->GetBinWidth(1);
    h1_EoverP_e->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));
    

    h1_EoverP_e->SetLineColor(kRed);
    h1_EoverP_mu->SetLineColor(kBlue);
    h1_EoverP_pi->SetLineColor(kGreen+2);

    h1_EoverP_e->Draw();
    h1_EoverP_mu->Draw("SAME");
    h1_EoverP_pi->Draw("SAME");

    TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(h1_EoverP_e,"Electrons","l");
    legend->AddEntry(h1_EoverP_mu,"Muons","l");
    legend->AddEntry(h1_EoverP_pi,"Pions","l");
    legend->Draw();

    canvas->SaveAs(savePath + "4_EoverP_particles.png");
}
// problem 3.0
{
   TCanvas* canvas = new TCanvas("canvas","canvas",1200,900);
    canvas->SetLeftMargin(0.15);
    double binWidth = h1_mJpsi_e->GetBinWidth(1);
    h1_mJpsi_e->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

    h1_mJpsi_e->Draw();
    // problem 3.2
    TF1* ffit = new TF1("ffit","crystalball",2.9,3.2);

    // good initial parameters
    ffit->SetParameters(
        20000,   // normalization
        3.097,   // mean
        0.02,    // sigma
        1.5,     // alpha
        2.0      // n
    );

    h1_mJpsi_e->Fit(ffit,"R");

    canvas->SaveAs(savePath + "Jpsi_mass_electron_CB.png");

}

//problem 3.1
//{
    //TCanvas* canvas = new TCanvas();

    //h1_mJpsi_mu->Draw();

    //TF1* ffit = new TF1("ffit","gaus",2.9,3.2);
    //ffit->SetParameters(100,3.097,0.02);

    //h1_mJpsi_mu->Fit(ffit,"R");

    //canvas->SaveAs(savePath + "6_Jpsi_mass_muon_fit.png");
//}
// problem 3.1
{
    TCanvas* canvas = new TCanvas("canvas","canvas",1200,900);
    canvas->SetLeftMargin(0.15);
    double binWidth = h1_mJpsi_mu->GetBinWidth(1);
    h1_mJpsi_mu->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

    h1_mJpsi_mu->Draw();

    TF1* ffit = new TF1(
    "ffit",
    "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + [3]*exp(-(x-[1])*(x-[1])/(2*[4]*[4]))",
    2.9,3.2);
    

    ffit->SetParameters(
    20000,   // amplitude narrow
    3.097,   // mean
    0.01,    // narrow sigma
    5000,    // amplitude broad
    0.04     // broad sigma
    );

    h1_mJpsi_mu->Fit(ffit,"R");

    // Narrow component
TF1* g1 = new TF1(
"g1",
"[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",
2.9,3.2);

g1->SetParameters(
ffit->GetParameter(0),
ffit->GetParameter(1),
ffit->GetParameter(2));

g1->SetLineColor(kGreen);
g1->SetLineStyle(2);
g1->Draw("same");


// Broad component
TF1* g2 = new TF1(
"g2",
"[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",
2.9,3.2);

g2->SetParameters(
ffit->GetParameter(3),
ffit->GetParameter(1),
ffit->GetParameter(4));

g2->SetLineColor(kMagenta);
g2->SetLineStyle(2);
g2->Draw("same");

    canvas->SaveAs(savePath + "6_Jpsi_mass_muon_doubleFit.png");
}
// problem 3.3
{
    TCanvas* canvas = new TCanvas("canvas","canvas",1200,900);
    canvas->SetLeftMargin(0.15);

    double binWidth = h1_mJpsi_recoil->GetBinWidth(1);
    h1_mJpsi_recoil->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

    h1_mJpsi_recoil->Draw();

    TF1* ffit_gaus = new TF1("ffit_gaus","gaus",2.95,3.15);

    ffit_gaus->SetParameters(
        20000,   // amplitude
        3.097,   // mean
        0.01     // sigma
    );

    h1_mJpsi_recoil->Fit(ffit_gaus,"R");

    canvas->SaveAs(savePath + "Jpsi_recoil_gaussian_fit.png");
}
{
    TCanvas* canvas = new TCanvas("canvas","canvas",1200,900);
    canvas->SetLeftMargin(0.15);

    h1_mJpsi_recoil->Draw();
    double binWidth = h1_mJpsi_recoil->GetBinWidth(1);
    double norm = h1_mJpsi_recoil->GetMaximum();

TF1* ffit_voigt = new TF1(
    "ffit_voigt",
    Form("[0]*TMath::Voigt(x-[1],[2],[3])*%f", binWidth),
    2.95,
    3.15
);

ffit_voigt->SetParameters(
    norm,     // normalization
    3.097,    // mean
    0.01,     // sigma (detector resolution)
    0.0001    // gamma (natural width)
);

ffit_voigt->SetLineColor(kRed);
ffit_voigt->SetLineWidth(2);

h1_mJpsi_recoil->Fit(ffit_voigt,"R");

    canvas->SaveAs(savePath + "Jpsi_recoil_voigt_fit.png");
}
        // h1_pTracks->Draw(); // Draw the histogram on the canvas

        // Save the canvas under the save path (set in the run.sh script)
        // canvas->SaveAs(savePath + "1_absoluteMomentum.png"); // You can use .png or .pdf or ...
    }
}
