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

// dalitz plot boundary
double lambda(double a, double b, double c)
{
    return a*a + b*b + c*c - 2*a*b - 2*a*c - 2*b*c;
}

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

    gStyle->SetStatX(0.3);
    gStyle->SetStatY(0.9);


    //=============================================================================
    // All histogram definitions should be placed here:
    // 1D parameters are: ("uniqueName", "title;x-axis-label;y-axis-label", bins, min, max)
    // 2D parameters are: ("uniqueName", "title;x-axis-label;y-axis-label;z-axis-label", x-bins, x-min, x-max, y-bins, y-min, y-max)
    //=============================================================================

   
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


// problem 3.6
TH2D* h2_dalitz = new TH2D(
"h2_dalitz",
"Dalitz plot; m^{2}_{#pi^{+}#pi^{-}} [GeV^{2}]; m^{2}_{J/#psi #pi} [GeV^{2}]",
300,0,1.5,
300,9,18
);
// problem 4.1
TH1D* h1_mJpsi_pipi = new TH1D(
    "h1_mJpsi_pipi",
    "Invariant mass of J/#psi #pi^{+}#pi^{-};M_{J/#psi#pi#pi} [GeV];Events",
    200, 3.0, 4.5
);

TH1D* h1_Egamma = new TH1D(
    "h1_Egamma",
    "ISR photon energy;E_{#gamma} [GeV];Events",
    200, 0, 2
);
TH1D* h1_Egamma_ISR = new TH1D(
    "h1_Egamma_ISR",
    "ISR photon energy from kinematics;E_{#gamma} [GeV];Events",
    200, 0, 2
);

TH1D* h1_Emiss = new TH1D(
    "h1_Emiss",
    "Missing energy;E_{miss} [GeV];Events",
    200, 0, 2
);// before cut for problem 4.2
TH1D* h1_mJpsi_pipi_before = new TH1D(
    "h1_mJpsi_pipi_before",
    "Invariant mass of J/#psi #pi^{+}#pi^{-};M_{J/#psi#pi#pi} [GeV];Events",
    200, 3.0, 4.5
);

TH1D* h1_Emiss_before = new TH1D(
    "h1_Emiss_before",
    "Missing energy;E_{miss} [GeV];Events",
    200, 0, 2
);
// problem 4.2
TH1D* h1_chi2_4C = new TH1D(
    "h1_chi2_4C",
    "4C kinematic fit #chi^{2};#chi^{2}_{4C};Events",
    200,0,200
);
//Problem 2.4
double m_jpsi = 3.0969;  // Gev
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
    for (Long64_t jentry = 0; jentry < nentries; jentry++) // event loop
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry); nbytes += nb;

        // ---------------- 4C KINEMATIC FIT FOUR-VECTORS ----------------

P4E pi_plus(
dblKinFit4CMomentaPx[0],
dblKinFit4CMomentaPy[0],
dblKinFit4CMomentaPz[0],
dblKinFit4CMomentaE[0]
);

P4E pi_minus(
dblKinFit4CMomentaPx[1],
dblKinFit4CMomentaPy[1],
dblKinFit4CMomentaPz[1],
dblKinFit4CMomentaE[1]
);

P4E lep_plus(
dblKinFit4CMomentaPx[2],
dblKinFit4CMomentaPy[2],
dblKinFit4CMomentaPz[2],
dblKinFit4CMomentaE[2]
);

P4E lep_minus(
dblKinFit4CMomentaPx[3],
dblKinFit4CMomentaPy[3],
dblKinFit4CMomentaPz[3],
dblKinFit4CMomentaE[3]
);

        //=============================================================================
        // Your selection for each event
        //=============================================================================

    
        // Loop over all charged tracks in the event (there are always 4 charged tracks: l+ l- pi+ pi-)

        // problem 3.5
        // Fill invariant mass histograms


bool electronEvent = (intNumberElectrons == 2);
bool muonEvent     = (intNumberMuons == 2);
// problem 3.3
// reconstruct J/psi from pion recoil

if (!(electronEvent || muonEvent))
continue;

P4E pions_recoil = pi_plus + pi_minus;
P4E pJpsi_recoil = pCMS - pions_recoil;

h1_mJpsi_recoil->Fill(pJpsi_recoil.M());
// problem 4.2
double chi2 = dblKinFit4CChiSq;

if (chi2 <= 0) continue;

h1_chi2_4C->Fill(chi2);

// problem 4.1
// reconstruct J/psi
P4E Jpsi = lep_plus + lep_minus;
// total final state 
// problem 4.2
P4E pions = pi_plus + pi_minus;
P4E total = Jpsi + pions;

// missing energy
double Emiss = pCMS.E() - total.E();

double s = sqrt_s * sqrt_s;
double M = total.M();

double Egamma_ISR = (s - M*M) / (2 * sqrt_s);

h1_Egamma_ISR->Fill(Egamma_ISR);


// ---------------- BEFORE χ² CUT ----------------
h1_mJpsi_pipi_before->Fill(total.M());
h1_Emiss_before->Fill(Emiss);



// ---------------- APPLY χ² CUT ----------------
if (chi2 > 40)
    continue;

selectedEvents++;   


// ---------------- AFTER χ² CUT ----------------
h1_mJpsi_pipi->Fill(total.M());
h1_Emiss->Fill(Emiss);


// detected ISR photons
for (int i = 0; i < intNumberGoodPhotons; i++)
{
    double Egamma = dblPhotonsE[i];
    if (Egamma < 0.05) continue;   // remove very soft clusters
    h1_Egamma->Fill(Egamma);
}
// Dalitz plot using 4C fit four-vectors

double m2_pipi = pions.M2();

double m2_jpsipi1 = (Jpsi + pi_plus).M2();
double m2_jpsipi2 = (Jpsi + pi_minus).M2();

h2_dalitz->Fill(m2_pipi, m2_jpsipi1);
h2_dalitz->Fill(m2_pipi, m2_jpsipi2);
 
    // // problem 3.0
   if (electronEvent)
{
    h1_mJpsi_e->Fill(Jpsi.M());
}

if (muonEvent)
{
    h1_mJpsi_mu->Fill(Jpsi.M());
}// end of muon loop
} // end of event loop
        //=============================================================================
        // End of selection
        //=============================================================================
    
    // problem 2.4
    std::cout << "Total events: " << nentries << std::endl;
    std::cout << "Selected events: " << selectedEvents << std::endl;


    //=============================================================================
    // All canvas definitions should be placed here
    //=============================================================================

    
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

// problem 3.6
{
    TCanvas* canvas = new TCanvas("canvas","canvas",1200,900);
    canvas->SetLeftMargin(0.15);

    //  draw Dalitz histogram
    h2_dalitz->Draw("COLZ");
    gPad->Update();

    double s = sqrt_s * sqrt_s;

    TF1* dalitz_upper = new TF1("dalitz_upper",
    [=](double *x,double*)
    {
        double m2_pipi = x[0];

        double term1 = (s - m_jpsi*m_jpsi - m2_pipi)*(m2_pipi - 2*m_pi*m_pi);
        double term2 = sqrt(lambda(m2_pipi,m_pi*m_pi,m_pi*m_pi) *
                            lambda(s,m2_pipi,m_jpsi*m_jpsi));

        return m_jpsi*m_jpsi + m_pi*m_pi + (term1 + term2)/(2*m2_pipi);

    },0.08,1.4,0);

   // quick diagnostic
    std::cout << "Dalitz upper at x=0.5 : "
              << dalitz_upper->Eval(0.5)
              << std::endl;
    //  define lower boundary
    TF1* dalitz_lower = new TF1("dalitz_lower",
[=](double *x,double*)
{
    double m2_pipi = x[0];

    double term1 = (s - m_jpsi*m_jpsi - m2_pipi)*(m2_pipi - 2*m_pi*m_pi);
    double term2 = sqrt(lambda(m2_pipi,m_pi*m_pi,m_pi*m_pi) *
                        lambda(s,m2_pipi,m_jpsi*m_jpsi));

    return m_jpsi*m_jpsi + m_pi*m_pi + (term1 - term2)/(2*m2_pipi);

},0.08,1.4,0);

    dalitz_upper->SetLineColor(kRed);
    dalitz_upper->SetLineWidth(3);
    dalitz_lower->SetLineColor(kRed);
    dalitz_lower->SetLineWidth(3);

    

    dalitz_upper->Draw("same");
    dalitz_lower->Draw("same");

    canvas->SaveAs(savePath + "Dalitz_plot.png");
}
// problem 4.1
{
    TCanvas* canvas = new TCanvas();
    double binWidth = h1_mJpsi_pipi_before->GetBinWidth(1);

    h1_mJpsi_pipi_before->GetYaxis()->SetTitle(
        Form("Events / %.3f GeV", binWidth)
    );

    h1_mJpsi_pipi_before->Draw();

    canvas->SaveAs(savePath + "m_Jpsi_pipi_beforeCut.png");
}
// problem 4.2
{
    TCanvas* canvas = new TCanvas();
    double binWidth = h1_mJpsi_pipi->GetBinWidth(1);

    h1_mJpsi_pipi->GetYaxis()->SetTitle(
        Form("Events / %.3f GeV", binWidth)
    );

    h1_mJpsi_pipi->Draw();

    canvas->SaveAs(savePath + "m_Jpsi_pipi_afterCut.png");
}

{
    TCanvas* canvas = new TCanvas();
    double binWidth = h1_Egamma->GetBinWidth(1);
    h1_Egamma->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));
    h1_Egamma->Draw();
   
    canvas->SaveAs(savePath + "ISR_photon_energy.png");
}
{
    TCanvas* canvas = new TCanvas();

    double binWidth = h1_Egamma_ISR->GetBinWidth(1);
    h1_Egamma_ISR->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

    h1_Egamma_ISR->Draw();

    canvas->SaveAs(savePath + "ISR_photon_energy_kinematic.png");
}
{
    TCanvas* canvas = new TCanvas();
    double binWidth = h1_Emiss_before->GetBinWidth(1);

    h1_Emiss_before->GetYaxis()->SetTitle(
        Form("Events / %.3f GeV", binWidth)
    );

    h1_Emiss_before->Draw();

    canvas->SaveAs(savePath + "Missing_energy_beforeCut.png");
}
// problem 4.2 
{
    TCanvas* canvas = new TCanvas();
    double binWidth = h1_Emiss->GetBinWidth(1);

    h1_Emiss->GetYaxis()->SetTitle(
        Form("Events / %.3f GeV", binWidth)
    );

    h1_Emiss->Draw();

    canvas->SaveAs(savePath + "Missing_energy_afterCut.png");
}
// problem 4.2
{
    TCanvas* canvas = new TCanvas("canvas_chi2","canvas_chi2",1200,900);
    canvas->SetLeftMargin(0.15);
    double binWidth = h1_chi2_4C->GetBinWidth(1);
    h1_chi2_4C->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));

    h1_chi2_4C->SetLineColor(kBlue);
    h1_chi2_4C->SetLineWidth(2);

    h1_chi2_4C->Draw();

    canvas->SaveAs(savePath + "chi2_4C_distribution.png");
}
      

        
        
    }

