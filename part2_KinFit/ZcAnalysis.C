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

// problem 4.4

TH1D* h1_p_e = new TH1D(
"h1_p_e",
"Electron momentum;|p_{e}| [GeV];Events",
200,0,3
);

TH1D* h1_p_mu = new TH1D(
"h1_p_mu",
"Muon momentum;|p_{#mu}| [GeV];Events",
200,0,3
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
TH1D* h1_mJpsi_pipi_afterBhabha = new TH1D(
    "h1_mJpsi_pipi_afterBhabha",
    "Invariant mass after Bhabha cut;M_{J/#psi#pi#pi} [GeV];Events",
    200,3.0,4.5
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
// problem 4.5
TH1D* h_cos_ee = new TH1D(
"cos_ee",
"cos(#theta) between e^{+}e^{-};cos(#theta);Events",
100,-1,1
);

TH1D* h_cos_epi_minus = new TH1D(
"cos_epi_minus",
"cos(#theta) between e^{+}#pi^{-};cos(#theta);Events",
100,-1,1
);

TH1D* h_cos_pie_minus = new TH1D(
"cos_pie_minus",
"cos(#theta) between #pi^{+}e^{-};cos(#theta);Events",
100,-1,1
);

TH1D* h_cos_pipi = new TH1D(
"cos_pipi",
"cos(#theta) between #pi^{+}#pi^{-};cos(#theta);Events",
100,-1,1
);
// problem 4.5
TH1D* h1_mJpsi_pipi_afterConversion = new TH1D(
"h1_mJpsi_pipi_afterConversion",
"Invariant mass after conversion cut;M_{J/#psi#pi#pi} [GeV];Events",
200,3.0,4.5
);
// J/psi mass spectra BEFORE cuts
TH1D* h_mJpsi_e_before = new TH1D("h_mJpsi_e_before",
"J/#psi mass e^{+}e^{-} before cuts;Mass [GeV];Events",100,2.9,3.2);

TH1D* h_mJpsi_mu_before = new TH1D("h_mJpsi_mu_before",
"J/#psi mass #mu^{+}#mu^{-} before cuts;Mass [GeV];Events",100,2.9,3.2);

TH1D* h_mJpsi_recoil_before = new TH1D("h_mJpsi_recoil_before",
"J/#psi recoil mass before cuts;Mass [GeV];Events",100,2.9,3.2);


// AFTER BHABHA CUT
TH1D* h_mJpsi_e_afterBhabha = new TH1D("h_mJpsi_e_afterBhabha",
"J/#psi mass e^{+}e^{-} after Bhabha cut;Mass [GeV];Events",100,2.9,3.2);

TH1D* h_mJpsi_mu_afterBhabha = new TH1D("h_mJpsi_mu_afterBhabha",
"J/#psi mass #mu^{+}#mu^{-} after Bhabha cut;Mass [GeV];Events",100,2.9,3.2);

TH1D* h_mJpsi_recoil_afterBhabha = new TH1D("h_mJpsi_recoil_afterBhabha",
"J/#psi recoil mass after Bhabha cut;Mass [GeV];Events",100,2.9,3.2);
// After conversion cut
TH1D* h_mJpsi_e_afterConversion = new TH1D(
"h_mJpsi_e_afterConversion",
"J/#psi mass e^{+}e^{-} after conversion cut;Mass [GeV];Events",
100,2.9,3.2
);

TH1D* h_mJpsi_mu_afterConversion = new TH1D(
"h_mJpsi_mu_afterConversion",
"J/#psi mass #mu^{+}#mu^{-} after conversion cut;Mass [GeV];Events",
100,2.9,3.2
);

TH1D* h_mJpsi_recoil_afterConversion = new TH1D(
"h_mJpsi_recoil_afterConversion",
"J/#psi recoil mass after conversion cut;Mass [GeV];Events",
100,2.9,3.2
);
// After chi^2 cut
TH1D* h_mJpsi_e_afterChi2 = new TH1D(
"h_mJpsi_e_afterChi2",
"J/#psi invariant mass after all selection cuts (e^{+}e^{-});Mass [GeV];Events",
100,2.9,3.2
);

TH1D* h_mJpsi_mu_afterChi2 = new TH1D(
"h_mJpsi_mu_afterChi2",
"J/#psi invariant mass after all selection cuts (#mu^{+}#mu^{-});Mass [GeV];Events",
100,2.9,3.2
);

TH1D* h_mJpsi_recoil_afterChi2 = new TH1D(
"h_mJpsi_recoil_afterChi2",
"J/#psi invariant mass after all selection cuts (recoil);Mass [GeV];Events",
100,2.9,3.2
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

// problem 4.6
int Ni_e = 0, Ni_mu = 0;

int N_bhabha_e = 0, N_bhabha_mu = 0;
int N_conv_e   = 0, N_conv_mu   = 0;
int N_chi2_e   = 0, N_chi2_mu   = 0;
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

// problem 4.4
double p_lp = lep_plus.P();
double p_lm = lep_minus.P();

// problem 4.5 photon conversion study (angular distributions)

double cos_ee = CosTheta(lep_plus, lep_minus);
double cos_epi_minus = CosTheta(lep_plus, pi_minus);
double cos_pie_minus = CosTheta(pi_plus, lep_minus);
double cos_pipi = CosTheta(pi_plus, pi_minus);

        //=============================================================================
        // Your selection for each event
        //=============================================================================

    
        // Loop over all charged tracks in the event (there are always 4 charged tracks: l+ l- pi+ pi-)

        // problem 3.5
        // Fill invariant mass histograms


bool electronEvent = (intNumberElectrons == 2);
bool muonEvent     = (intNumberMuons == 2);
// problem 4.6
if (electronEvent) Ni_e++;
if (muonEvent)     Ni_mu++;

// problem 4.5
// fill angular histograms
h_cos_ee->Fill(cos_ee);
h_cos_epi_minus->Fill(cos_epi_minus);
h_cos_pie_minus->Fill(cos_pie_minus);
h_cos_pipi->Fill(cos_pipi);

// reconstruct J/psi and total system BEFORE Bhabha cut
P4E Jpsi_temp  = lep_plus + lep_minus;
P4E pions_temp = pi_plus + pi_minus;
// BEFORE CUTS
if(electronEvent)
    h_mJpsi_e_before->Fill(Jpsi_temp.M());

if(muonEvent)
    h_mJpsi_mu_before->Fill(Jpsi_temp.M());

P4E recoil_temp = pCMS - pions_temp;
h_mJpsi_recoil_before->Fill(recoil_temp.M());
P4E total_temp = Jpsi_temp + pions_temp;

// fill invariant mass BEFORE Bhabha cut
h1_mJpsi_pipi_before->Fill(total_temp.M());
// problem 4.4
if (electronEvent)
{
    h1_p_e->Fill(p_lp);
    h1_p_e->Fill(p_lm);
}

if (muonEvent)
{
    h1_p_mu->Fill(p_lp);
    h1_p_mu->Fill(p_lm);
}
// ---------------- Bhabha background rejection ----------------

    if (p_lp > 2.0 || p_lm > 2.0)
        continue;
    // problem 4.6
    if (electronEvent) N_bhabha_e++;
    if (muonEvent)     N_bhabha_mu++;
// AFTER BHABHA CUT
if(electronEvent)
    h_mJpsi_e_afterBhabha->Fill(Jpsi_temp.M());

if(muonEvent)
    h_mJpsi_mu_afterBhabha->Fill(Jpsi_temp.M());

P4E recoil_bhabha = pCMS - pions_temp;

if (electronEvent || muonEvent)
    h_mJpsi_recoil_afterBhabha->Fill(recoil_bhabha.M());
// invariant mass after bhabha cut
h1_mJpsi_pipi_afterBhabha->Fill(total_temp.M());
// problem 3.3
// reconstruct J/psi from pion recoil


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
// problem 4.5
// ---- photon conversion rejection ----
bool passConversion = true;

if (electronEvent)
{
    if (cos_epi_minus > 0.98 || cos_pie_minus > 0.98)
        passConversion = false;
}

if (!passConversion) continue;

// problem 4.6
// electrons → affected by conversion cut
if (electronEvent) N_conv_e++;

// muons → unchanged from Bhabha stage
if (muonEvent) N_conv_mu++;
// AFTER CONVERSION CUT

if(electronEvent)
    h_mJpsi_e_afterConversion->Fill(Jpsi.M());

if(muonEvent)
    h_mJpsi_mu_afterConversion->Fill(Jpsi.M());

P4E recoil_conv = pCMS - pions;

if (electronEvent || muonEvent)
    h_mJpsi_recoil_afterConversion->Fill(recoil_conv.M());
// fill invariant mass after conversion cut
h1_mJpsi_pipi_afterConversion->Fill(total.M());

// missing energy
double Emiss = pCMS.E() - total.E();

double s = sqrt_s * sqrt_s;
double M = total.M();

double Egamma_ISR = (s - M*M) / (2 * sqrt_s);

h1_Egamma_ISR->Fill(Egamma_ISR);


// ---------------- BEFORE χ² CUT ----------------
h1_Emiss_before->Fill(Emiss);



// ---------------- APPLY χ² CUT ----------------
if (chi2 > 40)
    continue;
// problem 4.6
if (electronEvent) N_chi2_e++;
if (muonEvent)     N_chi2_mu++;
// AFTER χ² CUT

if(electronEvent)
    h_mJpsi_e_afterChi2->Fill(Jpsi.M());

if(muonEvent)
    h_mJpsi_mu_afterChi2->Fill(Jpsi.M());

P4E recoil_chi2 = pCMS - pions;

if (electronEvent || muonEvent)
    h_mJpsi_recoil_afterChi2->Fill(recoil_chi2.M());
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
// problem 5.1
TH1D* h_mJpsi_combined = (TH1D*)h_mJpsi_e_afterChi2->Clone("h_mJpsi_combined");
h_mJpsi_combined->Add(h_mJpsi_mu_afterChi2);
        //=============================================================================
        // End of selection
        //=============================================================================
    
    // problem 2.4
    std::cout << "Total events: " << nentries << std::endl;
    std::cout << "Selected events: " << selectedEvents << std::endl;
// problem 4.6
std::cout << "\n=== ELECTRONS ===" << std::endl;
std::cout << "Initial: " << Ni_e << std::endl;
std::cout << "After Bhabha: " << N_bhabha_e 
          << "  eff = " << (double)N_bhabha_e/Ni_e << std::endl;
std::cout << "After Conversion: " << N_conv_e 
          << "  eff = " << (double)N_conv_e/Ni_e << std::endl;
std::cout << "After Chi2: " << N_chi2_e 
          << "  eff = " << (double)N_chi2_e/Ni_e << std::endl;

std::cout << "\n=== MUONS ===" << std::endl;
std::cout << "Initial: " << Ni_mu << std::endl;
std::cout << "After Bhabha: " << N_bhabha_mu 
          << "  eff = " << (double)N_bhabha_mu/Ni_mu << std::endl;
std::cout << "After Conversion: " << N_conv_mu
          << "  eff = " << (double)N_bhabha_mu/Ni_mu << std::endl;
std::cout << "After Chi2: " << N_chi2_mu 
          << "  eff = " << (double)N_chi2_mu/Ni_mu << std::endl;
double Ni_gen_e = 500623.0;
double Ni_gen_mu = 499377.0;

// efficiencies
double eff_e = (double)N_chi2_e / Ni_gen_e;
double eff_mu = (double)N_chi2_mu / Ni_gen_mu;

// uncertainties (binomial)
double err_e = sqrt(eff_e * (1 - eff_e) / Ni_gen_e);
double err_mu = sqrt(eff_mu * (1 - eff_mu) / Ni_gen_mu);

// print
std::cout << "\n=== FINAL EFFICIENCY (MC truth) ===" << std::endl;

std::cout << "Electrons: " 
          << eff_e << " ± " << err_e << std::endl;

std::cout << "Muons: " 
          << eff_mu << " ± " << err_mu << std::endl;

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

    canvas->SaveAs(savePath + "m_Jpsi_pipi_after_chi^2_Cut.png");
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
// problem 4.4
{
TCanvas* canvas = new TCanvas();
double binWidth = h1_p_e->GetBinWidth(1);
h1_p_e->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));
h1_p_e->Draw();
canvas->SaveAs(savePath + "electron_momentum.png");
}

{
TCanvas* canvas = new TCanvas();
double binWidth = h1_p_mu->GetBinWidth(1);
h1_p_mu->GetYaxis()->SetTitle(Form("Events / %.3f GeV", binWidth));
h1_p_mu->Draw();
canvas->SaveAs(savePath + "muon_momentum.png");
}
{
    TCanvas* canvas = new TCanvas();
    double binWidth = h1_mJpsi_pipi_afterBhabha->GetBinWidth(1);

    h1_mJpsi_pipi_afterBhabha->GetYaxis()->SetTitle(
        Form("Events / %.3f GeV", binWidth)
    );

    h1_mJpsi_pipi_afterBhabha->Draw();

    canvas->SaveAs(savePath + "m_Jpsi_pipi_afterBhabhaCut.png");
}
// problem 4.5
// photon conversion angular distributions

{
TCanvas* canvas = new TCanvas();
h_cos_ee->Draw();
canvas->SaveAs(savePath + "cosTheta_ee.png");
}

{
TCanvas* canvas = new TCanvas();
h_cos_epi_minus->Draw();
canvas->SaveAs(savePath + "cosTheta_e_pi_minus.png");
}

{
TCanvas* canvas = new TCanvas();
h_cos_pie_minus->Draw();
canvas->SaveAs(savePath + "cosTheta_pi_plus_e_minus.png");
}

{
TCanvas* canvas = new TCanvas();
h_cos_pipi->Draw();
canvas->SaveAs(savePath + "cosTheta_pipi.png");
}
{
TCanvas* canvas = new TCanvas();
double binWidth = h1_mJpsi_pipi_afterConversion->GetBinWidth(1);

h1_mJpsi_pipi_afterConversion->GetYaxis()->SetTitle(
Form("Events / %.3f GeV", binWidth)
);

h1_mJpsi_pipi_afterConversion->Draw();

canvas->SaveAs(savePath + "m_Jpsi_pipi_afterConversionCut.png");
}
{
TCanvas* canvas = new TCanvas("canvas_mass","canvas_mass",1200,900);

h_mJpsi_e_afterBhabha->SetLineColor(kRed);
h_mJpsi_mu_afterBhabha->SetLineColor(kBlue);
h_mJpsi_recoil_afterBhabha->SetLineColor(kMagenta);

double max_e_b = h_mJpsi_e_afterBhabha->GetMaximum();
double max_mu_b = h_mJpsi_mu_afterBhabha->GetMaximum();
double max_recoil_b = h_mJpsi_recoil_afterBhabha->GetMaximum();

double max_all_b = std::max({max_e_b, max_mu_b, max_recoil_b});

// add margin
h_mJpsi_e_afterBhabha->SetMaximum(1.2 * max_all_b);

h_mJpsi_e_afterBhabha->Draw();
h_mJpsi_mu_afterBhabha->Draw("same");
h_mJpsi_recoil_afterBhabha->Draw("same");

TLegend* leg = new TLegend(0.15,0.7,0.35,0.85);
leg->AddEntry(h_mJpsi_e_afterBhabha,"e^{+}e^{-}","l");
leg->AddEntry(h_mJpsi_mu_afterBhabha,"#mu^{+}#mu^{-}","l");
leg->AddEntry(h_mJpsi_recoil_afterBhabha,"Recoil(#pi^{+}#pi^{-})","l");
leg->Draw();

canvas->SaveAs(savePath + "Jpsi_mass_spectrum_afterBhabha.png");
}
// after conversion cut
{
TCanvas* canvas = new TCanvas();

h_mJpsi_e_afterConversion->SetLineColor(kRed);
h_mJpsi_mu_afterConversion->SetLineColor(kBlue);
h_mJpsi_recoil_afterConversion->SetLineColor(kMagenta);

double max_e_c = h_mJpsi_e_afterConversion->GetMaximum();
double max_mu_c = h_mJpsi_mu_afterConversion->GetMaximum();
double max_recoil_c = h_mJpsi_recoil_afterConversion->GetMaximum();

double max_all_c = std::max({max_e_c, max_mu_c, max_recoil_c});

// add margin
h_mJpsi_e_afterConversion->SetMaximum(1.2 * max_all_c);

h_mJpsi_e_afterConversion->Draw();
h_mJpsi_mu_afterConversion->Draw("same");
h_mJpsi_recoil_afterConversion->Draw("same");

TLegend* leg = new TLegend(0.15,0.7,0.35,0.85);
leg->AddEntry(h_mJpsi_e_afterConversion,"e^{+}e^{-}","l");
leg->AddEntry(h_mJpsi_mu_afterConversion,"#mu^{+}#mu^{-}","l");
leg->AddEntry(h_mJpsi_recoil_afterConversion,"Recoil(#pi^{+}#pi^{-})","l");
leg->Draw();

canvas->SaveAs(savePath + "Jpsi_mass_afterConversion.png");
// after chi^2 cut
}
std::cout << "\n=== INTEGRALS AFTER CHI2 ===" << std::endl;

std::cout << "Electrons: "
          << h_mJpsi_e_afterChi2->Integral() << std::endl;

std::cout << "Muons: "
          << h_mJpsi_mu_afterChi2->Integral() << std::endl;

std::cout << "Recoil: "
          << h_mJpsi_recoil_afterChi2->Integral() << std::endl;
{

TCanvas* canvas_chi2 = new TCanvas("canvas_chi2_mass","canvas_chi2_mass",1200,900);
canvas_chi2->SetLeftMargin(0.15);

// set colors
h_mJpsi_e_afterChi2->SetLineColor(kRed);
h_mJpsi_mu_afterChi2->SetLineColor(kBlue);
h_mJpsi_recoil_afterChi2->SetLineColor(kMagenta);

h_mJpsi_e_afterChi2->SetLineWidth(2);
h_mJpsi_mu_afterChi2->SetLineWidth(2);
h_mJpsi_recoil_afterChi2->SetLineWidth(3);

// set y-axis title
double binWidth = h_mJpsi_e_afterChi2->GetBinWidth(1);
h_mJpsi_e_afterChi2->GetYaxis()->SetTitle(
Form("Events / %.3f GeV", binWidth)
);

// draw spectra
double max_e = h_mJpsi_e_afterChi2->GetMaximum();
double max_mu = h_mJpsi_mu_afterChi2->GetMaximum();
double max_recoil = h_mJpsi_recoil_afterChi2->GetMaximum();

double max_all = std::max({max_e, max_mu, max_recoil});

// add some margin (20%)
h_mJpsi_e_afterChi2->SetMaximum(1.2 * max_all);

h_mJpsi_e_afterChi2->Draw();
h_mJpsi_mu_afterChi2->Draw("same");
h_mJpsi_recoil_afterChi2->Draw("same");
// legend
TLegend* leg_chi2 = new TLegend(0.15,0.70,0.35,0.85);
leg_chi2->AddEntry(h_mJpsi_e_afterChi2,"e^{+}e^{-}","l");
leg_chi2->AddEntry(h_mJpsi_mu_afterChi2,"#mu^{+}#mu^{-}","l");
leg_chi2->AddEntry(h_mJpsi_recoil_afterChi2,"Recoil(#pi^{+}#pi^{-})","l");
leg_chi2->Draw();

canvas_chi2->SaveAs(savePath + "Jpsi_mass_after_all_cuts.png");
}

{
TCanvas* c = new TCanvas("c_jpsi_e_voigt","c_jpsi_e_voigt",1200,900);

double binWidth = h_mJpsi_e_afterChi2->GetBinWidth(1);

h_mJpsi_e_afterChi2->Draw();

TF1* ffit_e_voigt = new TF1(
    "ffit_e_voigt",
    Form("[0]*TMath::Voigt(x-[1],[2],[3])*%f + [4]", binWidth),
    3.05, 3.15
);

ffit_e_voigt->SetParameters(500, 3.097, 0.008, 0.000093, 10);
ffit_e_voigt->FixParameter(3, 0.000093);

h_mJpsi_e_afterChi2->Fit(ffit_e_voigt, "R");

ffit_e_voigt->SetRange(2.9, 3.2);
ffit_e_voigt->SetLineColor(kRed);
ffit_e_voigt->Draw("same");

c->SaveAs(savePath + "Jpsi_e_voigt.png");
}
{
TCanvas* c = new TCanvas("c_jpsi_mu_voigt","c_jpsi_mu_voigt",1200,900);

double binWidth = h_mJpsi_mu_afterChi2->GetBinWidth(1);

h_mJpsi_mu_afterChi2->Draw();

TF1* ffit_mu_voigt = new TF1(
    "ffit_mu_voigt",
    Form("[0]*TMath::Voigt(x-[1],[2],[3])*%f + [4]", binWidth),
    3.05, 3.15
);

ffit_mu_voigt->SetParameters(500, 3.097, 0.006, 0.000093, 10);
ffit_mu_voigt->FixParameter(3, 0.000093);

h_mJpsi_mu_afterChi2->Fit(ffit_mu_voigt, "R");

ffit_mu_voigt->SetRange(2.9, 3.2);
ffit_mu_voigt->SetLineColor(kBlue);
ffit_mu_voigt->Draw("same");

c->SaveAs(savePath + "Jpsi_mu_voigt.png");
}
{
TCanvas* c = new TCanvas("c_jpsi_recoil_voigt","c_jpsi_recoil_voigt",1200,900);

double binWidth = h_mJpsi_recoil_afterChi2->GetBinWidth(1);

h_mJpsi_recoil_afterChi2->Draw();

TF1* ffit_recoil_voigt = new TF1(
    "ffit_recoil_voigt",
    Form("[0]*TMath::Voigt(x-[1],[2],[3])*%f + [4]", binWidth),
    3.05, 3.15
);

ffit_recoil_voigt->SetParameters(500, 3.097, 0.010, 0.000093, 20);
ffit_recoil_voigt->FixParameter(3, 0.000093);

h_mJpsi_recoil_afterChi2->Fit(ffit_recoil_voigt, "R");

ffit_recoil_voigt->SetRange(2.9, 3.2);
ffit_recoil_voigt->SetLineColor(kMagenta);
ffit_recoil_voigt->Draw("same");

c->SaveAs(savePath + "Jpsi_recoil_voigt.png");
}
{
TCanvas* c = new TCanvas("c_jpsi_combined_voigt","c_jpsi_combined_voigt",1200,900);

double binWidth = h_mJpsi_combined->GetBinWidth(1);

h_mJpsi_combined->SetLineColor(kBlack);
h_mJpsi_combined->Draw();

TF1* ffit_combined_voigt = new TF1(
    "ffit_combined_voigt",
    Form("[0]*TMath::Voigt(x-[1],[2],[3])*%f + [4]", binWidth),
    3.05, 3.15
);

// slightly higher normalization since stats are larger
ffit_combined_voigt->SetParameters(1000, 3.097, 0.007, 0.000093, 10);
ffit_combined_voigt->FixParameter(3, 0.000093);

h_mJpsi_combined->Fit(ffit_combined_voigt, "R");

// extend drawing
ffit_combined_voigt->SetRange(2.9, 3.2);
ffit_combined_voigt->SetLineColor(kRed);
ffit_combined_voigt->SetLineWidth(2);
ffit_combined_voigt->Draw("same");

c->SaveAs(savePath + "Jpsi_combined_voigt.png");
}
{
int b1e = h_mJpsi_e_afterChi2->FindBin(3.05);
int b2e = h_mJpsi_e_afterChi2->FindBin(3.15);
double N_e = h_mJpsi_e_afterChi2->Integral(b1e, b2e);

int b1m = h_mJpsi_mu_afterChi2->FindBin(3.05);
int b2m = h_mJpsi_mu_afterChi2->FindBin(3.15);
double N_mu = h_mJpsi_mu_afterChi2->Integral(b1m, b2m);

int b1r = h_mJpsi_recoil_afterChi2->FindBin(3.05);
int b2r = h_mJpsi_recoil_afterChi2->FindBin(3.15);
double N_recoil = h_mJpsi_recoil_afterChi2->Integral(b1r, b2r);

int b1c = h_mJpsi_combined->FindBin(3.05);
int b2c = h_mJpsi_combined->FindBin(3.15);
double N_combined = h_mJpsi_combined->Integral(b1c, b2c);

std::cout << "\n=== FINAL EVENT YIELDS ===" << std::endl;
std::cout << "Electrons: " << N_e << std::endl;
std::cout << "Muons: " << N_mu << std::endl;
std::cout << "Recoil: " << N_recoil << std::endl;
std::cout << "Combined: " << N_combined << std::endl;
}
// problem 5.1
{
TCanvas* c = new TCanvas("c_jpsi_final","c_jpsi_final",1200,900);
c->SetLeftMargin(0.15);

// axis
double binWidth = h_mJpsi_combined->GetBinWidth(1);
h_mJpsi_combined->GetYaxis()->SetTitle(
Form("Events / %.3f GeV", binWidth)
);

// draw
h_mJpsi_combined->SetLineColor(kBlack);
h_mJpsi_combined->Draw();

// ===== FIT FUNCTION =====
// ===== VOIGT FIT FUNCTION =====
TF1* ffit_voigt = new TF1(
    "ffit_voigt",
    Form("[0]*TMath::Voigt(x-[1],[2],[3])*%f + [4]", binWidth),
    3.05, 3.15
);

// Initial parameters
ffit_voigt->SetParameters(
    1000,      // normalization (adjust if needed)
    3.097,     // mean (J/psi mass)
    0.008,     // sigma (~8 MeV resolution)
    0.000093,  // gamma (J/psi natural width)
    5          // background
);

// FIX gamma (very important)
ffit_voigt->FixParameter(3, 0.000093);

// Fit
h_mJpsi_combined->Fit(ffit_voigt, "R");

ffit_voigt->SetRange(2.9, 3.2);

// Draw
ffit_voigt->SetLineColor(kRed);
ffit_voigt->SetLineWidth(2);
ffit_voigt->Draw("same");

// ===== EXTRACT RESULTS =====
double mass  = ffit_voigt->GetParameter(1);
double sigma = ffit_voigt->GetParameter(2);

std::cout << "\n=== FINAL J/psi FIT ===" << std::endl;
std::cout << "Mass = " << mass << std::endl;
std::cout << "Width (sigma) = " << sigma << std::endl;

// save
c->SaveAs(savePath + "Jpsi_mass_final_combined_fit.png");
}
// ================= SIGNAL EXTRACTION FROM HISTOGRAM =================

int b1e = h_mJpsi_e_afterChi2->FindBin(3.05);
int b2e = h_mJpsi_e_afterChi2->FindBin(3.15);
double N_Jpsi_e = h_mJpsi_e_afterChi2->Integral(b1e, b2e);

int b1m = h_mJpsi_mu_afterChi2->FindBin(3.05);
int b2m = h_mJpsi_mu_afterChi2->FindBin(3.15);
double N_Jpsi_mu = h_mJpsi_mu_afterChi2->Integral(b1m, b2m);

std::cout << "N_Jpsi_e = " << N_Jpsi_e << std::endl;
std::cout << "N_Jpsi_mu = " << N_Jpsi_mu << std::endl;
// ================= CROSS SECTION =================

double L = 828.4; // pb^-1
double delta = 0.818;
double B_e = 0.0594;
double B_mu = 0.0594;

double sigma_e = N_Jpsi_e / (L * delta * eff_e * B_e);
double sigma_mu = N_Jpsi_mu / (L * delta * eff_mu * B_mu);

// ===== UNCERTAINTY (STATISTICAL) =====
double err_N_e = sqrt(N_Jpsi_e);
double err_N_mu = sqrt(N_Jpsi_mu);

double err_sigma_e = sigma_e * (err_N_e / N_Jpsi_e);
double err_sigma_mu = sigma_mu * (err_N_mu / N_Jpsi_mu);

// ===== PRINT RESULTS =====
std::cout << "\n=== CROSS SECTIONS ===" << std::endl;

std::cout << "sigma (e+e-) = " 
          << sigma_e << " ± " << err_sigma_e 
          << " pb" << std::endl;

std::cout << "sigma (mu+mu-) = " 
          << sigma_mu << " ± " << err_sigma_mu 
          << " pb" << std::endl;
} 

      

