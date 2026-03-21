






TH1D* g_phsp = nullptr;
TH1D* g_sig  = nullptr;

// ===== BACKGROUND FUNCTION =====
double bg_func(double *x, double *p)
{
    return p[0]*g_phsp->Interpolate(x[0]);
}

// ===== SIGNAL FUNCTION =====
double sig_func(double *x, double *p)
{
    return p[0]*g_sig->Interpolate(x[0]);
}
void combine()
{
    // ===== LOAD FILES =====
    TFile* f_data = new TFile("plots_data/output.root");
    TFile* f_phsp = new TFile("plots_pippimJpsi/output.root");
    TFile* f_sig  = new TFile("plots_pipmZc3900/output.root");

    // ===== GET HISTOGRAMS =====
    TH1D* h_data = (TH1D*)f_data->Get("h_mJpsipi_max_final")->Clone("h_data");
    TH1D* h_phsp = (TH1D*)f_phsp->Get("h_mJpsipi_max_final")->Clone("h_phsp");
    TH1D* h_sig  = (TH1D*)f_sig->Get("h_mJpsipi_max_final")->Clone("h_sig");

    h_data->Rebin(2);
    h_phsp->Rebin(2);
    h_sig->Rebin(2);

    h_data->SetDirectory(0);
    h_phsp->SetDirectory(0);
    h_sig->SetDirectory(0);

    // assign globals (IMPORTANT)
    g_phsp = h_phsp;
    g_sig  = h_sig;

    // ===== SMOOTH MC =====
    h_phsp->Smooth(1);
    h_sig->Smooth(1);

    // ===== NORMALIZE MC =====
   // Normalize to unit area ONLY
    h_phsp->Scale(1.0 / h_phsp->Integral());
    h_sig->Scale(1.0 / h_sig->Integral());

    // ===== FIT FUNCTION =====
    TF1* fit = new TF1("fit",
    [](double *x, double *p)
    {
        double phsp = g_phsp->Interpolate(x[0]);
        double sig  = g_sig->Interpolate(x[0]);

        return p[0]*phsp + p[1]*sig;
    },
    3.6, 4.0, 4
);

    fit->SetParameters(1000, 500);
    fit->SetParLimits(0, 0, 1e6);
    fit->SetParLimits(1, 0, 1e6);
    // ===== FIT =====
    h_data->Fit(fit, "R");

    double a = fit->GetParameter(0);
    double b = fit->GetParameter(1);

    std::cout << "\n=== FIT RESULTS ===" << std::endl;
    std::cout << "a (PHSP) = " << a << std::endl;
    std::cout << "b (Zc)   = " << b << std::endl;


    // ===== COMPONENT CURVES =====
    TF1* f_bg = new TF1("f_bg", bg_func, 3.65, 4.0, 1);

    f_bg->SetParameter(0, a);
    f_bg->SetLineColor(kBlue);
    f_bg->SetLineStyle(2);
    f_bg->SetLineWidth(3);

    TF1* f_sig_func = new TF1("f_sig_func", sig_func, 3.65, 4.0, 1);
    f_sig_func->SetParameter(0, b);
    f_sig_func->SetLineColor(kRed);
    f_sig_func->SetLineStyle(3);
    f_sig_func->SetLineWidth(2);

    // ===== PLOT =====
    TCanvas* c = new TCanvas("c","c",1200,900);

    h_data->SetMaximum(1.4 * h_data->GetMaximum());
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.2);

    h_data->Draw("E");

    f_bg->Draw("SAME");

    fit->SetLineColor(kRed);
    fit->SetLineWidth(3);
    fit->Draw("SAME");

    f_sig_func->Draw("SAME");

    // ===== LEGEND =====
    TLegend* leg = new TLegend(0.55,0.65,0.85,0.85);
    leg->AddEntry(h_data,"Data","lep");
    leg->AddEntry(fit,"Total fit","l");
    leg->AddEntry(f_bg,"Background (PHSP)","l");
    leg->AddEntry(f_sig_func,"Zc signal","l");
    leg->Draw();

    c->SaveAs("Zc_combined_fit.png");

    // ===== Zc YIELD =====
    double N_Zc = b;

    std::cout << "\n=== Zc EVENTS FROM FIT ===" << std::endl;
    std::cout << "N_Zc = " << N_Zc << std::endl;
}