






TH1D* g_phsp = nullptr;
TH1D* g_sig  = nullptr;

// ===== BACKGROUND FUNCTION =====
double bg_func(double *x, double *p)
{
    int bin = g_phsp->FindBin(x[0]);
    double val = g_phsp->GetBinContent(bin);
    double dx = x[0] - 3.8;
    return p[0]*val*(1 + p[1]*dx + p[2]*dx*dx);
}

double sig_func(double *x, double *p)
{
    int bin = g_sig->FindBin(x[0]);
    double val = g_sig->GetBinContent(bin);
    return p[0]*val;
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
    TH1D* h_sb = (TH1D*)f_data->Get("h_mJpsipi_max_sb")->Clone("h_sb");

    h_data->Rebin(2);
    h_phsp->Rebin(2);
    h_sig->Rebin(2);
    h_sb->Rebin(2);

    h_data->SetDirectory(0);
    h_phsp->SetDirectory(0);
    h_sig->SetDirectory(0);
    h_sb->SetDirectory(0);

    // assign globals (IMPORTANT)
    g_phsp = h_phsp;
    g_sig  = h_sig;

    // ===== SMOOTH MC =====
    h_phsp->Smooth(1);
    h_sig->Smooth(1);

    double data_int = h_data->Integral();

    h_phsp->Scale(data_int / h_phsp->Integral());
    h_sig->Scale(data_int / h_sig->Integral());

    // ===== NORMALIZE MC =====
   // Normalize to unit area ONLY
    //h_phsp->Scale(1.0 / h_phsp->Integral());
    //h_sig->Scale(1.0 / h_sig->Integral());
    
    // ===== FIT FUNCTION =====
   TF1* fit = new TF1("fit",
[](double *x, double *p)
{
    int bin_phsp = g_phsp->FindBin(x[0]);
    double phsp = g_phsp->GetBinContent(bin_phsp);

    int bin_sig = g_sig->FindBin(x[0]);
    double sig  = g_sig->GetBinContent(bin_sig);

    // linear tilt added here
    double dx = x[0] - 3.8;
    return p[0]*phsp*(1 + p[2]*dx + p[3]*dx*dx) + p[1]*sig;
},
3.7, 4.0, 4
);

    fit->SetParameters(1, 0.2, 0.0, 0.0);
    fit->SetParLimits(2, -5, 5);
    fit->SetParLimits(3, -5, 5);
    // ===== FIT =====
    h_data->Fit(fit, "R");

    double a = fit->GetParameter(0);
    double b = fit->GetParameter(1);

    double sb_scale = 0.04 / 0.12;
    h_sb->Scale(sb_scale);


    // style
    h_phsp->SetLineColor(kRed);
    h_phsp->SetLineStyle(2);
    h_phsp->SetLineWidth(2);

    h_sb->SetFillColor(kGreen+1);
    h_sb->SetFillStyle(3001);
    h_sb->SetLineColor(kGreen+3);

    std::cout << "\n=== FIT RESULTS ===" << std::endl;
    std::cout << "a (PHSP) = " << a << std::endl;
    std::cout << "b (Zc)   = " << b << std::endl;


    // ===== COMPONENT CURVES =====
    TF1* f_bg = new TF1("f_bg", bg_func, 3.65, 4.0, 3);

    f_bg->SetParameters(a, fit->GetParameter(2), fit->GetParameter(3));
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

    // subtle grid (publication style)
    c->SetGridx();
    c->SetGridy();

    // axis formatting
    h_data->GetXaxis()->SetTitle("M_{max}(J/#psi#pi) [GeV/c^{2}]");
    h_data->GetYaxis()->SetTitle("Events / bin");

    h_data->GetXaxis()->SetTitleSize(0.045);
    h_data->GetYaxis()->SetTitleSize(0.045);

    h_data->GetXaxis()->SetLabelSize(0.04);
    h_data->GetYaxis()->SetLabelSize(0.04);

    h_data->GetYaxis()->SetTitleOffset(1.2);

    // cleaner scaling
    h_data->SetMaximum(1.25 * h_data->GetMaximum());

    h_data->GetXaxis()->SetRangeUser(3.6, 4.1);

    // data style
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.0);
    h_data->SetLineWidth(2);

    // sideband (lighter)
    h_sb->SetFillColorAlpha(kGreen+2, 0.35);
    h_sb->SetLineColor(kGreen+3);

    // ===== DRAW ORDER =====
    h_data->Draw("E");

    h_sb->Draw("HIST SAME");

    // ===== PHSP FITTED =====
    TH1D* h_phsp_fit = (TH1D*)h_phsp->Clone("h_phsp_fit");

    for (int i = 1; i <= h_phsp_fit->GetNbinsX(); i++)
    {
        double x = h_phsp_fit->GetBinCenter(i);
        double dx = x - 3.8;

        double val = g_phsp->GetBinContent(i);
        double fitted = a * val * (1 + fit->GetParameter(2)*dx + fit->GetParameter(3)*dx*dx);

        h_phsp_fit->SetBinContent(i, fitted);
    }
    // background
    f_bg->SetLineColor(kRed);
    f_bg->SetLineStyle(2);
    f_bg->SetLineWidth(2);
    f_bg->Draw("SAME");

    // PHSP = subtle reference
    h_phsp_fit->SetLineColor(kGray+1);
    h_phsp_fit->SetLineStyle(2);   // simple dashed
    h_phsp_fit->SetLineWidth(2);
    h_phsp_fit->Draw("HIST SAME");


    // signal
    f_sig_func->SetLineColor(kMagenta+2);
    f_sig_func->SetLineStyle(3);
    f_sig_func->SetLineWidth(2);
    f_sig_func->Draw("SAME");

    // total fit (dominant)
    fit->SetLineColor(kRed);
    fit->SetLineWidth(3);
    fit->Draw("SAME");

    h_data->Draw("E SAME");

    // ===== LEGEND =====
    TLegend* leg = new TLegend(0.65,0.55,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    leg->AddEntry(h_data,"Data","lep");
    leg->AddEntry(fit,"Total fit","l");
    leg->AddEntry(f_bg,"Background fit","l");
    leg->AddEntry(f_sig_func,"Signal","l");
    leg->AddEntry(h_phsp_fit,"PHSP MC (fitted)","l");
    leg->AddEntry(h_sb,"Sideband","f");

    leg->Draw();

    c->SaveAs("Zc_combined_fit.png");

    // ===== Zc YIELD =====
    
    double N_Zc = b * h_sig->Integral();

    std::cout << "\n=== Zc EVENTS FROM FIT ===" << std::endl;
    std::cout << "N_Zc = " << N_Zc << std::endl;
}