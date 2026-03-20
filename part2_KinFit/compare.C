void compare()
{
    TFile* f_data = new TFile("plots_data/output.root");
    TFile* f_mc   = new TFile("plots_pippimJpsi/output.root");

    TH1D* h_data = (TH1D*)f_data->Get("h_mJpsipi_max_final");
    TH1D* h_mc   = (TH1D*)f_mc->Get("h_mJpsipi_max_final");

    h_data->Scale(1.0 / h_data->Integral());
    h_mc->Scale(1.0 / h_mc->Integral());

    h_data->SetLineColor(kBlack);
    h_mc->SetLineColor(kRed);

    h_data->SetLineWidth(2);
    h_mc->SetLineWidth(2);

    TCanvas* c = new TCanvas("c","c",1200,900);

    h_data->Draw("HIST");
    h_mc->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.6,0.7,0.85,0.85);
    leg->AddEntry(h_data,"Data","l");
    leg->AddEntry(h_mc,"MC","l");
    leg->Draw();

    c->SaveAs("Zc_data_vs_mc.png");
}