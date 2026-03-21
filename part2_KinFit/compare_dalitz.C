void compare_dalitz()
{
    TFile* f_data = new TFile("plots_data/output.root");
    TFile* f_mc   = new TFile("plots_pipmZc3900/output.root");

    TH2D* h_data = (TH2D*)f_data->Get("h2_dalitz_final");
    TH2D* h_mc   = (TH2D*)f_mc->Get("h2_dalitz_final");

    TCanvas* c = new TCanvas("c","c",1400,600);

    c->Divide(2,1);

    c->cd(1);
    h_data->Draw("COLZ");
    gPad->SetRightMargin(0.15);

    c->cd(2);
    h_mc->Draw("COLZ");
    gPad->SetRightMargin(0.15);

    c->SaveAs("Dalitz_data_vs_mc.png");
}