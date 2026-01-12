#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>

void draw_syst() {
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open("RSG_W_SYST/signal/output_RSGravitonToGG_M700_kMpl001_13TeV_pythia8.root");
    TTree *tree_nom = (TTree*)f->Get("DiphotonTree/rsg_700_001_13TeV_rsg_cut");
    TTree *tree_up = (TTree*)f->Get("DiphotonTree/rsg_700_001_13TeV_rsg_cut_ScaleUp01sigma");

    // --- Create canvas ---
    TCanvas *c_all = new TCanvas("c_all", "CMS_hgg_mass PDFs", 800, 600);

    // --- Nominal histogram ---
    TH1F *h_nom = new TH1F("h_nom", "Systematics - Scale up 1 sigma;CMS_hgg_mass [GeV]; Events (norm.)",
                           100, 650, 750);
    h_nom->Sumw2();

    tree_nom->Draw("CMS_hgg_mass >> h_nom", "weight", "goff");

    double integral_nom = h_nom->Integral("width"); // PDF normalization
    double sumW_nom = h_nom->GetSumOfWeights();
    if (integral_nom > 0) h_nom->Scale(1.0 / sumW_nom);

    // --- Up variation histogram ---
    TH1F *h_up = new TH1F("h_up", "",
                          100, 650, 750);
    h_up->Sumw2();

    tree_up->Draw("CMS_hgg_mass >> h_up", "weight", "goff");

    double integral_up = h_up->Integral("width");
    double sumW_up = h_up->GetSumOfWeights();
    if (integral_up > 0) h_up->Scale(1.0 / sumW_up);

    // --- Styling ---
    h_nom->SetLineColor(kBlack);
    h_nom->SetLineWidth(2);

    h_up->SetLineColor(kRed);
    h_up->SetLineWidth(2);

    // --- Set common y-axis range ---
    double ymax = std::max(h_nom->GetMaximum(), h_up->GetMaximum());
    h_nom->SetMaximum(1.1 * ymax);
    h_nom->SetMinimum(0);

    // --- Draw histograms ---
    h_nom->Draw("HIST");
    h_up->Draw("HIST SAME");

    // ---------- Load RooFit PDF ----------
    TFile *f_sig = TFile::Open("flashggFinalFit/Signal/outdir_RSG_500-1000_postEE_syst/signalFit/output/CMS-HGG_sigfit_RSG_500-1000_postEE_syst_rsg_2022postEE_rsg_cut.root");
    RooWorkspace *w_sig = (RooWorkspace*)f_sig->Get("wsig_13TeV");
    RooAbsPdf *pdf = w_sig->pdf("hggpdfsmrel_rsg_2022postEE_rsg_cut_13TeV");

    // ---------- Set variables ----------
    RooRealVar *MH = w_sig->var("MH");
    RooRealVar *G0 = w_sig->var("G0");
    RooRealVar *nuis = w_sig->var("CMS_hgg_nuisance_Scale_2022postEE");

    MH->setVal(700);
    G0->setVal(0.099);

    // ---------- Plot PDF ----------
    RooRealVar *x = w_sig->var("CMS_hgg_mass");
    RooBinning b(100, 650, 750);
    x->setBinning(b);
    RooPlot *frame = x->frame();

    nuis->setVal(0);
    pdf->plotOn(frame, RooFit::Normalization(1, RooAbsReal::Relative), RooFit::LineColor(kBlue), RooFit::LineWidth(2), RooFit::Name("pdf_nom"));

    nuis->setVal(1);
    pdf->plotOn(frame, RooFit::Normalization(1, RooAbsReal::Relative), RooFit::LineColor(kViolet), RooFit::LineWidth(2), RooFit::Name("pdf_up"));

    frame->Draw("SAME"); // draw PDF on same canvas

    // --- Legend ---
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(h_nom, "Histo nom", "l");
    leg->AddEntry(h_up, "Histo up", "l");
    leg->Draw();
    leg->AddEntry("pdf_nom","PDF nom.","LP");
    leg->AddEntry("pdf_up","PDF up","LP");
    leg->Draw();

    // --- Update canvas and save ---
    c_all->Update();
    c_all->SaveAs("PDF_syst_check_Scale.png");  // save as PNG
    c_all->SaveAs("PDF_syst_check_Scale.root"); // ROOT file
}

int main() {
    draw_syst();
    return 0;
}