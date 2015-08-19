#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TF1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
   
void ntuple_pseudo()
{

  filltuple("E1E1_fit:generated_E1E1:E1E1_error",
            "c_E1E1",
            "E1E1 reconstructed with B#chiPT",
            "results/Pseudo_E1E1_lowE.pdf",
            -1000, -7, 7);

  filltuple("M1M1_fit:generated_E1E1:M1M1_error",
            "c_M1M1",
            "M1M1 reconstructed with B#chiPT",
            "results/Pseudo_M1M1_lowE.pdf",
            2.9, -7, 7);

  filltuple("E1M2_fit:generated_E1E1:E1M2_error",
            "c_E1M2",
            "E1M2 reconstructed with B#chiPT",
            "results/Pseudo_E1M2_lowE.pdf",
            -0.02, -7, 7);

  filltuple("M1E2_fit:generated_E1E1:M1E2_error",
            "c_M1E2",
            "M1E2 reconstructed with B#chiPT",
            "results/Pseudo_M1E2_lowE.pdf",
            2.2,  -7, 7);

  filltuple("alpha_fit:generated_E1E1:alpha_error",
            "c_alpha",
            "#alpha reconstructed with B#chiPT",
            "results/Pseudo_alpha_lowE.pdf",
            11.2, 7, 12);

  filltuple("beta_fit:generated_E1E1:beta_error",
            "c_beta",
            "beta reconstructed with B#chiPT",
            "results/Pseudo_beta_lowE.pdf",
            2.5,  0, 9);

}

void filltuple(Char_t* name, Char_t* namec, Char_t* namey, Char_t* namep, double nom, double rmin, double rmax)
{
   double erx[13] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   TFile *f1   = TFile::Open("results/save_root/Pseudo_fitter_search_gen_with_pasquini_lowE.root");
   TNtuple *r1 = (TNtuple*)f1->Get("result");

   TFile *f2 = TFile::Open("results/save_root/Pseudo_fitter_search_gen_with_pascalutsa_lowE.root");
   TNtuple *r2 = (TNtuple*)f2->Get("result");

   TCanvas* c1 = new TCanvas(namec,namec);
   c1->SetGrid(1,1);

   r1->Draw(name);
   TGraphErrors *hist = new TGraphErrors(r1->GetSelectedRows(),r1->GetV2(), r1->GetV1(), erx, r1->GetV3());

   r2->Draw(name);
   TGraphErrors *hist2 = new TGraphErrors(r2->GetSelectedRows(),r2->GetV2(), r2->GetV1(), erx, r2->GetV3());

   hist->SetMarkerStyle(21);
   hist->SetMarkerSize(0.7);
   hist->SetLineWidth(2);
   hist2->SetLineWidth(2);
   hist2->SetLineColor(kBlue);
   hist2->SetMarkerStyle(21);
   hist2->SetMarkerSize(0.7);
   hist2->SetMarkerColor(kBlue);   

   hist->Draw("AP");
   hist->GetYaxis()->SetRangeUser(rmin,rmax);
   hist2->Draw("P");

   hist->GetXaxis()->SetTitle("E1E1 generated with HDPV (black) - B#chiPT (blue)");
   hist->GetYaxis()->SetTitle(namey);
   hist->GetXaxis()->SetTitleSize(0.045);
   hist->GetYaxis()->SetTitleSize(0.045);
   hist->SetTitle("");

   if (nom == -1000)
   {
   TLine *line = new TLine(-7,-7,7,7);
   line->SetLineColor(kRed);
   line->SetLineWidth(4);
   line->Draw("same");
   }
   else
   {
   TLine *line = new TLine(-7,nom,7,nom);
   line->SetLineColor(kGreen+2);
   line->SetLineWidth(4);
   line->Draw("same");
   }

   c1->SaveAs(namep);
}

void ntuple1()
{
//    params fitparam(11.2, 2.5, -4.3, 2.9, -0.02, 2.2)
//
    fill("E1E1_pascalutsa:E1E1_pasquini>>h1(1,-7,7,1,-7,7)",
         "E1E1_pascalutsa:E1E1_pasquini>>h2(1,-7,7,1,-7,7)",
         "c_E1E1",
         "E1E1 reconstructed with B#chiPT",
         "results/Pseudo_E1E1_lowE.pdf",
         -4.3);

    fill("M1M1_pascalutsa:E1E1_pasquini>>h1(1,-7,7,1,-1,4)",
         "M1M1_pascalutsa:E1E1_pasquini>>h2(1,-7,7,1,-1,4)",
         "c_M1M1",
         "M1M1 reconstructed with B#chiPT",
         "results/Pseudo_M1M1_lowE.pdf",
         2.9);

    fill("E1M2_pascalutsa:E1E1_pasquini>>h1(1,-7,7,1,-1,6)",
         "E1M2_pascalutsa:E1E1_pasquini>>h2(1,-7,7,1,-1,6)",
         "c_E1M2",
         "E1M2 reconstructed with B#chiPT",
         "results/Pseudo_E1M2_lowE.pdf",
         -0.02);

    fill("M1E2_pascalutsa:E1E1_pasquini>>h1(1,-7,7,1,-0.5,2.5)",
         "M1E2_pascalutsa:E1E1_pasquini>>h2(1,-7,7,1,-0.5,2.5)",
         "c_M1E2",
         "M1E2 reconstructed with B#chiPT",
         "results/Pseudo_M1E2_lowE.pdf",
         2.2);

    fill("alpha_pascalutsa:E1E1_pasquini>>h1(1,-7,7,1,7,12)",
         "alpha_pascalutsa:E1E1_pasquini>>h2(1,-7,7,1,7,12)",
         "c_alpha",
         "#alpha reconstructed with B#chiPT",
         "results/Pseudo_alpha_lowE.pdf",
         11.2);

    fill("beta_pascalutsa:E1E1_pasquini>>h1(1,-7,7,1,2,8)",
         "beta_pascalutsa:E1E1_pasquini>>h2(1,-7,7,1,2,8)",
         "c_beta",
         "#beta reconstructed with B#chiPT",
         "results/Pseudo_beta_lowE.pdf",
         2.5);


}

void fill(Char_t* name1, Char_t* name2, Char_t* namec, Char_t* namey, Char_t* namep, double nom) 
{
// TFile *f1   = TFile::Open("results/Pseudo_results_pasquini_gen.root");
   TFile *f1   = TFile::Open("results/Pseudo_fitter_search_lowE_gen_with_pasquini.root");
   TNtuple *r1 = (TNtuple*)f1->Get("result");
   r1->SetMarkerStyle(21);

//   TFile *f2 = TFile::Open("results/Pseudo_results_pascalutsa_gen.root");
   TFile *f2 = TFile::Open("results/Pseudo_fitter_search_lowE_gen_with_pascalutsa.root");
   TNtuple *r2 = (TNtuple*)f2->Get("result");
   r2->SetMarkerColor(kBlue);
   r2->SetMarkerStyle(21);

   TCanvas* c1 = new TCanvas(namec,namec);
   c1->SetGrid(1,1);

   r1->Draw(name1);
   r2->Draw(name2,"","same");

   TH1* hist = (TH1*) gPad->GetListOfPrimitives()->FindObject("h1");
   hist->GetXaxis()->SetTitle("E1E1 generated with HDPV (black) - B#chiPT (blue)");
   hist->GetYaxis()->SetTitle(namey);
   hist->GetXaxis()->SetTitleSize(0.045);
   hist->GetYaxis()->SetTitleSize(0.045);
   hist->SetTitle("");
   TLine *line = new TLine(-7,nom,7,nom);
   line->SetLineColor(kGreen+2);
   line->SetLineWidth(4);
   line->Draw("same");

   c1->SaveAs(namep);
}



