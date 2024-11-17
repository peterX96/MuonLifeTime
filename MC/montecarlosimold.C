#include <iostream>
#include <cmath>

#include "TRandom3.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"

#include "generator.h"
#include "constants.h"

int main(int argc, char** argv){
	Int_t N=1000; // Iteration
	Int_t B=2; // Baseline
	double fit_bin, fit_err, avg_bin;
	TH1D* base_hist; // Array of histograms
	double min_bin=B*0.9, max_bin=B*1.1, fit_width = std::stof(argv[1]);
	
	TH1D* fit_hist = new TH1D("FitHist", ("Baseline " + std::to_string(B) + " estimation from part of the hist").c_str(), 50, min_bin, max_bin);
	TH1D* fit_err_hist = new TH1D("FitErrors", ("Fit errors for B = " + std::to_string(B)).c_str(), 50, 0.02, 0.07);
	TH1D* avg_hist = new TH1D("AvgHist", ("Baseline " + std::to_string(B) + " average of part of the hist").c_str(), 50, min_bin, max_bin);
	
	TF1 *fit0 = new TF1("pol0_left", "pol0", min_edge, min_edge+fit_width);
	TF1 *fit1 = new TF1("pol0_right", "pol0", (min_edge+max_edge)/2, (min_edge+max_edge)/2 + fit_width);
	
	for (Int_t i = 0; i < N; i++) {
		base_hist = generator(B, i);
		std::cout << "LOADING : " << i << " / " << N << std::endl;
		base_hist->Fit(fit0, "LR");
		base_hist->Fit(fit1, "LR");
		
		fit0->GetParameters(&fit_bin);
		fit_hist->Fill(fit_bin);
		fit_err = fit0->GetParError(0);
		fit_err_hist->Fill(fit_err);
		
		fit1->GetParameters(&fit_bin);
		fit_hist->Fill(fit_bin);
		fit_err = fit1->GetParError(0);
		fit_err_hist->Fill(fit_err);
		
		avg_bin = base_hist->Integral(min_edge, min_edge+fit_width)/fit_width;
		avg_hist->Fill(avg_bin);
		
		avg_bin = base_hist->Integral((min_edge+max_edge)/2, (min_edge+max_edge)/2 + fit_width)/fit_width;
		avg_hist->Fill(avg_bin);
		
	}
	
	TCanvas *canvas = new TCanvas("canvas", "Histogram and Fit", 1000, 1000);
	canvas->cd();
	fit_hist->SetLineColor(kRed);
	fit_hist->Draw("same");
	avg_hist->SetLineColor(kBlack);
	avg_hist->Draw("same");
	canvas->Update();
	canvas->SaveAs(("example_"+std::to_string(fit_width)+".png").c_str());
	
	TFile* res_file = new TFile(("baseline_"+std::to_string(fit_width)+".root").c_str(), "RECREATE");
	res_file->cd();
	fit_hist->Write();
	fit_err_hist->Write();
	avg_hist->Write();
	res_file->Close();
	return 0;
}


