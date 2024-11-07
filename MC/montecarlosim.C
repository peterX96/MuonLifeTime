#include <iostream>
#include <cmath>

#include "TRandom3.h"
#include "TH1D.h"

const int Nchannels = 2048;
const double min_edge = 0., max_edge = 2048.;

TH1D* generator(int B=2){
	TRandom3 *rg = new TRandom3(static_cast<int>(time(0))); // I MADE ONLY THIS CHANGE
	float channel;
	// Bin 0 - overflow; Bin Nchannels+1 - overflow; Bins 1 and Nchannels are bounds
	TH1D* baseline_hist = new TH1D("Base", "Baseline", Nchannels, min_edge, max_edge);
	
	for (int i = 0; i < B*Nchannels; i ++){
		channel = (int)std::floor((rg->Uniform(min_edge, max_edge)))+0.5;
		baseline_hist->Fill(channel);
	}
	
	baseline_hist->Draw();
	std::cout << baseline_hist->GetBinContent(0) << ' ' << baseline_hist->GetBinContent(1) << ' ' 
		<< baseline_hist->GetBinContent(Nchannels) << ' ' << baseline_hist->GetBinContent(Nchannels+1) << std::endl;
	return baseline_hist;
}


void montecarlosim(){

TH1D* baseline_sum = new TH1D("Base", "Baseline average", 2048,0.,2048.);

Int_t N=500; // Iteration
Int_t B=2; // Baseline
float scale = 1.0f / N; // Normalization term

TH1D* histograms[N]; // Array of histograms

for (Int_t i = 0; i < N+1; i++) {

histograms[i] = generator();

baseline_sum->Add(histograms[i]);

cout << "LOADING : " << i << " / " << N << endl;
}

TH1F* result = (TH1F*)baseline_sum->Clone("result");

result ->Scale(scale);

TCanvas *canvas = new TCanvas("canvas", "Histogram and Fit", 800, 600);

result->Draw("HIST");

// fit with pol0 & pol1

	// Fit a linear (pol1) to the histogram
        TF1 *fit_pol0 = new TF1("fit_pol0", "pol0", 0, 2048);
        result->Fit(fit_pol0, "R");  

        // Fit a linear (pol1) to the histogram
        TF1 *fit_pol1 = new TF1("fit_pol1", "pol1", 0, 2048);
        result->Fit(fit_pol1, "R");
	
canvas->Update();

}


