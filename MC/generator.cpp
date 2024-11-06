#include <iostream>
#include <cmath>

#include "TRandom3.h"
#include "TH1D.h"

const int Nchannels = 2048;
const double min_edge = 0., max_edge = 2048.;

TH1D* generator(int B=2){
	TRandom3 *rg = new TRandom3();
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
