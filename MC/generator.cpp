#include <iostream>
#include <cmath>
#include <ctime>

#include "TRandom3.h"
#include "TH1D.h"

#include "generator.h"
#include "constants.h"

TH1D* generator(int B, int seed){
	TRandom3 *rg = new TRandom3(seed);//static_cast<int>(std::time(0)));
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

TH1D* generator_exp(float tau, int integral, int seed){
	//tau in \mu s, integral is a total number of events

	TRandom3 *rg = new TRandom3(seed);//static_cast<int>(std::time(0)));

	//TRandom3 *rg = new TRandom3(static_cast<int>(std::time(0)));


	float tau_ch = (tau-cal_par0)/cal_par1; //reevaluating tau in channels
	TH1D* exp_hist = new TH1D("Exp", ("Exponent with tau = "+std::to_string(tau)).c_str(), Nchannels, min_edge, max_edge);
	float entry;
	
	for (int i = 0; i < integral; i++){
		entry = rg->Exp(tau_ch);
		while (entry < min_edge or entry > max_edge){
		entry = rg->Exp(tau_ch);
		}
		exp_hist->Fill(entry);
	}
	return exp_hist;
}
