#include <iostream>
#include <cmath>
#include <ctime>

#include "TRandom3.h"
#include "TH1D.h"

#include "generator.h"
//#include "constants.h"

TH1D* generator(int B, int seed, float min_edge_hist, float max_edge_hist){
	TRandom3 *rg = new TRandom3(seed);//static_cast<int>(std::time(0)));
	float channel;
	// Bin 0 - overflow; Bin Nchannels+1 - overflow; Bins 1 and Nchannels are bounds
	TH1D* baseline_hist = new TH1D(("Base"+std::to_string(seed)).c_str(), "Baseline", Nchannels, min_edge_hist, max_edge_hist);
	
	for (int i = 0; i < B*Nchannels; i ++){
		channel = (int)std::floor((rg->Uniform(min_edge, max_edge)))+0.5;
		baseline_hist->Fill(channel);
	}
	
	baseline_hist->Draw();
	return baseline_hist;
}

TH1D* generator_exp(float tau, int integral, int seed, float min_edge_hist, float max_edge_hist){
	//tau in \mu s, integral is a total number of events
	TRandom3 *rg = new TRandom3(seed);//static_cast<int>(std::time(0)));
	float tau_ch = (tau-cal_par0)/cal_par1; //reevaluating tau in channels
	TH1D* exp_hist = new TH1D(("Exp"+std::to_string(seed)).c_str(), ("Exponent with tau = "+std::to_string(tau)).c_str(), Nchannels, min_edge_hist, max_edge_hist);
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

template <typename T>
T lin_interpol(T x, T* x_arr, T* y_arr, size_t size){
	if(x > x_arr[size-1] or x < x_arr[0]){
		throw std::out_of_range("x out of interpolation range");
	}
	size_t init_guess = static_cast<size_t>(std::round((x-x_arr[0]) / (x_arr[size-1]-x_arr[0]) * size )) - 10;
	for(size_t i = 0; i < size-1; i++){
		if (x > x_arr[i] and x <= x_arr[i+1]){
			return y_arr[i] + (x-x_arr[i]) * (y_arr[i+1]-y_arr[i]) / (x_arr[i+1]-x_arr[i]);
		}
	}
}

double generator_theta(TRandom3 *rg, double *theta_arr, double *cdf_arr, size_t arr_size){
	double cdf_value = rg->Uniform();
	return lin_interpol(cdf_value, cdf_arr, theta_arr, arr_size);
}
