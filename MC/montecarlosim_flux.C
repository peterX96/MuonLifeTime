#include <iostream>
#include <cmath>
#include <TTree.h>

#include "TRandom3.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"

#include "generator.h"
#include "constants.h"

template <typename T>
T* linspace_pointer(T start, T end, size_t num) {
    if (num <= 0) return nullptr;
    T* result = new T[num];
    T step = (end - start) / static_cast<T>(num - 1);

    for (size_t i = 0; i < num; ++i) {
        result[i] = start + step * i; 
    }
    
    return result;
}

double cdf(double theta){
	return 2/TMath::Pi()*(theta-TMath::Sin(4*theta)/4);
}

int main(int argc, char** argv){
	//some constants for current run
	int Nweeks = std::atoi(argv[1]); //we pass Nweeks from command line
	int Nevents = 11000*7*Nweeks; // day*7
	size_t arr_size = 10000; //for interpolation
	
	double *theta_arr = linspace_pointer(0., TMath::Pi(), arr_size);
	double *cdf_arr;
	cdf_arr = new double[arr_size];
	for (int i = 0; i < arr_size; i++){
		cdf_arr[i] = cdf(theta_arr[i]);
	}
	
	TTree *res_tree;
	double gen_x, gen_y, gen_phi, gen_theta;
	res_tree = new TTree("Tree_muons", "Tree with fir results");
	res_tree->Branch("x", &gen_x);
	res_tree->Branch("y", &gen_y);
	res_tree->Branch("phi", &gen_phi);
	res_tree->Branch("theta", &gen_theta);

	TFile* out_file = new TFile(("Flux_"+std::to_string(Nweeks)+"_weeks.root").c_str(), "RECREATE");
	out_file->cd();

	TRandom3 *rg = new TRandom3();

	for (int i = 0; i < Nevents; i++){
		gen_x = rg->Uniform()*slab_lenght - slab_length/2;
		gen_y = rg->Uniform()*slab_width - slab_width/2;
		gen_phi = rg->Uniform()*2*TMath::Pi();
		gen_theta = generator_theta(rg, theta_arr, cdf_arr, arr_size);
		
		res_tree->Fill();
	}
	
	res_tree->Write();
	out_file->Close();
	return 0;
}
