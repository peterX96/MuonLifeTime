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

int main(int argc, char** argv){
	//some constants for current run
	int Nevents = 11300*7; // day*7
	int B = 2;
	int used_ch = 1880;
	
	double tauplus_ch = (tauplus-cal_par0)/cal_par1, tauminus_ch = (tauminus-cal_par0)/cal_par1;
	int exp_int = Nevents-B*used_ch; //how many events correspond to the exp distribution
	double b_meas, b_err, tau_meas, tau_err;
	
	TH1D *exp_hist, *base_hist;
	TTree *res_tree;
	TF1 *base_func = new TF1("fit_pol0", "pol0", min_edge, used_ch);
	TF1 *exp_func = new TF1("fit_exp", "expo", min_edge, used_ch);
	exp_func->SetParameter(0, exp_int/500);
	exp_func->SetParameter(1, 0.8/tauplus_ch);
	
	TFile* out_file = new TFile(("exp"+std::to_string(tauplus)+"_b"+std::to_string(B)+"_fit_results.root").c_str(), "RECREATE");
	out_file->cd();
	
	int Niter = 10;
	for (int fit_ch = 100; fit_ch <= used_ch-5*tauminus_ch; fit_ch += 200){
		res_tree = new TTree(("Tree_exp_b"+std::to_string(fit_ch)).c_str(), "Tree with fir results");
		res_tree->Branch("baseline", &b_meas);
		res_tree->Branch("baseline_err", &b_err);
		res_tree->Branch("tau", &tau_meas);
		res_tree->Branch("tau_err", &tau_err);
		for (int i = 0; i < Niter; i++){
			base_hist = generator(B, i, min_edge, used_ch);
			exp_hist = generator_exp(tauplus, exp_int, i, min_edge, used_ch); //exp with tauplus
			exp_hist->Add(base_hist); //now exp_hist contain sum of exp and baseline
			
			exp_hist->Fit(base_func, "LRQ", "", used_ch-fit_ch, used_ch); //fit baseline in last fit_ch channels
			b_meas = base_func->GetParameter(0);
			b_err = base_func->GetParError(0);
			
			exp_hist->Add(base_func, -1.); //we subtract obtained baseline
			exp_hist->Fit(exp_func, "LRQ", "", min_edge, used_ch-fit_ch);
			tau_meas = 1/exp_func->GetParameter(1);
			tau_err = exp_func->GetParError(1)/(tau_meas*tau_meas);
			tau_meas = -tau_meas*cal_par1 + cal_par0;
			tau_err = tau_err*cal_par1 + cal_par0;
			
			res_tree->Fill();
			delete base_hist, exp_hist;
		}
		res_tree->Write();
		delete res_tree;
	}
	
	out_file->Close();
	return 0;
}
