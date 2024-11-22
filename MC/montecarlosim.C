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
// variable exponential
	Int_t Nevents = 11300;
	Double_t tauplus = 2.2;
	
	double a1,a2,a3,a4,a5,a6;
	double b1,b2,b3,b4;

	Int_t N = 500; // Iteration	
	TH1D* exp_hist; // Array of histograms
	
	TH1D* fit_bin = new TH1D("MeanDistr", "", 1000, 2, 2.4);
	TH1D* fit_err = new TH1D("SigmaDistr","" , 1000, 0, 1);
	TH1D* dev = new TH1D("Pull", "", 1000, -5, 5);
	
	TF1 *fit0 = new TF1("f1", "expo", min_edge, max_edge);
	TF1 *fit_custom = new TF1("f_custom", "[0]*exp(-x/[1])", min_edge, max_edge);
	

	fit0->SetParameter(0, 30); 
    fit0->SetParameter(1, 3); 

	fit_custom->SetParameter(0, 30);
	fit_custom->SetParameter(1, 1);

	// Storing data
	
	TFile *file = new TFile("ExponentialFit.root", "UPDATE"); 
	double tau_ch = (tauplus-cal_par0)/cal_par1;
	for (int sigmas = 0; sigmas <= 7; sigmas ++){
		TString treeName = Form(("Data"+std::to_string(sigmas)).c_str());
		
		TTree *tree = new TTree(treeName, treeName); 

		tree->Branch("fit_bin_tau", &a1);
		tree->Branch("fit_err_tau", &a2);
		tree->Branch("dev_tau", &a3);  
		tree->Branch("fit_bin_Amp", &a4);
		
		tree->Branch("fit_bin_tau_c", &b1);
		tree->Branch("fit_err_tau_c", &b2);
		tree->Branch("dev_tau_c", &b3);
		
		tree->Branch("fit_err_Amp", &a5);
		tree->Branch("dev_Amp", &a6);  

		for (Int_t i = 0; i < N; i++) {

			std::cout << "LOADING : " << i << " / " << N << std::endl;

			exp_hist = generator_exp(tauplus,Nevents,i);
			exp_hist->Fit(fit0, "LR", "", tau_ch*sigmas, max_edge);
			
			/*exp_hist->Fit(fit_custom, "LR");
			b1 = fit_custom->GetParameter(1);
			b2 = fit_custom->GetParError(1);
			b2 = b2*cal_par1;           //propagation
			b1 = b1*cal_par1+cal_par0;   // from ch to time
			b3 = (b1-tauplus)/b2;*/
			
			a1 = fit0->GetParameter(1);
			a2 = fit0->GetParError(1);
			a2 = a2/(a1*a1);
			a2 = a2*cal_par1;                   //propagation
			a1 = (-1/a1)*cal_par1+cal_par0;   // from ch to time

			a3 = (a1-tauplus)/a2;


			a4 = fit0->GetParameter(0);
			a5 = fit0->GetParError(0);
			a5 = exp(a4)*a5;
			a4 = exp(a4);
			//a6 = (a4-)/a6;

			tree->Fill();
		}

		file->Write();
	}
    file->Close();


}


