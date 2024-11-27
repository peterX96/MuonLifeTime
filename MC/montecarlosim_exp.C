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

	Int_t Nevents = 11000*7; // day*7
	double used_ch = 1880. ;

	double tauplus_ch = (tauplus-cal_par0)/cal_par1, tauminus_ch = (tauminus-cal_par0)/cal_par1;
	double a1,a2,a3,a4,a5,a6,c1;

	Int_t N = 500; // Iteration	
	TH1D* exp_hist; // Array of histograms
	
	TH1D* fit_bin = new TH1D("MeanDistr", "", 1000, 2, 2.4);
	TH1D* fit_err = new TH1D("SigmaDistr","" , 1000, 0, 1);
	TH1D* dev = new TH1D("Pull", "", 1000, -5, 5);
	
	TF1 *fit0 = new TF1("f1", "expo", min_edge, used_ch);
	

	fit0->SetParameter(0, 30); 
        fit0->SetParameter(1, 3); 

	
	// Storing data

	TFile *file = new TFile("ExponentialFit_week_tauminus.root", "UPDATE"); 

	double tau_ch = (tauplus-cal_par0)/cal_par1;

	// WEEK LOOP

	int sigmas;

	int sigma[4] = {0,1,3,5};

	for (int k = 0; k < 4; k++){

	sigmas = sigma[k];

	for (int j = 1; j < 9; j++){

 	TString treeName = Form("Data_week_%d_tauminus_%d",j,sigmas);

		TTree *tree = new TTree(treeName, treeName); 

		tree->Branch("fit_bin_tau", &a1);
		tree->Branch("fit_err_tau", &a2);
		tree->Branch("dev_tau", &a3);  
		tree->Branch("entries_cut", &c1);  

		/*
		tree->Branch("fit_bin_Amp", &a4);
		tree->Branch("fit_bin_tau_c", &b1);
		tree->Branch("fit_err_tau_c", &b2);
		tree->Branch("dev_tau_c", &b3);
		tree->Branch("fit_err_Amp", &a5);
		tree->Branch("dev_Amp", &a6);  
		*/

		for (Int_t i = 0; i < N; i++) {

			std::cout << "LOADING : "<< k << " / " << "3"  << "\t "  << j << " / " << "8"  << "\t " << i << " / " << N << std::endl;

			
			exp_hist = generator_exp(tauplus,Nevents*j,i,min_edge,used_ch);
			exp_hist->Fit(fit0, "LRQ", "", tauminus_ch*sigmas, used_ch);
			
			a1 = fit0->GetParameter(1);
			a2 = fit0->GetParError(1);
			a2 = a2/(a1*a1);
			a2 = a2*cal_par1;                   //propagation
			a1 = (-1/a1)*cal_par1+cal_par0;   // from ch to time
			a3 = (a1-tauplus)/a2;
			c1 = exp_hist-> Integral(tauminus_ch*sigmas, used_ch);

			/*
			a4 = fit0->GetParameter(0);
			a5 = fit0->GetParError(0);
			a5 = exp(a4)*a5;
			a4 = exp(a4);
			a6 = (a4-)/a6;
			*/
				
			tree->Fill();
	
			delete exp_hist; 
		}
		
	//file->Write(); 	
	}

// end loop

	file->Write();
        

}
file->Close();
}


