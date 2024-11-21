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
	Int_t N = 2000; // Iteration	
	TH1D* exp_hist; // Array of histograms
	
	TH1D* fit_bin = new TH1D("MeanDistr", "", 1000, 2, 2.4);
	TH1D* fit_err = new TH1D("SigmaDistr","" , 1000, 0, 1);
	TH1D* dev = new TH1D("Pull", "", 1000, -5, 5);
	
	TF1 *fit0 = new TF1("f1", "expo", min_edge, max_edge);

	fit0->SetParameter(0, 30); 
        fit0->SetParameter(1, 3); 

	// Storing data
	
	TFile *file = new TFile("ExponentialFit.root", "UPDATE"); 

	TString treeName = Form("Data"); 

        TTree *tree = new TTree(treeName, treeName); 

        double a1,a2,a3,a4,a5,a6;

        tree->Branch("fit_bin_tau", &a1);
	tree->Branch("fit_err_tau", &a2);
	tree->Branch("dev_tau", &a3);  
	tree->Branch("fit_bin_Amp", &a4);
	tree->Branch("fit_err_Amp", &a5);
	tree->Branch("dev_Amp", &a6);  

	for (Int_t i = 0; i < N; i++) {

		std::cout << "LOADING : " << i << " / " << N << std::endl;

		exp_hist = generator_exp (tauplus,Nevents,i);
		exp_hist->Fit(fit0,"LR");
		
		a1 = fit0->GetParameter(1);
		a2 = fit0->GetParError(1);
		a2 = a2/(a1*a1);
		a2 = a2*0.01182;                   //propagation
		a1 = (-1/a1)*0.01182+0.02727;   // from ch to time

		a3 = (a1-tauplus)/a2;


		a4 = fit0->GetParameter(0);
		a5 = fit0->GetParError(0);
		a5 = exp(a4)*a5;
		a4 = exp(a4);
		//a6 = (a4-)/a6;

		tree->Fill();
	}

    file->Write();
    file->Close();


}


