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
	double Bdata = 1.5;
	int indices[5] = {2, 2, 1, 2, 1}; //to be passed through command line
	double tauplus = tauplus_values[indices[0]];
	double tauminus = tauminus_values[indices[1]];
	double R = R_values[indices[2]];
	double Nweeks = Nweeks_ratio[indices[3]];
	double B = B_values[indices[4]]*Bdata;
	
	int Nevents = 10775*7*Nweeks; // day*7*Nweeks
	int used_ch = 1880;
	
	double tauplus_ch = (tauplus-cal_par0)/cal_par1, tauminus_ch = (tauminus-cal_par0)/cal_par1;
	//for t+ = 2.2 t+_ch = 184 and for t- = 0.9 t-_ch = 74
	
	double Aplus = R*(Nevents-B*used_ch*Nweeks)/(tauminus_ch + R*tauplus_ch);
	double Aminus = Aplus/R;
	int exp_plus_int = std::round(Aplus*tauplus_ch), exp_minus_int = std::round(Aminus*tauminus_ch);
	
	std::cout << tauplus << ' ' << tauminus << ' ' << R << ' ' << Nweeks << ' ' << B << ' ' << Aplus << ' ' << Aminus << std::endl;
	
	double b_meas, b_meas_new, b_err, b_err_new, b_chi, tau_p_meas, tau_p_err, tau_m_meas, tau_m_err, tau_chi, Aplus_meas, Aplus_err, Aminus_meas, Aminus_err, r_meas, r_err;
	
	TH1D *exp1_hist, *exp2_hist, *base_hist;
	TTree *res_tree;
	TF1 *base_func = new TF1("fit_pol0", "pol0", min_edge, used_ch);
	//									[0]: B, [1]: A+, [2]: tau+ ; taus are in channels
	TF1 *exp_b_func = new TF1("fit_1exp", "[0]+[1]*exp(-x/[2])", min_edge, used_ch);
	//									[0]: B, [1]: A-, [2]: tau-, [3]: A+, [4]: tau+
	TF1 *exp2_b_func = new TF1("fit_2exp", "[0]+[1]*exp(-x/[2]) + [3]*exp(-x/[4])", min_edge, used_ch);
	//									[0]: B, [1]: A-, [2]: tau-, [3]: R, [4]: tau+	
	TF1 *exp2_b_func_r = new TF1("fit_2exp", "[0]+[1]*(exp(-x/[2]) + [3]*exp(-x/[4]))", min_edge, used_ch);
	TRandom3 *rg_seed = new TRandom3();
	
	exp_b_func->SetParameter(1, 1.1*Aplus);
	exp_b_func->SetParameter(2, tauplus_ch+1);
	
	exp2_b_func->SetParameter(1, 1.1*Aminus);
	exp2_b_func->SetParameter(2, tauminus_ch - 23);
	exp2_b_func->SetParameter(3, 0.8*Aplus);
	exp2_b_func->SetParameter(4, tauplus_ch + 12);
	
	exp2_b_func->SetParameter(1, 1.1*Aminus);
	exp2_b_func_r->SetParameter(2, tauminus_ch + 12);
	exp2_b_func->SetParameter(3, 0.8*Aplus);
	exp2_b_func_r->SetParameter(4, tauplus_ch - 20);
	
	//do Tstring
	TFile* out_file = new TFile("2exp_b_bounds.root", "RECREATE");
	out_file->cd();
	
	int Niter = 1000;
	for (int b_ch = 5; b_ch < 10; b_ch += 1){     //b_ch - start for baseline study (in tau+ units)
		for (int t_ch = 5; t_ch < 10; t_ch += 1){ //t_ch - start for tau+     study (in tau- units), 10*tau- <= 5*tau+ so they dont overlap
			res_tree = new TTree(("Tree_b"+std::to_string(b_ch)+"_t"+std::to_string(t_ch)).c_str(), "Tree with fir results");
			res_tree->Branch("baseline", &b_meas);
			res_tree->Branch("baseline_new", &b_meas_new);
			res_tree->Branch("baseline_err", &b_err);
			res_tree->Branch("baseline_err_new", &b_err_new);
			res_tree->Branch("baseline_chi", &b_chi);
			res_tree->Branch("tau_p", &tau_p_meas);
			res_tree->Branch("tau_p_err", &tau_p_err);
			res_tree->Branch("tau_m", &tau_m_meas);
			res_tree->Branch("tau_m_err", &tau_m_err);
			res_tree->Branch("tau_chi", &tau_chi);
			res_tree->Branch("Aplus", &Aplus_meas);
			res_tree->Branch("Aplus_err", &Aplus_err);
			res_tree->Branch("Aminus", &Aminus_meas);
			res_tree->Branch("Aminus_err", &Aminus_err);
			res_tree->Branch("R", &r_meas);
			res_tree->Branch("R_err", &r_err);
			for (int i = 0; i < Niter; i++){
				base_hist = generator(B*Nweeks, rg_seed->Uniform()*100000, min_edge, used_ch);
				exp1_hist = generator_exp(tauplus, exp_plus_int, rg_seed->Uniform()*100000, min_edge, used_ch); //exp with tauplus
				exp2_hist = generator_exp(tauminus, exp_minus_int, rg_seed->Uniform()*100000, min_edge, used_ch); //exp with tauminus
				
				exp2_hist->Add(base_hist); 
				exp2_hist->Add(exp1_hist); //now exp2_hist contains sum of two exp-s and baseline
				//1st fit - baseline
				exp2_hist->Fit(base_func, "LCQ", "", b_ch*tauplus_ch, used_ch); //fit baseline in last fit_ch channels
				b_meas = base_func->GetParameter(0)/Nweeks;
				b_err = base_func->GetParError(0)/Nweeks;
				exp_b_func->SetParameter(0, b_meas*Nweeks);
				exp_b_func->SetParLimits(0, (b_meas-3*b_err)*Nweeks, (b_meas+3*b_err)*Nweeks);
				
				//2nd fit - tauplus
				exp2_hist->Fit(exp_b_func, "LCQ", "", t_ch*tauminus_ch, used_ch);//-fit_ch);
				b_meas_new = exp_b_func->GetParameter(0)/Nweeks;
				b_err_new = exp_b_func->GetParError(0)/Nweeks;
				tau_p_meas = exp_b_func->GetParameter(2);
				tau_p_err = exp_b_func->GetParError(2);
				exp2_b_func->SetParameter(0, b_meas_new*Nweeks); //baseline
				exp2_b_func->SetParLimits(0, (b_meas_new-3*b_err_new)*Nweeks, (b_meas_new+3*b_err_new)*Nweeks);
				exp2_b_func->SetParameter(4, tau_p_meas); //tauplus
				exp2_b_func->SetParLimits(4, tau_p_meas-3*tau_p_err,  tau_p_meas+3*tau_p_err);
				
				//3rd fit - tauminus
				exp2_hist->Fit(exp2_b_func, "LQ", "", min_edge, used_ch);
				tau_chi = exp2_b_func->GetChisquare();
				b_meas_new = exp2_b_func->GetParameter(0)/Nweeks;
				b_err_new = exp2_b_func->GetParError(0)/Nweeks;
				
				tau_p_meas = exp2_b_func->GetParameter(4);
				tau_p_err = exp2_b_func->GetParError(4);
				tau_p_meas = tau_p_meas*cal_par1 + cal_par0;
				tau_p_err = tau_p_err*cal_par1; // + cal_par0;
				
				tau_m_meas = exp2_b_func->GetParameter(2);
				tau_m_err = exp2_b_func->GetParError(2);
				tau_m_meas = tau_m_meas*cal_par1 + cal_par0;
				tau_m_err = tau_m_err*cal_par1; // + cal_par0;
				
				Aplus_meas = exp2_b_func->GetParameter(3);
				Aplus_err = exp2_b_func->GetParError(3);
				
				Aminus_meas = exp2_b_func->GetParameter(1);
				Aminus_err = exp2_b_func->GetParError(1);
				
				r_meas = Aplus_meas/Aminus_meas;
				r_err = r_meas*(TMath::Sqrt( std::pow(Aplus_err/Aplus, 2) + std::pow(Aminus_err/Aminus, 2) )); //uncorrelated case
				
				res_tree->Fill();
				delete base_hist, exp1_hist, exp2_hist;
				if (b_ch==5 and i==0){
					TCanvas *c1 = new TCanvas("", "", 900, 900);
					c1->cd();
					exp2_hist->Draw();
					c1->SaveAs("aboba.png");
					std::cout << tau_p_meas << ' ' << tau_m_meas << ' ' << b_meas << ' ' << b_meas_new << ' ' << r_meas << ' ' << Aplus_meas << ' ' << Aminus_meas << std::endl;
					//exp2_hist->Write();
					//return 0;
				}
			}
			std::cout << "Simulated baseline from " << b_ch << " to " << used_ch << ", and tau+ from " << t_ch << " to " << used_ch << std::endl;
			res_tree->Write();
			delete res_tree;
		}
	}
	
	out_file->Close();
	return 0;
}
