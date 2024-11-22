
void PlotComp() {

TFile *inputFile = new TFile("ExponentialFitWeek.root");

Double_t fit_bin,fit_err,dev;

Double_t mean,sigma,meanerr,out3sig,entries;

cout << "Number-"<< "Mean-" << "Sigma-" << "MeanError-" << "Outpoints-" << "Entries-" << endl ;
for (Int_t j = 1; j < 9; j++){

TH1F *h1 = new TH1F("h1", "", 10000, 1, 5);
 
TH1F *h2 = new TH1F("h2", "", 10000, 0, 1);

//TH1F *h3 = new TH1F("h3", "", 1000, 1.7, 2.3);

TH1F *h4 = new TH1F("h4", "", 1000,-5, 5);

  TString treeStartName = Form("Data_week_%d",j); 
  
  TTree *tree = (TTree*)inputFile->Get(treeStartName);

  tree->SetBranchAddress("fit_bin_tau",&fit_bin);
  tree->SetBranchAddress("fit_err_tau",&fit_err);
  tree->SetBranchAddress("dev_tau",&dev);

for (Int_t i = 0; i < tree->GetEntries(); i++) {
     
        tree->GetEntry(i);   // load current entry

	h1->Fill(fit_bin);
	h2->Fill(fit_err);

	if (dev >3 || dev < -3){
	h4->Fill(dev);
	}

    }

TF1 *fit = new TF1("fit","gaus",1,5);

h1->Fit("fit","LQ");

mean = fit->GetParameter(1);

entries = h1->GetEntries();

sigma = fit->GetParameter(2);

meanerr = h2->GetMean();

out3sig = h4->GetEntries();

cout << j << "\t " << mean << "\t " << sigma <<  "\t "<< meanerr << "\t "<< out3sig <<"\t "<< entries << endl ;

h1->Draw();

delete h1;
delete h2;
//delete h3;
delete h4;


}

}