
void PlotComp1() {

TFile *inputFile = new TFile("ExponentialFitRange.root");

Double_t fit_bin,fit_err,avg_bin,dev,entries;

Double_t mean,sigma,meanerr,out3sig,entriescut;


cout << "tauminus" << "\t " << "weeks" << "\t " << "mean" << "\t " << "sigma" <<  "\t "<< "meanerr" << "\t "<< "out3sig "<< "\t " << "entries"<< endl;

	for (int j = 0; j < 8; j++){

 	TString treeStartName = Form("Data_week1_tauminus%d",j);


TH1F *h1 = new TH1F("h1", "", 1000, 1.7, 2.3);
 
TH1F *h2 = new TH1F("h2", "", 1000, 0, 1);

TH1F *h4 = new TH1F("h4", "", 1000,-5, 5);

TH1F *h3 = new TH1F("h3", "", 10000,0, 7700000);
  
  TTree *tree = (TTree*)inputFile->Get(treeStartName);

  tree->SetBranchAddress("fit_bin_tau",&fit_bin);
  tree->SetBranchAddress("fit_err_tau",&fit_err);
  tree->SetBranchAddress("dev_tau",&dev);
  //tree->SetBranchAddress("entries_cut",&entries);


for (Int_t i = 0; i < tree->GetEntries(); i++) {
     
    tree->GetEntry(i);   // load current entry

h1->Fill(fit_bin);
h2->Fill(fit_err);
//h3->Fill(entries);

	if (dev >3 || dev < -3){
	h4->Fill(dev);
	}

    }

TF1 *fit = new TF1("fit","gaus",1.6,2.3);

h1->Fit("fit","LQ");


mean = fit->GetParameter(1);

sigma = fit->GetParameter(2);

meanerr = h2->GetMean();

out3sig = h4->GetEntries();

//entriescut = h3->GetMean();

cout << j << "\t " << mean << "\t " << sigma <<  "\t "<< meanerr << "\t "<< out3sig  << endl;


delete h1;
delete h2;
delete h3;
delete h4;



}
}