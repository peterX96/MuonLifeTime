
void PlotComp() {

TFile *inputFile = new TFile("PartialFit.root");

Double_t fit_bin,fit_err,avg_bin,dev;

Double_t mean,sigma,meanerr,out3sig;

Int_t j;

for (Int_t j = 300; j < 850; j+= 50){

TH1F *h1 = new TH1F("h1", "", 1000, 1.7, 2.3);
 
TH1F *h2 = new TH1F("h2", "", 1000, 0.05, 0.09);

TH1F *h3 = new TH1F("h3", "", 1000, 1.7, 2.3);

TH1F *h4 = new TH1F("h4", "", 1000,-5, 5);

  TString treeStartName = Form("Data_%d",j); 
  
  TTree *tree = (TTree*)inputFile->Get(treeStartName);

  tree->SetBranchAddress("fit_bin",&fit_bin);
  tree->SetBranchAddress("fit_err",&fit_err);
  tree->SetBranchAddress("avg_bin",&avg_bin);
  tree->SetBranchAddress("dev",&dev);


for (Int_t i = 0; i < tree->GetEntries(); i++) {
     
    tree->GetEntry(i);   // load current entry

h1->Fill(fit_bin);
h2->Fill(fit_err);
h3->Fill(avg_bin);

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

cout << j << "\t " << mean << "\t " << sigma <<  "\t "<< meanerr << "\t "<< out3sig << endl;

h1->Draw();

delete h1;
delete h2;
delete h3;
delete h4;


}

}