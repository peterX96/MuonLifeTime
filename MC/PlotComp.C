
void PlotComp() {

  TFile *inputFile = new TFile("PartialFit.root");


Double_t fit_bin,fit_err,avg_bin,dev;

	 TH1F *h1 = new TH1F("h1", "", 100, 0.07, 0.09);
  
  int j=300;
  TString treeStartName = Form("Data_%d",j); 
  TString Leaf1name = Form("fit_bin_%d",j); 
  TString Leaf2name = Form("fit_err_%d",j); 
  TString Leaf3name = Form("avg_bin_%d",j); 
  TString Leaf4name = Form("dev_%d",j); 
  
  TTree *tree = (TTree*)inputFile->Get(treeStartName);

  tree->SetBranchAddress("fit_bin",&fit_bin);
  tree->SetBranchAddress("fit_err",&fit_err);
  tree->SetBranchAddress("avg_bin",&avg_bin);
  tree->SetBranchAddress("dev",&dev);


for (Int_t i = 0; i < tree->GetEntries(); i++) {
     
    tree->GetEntry(i);   // load current entry

    h1->Fill(fit_err);  
    
    }

h1->Draw();





}