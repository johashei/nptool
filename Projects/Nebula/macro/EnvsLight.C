void EnvsLight(){
  TFile* f = TFile::Open("root/simulation/SimulatedTree.root");
  TTree* t = (TTree*)f->Get("SimulatedTree");
  TCanvas *c1 = new TCanvas("c1","Energy vs Light",1080,720); 

  double En, Light;

  TBranch* branch_En;
  TBranch* branch_Light;

  TH1F *h_En = new TH1F("h_En","Energy",10,0.1,2);
  TH1F *h_Light = new TH1F("h_Light","Light",10,0.1,2);
  t->SetBranchAddress("Energy", &En, &branch_En);
  t->SetBranchAddress("Light", &Light, &branch_Light);
  
  int nentries = t->GetEntries();
  for(int i=0 ; i<nentries ; i++){
    t->GetEntry(i);
    h_En->Fill(En);
    h_Light->Fill(Light);
  }

  h_Light->SetLineColor(2);
  h_En->Draw();
  h_Light->Draw("SAME");
}
