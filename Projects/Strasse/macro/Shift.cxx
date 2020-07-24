void Shift(){
  TFile* file_ok = TFile::Open("../../Outputs/Analysis/strasse_ok.root");
  TFile* file_shifted = TFile::Open("../../Outputs/Analysis/strasse_shifted.root");
  TTree* ok= (TTree*) file_ok->FindObjectAny("PhysicsTree");
  TTree* shifted= (TTree*) file_shifted->FindObjectAny("PhysicsTree");
  TCanvas* ctheta= new TCanvas("ControlTheta","ControlTheta",1000,1000);
  string cond = "Theta12!=-1000";
  ok->Draw("Theta12>>ht(5000)",cond.c_str(),"") ; 
  shifted->Draw("Theta12>>hts(5000)",cond.c_str(),"same") ; 
  TH1* hts= (TH1*) gDirectory->FindObjectAny("hts");
  hts->SetFillColor(kOrange-3);
  hts->SetLineColor(kOrange-3);

}
