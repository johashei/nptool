void Shift(){
  TFile* file_ok = TFile::Open("../../Outputs/Analysis/strasse_ok.root");
  TFile* file_shifted = TFile::Open("../../Outputs/Analysis/strasse_shifted.root");
  TTree* ok= (TTree*) file_ok->FindObjectAny("PhysicsTree");
  TTree* shifted= (TTree*) file_shifted->FindObjectAny("PhysicsTree");
  TCanvas* ctheta= new TCanvas("ControlTheta","ControlTheta",2000,1000);
  ctheta->Divide(2,1);
  ctheta->cd(1);
  string cond = "Theta12!=-1000";
  ok->Draw("Theta12>>ht(5000)",cond.c_str(),"") ; 
  shifted->Draw("Theta12>>hts(5000)",cond.c_str(),"same") ; 
  TH1* hts= (TH1*) gDirectory->FindObjectAny("hts");
  hts->SetFillColor(kOrange-3);
  hts->SetLineColor(kOrange-3);
  ctheta->cd(2);
  cond = "deltaPhi!=-1000";
  ok->Draw("deltaPhi>>hp(5000)",cond.c_str(),"") ; 
  shifted->Draw("deltaPhi>>hps(5000)",cond.c_str(),"same") ; 
  TH1* hps= (TH1*) gDirectory->FindObjectAny("hps");
  hps->SetFillColor(kOrange-3);
  hps->SetLineColor(kOrange-3);


}
