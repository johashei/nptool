#include "/local/lemair/nptool/NPLib/Detectors/Nebula/TNebulaPhysics.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::find
#include <vector>   


void EnvsLight(){
  TFile* f = TFile::Open("root/simulation/SimulatedTree.root");
  TTree* t = (TTree*)f->Get("SimulatedTree");
  TCanvas *c1 = new TCanvas("c1","Position Z",1920,1200);

  double Energy, Light;
  TBranch* branch_En;
  TBranch* branch_Li;
  
  TH1F *h_en = new TH1F("h_en","Deposited Energy",100,1,200);
  TH1F *h_li = new TH1F("h_li","Deposited Light",100,1,200);
  t->SetBranchAddress("PlasticBar_Energy", &Energy, &branch_En);
  t->SetBranchAddress("PlasticBar_Light", &Light, &branch_Li);
  
  int nentries = t->GetEntries();
  for(int i=0 ; i<nentries ; i++){
    t->GetEntry(i);
    h_en->Fill(Energy);
    h_li->Fill(Light);
  }
  
  h_en->SetLineColor(kBlue-7);
  h_en->Draw();
  h_li->SetLineColor(kRed-7);
  h_li->Draw("SAME");
}
 
