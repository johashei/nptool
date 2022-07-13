#include "/local/lemair/nptool/NPLib/Detectors/Nebula/TNebulaPhysics.h"

void PlasticBarTimeAndPosition(){
  TFile* f = TFile::Open("root/simulation/SimulatedTree.root");
  TTree* t = (TTree*)f->Get("SimulatedTree");

  double Time;
  //int Eff;

  //TBranch* branch_Eff;
  TBranch* branch;
  //TH1F *h_Eff = new TH1F("h_Eff","test");
  TH1F *h_Time = new TH1F("h_Time","PlasticBar_Time",500,55,65);
  //t->SetBranchAddress("Efficiency",&Efficiency,&branch_Eff);
  t->SetBranchAddress("PlasticBar_Time",&Time,&branch);
  //int nentries_Eff = branch_Eff->GetEntries();
  int nentries = branch->GetEntries();
  
  for(int i = 0; i < nentries ; i++){
    t->GetEntry(i);
    h_Time->Fill(Time);
  }

  TCanvas *c2 = new TCanvas("c2","Neutron Energy",1080,720); 
  h_Time->Draw("");
}
