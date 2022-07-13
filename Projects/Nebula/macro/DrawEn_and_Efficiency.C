#include "TNebulaPhysics.h"

void DrawEn_and_Efficiency(){
  TFile* f = TFile::Open("root/analysis/PhysicsTree.root");
  TTree* t = (TTree*)f->Get("PhysicsTree");
  TCanvas *c1 = new TCanvas("c1","Efficiency",1080,720); 
  t->Draw("Efficiency");

  TNebulaPhysics *Neb_Phys = new TNebulaPhysics();
  std::vector<double> NeutronEnergy;
  //std::vector<int>* Efficiency = 0;
  double En, counter;
  //int Eff;

  //TBranch* branch_Eff;
  TBranch* branch_Nebula;
  //TH1F *h_Eff = new TH1F("h_Eff","test");
  TH1F *h_En = new TH1F("h_En","Neutrons Retrieved Energy",140,140,210);
  //t->SetBranchAddress("Efficiency",&Efficiency,&branch_Eff);
  t->SetBranchAddress("Nebula",&Neb_Phys,&branch_Nebula);
  //int nentries_Eff = branch_Eff->GetEntries();
  int nentries_Neb = branch_Nebula->GetEntries();
  
  for(int i = 0; i < nentries_Neb ; i++){
    t->GetEntry(i);
    NeutronEnergy = Neb_Phys->NeutronEn;
    int vec_size = NeutronEnergy.size();
    for(int j = 0; j < vec_size; j++){
      En = NeutronEnergy[j];
      h_En->Fill(En);
    }
  }

  TCanvas *c2 = new TCanvas("c2","Neutron Energy",1080,720); 
  h_En->SetLineWidth(6);
  h_En->SetLineColor(6);
  c2->SetFrameLineColor(0);
  h_En->Draw("C9");
}
