#include "/local/lemair/nptool/NPLib/Detectors/Nebula/TNebulaPhysics.h"

void ChargeDraw(){
  TFile* f = TFile::Open("root/analysis/PhysicsTree.root");
  TTree* t = (TTree*)f->Get("PhysicsTree");

  TNebulaPhysics *Neb_Phys = new TNebulaPhysics();
  std::vector<double> Charge;
  //std::vector<int>* Efficiency = 0;
  double Ch, counter;
  //int Eff;

  //TBranch* branch_Eff;
  TBranch* branch ;
  //TH1F *h_Eff = new TH1F("h_Eff","test");
  TH1F *h_En = new TH1F("h_En","Charge",160,0,160);
  //t->SetBranchAddress("Efficiency",&Efficiency,&branch_Eff);
  t->SetBranchAddress("Nebula",&Neb_Phys,&branch);
  //int nentries_Eff = branch_Eff->GetEntries();
  int nentries_Neb = branch->GetEntries();
  
  for(int i = 0; i < nentries_Neb ; i++){
    t->GetEntry(i);
    Charge = Neb_Phys->Charge;
    int vec_size = Charge.size();
    for(int j = 0; j < vec_size; j++){
      Ch = Charge[j];
      if(Ch>160) cout << "Ch " << Ch << endl;
      h_En->Fill(Ch);
    }
  }

  TCanvas *c2 = new TCanvas("c2","Charge",1080,720); 
  h_En->Draw("");
}
