#include "/local/lemair/nptool/NPLib/Detectors/Nebula/TNebulaPhysics.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::find
#include <vector>   


void PositionZ(){
  TFile* f = TFile::Open("root/simulation/SimulatedTree.root");
  TTree* t = (TTree*)f->Get("SimulatedTree");
  TCanvas *c1 = new TCanvas("c1","Position Z",1080,720);

  vector<double>* PosZ;
  int vec_size;
  TBranch* branch;
  
  TH1F *h = new TH1F("h","Position Z",2400,10800,12200);
  t->SetBranchAddress("Pos_Z", &PosZ, &branch);
  
  int nentries = t->GetEntries();
  for(int i=0 ; i<nentries ; i++){
    t->GetEntry(i);
    vec_size = (*PosZ).size();
    for(int j = 0 ; j<vec_size ; j++)
      h->Fill((*PosZ)[j]);
  }
  
  h->Draw("SAME");
}
 
