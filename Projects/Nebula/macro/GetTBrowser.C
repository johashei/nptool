void GetTBrowser(){
  
  TH1F *h1 = new TH1F("h1","Somme des Temps",100,-4,4);

  TFile* f = TFile::Open("root/simulation/SimulatedTree.root");
  TTree *t = (TTree*)f->Get("SimulatedTree");
  TBranch *b = (TBranch*)t->GetBranch("Nebula");
  TLeaf *tu = b->GetLeaf("fNebula_Tu_Time");
  //Float_t td;
  //b->SetAddress("fNebula_Td_Time", &td);
  Int_t len = b->GetEntries();
  Float_t sum;
  for(int i = 0; i < len; i++){
    sum = 0;
    sum += tu->GetValue(i);
    //sum += td;
    h1->Fill(sum);
  }
}
