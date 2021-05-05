void check(){

  auto file = new TFile("../../Outputs/Analysis/test582.root","READ");
  auto tree = (TTree*) file->FindObjectAny("PhysicsTree");
  tree->Print();
  auto fdc2 = new TSamuraiFDC2Physics();

  tree->SetBranchAddress("SamuraiFDC2", &fdc2);

  auto c = new TCanvas();
  unsigned int entries = tree->GetEntries();
  for(unsigned int e = 0 ; e < entries ; e++){
    tree->GetEntry(e);
    unsigned int mult = fdc2->DriftLength.size();      
    for(unsigned int i = 0 ; i < mult ; i++){
      double r = fdc2->DriftLength[i]
      double r = fdc2->DriftLength[i]
      double r = fdc2->DriftLength[i]
    }
  }


}
