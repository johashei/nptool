void MakeHisto(){

  auto file = new TFile("../../../Outputs/Analysis/s034_test.root","READ");
  auto tree = (TTree*) file->FindObjectAny("PhysicsTree");
  auto out = new TFile("TPad_distrib_new.root","RECREATE");
  for(unsigned int i = 1 ; i < 2 ; i++){
    cout << i << endl;
   tree->Draw(Form("T_Pad>>TPad_ring%d(2000,0,400)",i+1),Form("Ring_Pad==%d",i+1));
  // auto h = (TH1F*) gDirectory->FindObjectAny(Form("TPad_ring%d",i+1));
  }
  out->Write();
  out->Close();
}
