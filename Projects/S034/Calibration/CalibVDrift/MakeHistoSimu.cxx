
void MakeHistoSimu(){
  auto file  = new TFile("ZPad_distrib_d5mm.root","RECREATE");
  auto chain = new TChain("PhysicsTree");
  chain->Add("../../root/analysis/simu_tpad_d5mm*.root");
  
  for(unsigned int i = 1 ; i < 19 ; i++){
    cout << "\r " << i << endl;
    chain->Draw(Form("Z_Pad>>ZPad_Ring%d(1000,0,600)",i),Form("Ring_Pad==%d&&Dali.Energy@.size()>0",i));
    }
  file->Write();
}
