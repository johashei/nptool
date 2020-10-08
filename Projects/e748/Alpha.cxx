void Alpha(){
  
  // Load the Main Tree
  TFile* file = new TFile("../../Outputs/Analysis/e748_3A.root");
  TTree* tree = (TTree*) file->FindObjectAny("PhysicsTree");

  // TOF per telescope
  for(unsigned int i = 1 ; i < 5 ; i++){
      new TCanvas();
      tree->Draw(Form("Si_EY>>hTOFc%i(1000,0,8)",i),Form(" TelescopeNumber==%i",i),"colz");
    }

}
