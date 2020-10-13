void Target(){
  
  // Load the Main Tree
  TFile* file = new TFile("../../Outputs/Analysis/e748_Mask.root");
  TTree* tree = (TTree*) file->FindObjectAny("PhysicsTree");


  new TCanvas();
   tree->Draw("PositionOnTargetY:PositionOnTargetX>>hT1(1000,-30,30,1000,-30,30)","","colz");  
   TMarker* star = new TMarker (-1.4,-0.1,29);
   star->SetMarkerSize(2);
   star->SetMarkerColor(kRed);
   star->Draw();

}
