void CheckFDC(){
   auto tree = new TChain("PhysicsTree");
   tree->Add("../../Outputs/Analysis/PhysicsTree.root");
   
   TString Cond="SamuraiFDC0.PosX!=-10000 && SamuraiFDC0.PosY!=-10000 && SamuraiFDC0.PileUp==0 && SamuraiFDC2.PosX!=-10000 && SamuraiFDC2.PosY!=-10000 && SamuraiFDC2.PileUp==0";

   auto c = new TCanvas();
   tree->Draw("SamuraiFDC0.PosX:SamuraiFDC2.PosX>>hX(500,-1500,1500,500,-100,100)",Cond,"col");

   c = new TCanvas();
   tree->Draw("SamuraiFDC0.PosY:SamuraiFDC2.PosY>>hY(500,-1500,1500,500,-100,100)",Cond,"col");
 
   new TCanvas();
   tree->Draw("SamuraiFDC0.ThetaX*180./3.14159:SamuraiFDC2.ThetaX*180./3.14159>>hT(500,-50,50,500,-50,50)",Cond,"col");

   new TCanvas();
   tree->Draw("SamuraiFDC2.ThetaX*180./3.14159:SamuraiFDC2.PosX>>hXT(500,-1500,1500,500,-50,50)",Cond,"col");
 
  }
