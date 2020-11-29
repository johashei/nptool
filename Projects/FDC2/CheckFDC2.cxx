void CheckFDC2(){
   auto tree = new TChain("PhysicsTree");
   tree->Add("../../Outputs/Analysis/PhysicsTree.root");
   
   TString Cond="SamuraiFDC2.PosX!=-10000 && SamuraiFDC2.PosY!=-10000 && SamuraiFDC2.PileUp==0";

   //TString Cond="SamuraiFDC2.PosX!=-10000 && SamuraiFDC2.PosY!=-10000";

   auto c = new TCanvas();
   tree->Draw("SamuraiFDC2.PosY:SamuraiFDC2.PosX>>hXY(500,-1500,1500,500,-500,500)",Cond,"col");
   //c->SaveAs("~/XY.pdf");
   new TCanvas();
   tree->Draw("SamuraiFDC2.ThetaX*180./3.14159:SamuraiFDC2.PosX>>hXT(500,-1500,1500,500,-50,50)",Cond,"col");
   //c->SaveAs("~/ThetaX.pdf");
   //
   new TCanvas();
   tree->Draw("SamuraiFDC2.PhiY*180./3.14159:SamuraiFDC2.PosY>>hYP(500,-500,500,500,-50,50)",Cond,"col");
   //c->SaveAs("~/PhiY.pdf)");
   //
   new TCanvas();
   tree->Draw("SamuraiFDC2.PosX>>hX(500,-1500,1500)",Cond,"col");
   //c->SaveAs("~/X.pdf");
   new TCanvas();
   tree->Draw("SamuraiFDC2.PosY>>hY(500,-500,500)",Cond,"col");
   //c->SaveAs("~/X.pdf");


  }
