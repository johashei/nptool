void CheckFDC0(){
   auto tree = new TChain("PhysicsTree");
   tree->Add("../../Outputs/Analysis/PhysicsTree.root");
   
   TString Cond="SamuraiFDC0.PosX!=-10000 && SamuraiFDC0.PosY!=-10000 && SamuraiFDC0.PileUp==0";

   //TString Cond="SamuraiFDC0.PosX!=-10000 && SamuraiFDC0.PosY!=-10000";

   auto c = new TCanvas();
   tree->Draw("SamuraiFDC0.PosY:SamuraiFDC0.PosX>>hXY(500,-100,100,500,-100,100)",Cond,"col");
   //c->SaveAs("~/XY.pdf");
   new TCanvas();
   tree->Draw("SamuraiFDC0.ThetaX*180./3.14159:SamuraiFDC0.PosX>>hXT(500,-100,100,500,-50,50)",Cond,"col");
   //c->SaveAs("~/ThetaX.pdf");
   //
   new TCanvas();
   tree->Draw("SamuraiFDC0.PhiY*180./3.14159:SamuraiFDC0.PosY>>hYP(500,-100,100,500,-50,50)",Cond,"col");
   //c->SaveAs("~/PhiY.pdf)");
   //
   new TCanvas();
   tree->Draw("SamuraiFDC0.PosX>>hX(500,-100,100)",Cond,"col");
   //c->SaveAs("~/X.pdf");
   new TCanvas();
   tree->Draw("SamuraiFDC0.PosY>>hY(500,-100,100)",Cond,"col");
   //c->SaveAs("~/X.pdf");


  }
