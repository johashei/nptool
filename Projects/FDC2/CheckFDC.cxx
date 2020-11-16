void CheckFDC(){
   auto tree = new TChain("PhysicsTree");
   tree->Add("../../Outputs/Analysis/PhysicsTree.root");
   
   TString Cond="PosX!=-10000&&PosY!=-10000&&PileUp==0&&devX<10&&devY<20";

   auto c = new TCanvas();
   tree->Draw("PosY:PosX>>h(1500,-1500,1500,500,-500,500)",Cond,"col");
   //c->SaveAs("~/XY.pdf");
   new TCanvas();
   tree->Draw("ThetaX*180./3.14159:PosX>>hXT(500,-1200,1200,500,-50,50)",Cond,"col");
   //c->SaveAs("~/ThetaX.pdf");
   //
   new TCanvas();
   tree->Draw("PhiY*180./3.14159:PosY>>hY(500,-1200,1200,500,-50,50)",Cond,"col");
   //c->SaveAs("~/PhiY.pdf)");
   //
   new TCanvas();
   tree->Draw("PosX>>hX(1500,-1500,1500)",Cond,"col");
   //c->SaveAs("~/X.pdf");

  }
