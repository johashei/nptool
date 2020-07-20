void Control(){
  TFile* file = TFile::Open("../../Outputs/Analysis/PhysicsTree.root");
  TTree* PhysicsTree= (TTree*) file->FindObjectAny("PhysicsTree");
  string cond = "InnerPosX!=-1000";
/*
  TCanvas* cInner = new TCanvas("ControlInner","ControlInner",1000,1000);
  cInner->Divide(2,2);
  cInner->cd(1);
  PhysicsTree->Draw("InnerPosX:fDetected_Position_X",cond.c_str(),"col") ; 
  cInner->cd(2);
  PhysicsTree->Draw("InnerPosY:fDetected_Position_Y",cond.c_str(),"col") ; 
  cInner->cd(3);
  PhysicsTree->Draw("InnerPosZ:fDetected_Position_Z",cond.c_str(),"col") ; 
  cInner->cd(4);
  PhysicsTree->Draw("InnerPosY:InnerPosX:InnerPosZ", cond.c_str(),"");            
*/
  TCanvas* cOuter = new TCanvas("ControlOuter","ControlOuter",1000,1000);
  cOuter->Divide(2,2);
  cond = "OuterPosX!=-1000";
  cOuter->cd(1);
  PhysicsTree->Draw("OuterPosX:fDetected_Position_X[3]",cond.c_str(),"col") ; 
  cOuter->cd(2);
  PhysicsTree->Draw("OuterPosY:fDetected_Position_Y[3]",cond.c_str(),"col") ; 
  cOuter->cd(3);
  PhysicsTree->Draw("OuterPosZ:fDetected_Position_Z[3]",cond.c_str(),"col") ; 
  cOuter->cd(4);
  PhysicsTree->Draw("OuterPosY:OuterPosX:OuterPosZ", cond.c_str(),"");            


  TCanvas* c2 = new TCanvas("Control 2", "Control2",500,500,2000,1000);  
  c2->Divide(2,1);
  c2->cd(1);
  PhysicsTree->Draw("OuterPosY:OuterPosX",cond.c_str());            
  PhysicsTree->Draw("InnerPosY:InnerPosX",cond.c_str(),"same");            

  c2->cd(2);
  PhysicsTree->Draw("OuterPosY:OuterPosZ",cond.c_str());            
  PhysicsTree->Draw("InnerPosY:InnerPosZ",cond.c_str(),"same");            
}
