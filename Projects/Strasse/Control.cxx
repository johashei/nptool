void Control(){
  TFile* file = TFile::Open("../../Outputs/Analysis/PhysicsTree.root");
  TTree* PhysicsTree= (TTree*) file->FindObjectAny("PhysicsTree");
  TCanvas* c = new TCanvas("Control","Control",1000,1000);
  c->Divide(2,2);
  //string cond = "InnerPosX!=-60 && InnerStripL[0]==64  && InnerStripT[0]==152";

  string cond = "InnerPosX!=-60 && (DetectorNumber==3 || DetectorNumber==1)";
  c->cd(1);
//  PhysicsTree->Draw("InnerPosX:fDetected_Position_X",cond.c_str(),"col") ; 
//
  PhysicsTree->Draw("InnerPosX-fDetected_Position_X>>h(1000,-2,2)",cond.c_str(),"col") ; 
  c->cd(2);
  //PhysicsTree->Draw("InnerPosY:fDetected_Position_Y",cond.c_str(),"col") ; 

  PhysicsTree->Draw("InnerPosY-fDetected_Position_Y>>h2(1000,-2,2)",cond.c_str(),"col") ; 
  c->cd(3);
  //PhysicsTree->Draw("InnerPosZ:fDetected_Position_Z",cond.c_str(),"col") ; 

  PhysicsTree->Draw("InnerPosZ-fDetected_Position_Z>>h3(1000,-2,2)",cond.c_str(),"col") ; 
  c->cd(3);
  c->cd(4);
  PhysicsTree->Draw("InnerPosY:InnerPosX:InnerPosZ", cond.c_str(),"");            

 /* TCanvas* c2 = new TCanvas("Control 2", "Control2",500,500,2000,1000);  
  c2->Divide(2,1);
  c2->cd(1);
  PhysicsTree->Draw("InnerPosY:InnerPosX","InnerPosX!=-60");            

  c2->cd(2);
  PhysicsTree->Draw("fDetected_Position_Y:fDetected_Position_X","InnerPosX!=-60","");            
*/
}
