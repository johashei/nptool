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

  TCanvas* cVertex = new TCanvas("ControlVertex","ControlVertex",1000,1000);
  cVertex->Divide(2,2);
  cond = "VertexX!=-1000";
  cVertex->cd(1);
  PhysicsTree->Draw("VertexX-fRC_Vertex_Position_X>>hx(500,-2,2)",cond.c_str(),"col") ; 
  cVertex->cd(2);
  PhysicsTree->Draw("VertexY-fRC_Vertex_Position_Y>>hy(500,-2,2)",cond.c_str(),"col") ; 
  cVertex->cd(3);
  PhysicsTree->Draw("VertexZ-fRC_Vertex_Position_Z>>hz(5000,-20,20)",cond.c_str(),"col") ; 
  cVertex->cd(4);

  PhysicsTree->Draw("Distance>>hd(500,0,80)", cond.c_str(),"");            
  //PhysicsTree->Draw("VertexY:VertexX:VertexZ", cond.c_str(),"");            

  TCanvas* cdelta= new TCanvas("ControlDelta","ControlDelta",1000,1000);
  cdelta->Divide(2,2);
  cond = "deltaX!=-1000";
  cdelta->cd(1);
  PhysicsTree->Draw("deltaX>>dx(500,-0.01,0.01)",cond.c_str(),"col") ; 
  cdelta->cd(2);
  PhysicsTree->Draw("deltaY>>dy(500,-0.01,0.01)",cond.c_str(),"col") ; 
  cdelta->cd(3);
  PhysicsTree->Draw("deltaZ>>dz(500,-0.01,0.01)",cond.c_str(),"col") ; 
  cdelta->cd(4);
  PhysicsTree->Draw("Distance>>hd(500,0,0.01)", cond.c_str(),"");           
*/
  TCanvas* ctheta= new TCanvas("ControlTheta","ControlTheta",1000,1000);
  cond = "Theta12!=-1000";
  PhysicsTree->Draw("Theta12>>ht",cond.c_str(),"col") ; 
 

/*
  TCanvas* c2 = new TCanvas("Control 2", "Control2",500,500,2000,1000);  
  c2->Divide(2,1);
  c2->cd(1);
  PhysicsTree->Draw("VertexY:VertexX:VertexZ",cond.c_str());            
  c2->cd(2);

*/
}
