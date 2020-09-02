void Shift(){
  
  TFile* file_ok = TFile::Open("../../Outputs/Analysis/strasse_ok.root");
  TFile* file_shifted = TFile::Open("../../Outputs/Analysis/strasse_shifted.root");
  TTree* ok= (TTree*) file_ok->FindObjectAny("PhysicsTree");
  TTree* shifted= (TTree*) file_shifted->FindObjectAny("PhysicsTree");
 
  string cond = "sqrt(VertexX*VertexX+VertexY*VertexY)<20 && VertexZ>-80 && VertexZ<80";
  // Vertex 
  TCanvas* cVertex = new TCanvas("ControlVertex","ControlVertex",1000,1000);
  cVertex->Divide(2,2);
  cond = "VertexX!=-1000";
  cVertex->cd(1);
  ok->Draw("VertexX-fRC_Vertex_Position_X>>hxok(500,-2,2)",cond.c_str(),"") ; 
  shifted->Draw("VertexX-fRC_Vertex_Position_X>>hx(500,-2,2)",cond.c_str(),"same") ; 
  TH1* hx= (TH1*) gDirectory->FindObjectAny("hx");
  hx->SetLineColor(kOrange-3);
  hx->SetFillColor(kOrange-3); 

  TH1* hxok= (TH1*) gDirectory->FindObjectAny("hxok");
  hxok->GetXaxis()->SetTitle("X_{reconstructed}-X_{real} (mm)");

  cVertex->cd(2);

  ok->Draw("VertexY-fRC_Vertex_Position_Y>>hyok(500,-2,2)",cond.c_str(),"") ; 
  shifted->Draw("VertexY-fRC_Vertex_Position_Y>>hy(500,-2,2)",cond.c_str(),"same") ; 

  TH1* hy= (TH1*) gDirectory->FindObjectAny("hy");
  hy->SetLineColor(kOrange-3);
  hy->SetFillColor(kOrange-3); 
  TH1* hyok= (TH1*) gDirectory->FindObjectAny("hyok");
  hyok->GetXaxis()->SetTitle("Y_{reconstructed}-Y_{real} (mm)");


  cVertex->cd(3);

  ok->Draw("VertexZ-fRC_Vertex_Position_Z>>hzok(2000,-20,20)",cond.c_str(),"") ; 
  shifted->Draw("VertexZ-fRC_Vertex_Position_Z>>hz(2000,-20,20)",cond.c_str(),"same") ; 

  TH1* hz= (TH1*) gDirectory->FindObjectAny("hz");
  hz->SetLineColor(kOrange-3);
  hz->SetFillColor(kOrange-3); 
  TH1* hzok= (TH1*) gDirectory->FindObjectAny("hzok");
  hzok->GetXaxis()->SetTitle("Z_{reconstructed}-Z_{real} (mm)");


  cVertex->cd(4);
  ok->Draw("Distance>>hdok(500,0,2)", cond.c_str(),"");            
  shifted->Draw("Distance>>hd(500,0,2)", cond.c_str(),"same");            
  TH1* hd= (TH1*) gDirectory->FindObjectAny("hd");
  hd->SetLineColor(kOrange-3);
  hd->SetFillColor(kOrange-3); 
  TH1* hdok= (TH1*) gDirectory->FindObjectAny("hdok");
  hdok->GetXaxis()->SetTitle("Minimum track distance (mm)");



  // Angle
  TCanvas* ctheta= new TCanvas("ControlTheta","ControlTheta",2000,1000);
  ctheta->Divide(2,1);
  ctheta->cd(1);
  ok->Draw("Theta12>>ht(2000)",cond.c_str(),"") ; 
  shifted->Draw("Theta12>>hts(2000)",cond.c_str(),"same") ; 
  TH1* hts= (TH1*) gDirectory->FindObjectAny("hts");
  hts->SetFillColor(kOrange-3);
  hts->SetLineColor(kOrange-3);
  ctheta->cd(2);
  cond = "deltaPhi!=-1000";
  ok->Draw("deltaPhi>>hp(2000)",cond.c_str(),"") ; 
  shifted->Draw("deltaPhi>>hps(2000)",cond.c_str(),"same") ; 
  TH1* hps= (TH1*) gDirectory->FindObjectAny("hps");
  hps->SetFillColor(kOrange-3);
  hps->SetLineColor(kOrange-3);


}
