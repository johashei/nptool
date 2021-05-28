
SamuraiFieldMap field;

void DrawT(TVector3 pos, TVector3 dir, double Brho){
std::vector< TVector3 > track = field.Propagate(3000,Brho,pos,dir);
  auto g = new TGraph();
  unsigned int size = track.size();
  g->Set(size);
  for(unsigned int i = 0 ; i < size ; i++){
    g->SetPoint(i,-track[i].X(),track[i].Z());
  }
  g->SetLineWidth(6);
  g->Draw("l");
}
////////////////////////////////////////////////////////////////////////////////
void testB(){
  auto c = new TCanvas("trajectory","trajectory",1000,1000);
  double angle = 30*deg;
  field.LoadMap(angle,"field_map/180702-2,40T-3000.table.bin",10);
  field.SetFDC2Angle((59.930-90.0)*deg);
  field.SetFDC2R(3686.77 + 880.745/2.);
  auto scan  = field.BrhoScan(2.5,10,0.1);
  
  unsigned int size = 1000;
  vector<float> pos = {0,0,0};
  auto h = new TH2F("h","h", size,-4000,4000,size,-4000,4000);
  for(int x = -3000 ; x < 3000 ; x+=10){
    for(int z = -3000 ; z < 3000 ; z+=10){
      pos[0]=x;pos[1]=0;pos[2]=z; 
      TVector3 p(x,0,z);
      //float b = 1000*field.GetB(pos)[2];
      float b = 1000*field.InterpolateB(p)[1];
      p.RotateY(angle);
      if(b!=0)
        h->Fill(p.X(),p.Z(),b);
    }
  }
  h->Draw("colz");
  h->GetZaxis()->SetRangeUser(-1,3);
  DrawT(TVector3(0,0,-3500),TVector3(0,0,1),5.48);
  DrawT(TVector3(0,0,-3500),TVector3(0,0,1),3.62);

  TVector3 p(1,1,-3500); TVector3 d(0.01,-0.01,1); double b = 5.481923;
  DrawT(p,d,b);
  std::vector< TVector3 > track = field.Propagate(3000,b,p,d);
  cout << field.FindBrho(p,d,track.back(),d)<< endl;;

  double rFDC2 = 3686.77 + 880.745/2.;
  double phiFDC2 = (59.930-90.0-angle/deg)*deg;
  TVector3 C(-252 ,0, rFDC2 );
  TVector3 L(-252+1308 ,0, rFDC2 );
  TVector3 R(-252-1308, 0, rFDC2 );
  C.RotateY(phiFDC2);
  L.RotateY(phiFDC2);
  R.RotateY(phiFDC2);
  auto FDC2 = new TLine(-L.X(),L.Z(),-R.X(),R.Z());
  FDC2->SetLineWidth(4);
  FDC2->SetLineColor(kRed);
  FDC2->Draw();

  auto  CFDC2= new TMarker(-C.X(),C.Z(),20);
  C.Print();
  CFDC2->SetMarkerColor(kRed);
  CFDC2->Draw();

  auto mag = new TEllipse (0,0, 1000, 0);
  mag->SetLineColor(kBlue);
  mag->SetLineWidth(4);
  mag->SetFillStyle(0);
  mag->Draw();

  new TCanvas();
  scan->Draw("ac");
  cout << scan->Eval(700-252) << endl;
  
  

}
