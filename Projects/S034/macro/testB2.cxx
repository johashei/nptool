void testB2(){

  SamuraiFieldMap field;
  double angle = 30*deg;
  field.LoadMap(angle,"field_map/180702-2,40T-3000.table.bin",10);
  field.SetFDC2Angle((59.930-90.0)*deg);
  field.SetFDC2R(3686.77 + 880.745/2.);
  auto scan  = field.BrhoScan(2.5,10,0.1);
  auto h=new TH1D("dB","dB",1000,-2e-3,2e-3); 
  double r = gRandom->Uniform(-1,1);
  for(unsigned int i = 0 ; i < 1000 ; i++){
  TVector3 p(gRandom->Uniform(-5,5),gRandom->Uniform(-5,5),-3500); TVector3 d(gRandom->Uniform(-0.05,0.05),gRandom->Uniform(-0.05,0.05),1); 
  double b = gRandom->Uniform(3,7);
  std::vector< TVector3 > track = field.Propagate(3000,b,p,d);
  h->Fill(field.FindBrho(p,d,track.back(),d)-b);
  }
  
  h->Draw();
} 
