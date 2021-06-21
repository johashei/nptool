TChain* MakeChain1();
TChain* MakeChain2();
TChain* MakeChain();
TH2F*   MakeTH2();
TH2F*   GetTH2();
TGraph* graph = new TGraph(200);
double off;
double c_light=299.792458;//mm/ns
auto chain = MakeChain();
void process1bar(int b);
auto h = new TH2F("h","h",200,0,200,500,0,1000);
auto r = new TH2F("r","r",200,0,200,10000,17000,20000);
ofstream output("Calibration/Nebula/offset_gamma.txt");
////////////////////////////////////////////////////////////////////////////////
void gamma(){
  double beta=0.5504;
  // Distance from SBT(-7377.56) to target(3774.7 for FP + a value for each target)
  double D1=7377.56-3774.7+2.795;
  double D2=7377.56-3774.7+6.972;
  // time offset for each case
  double off1 = D1/(beta*c_light);
  double off2 = D2/(beta*c_light);
  off = (off1+off2)*0.5;
  cout << "Offset1 : " << off1 << endl;
  cout << "Offset2 : " << off2 << endl;
  cout << "Mean : " << off << endl;
  auto c = new TCanvas("tof","tof");

  chain->SetAlias("R","sqrt(Nebula.PosX*Nebula.PosX+Nebula.PosY*Nebula.PosY+(Nebula.PosZ+3774.7)*(Nebula.PosZ+3774.7))");
  h=GetTH2();
  h->Draw("colz");
  new TCanvas();
  for(unsigned int i = 0 ; i < 160 ; i++){
    if(i!=60)
    process1bar(i); 
  }

  process1bar(60); 
  output.close();
  new TCanvas();
  graph->Draw("ap");
}

////////////////////////////////////////////////////////////////////////////////
void process1bar(int b){
  //new TCanvas();
  // Get the Radius for the distance to this barre
  auto r1 = r->ProjectionY(Form("h%d",b),b,b+1);  
  double R =  r1->GetBinCenter(r1->GetMaximumBin());

  auto h1 = h->ProjectionY(Form("h%d",b),b,b+1);
  //h1->Rebin(4);
  double max = h1->GetBinCenter(h1->GetMaximumBin());
  //h1->Draw();
  auto f = new TF1("f","crystalball(0)+pol1(5)",max-50,max+50);
  f->SetParameter(0,h1->GetMaximum());
  f->SetParameter(1,max);
  f->SetParameter(2,35);
  f->SetParameter(3,0.1);
  f->SetParameter(4,1);

  h1->Fit(f,"W");
  
    // Vbad = R/(TOF) -> TOF/R = 1/Vbad
    // c= R/(TOF+X)   -> (TOF+X)/R = 1/c -> TOF/R+X/R=1/c
    // 1/Vbad +X/R = 1/c
    // X=R*(1/c-1/Vbad)
    
  double offset=R*(1/c_light-1/f->GetParameter(1)) ;
 // double offset=R*(1/c_light-1/max) ;
  
  cout <<f->GetParameter(1) << " " <<  offset << " " << R/(f->GetParameter(1)-offset) << endl;
  if(offset>0){
    output << "NEBULA_T_ID"  << b << " " << offset << endl; 
    graph->SetPoint(b,b,offset);
    }
}
////////////////////////////////////////////////////////////////////////////////
TH2F* GetTH2(){
  auto File= new TFile("Calibration/Nebula/hist_v.root");
  h = (TH2F*) File->FindObjectAny("h");

  File= new TFile("Calibration/Nebula/hist_r.root");
  r = (TH2F*) File->FindObjectAny("r");
  return h;
  }

////////////////////////////////////////////////////////////////////////////////
TH2F* MakeTH2(){
  TString cond=Form("(Nebula.TOF-%f)>20&&(Nebula.TOF-%f)<38",off,off);
  TString draw=Form("R/(Nebula.TOF-%f):Nebula.DetectorNumber>>h",off); 
  chain->Draw(draw,cond,"colz");
  h->SaveAs("Calibration/Nebula/hist_v.root");

  chain->Draw("R:Nebula.DetectorNumber>>r","","colz");
  r->SaveAs("Calibration/Nebula/hist_r.root");

  return h;
  }

////////////////////////////////////////////////////////////////////////////////
TChain* MakeChain(){
  auto chain = new TChain("PhysicsTree");
  chain->Add("root/analysis/gamma/run*.root");
  return chain;
}

////////////////////////////////////////////////////////////////////////////////
TChain* MakeChain1(){
  auto chain = new TChain("PhysicsTree");
  for(unsigned int i = 418 ; i < 442  ; i++){
    chain->Add(Form("root/analysis/gamma/run%d.root",i));
  }
  return chain;
}

////////////////////////////////////////////////////////////////////////////////
TChain* MakeChain2(){
  auto chain = new TChain("PhysicsTree");
  for(unsigned int i = 488 ; i < 498  ; i++){
    chain->Add(Form("root/analysis/gamma/run%d.root",i));
  }
  return chain;
}
