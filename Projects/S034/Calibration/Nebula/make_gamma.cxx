
TChain* MakeChain();

double c_light=299.792458;//mm/ns
void make_gamma(){

  auto file = new TFile("Calibration/Nebula/hist_gamma.root","RECREATE");
  // Create histograms
  vector<TH1F*> h,r;
  for(unsigned int i = 0 ; i < 150 ; i++){
    TString name = Form("h%d",i);
    h.push_back(new TH1F(name,name,500,0,1000));
    name = Form("r%d",i);
    r.push_back(new TH1F(name,name,10000,17000,20000));
  }

  auto chain = MakeChain();
  auto neb = new TNebulaPhysics();

  double beta=0.5504;
  // Distance from SBT(-7377.56) to target(3774.7 for FP + a value for each target)
  double D1=7377.56-3774.7+2.795;
  double D2=7377.56-3774.7+6.972;
  // time offset for each case
  double off1 = D1/(beta*c_light);
  double off2 = D2/(beta*c_light);
  double off = (off1+off2)*0.5;

  cout << "Offset1 : " << off1 << endl;
  cout << "Offset2 : " << off2 << endl;
  cout << "Mean : " << off << endl;

  chain->SetBranchAddress("Nebula",&neb);
  unsigned int size = chain->GetEntries();

 // unsigned int size = 10000;
  for(unsigned int i = 0 ; i < size ; i++){
    if(i%1000==0)
      cout << "\r entry " << i << " / " << size << flush; 
    chain->GetEntry(i); 
    if(neb->DetectorNumber.size()==1){
      double R = sqrt(neb->PosX[0]*neb->PosX[0]+neb->PosY[0]*neb->PosY[0]+(neb->PosZ[0]+3774.7)*(neb->PosZ[0]+3774.7));
      double v = R/(neb->TOF[0]-off); 
      h[neb->DetectorNumber[0]]->Fill(v);
      r[neb->DetectorNumber[0]]->Fill(R);
    }
  }
 
 for(auto it = h.begin(); it !=h.end(); it++)
   (*it)->Write();
  for(auto it = r.begin(); it !=r.end(); it++)
   (*it)->Write();
  
  file->Close();
cout << endl; 
}
////////////////////////////////////////////////////////////////////////////////
TChain* MakeChain(){
  auto chain = new TChain("PhysicsTree");
  chain->Add("root/analysis/gamma/run*.root");
  return chain;
}
