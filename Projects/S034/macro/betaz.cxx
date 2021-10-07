void betaz(){

  auto fz = new TFile("root/zaihong/run0582_RIG20210424_6He.root");
  auto tz = (TTree*) fz->FindObjectAny("rig");
  auto fl = new TFile("root/analysis/Result582.root");
  auto tl = (TTree*) fl->FindObjectAny("PhysicsTree");
  tl->AddFriend(tz); 
  double Brho,betaZ,beta;
  int    FragID;
  tl->SetBranchAddress("Brho",&Brho);
  tl->SetBranchAddress("fBeta",&betaZ);
  tl->SetBranchAddress("FragID",&FragID);
  NPL::Particle H2("2H");
  NPL::Particle H3("3H");
  NPL::Particle He4("4He");
  NPL::Particle He6("6He");

  auto z = new TH1D("betaZ","betaZ",100,0.45,0.6);
  auto b = new TH1D("betaF","betaZ",100,0.45,0.6);
  auto d = new TH1D("betaF","betaZ",1000,-0.01,0.01);
  unsigned int entries = tl->GetEntries();

  for(unsigned int i = 0 ; i < entries; i++){
    tl->GetEntry(i);
    if(FragID>0 && FragID<27){
      // compute Brho based on beta and FragID
      double rig ;
      if(FragID==26){
        He6.SetBrho(Brho);
        beta= He6.GetBeta();
        b->Fill(beta);
        z->Fill(betaZ);
        d->Fill(beta-betaZ);
      }
    }
  }
  //  h->Scale(1./h->Integral());
  //  h->Draw(); h->SetLineColor(kBlack);
  b->Draw(""); 
  b->Scale(1./b->Integral());
  b->SetLineColor(kAzure+7);b->SetLineWidth(2);
  z->Scale(1./z->Integral());
  z->SetLineColor(kOrange+7);z->SetLineWidth(2); 
  z->Draw("same");
  new TCanvas();
  d->Draw();
  d->SetLineColor(kAzure+7);b->SetLineWidth(2);
}

