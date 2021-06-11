void gamma(){

  auto chain = new TChain("PhysicsTree");
  chain->Add("root/analysis/gamma/run418.root");

  chain->Draw("Nebula.Time:Nebula.DetectorNumber>>h(150,0,150,1000,0,150)","","colz");

  }
