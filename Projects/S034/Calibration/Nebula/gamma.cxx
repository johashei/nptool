void gamma(){

  auto chain = new TChain("PhysicsTree");
  chain->Add("root/analysis/gamma/run418.root");

  chain->Draw("(Nebula.PosZ+4650)/Nebula.TOF:Nebula.DetectorNumber>>h(150,0,150,2000,0,500)","","colz");

  }
