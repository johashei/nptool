TChain* chain=NULL;

//////////////////////////////////////////////////////////////////
void LoadRootFile(){
  chain = new TChain("SimulatedTree");

  chain->Add("../../../Outputs/Simulation/sofia_simu.root");
  
  //chain->Add("../../../Outputs/Simulation/sofia_simu_1.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_2.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_3.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_4.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_5.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_6.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_7.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_8.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_9.root");
  //chain->Add("../../../Outputs/Simulation/sofia_simu_10.root");
}

//////////////////////////////////////////////////////////////////
void ShowResultSimu()
{
  LoadRootFile();

  chain->Draw("fDetected_Position_Y:fDetected_Position_X>>hpos(500,-1900,-900,500,-350,350)","fTOF_Energy@.size()==2","colz");
  TH2F* hpos = (TH2F*)gDirectory->FindObjectAny("hpos");

  chain->Draw("fFC_Fragment_Z>>h1(35,30,65)");
  TH1F* h1 = (TH1F*)gDirectory->FindObjectAny("h1");

  chain->Draw("fDetected_Z>>h2(35,30,65)","fTOF_Energy@.size()==2");
  TH1F* h2 = (TH1F*)gDirectory->FindObjectAny("h2");

  h1->Sumw2();
  h2->Sumw2();

  TH1F* hratio = (TH1F*)h2->Clone();
  hratio->Divide(h1);

  TCanvas* c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(2,2);

  c1->cd(1);
  hpos->Draw("colz");
  c1->cd(2);
  h1->Draw();
  c1->cd(3);
  h2->Draw();
  c1->cd(4);
  hratio->Draw("E1");


}
