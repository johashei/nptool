TChain* chain=NULL;

//////////////////////////////////////////////////////////////////
void LoadRootFile(string nucleus){
  chain = new TChain("SimulatedTree");

  string rootfilename = "../../../Outputs/Simulation/sofia_simu_" + nucleus + "_*";

  chain->Add(rootfilename.c_str());
}

//////////////////////////////////////////////////////////////////
void ShowResultSimu(string nucleus="182Hg")
{
  LoadRootFile(nucleus);

  gROOT->SetStyle("pierre_style");

  chain->Draw("fDetected_Position_Y:fDetected_Position_X>>hpos(200,-1900,-900,200,-400,400)","fTwim_AnodeNbr@.size()==32 && fTOF_Energy@.size()>=0","colz");
  TH2F* hpos = (TH2F*)gDirectory->FindObjectAny("hpos");
  hpos->GetXaxis()->SetTitle("X_{MWPC3} (mm)");
  hpos->GetYaxis()->SetTitle("Y_{MWPC3} (mm)");
  hpos->GetXaxis()->CenterTitle();
  hpos->GetYaxis()->CenterTitle();


  chain->Draw("fFC_Fragment_Z>>h1(45,20,65)");
  TH1F* h1 = (TH1F*)gDirectory->FindObjectAny("h1");
  h1->GetXaxis()->SetTitle("Simulated Charge Z");
  h1->GetXaxis()->CenterTitle();
    

  //chain->Draw("fDetected_Z>>h2(40,25,65)","fTOF_Energy@.size()==2");
  chain->Draw("fDetected_Z>>h2(45,20,65)","fTOF_Energy@.size()==2 && fTwim_AnodeNbr@.size()==32");
  TH1F* h2 = (TH1F*)gDirectory->FindObjectAny("h2");
  h2->GetXaxis()->SetTitle("Detected Charge Z");
  h2->GetXaxis()->CenterTitle();

  h1->Sumw2();
  h2->Sumw2();

  TH1F* hratio = (TH1F*)h2->Clone();
  hratio->Divide(h1); 
  hratio->GetXaxis()->SetTitle("Charge Z");
  hratio->GetYaxis()->SetTitle("Detection efficiency");
  hratio->GetXaxis()->CenterTitle();
  hratio->GetYaxis()->CenterTitle();

  TCanvas* c0 = new TCanvas("c0","c0",600,600);
  TCanvas* c1 = new TCanvas("c1","c1",1800,600);
  c1->Divide(3,1);

  c0->cd();
  hpos->Draw("colz");
  
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();
  c1->cd(3);
  hratio->Draw("E1");

  string histo_name = "./efficiency/eff_" + nucleus + ".root"; 
  hratio->SaveAs(histo_name.c_str());

}
