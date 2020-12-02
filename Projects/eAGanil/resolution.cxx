void AddNuclei(string,double);
void resolution(){

  auto f = new TFile("../../Outputs/Analysis/PhysicsTree.root");
  auto tree= (TTree*) f->FindObjectAny("PhysicsTree");

  vector<double> Resolution={5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6};
  unsigned int sizeR = Resolution.size();
  auto c = new TCanvas();
  c->Divide(3,2);
  vector<double> R ;
  for(unsigned int i = 0 ; i < sizeR ; i++){
    c->cd(i+1);
    if(i<3)
      tree->Draw(Form("Ex>>h%d(1000,-10,10)",i),Form("Resolution==%f",Resolution[i]),"");
    else
      tree->Draw(Form("Ex>>h%d(1000,-0.5,0.5)",i),Form("Resolution==%f",Resolution[i]),"");

    auto h = (TH1*) gDirectory->FindObjectAny(Form("h%d",i)); 
    h->Fit("gaus");
    R.push_back(2*sqrt(2*log(2))*h->GetFunction("gaus")->GetParameter(2)*1000.);
    cout << R[i] << endl;
  }
   
  auto g = new TGraph(R.size(),&Resolution[0],&R[0]);
  new TCanvas();
  g->Draw("apc");
  gPad->SetLogx();

  g->GetYaxis()->SetRangeUser(-100,1400);
  g->GetXaxis()->SetTitle("E/#deltaE");
  g->GetYaxis()->SetTitle("FWHM (keV)");


  AddNuclei("134Sn",725.6);
  AddNuclei("132Sn",64.4);
  AddNuclei("133Sn",853.7);
  AddNuclei("70Ni",183.11);
  AddNuclei("71Ni",252.2);
  AddNuclei("72Ni",455);
  //AddNuclei("73Ni",239.2);
  AddNuclei("74Ni",739);
}

void AddNuclei(string name,double keV){
  static bool alt = true;
  static double off =1e-6;
  auto L = new TLine(1e-6,keV,5e-2,keV);
  L->SetLineStyle(2);
  L->SetLineColor(kAzure+7);
  L->Draw("same");
  TLatex* t;
  if(alt){
    t = new TLatex(3e-6-off,keV-30,name.c_str());
    alt=false;
  }
  else{
    t = new TLatex(3e-6+off,keV+10,name.c_str());
    alt=true;
  }
   t->SetTextSize(0.02);
   t->Draw(); 
}
