void ProcVsKinEn2D(bool DoPads = 1, const char* PartPara = "all"){
  
  TFile *f = TFile::Open("root/simulation/SimulatedTree.root");
  TTree *t = (TTree*)f->Get("SimulatedTree");

  std::vector<double>* v_kinen = 0;
  std::vector<std::string>* v_partname = 0;
  std::vector<std::string>* v_procname = 0;
  const char* procname;
  const char* partname;
  double kinEn, prevkinEn;
  double kinDiff;

  TBranch* branch_kin;
  TBranch* branch_proc;
  TBranch* branch_part;

  TCanvas* c = new TCanvas("c", "Processes by Kinetic Energy", 1920, 1200);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.1);
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.1);

  int NbBin = 100;
  int Xmax = 200;
  int Xmin = 0;
  int Ymin = 1;
  int Ymax = 1e7; 

  THStack *hs = new THStack("hs","Process CS according to Kinetic Energy");

  TH1F *h_neutIn = new TH1F("h_neutIn","Neutron Inelastic",NbBin,Xmin,Xmax);
  h_neutIn->GetYaxis()->SetRangeUser(Ymin, Ymax);
  TH1F *h_hadEl = new TH1F("h_hadEl","Hadron Elastic",NbBin,Xmin,Xmax);
  h_hadEl->GetYaxis()->SetRangeUser(Ymin, Ymax);
  TH1F *h_ionIoni = new TH1F("h_ionIoni","Ion Ionization",NbBin,Xmin,Xmax);
  h_ionIoni->GetYaxis()->SetRangeUser(Ymin, Ymax);
  TH1F *h_hIoni = new TH1F("h_hIoni","Hadron Ionization",NbBin,Xmin,Xmax);
  h_hIoni->GetYaxis()->SetRangeUser(Ymin, Ymax);

  t->SetBranchAddress("KinEn", &v_kinen, &branch_kin);
  t->SetBranchAddress("Process", &v_procname, &branch_proc);
  t->SetBranchAddress("Particle", &v_partname, &branch_part);
  int nentries = branch_kin->GetEntries();

  for(int i = 0; i < nentries ; i++){
    kinDiff = 0;
    t->GetEntry(i);
    int vec_size = (*v_kinen).size();

    for(int j = 1; j < vec_size; j++){
      procname = ((*v_procname)[j]).c_str();;
      partname = ((*v_partname)[j]).c_str();;
      prevkinEn = (*v_kinen)[j-1];
      kinDiff = prevkinEn - kinEn;
      if((strcmp(PartPara, "all")==0 || strcmp(partname, PartPara)==0) && 
         strcmp(procname, "Transportation")!=0 &&
         kinDiff>0){
           if(strcmp(procname, "neutronInelastic")==0) h_neutIn->Fill(kinDiff);
           else if(strcmp(procname, "hadElastic")==0) h_hadEl->Fill(kinDiff);
           else if(strcmp(procname, "ionIoni")==0) h_ionIoni->Fill(kinDiff);
           else if(strcmp(procname, "hIoni")==0) h_hIoni->Fill(kinDiff);
      }
    }
  }
  h_neutIn->GetYaxis()->SetTitle("Nb of Hits");
  h_neutIn->GetXaxis()->SetTitle("Kinetic Energy (MeV)");

  h_hadEl->GetYaxis()->SetTitle("Nb of Hits");
  h_hadEl->GetXaxis()->SetTitle("Kinetic Energy (MeV)");

  h_hIoni->GetYaxis()->SetTitle("Nb of Hits");
  h_hIoni->GetXaxis()->SetTitle("Kinetic Energy (MeV)");

  h_ionIoni->GetYaxis()->SetTitle("Nb of Hits");
  h_ionIoni->GetXaxis()->SetTitle("Kinetic Energy (MeV)");

  hs->Add(h_neutIn);
  hs->Add(h_hadEl);
  hs->Add(h_ionIoni);
  hs->Add(h_hIoni);

  if(DoPads){
    h_neutIn->SetFillColor(kBlue-7);
    h_hadEl->SetFillColor(kRed-7);
    h_hIoni->SetFillColor(kYellow-7);
    h_ionIoni->SetFillColor(kGreen-7);

    hs->Draw("pads");

    gPad->Update();
    c->cd(1);
      gPad->SetLogy(1);

      gStyle->SetHistTopMargin(0.);
      gStyle->SetOptStat(0);

    c->cd(2);
      gPad->SetLogy(1);

      gStyle->SetHistTopMargin(0.);
      gStyle->SetOptStat(0);

    c->cd(3);
      gPad->SetLogy(1);

      gStyle->SetHistTopMargin(0.);
      gStyle->SetOptStat(0);

    c->cd(4);
      gPad->SetLogy(1);

      gStyle->SetHistTopMargin(0.);
      gStyle->SetOptStat(0);
  }
  else{
    gStyle->SetPalette(60);
    hs->Draw("plc");
    gPad->SetLogy(1);
    hs->SetMinimum(1);
    hs->SetMaximum(1e6);
    hs->GetYaxis()->SetLimits(1,1e6);
    //hs->Draw("plc");
  }
}
