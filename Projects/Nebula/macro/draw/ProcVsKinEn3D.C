void ProcVsKinEn3D(){

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

  TH2F *h_process_kin = new TH2F("h_process_kin","Processes by Kinetic Energy",80,0,200,1,0,1);
  //h_process_kin->SetCanExtend(TH1::kAllAxes);
  h_process_kin->SetStats(0);

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
      if(/*strcmp(partname, "neutron")==0 && */
        strcmp(procname, "Transportation")!=0 &&
        kinDiff>0){
        h_process_kin->Fill(kinDiff,procname,1);
      }
    }
  }
  //h_process_kin->SetMinimum(-1);
  //h_process_kin->SetMaximum();
  h_process_kin->GetYaxis()->SetTitle("Process Name");
  h_process_kin->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  //gStyle->SetPalette(60);
  //gStyle->SetPalette(53);
  c->SetFrameFillColor(kAzure-6);
  h_process_kin->Draw("COLZ");
  c->SetLogz();
}
