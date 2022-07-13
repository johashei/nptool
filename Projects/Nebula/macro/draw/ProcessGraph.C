void DrawProcessGraph(){
  TFile* f = TFile::Open("root/simulation/SimulatedTree.root");
  TTree *t = (TTree*)f->Get("SimulatedTree");
  //new TCanvas;
  //t->Draw("Process");
  //new TCanvas;
  //t->Draw("El_is_1_and_InEl_is_2");
  //new TCanvas;
  //t->Draw("Efficiency");


  std::vector<std::string>* processname = 0;
  const char* name;
  std::vector<double>* pos_z_vec = 0;
  double pos_z;

  TBranch* branch;
  TCanvas* c = new TCanvas();
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.13);
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.08);
  TH2F *h = new TH2F("h","test",500,9900,10100,2,0,2);
  h->SetCanExtend(TH1::kAllAxes);
  h->SetStats(0);
  t->SetBranchAddress("Pos_Z",&pos_z_vec);
  t->SetBranchAddress("Process",&processname,&branch);
  int nentries = branch->GetEntries();
    
  for(int i = 0; i < nentries ; i++){
    
    t->GetEntry(i);
    int vec_size = (*processname).size();
    
    for(int j = 0; j < vec_size; j++){
      name = ((*processname)[j]).c_str();
      pos_z = (*pos_z_vec)[j];
      if(strcmp(name,"Transportation")!=0){
        h->Fill(pos_z,name,1);
      }
    }
  }

  //h->SetMaximum(5000);
  h->SetMinimum(-1);
  h->LabelsDeflate("X");
  h->LabelsDeflate("Y");
  h->LabelsOption("v");
  //h->Draw("LEGO2");
  h->Draw("SCAT");
}
