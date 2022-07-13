void Colorbyprocess(){
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
  TH1F *Inel = new TH1F("Inel","test",3000,10500,12500);
  TH1F *El = new TH1F("El","test",3000,10500,12500);
  Inel->SetStats(0);
  t->SetBranchAddress("Pos_Z",&pos_z_vec);
  t->SetBranchAddress("Process",&processname,&branch);
  int nentries = branch->GetEntries();
    
  for(int i = 0; i < nentries ; i++){
    
    t->GetEntry(i);
    int vec_size = (*processname).size();
    
    for(int j = 0; j < vec_size; j++){
      name = ((*processname)[j]).c_str();
      pos_z = (*pos_z_vec)[j];
      if(strcmp(name,"hadElastic")==0){
        El->Fill(pos_z);
      }
      else if(strcmp(name,"neutronInelastic")==0){
        Inel->Fill(pos_z);
      }
    }
  }


  TCanvas* c = new TCanvas();

  Inel->SetLineColor(kRed);
  El->SetLineColor(kBlue);

  El->Draw();
  Inel->Draw("SAME");
}
