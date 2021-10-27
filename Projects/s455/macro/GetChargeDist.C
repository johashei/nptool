TChain* chain=NULL;

TCutG* cut_Pb1[13];
TCutG* cut_Pb2[13];
TCutG* cut_C[13];

int A;
int Z;
NPL::Particle* iso;

TH2F* hbid;
TH1F* hsumpb;
TH1F* hsumc;
TH1F* hzpb;
TH1F* hzc;
TH1F* hzcscaled;
TH1F* hsumcscaled;
TH1F* hzsub;

void LoadRootFiles(string iso);


//////////////////////////////////////
void LoadCuts(){

  TString element[13] = {"180Hg_1", "180Hg_2", "182Hg", "184Hg", "187Pb", "189Pb", "175Pt", "204Fr", "207Fr", "199At_1", "199At_2", "197At", "216Th"};


  TFile* file;
  TString filename1[13];
  TString filename2[13];
  TString filename3[13];
  for(int i=0; i<13; i++){
    filename1[i]= "cuts/"+element[i]+"/cut_Pb1.root";
    filename2[i]= "cuts/"+element[i]+"/cut_Pb2.root";
    filename3[i]= "cuts/"+element[i]+"/cut_C.root";

    file = new TFile(filename1[i],"read");
    cut_Pb1[i] = (TCutG*) file->FindObjectAny("cut_Pb1");
    cut_Pb1[i]->SetName(element[i]+"_Pb1");

    file = new TFile(filename2[i],"read");
    cut_Pb2[i] = (TCutG*) file->FindObjectAny("cut_Pb2");
    cut_Pb2[i]->SetName(element[i]+"_Pb2");
    
    file = new TFile(filename3[i],"read");
    cut_C[i] = (TCutG*) file->FindObjectAny("cut_C");
    cut_C[i]->SetName(element[i]+"_C");
  }
}

//////////////////////////////////////
void GetChargeDist(string nucleus_name){

  iso = new NPL::Particle(nucleus_name);
  A = iso->GetA();
  Z = iso->GetZ();
  double AoQ=(double)A/(double)Z;

  cout << "*******************************" << endl;
  cout << "Nucleus name = " << nucleus_name << endl;
  cout << "A= " << A << endl;
  cout << "Z= " << Z << endl;
  cout << "AoQ= " << AoQ << endl;

  LoadRootFiles(nucleus_name);
  LoadCuts();


  //*** opening output root file ***//
  TString output_filename = "root/charge_dist_" + nucleus_name + ".root";
  TFile* ofile = new TFile(output_filename, "recreate");
  cout << "Creating file: " << output_filename << endl;
  ofile->cd();

  hbid        = new TH2F("hbid","hbid",1000,2.2,2.5,1000,70,95);
  hsumpb      = new TH1F("hsumpb","hsumpb",1000,50,95);
  hsumc       = new TH1F("hsumc","hsumc",1000,50,95);
  hzpb        = new TH1F("hzpb","hzpb",1000,20,60);
  hzc         = new TH1F("hzc","hzc",1000,20,60);
  hzcscaled   = new TH1F("hzcscaled","hzcscaled",1000,20,60);
  hsumcscaled = new TH1F("hsumcscaled","hsumcscaled",1000,50,95);
  hzsub       = new TH1F("hzsub","hzsub",1000,20,60);


  //*** Defining beam cut ***//
  TCutG* cut_beam = new TCutG("cut_beam");
  cut_beam->SetName("cut_beam");
  cut_beam->SetVarX("fBeam_AoQ");
  cut_beam->SetVarY("fBeam_Z");
  double Rx=0.006;
  double Ry=0.5;
  AoQ = AoQ;// + 0.001;
  double Zcenter = Z;// - 0.5;
  for(int i=0; i<360; i++){
    double x = AoQ+Rx*cos(i*TMath::Pi()/180);
    double y = Zcenter+Ry*sin(i*TMath::Pi()/180);
    cut_beam->SetPoint(i,x,y);
  }
  cut_beam->SetLineWidth(2);
  cut_beam->SetLineColor(2);


   
  int nentries = chain->GetEntries();

  //*** Init input branch ***//
  int RunID; 
  TSofBeamID* SofBeam        = new TSofBeamID();
  TSofAtPhysics* SofAt       = new TSofAtPhysics();
  TSofFissionFragment* SofFF = new TSofFissionFragment();
    
  chain->SetBranchAddress("RunID",&RunID);
  chain->SetBranchAddress("SofBeamID",&SofBeam);
  chain->SetBranchAddress("SofAt",&SofAt);
  chain->SetBranchAddress("SofFissionFragment",&SofFF);

  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);
   
    if(i%100000==0){
      cout << "\033[34m\r Processing tree... " << (double)i/nentries*100 << "\% done" << flush;
    }

    if(SofAt->Energy.size()==3 || SofAt->Energy.size()==4){    
       
      double Zbeam = SofBeam->GetZbeam();
      double AoQ   = SofBeam->GetAoQ();
        
      double Anode1 = SofAt->Energy[0];
      double Anode2 = SofAt->Energy[1];
      double Anode3 = SofAt->Energy[2];
      
      int k=RunID-1;
      if(cut_Pb1[k]->IsInside(Anode1,Anode2) || cut_Pb2[k]->IsInside(Anode2,Anode3)){
        hbid->Fill(AoQ,Zbeam);
      }

      if(cut_beam->IsInside(AoQ,Zbeam)){
        //*** Fission Fragment information ***//
        double Zsum = SofFF->GetZsum();
        int iZsum   = SofFF->GetIntZsum();
        double Z1 = -1;
        double Z2 = -1;
        if(SofFF->GetMult()>0){
          Z1   = SofFF->GetZ(0);
          Z2   = SofFF->GetZ(1);
        }

        //*** Taking event from fission in Lead cathod ***//
        if(cut_Pb1[k]->IsInside(Anode1,Anode2) || cut_Pb2[k]->IsInside(Anode2,Anode3)){
          hsumpb->Fill(Zsum);
          if(iZsum==Z){
            hzpb->Fill(Z1);
            hzpb->Fill(Z2);
          
            hzsub->Fill(Z1);
            hzsub->Fill(Z2);
          }
        }

        //*** Taking event from fission in Carbon cathod ***//
        if(cut_C[k]->IsInside(Anode2,Anode3)){
          hsumc->Fill(Zsum);
          hsumcscaled->Fill(Zsum);
          if(iZsum==Z){
            hzc->Fill(Z1);
            hzc->Fill(Z2);

            hzcscaled->Fill(Z1);
            hzcscaled->Fill(Z2);
          }
        }
          }
    }
  } // end of loop ovent nentries

  double integral_pb = hsumpb->Integral(hsumpb->GetXaxis()->FindBin(65.5), hsumpb->GetXaxis()->FindBin(Z-1.5));
  double integral_c  = hsumc->Integral(hsumc->GetXaxis()->FindBin(65.5), hsumc->GetXaxis()->FindBin(Z-1.5));
  double ratio = integral_pb/integral_c;
  if(ratio>0){
    cout << "*** ratio= " << ratio << endl; 
    hsumcscaled->Scale(ratio);
    hzcscaled->Scale(ratio);

    hzsub->Add(hzcscaled,-1);
  }
  cout << endl;
  hbid->Write();
  hsumpb->Write();
  hsumc->Write();
  hzpb->Write();
  hzc->Write();
  hzcscaled->Write();
  hsumcscaled->Write();
  hzsub->Write();
  cut_beam->Write();
}


//////////////////////////////////////
void LoadRootFiles(string iso){
  chain = new TChain("PhysicsTree");

  //*** Th isotopes ***//
  if(iso=="217Th" || iso=="216Th" || iso=="209Ra" || iso=="210Ra" || iso=="211Ra" || iso=="212Ra" || iso=="213Ra" || iso=="214Ra"){
    chain->Add("../../../Outputs/Analysis/216Th_1.root");
  }

  //*** Fr isotopes ***//
  if(iso=="203Fr" || iso=="204Fr" || iso=="205Fr" || iso=="206Fr"){
    chain->Add("../../../Outputs/Analysis/204Fr.root");
  }
  if(iso=="205Fr" || iso=="206Fr" || iso=="207Fr" || iso=="208Fr"){
    chain->Add("../../../Outputs/Analysis/207Fr_2.root");
  }
  if(iso=="208Fr" || iso=="209Fr" || iso=="210Fr" || iso=="211Fr"){
    chain->Add("../../../Outputs/Analysis/216Th_1.root");
  }

  //*** Rn isotopes ***//
  if(iso=="200Rn" || iso=="201Rn" || iso=="202Rn" || iso=="203Rn" || iso=="204Rn" || iso=="205Rn"){
    chain->Add("../../../Outputs/Analysis/204Fr.root");
  }
  if(iso=="203Rn" || iso=="204Rn" || iso=="205Rn" || iso=="206Rn"){
    chain->Add("../../../Outputs/Analysis/207Fr.root");
  }

  //*** At isotopes ***//
  if(iso=="197At" || iso=="198At" || iso=="199At" || iso=="200At" || iso=="201At" || iso=="202At"){
    chain->Add("../../../Outputs/Analysis/204Fr.root");
  }
  if(iso=="200At" || iso=="201At" || iso=="202At" || iso=="203At"){
    chain->Add("../../../Outputs/Analysis/207Fr.root");
  }
  if(iso=="198At" || iso=="199At" || iso=="200At" || iso=="201At"){
    chain->Add("../../../Outputs/Analysis/199At_1.root");
    chain->Add("../../../Outputs/Analysis/199At_2.root");
  }
  if(iso=="196At" || iso=="197At" || iso=="198At" || iso=="199At"){
    chain->Add("../../../Outputs/Analysis/197At.root");
  }

  //*** Po isotopes ***//
  if(iso=="192Po" || iso=="193Po" || iso=="194Po" || iso=="195Po" || iso=="196Po" || iso=="197Po"){
    chain->Add("../../../Outputs/Analysis/197At.root");
  }
  if(iso=="195Po" || iso=="196Po" || iso=="197Po" || iso=="198Po" || iso=="199Po" || iso=="200Po"){
    chain->Add("../../../Outputs/Analysis/199At_1.root");
    chain->Add("../../../Outputs/Analysis/199At_2.root");
  }

  //*** Bi isotopes ***//
  if(iso=="189Bi" || iso=="190Bi" || iso=="191Bi" || iso=="192Bi" || iso=="193Bi" || iso=="194Bi" || iso=="195Bi"){
    chain->Add("../../../Outputs/Analysis/197At.root");
  }
  if(iso=="192Bi" || iso=="193Bi" || iso=="194Bi" || iso=="195Bi" || iso=="196Bi" || iso=="197Bi"){
    chain->Add("../../../Outputs/Analysis/199At_1.root");
    chain->Add("../../../Outputs/Analysis/199At_2.root");
  }

  //*** Pb isotopes ***//
  if(iso=="186Pb" || iso=="187Pb" || iso=="188Pb" || iso=="189Pb" || iso=="190Pb" || iso=="191Pb" || iso=="192Pb"){
    chain->Add("../../../Outputs/Analysis/197At.root");
  }
  if(iso=="191Pb" || iso=="192Pb" || iso=="193Pb" || iso=="194Pb"){
    chain->Add("../../../Outputs/Analysis/199At_1.root");
    chain->Add("../../../Outputs/Analysis/199At_2.root");
  }
  if(iso=="188Pb" || iso=="189Pb" || iso=="190Pb" || iso=="191Pb"){
    chain->Add("../../../Outputs/Analysis/189Pb.root");
  }
  if(iso=="186Pb" || iso=="187Pb" || iso=="188Pb" || iso=="189Pb"){
    chain->Add("../../../Outputs/Analysis/187Pb.root");
  }

  //*** Tl isotopes ***//
  if(iso=="186Tl" || iso=="187Tl" || iso=="188Tl" || iso=="189Tl" || iso=="190Tl"){
    chain->Add("../../../Outputs/Analysis/197At.root");
  }
  if(iso=="185Tl" || iso=="186Tl" || iso=="187Tl" || iso=="188Tl" || iso=="189Tl"){
    chain->Add("../../../Outputs/Analysis/189Pb.root");
  }
  if(iso=="183Tl" || iso=="184Tl" || iso=="185Tl" || iso=="186Tl" || iso=="187Tl"){
    chain->Add("../../../Outputs/Analysis/187Pb.root");
  }
  if(iso=="185Tl" || iso=="186Tl" || iso=="187Tl"){
    chain->Add("../../../Outputs/Analysis/184Hg.root");
  }
  if(iso=="183Tl" || iso=="184Tl" || iso=="185Tl"){
    chain->Add("../../../Outputs/Analysis/182Hg.root");
  }

  //*** Hg isotopes ***//
  if(iso=="182Hg" || iso=="183Hg" || iso=="184Hg" || iso=="185Hg" || iso=="186Hg" || iso=="187Hg"){
    chain->Add("../../../Outputs/Analysis/189Pb.root");
  }
  if(iso=="180Hg" || iso=="181Hg" || iso=="182Hg" || iso=="183Hg" || iso=="184Hg" || iso=="185Hg"){
    chain->Add("../../../Outputs/Analysis/187Pb.root");
  }
  if(iso=="183Hg" || iso=="184Hg" || iso=="185Hg" || iso=="186Hg"){
    chain->Add("../../../Outputs/Analysis/184Hg.root");
  }
  if(iso=="181Hg" || iso=="182Hg" || iso=="183Hg" || iso=="184Hg"){
    chain->Add("../../../Outputs/Analysis/182Hg.root");
  }
  if(iso=="178Hg" || iso=="179Hg" || iso=="180Hg" || iso=="181Hg"){
    chain->Add("../../../Outputs/Analysis/180Hg_1.root");
    chain->Add("../../../Outputs/Analysis/180Hg_2.root");
    chain->Add("../../../Outputs/Analysis/175Pt.root");
  }

  //*** Au isotopes ***//
  if(iso=="177Au" || iso=="178Au" || iso=="179Au"){
    chain->Add("../../../Outputs/Analysis/175Pt.root");
  }
  if(iso=="180Au" || iso=="181Au" || iso=="182Au"){
    chain->Add("../../../Outputs/Analysis/189Pb.root");
  }
  if(iso=="179Au" || iso=="180Au" || iso=="181Au" || iso=="182Au"){
    chain->Add("../../../Outputs/Analysis/187Pb.root");
  }
  if(iso=="180Au" || iso=="181Au" || iso=="182Au"){
    chain->Add("../../../Outputs/Analysis/184Hg.root");
  }
  if(iso=="178Au" || iso=="179Au" || iso=="180Au" || iso=="181Au" || iso=="182Au"){
    chain->Add("../../../Outputs/Analysis/182Hg.root");
  }
  if(iso=="177Au" || iso=="178Au" || iso=="179Au" || iso=="180Au"){
    chain->Add("../../../Outputs/Analysis/180Hg_1.root");
    chain->Add("../../../Outputs/Analysis/180Hg_2.root");
  }

  //*** Pt isotopes ***//
  if(iso=="173Pt" || iso=="174Pt" || iso=="175Pt" || iso=="176Pt" || iso=="177Pt" || iso=="178Pt"){
    chain->Add("../../../Outputs/Analysis/175Pt.root");
  }
  if(iso=="174Pt" || iso=="175Pt" || iso=="176Pt" || iso=="177Pt" || iso=="178Pt"){
    chain->Add("../../../Outputs/Analysis/180Hg_1.root");
    chain->Add("../../../Outputs/Analysis/180Hg_2.root");
  }

  //*** Ir isotopes ***//
  if(iso=="172Ir" || iso=="173Ir" || iso=="174Ir" || iso=="175Ir" || iso=="176It"){
    chain->Add("../../../Outputs/Analysis/175Pt.root");
    chain->Add("../../../Outputs/Analysis/180Hg_1.root");
    chain->Add("../../../Outputs/Analysis/180Hg_2.root");
  }

  cout << "Number of entries: " << chain->GetEntries() << endl;
}


