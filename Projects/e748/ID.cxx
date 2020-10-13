void ID(){

  // Load the Main Tree
  TFile* file = new TFile("../../Outputs/Analysis/e748_Physics_12Be.root");
  TTree* tree = (TTree*) file->FindObjectAny("PhysicsTree");

  TFile* fileR = new TFile("../../Outputs/Analysis/e748_12Be.root");
  TTree* treeR = (TTree*) fileR->FindObjectAny("ResultTree");
  tree->AddFriend(treeR);

  // Load the IC chain
  TChain* IC = new TChain("Numexo2");
  IC->Add("/data/Transfert/e748/Merged/run_0315.root");
  IC->Add("/data/Transfert/e748/Merged/run_0316.root");
  IC->Add("/data/Transfert/e748/Merged/run_0317.root");
  IC->Add("/data/Transfert/e748/Merged/run_0318.root");
  IC->Add("/data/Transfert/e748/Merged/run_0320.root");
  IC->Add("/data/Transfert/e748/Merged/run_0321.root");
  IC->Add("/data/Transfert/e748/Merged/run_0323.root");
  IC->Add("/data/Transfert/e748/Merged/run_0325.root");
  IC->Add("/data/Transfert/e748/Merged/run_0326.root");
  IC->Add("/data/Transfert/e748/Merged/run_0327.root");
  IC->Add("/data/Transfert/e748/Merged/run_0328.root");
  IC->Add("/data/Transfert/e748/Merged/run_0329.root");
  IC->Add("/data/Transfert/e748/Merged/run_0330.root");
  IC->Add("/data/Transfert/e748/Merged/run_0331.root");
  IC->Add("/data/Transfert/e748/Merged/run_0339.root");
  IC->Add("/data/Transfert/e748/Merged/run_0341.root");
  IC->Add("/data/Transfert/e748/Merged/run_0342.root");
  IC->Add("/data/Transfert/e748/Merged/run_0346.root");
  IC->Add("/data/Transfert/e748/Merged/run_0347.root");
  IC->Add("/data/Transfert/e748/Merged/run_0348.root");
  /*
  // Check that the number of entry is the same
  if( tree->GetEntries() != IC->GetEntries() ){
  cout << "ERROR : number of entries different \n " <<  tree->GetEntries() << " " << IC->GetEntries() << endl;
  exit(1);
  }
  */
  // Friend the two trees
  tree->AddFriend(IC);

  new TCanvas;
  tree->Draw("T_CATS1_CAV:BeamEnergy","T_CATS1_CAV>0","colz");  
  gPad->SetLogy();

  // Limit Entries to Qualifying events
  TFile* elf= new TFile("EntryList.root");
  TEntryList* elist = (TEntryList*) elf->FindObjectAny("elist");
  if(!elist){
    tree->Draw(">>elist", "PositionOnTargetY > -10 && PositionOnTargetY < 10 && PositionOnTargetX>-20 && PositionOnTargetX< 20 &&TelescopeNumber<5 && CsI_E<0 && IC_E > 0 && QPlast>0 && T_CATS1_CAV>160 && T_CATS1_CAV<235 ", "entrylist");
    elist = (TEntryList*) gDirectory->Get("elist");
    elist->SaveAs("EntryList.root");
  }
  tree->SetEntryList(elist);
  // Load Cut
  TFile* cutf = new TFile("Cut.root");
  TCutG* cutg = (TCutG*) cutf->FindObjectAny("Cut_IC");
  // IC for MUST2 qualifying events
  new TCanvas();
  tree->Draw("IC_E:QPlast>>hIC(4000,0,16000,900,0,900)","","colz");
   cutg->Draw("same");

  // TOF per telescope

  new TCanvas();
  TFile* cutTOF = new TFile("TOF.root");
  TCutG* cutt = (TCutG*) cutTOF->FindObjectAny("TOF");

  tree->Draw("Si_T+TimeCorr+((TelescopeNumber==1)*-2.521+(TelescopeNumber==2)*0.148+(TelescopeNumber==3)*-1.922+(TelescopeNumber==4)*-7.176):Si_E>>hTOFc(1000,0,30,1000,460,580)","","colz");
  cutt->Draw("same");

  /*new TCanvas();
    for(unsigned int i = 1 ; i < 5 ; i++){
  //  new TCanvas();
  // TFile* cutTOF = new TFile(Form("TOF%i.root",i));
  //TCutG* cutt = (TCutG*) cutTOF->FindObjectAny(Form("TOF%i",i));
  tree->Draw(Form("Si_T+TimeCorr+((TelescopeNumber==1)*-2.521+(TelescopeNumber==2)*0.148+(TelescopeNumber==3)*-1.922+(TelescopeNumber==4)*-7.176):Si_E>>hTOFc%i(1000,0,30,1000,460,580)",i),Form(" TelescopeNumber==%i",i),"colz");
  else
  }

*/



  new TCanvas();
  TString cond = "Cut_IC && TOF";

  tree->Draw("ELab:ThetaLab>>hKin(1000,0,60,1000,0,30)",cond,""); 
  TH1* hKin=(TH1*) gDirectory->FindObjectAny("hKin");
  hKin->SetMarkerColor(kBlack);hKin->SetMarkerStyle(21);
  hKin->SetMarkerSize(0.5);
  NPL::Reaction r("12Be(d,3He)11Li@360");
  r.GetKinematicLine3()->Draw("c");

  new TCanvas();
  tree->Draw("Ex>>hEx(100,-10,10)",cond,""); 
  TH1* hEx=(TH1*) gDirectory->FindObjectAny("hEx");
  double max = 22;
  hEx->GetYaxis()->SetRangeUser(0,max);
  hEx->GetXaxis()->SetRangeUser(-2,10);

  TLine* line1 = new TLine(0,max,0,0);
  line1->SetLineWidth(2); //line1->SetLineColor(kOrange+3);
  TLine* line2 = new TLine(2.474,max,2.474,0);
  line2->SetLineWidth(2); //line2->SetLineColor(kOrange+3);
  TLine* line3 = new TLine(4.86,max,4.86,0);
  line3->SetLineWidth(2); //line3->SetLineColor(kOrange+3);
  TLine* linen = new TLine(0.396,max,0.396,0);
  linen->SetLineWidth(2);linen->SetLineStyle(2); linen->SetLineColor(kBlack);

  linen->Draw("same");
  line1->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
}
