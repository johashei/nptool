void MT(){
using namespace ROOT::Experimental; // TDataFrame's namespace
//ROOT::EnableImplicitMT();
  // Load the Main Tree
  TChain chain("PhysicsTree");
  chain.Add("../../Outputs/Analysis/e748_12Be.root");

  // Load the IC chain
  TChain IC("Numexo2");
  IC.Add("/data/Transfert/e748/Merged/run_0315.root");
  IC.Add("/data/Transfert/e748/Merged/run_0316.root");
  IC.Add("/data/Transfert/e748/Merged/run_0317.root");
  IC.Add("/data/Transfert/e748/Merged/run_0318.root");
  IC.Add("/data/Transfert/e748/Merged/run_0320.root");
  IC.Add("/data/Transfert/e748/Merged/run_0321.root");
  IC.Add("/data/Transfert/e748/Merged/run_0323.root");
  IC.Add("/data/Transfert/e748/Merged/run_0325.root");
  IC.Add("/data/Transfert/e748/Merged/run_0326.root");
  IC.Add("/data/Transfert/e748/Merged/run_0327.root");
  IC.Add("/data/Transfert/e748/Merged/run_0328.root");
  IC.Add("/data/Transfert/e748/Merged/run_0329.root");
  IC.Add("/data/Transfert/e748/Merged/run_0330.root");
  IC.Add("/data/Transfert/e748/Merged/run_0331.root");
  IC.Add("/data/Transfert/e748/Merged/run_0339.root");
  IC.Add("/data/Transfert/e748/Merged/run_0341.root");
  IC.Add("/data/Transfert/e748/Merged/run_0342.root");
  IC.Add("/data/Transfert/e748/Merged/run_0346.root");
  IC.Add("/data/Transfert/e748/Merged/run_0347.root");
  IC.Add("/data/Transfert/e748/Merged/run_0348.root");
/*
  // Check that the number of entry is the same
  if( tree->GetEntries() != IC->GetEntries() ){
    cout << "ERROR : number of entries different \n " <<  tree->GetEntries() << " " << IC->GetEntries() << endl;
    exit(1);
    }
*/
  // Friend the two trees
  chain.AddFriend(&IC);

  // Limit Entries to Qualifying events
  TFile* elf= new TFile("EntryList.root");
  TEntryList* elist = (TEntryList*) elf->FindObjectAny("elist");
  if(!elist){
    chain.Draw(">>elist", "PositionOnTargetY > -10 && PositionOnTargetY < 10 && PositionOnTargetX>-20 && PositionOnTargetX< 20 &&TelescopeNumber<5 && IC_E > 0 && QPlast>0 && T_CATS1_CAV>160 && T_CATS1_CAV<235 ", "entrylist");
    elist = (TEntryList*) gDirectory->Get("elist");
    elist->SaveAs("EntryList.root");
  }
  chain.SetEntryList(elist);

  TDataFrame DF(chain);


  // TOF per telescope
      auto h1 = DF.Filter("TelescopeNumber==1").Histo1D("Si_T+TimeCorr:Si_E>>hTOFc1(1000,0,30,1000,460,580)");
  

      new TCanvas();
      h1->Draw("colz");
      

}
