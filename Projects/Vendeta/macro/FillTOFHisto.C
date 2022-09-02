TChain* chain;

int NumberOfDetectors= 72;
int NumberOfAnodes= 1;
int nentries=1e6;

/////////////////////////////////////
void LoadRootFile(){
	chain = new TChain("PhysicsTree");
	chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_1.root");
	//chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_2.root");
	//chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_3.root");
	//chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_4.root");
}

/////////////////////////////////////
void FillTOFHisto(){

	LoadRootFile();
	nentries = chain->GetEntries();
	cout << "Number of entries: " << nentries << endl;

	TFile* ofile = new TFile("histo_tof_file.root","recreate");
	TH1F* hLG[792];
	TH1F* hHG[792];
	
	vector<double>* FC_Q1 = new vector<double>();
 	
	vector<double>* LG_Tof = new vector<double>();
	vector<int>* LG_ID = new vector<int>();
	vector<int>* LG_Anode_ID = new vector<int>();
	
	vector<double>* HG_Tof = new vector<double>();
	vector<int>* HG_ID = new vector<int>();
	vector<int>* HG_Anode_ID = new vector<int>();

	TFissionChamberPhysics* FC = new TFissionChamberPhysics();
	chain->SetBranchAddress("FissionChamber",&FC);

	chain->SetBranchAddress("FC_Q1",&FC_Q1);
	chain->SetBranchAddress("LG_Tof",&LG_Tof);
	chain->SetBranchAddress("LG_ID",&LG_ID);
	chain->SetBranchAddress("LG_Anode_ID",&LG_Anode_ID);
	chain->SetBranchAddress("HG_Tof",&HG_Tof);
	chain->SetBranchAddress("HG_ID",&HG_ID);
	chain->SetBranchAddress("HG_Anode_ID",&HG_Anode_ID);
	
for(int i=0; i<NumberOfDetectors; i++){
			for(int j=0; j<NumberOfAnodes; j++){
					int anode = 6;
					//int index = (i+1) * (j+1);
					int index = (i+1) * anode;
					//TString histo_name = Form("hLG_Det%i_Anode%i",i+1,j+1);
					TString histo_name = Form("hLG_Det%i_Anode%i",i+1,anode);
					hLG[index-1] = new TH1F(histo_name,histo_name,4000,-100,300);
					//histo_name = Form("hHG_Det%i_Anode%i",i+1,j+1);
					histo_name = Form("hHG_Det%i_Anode%i",i+1,anode);
					hHG[index-1] = new TH1F(histo_name,histo_name,3000,0,300);
			}
	}
	
	for(int i=0; i<nentries; i++){
			chain->GetEntry(i);

			if(i%100000==0){
				cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
			}
			//int FC_mult = FC->AnodeNumber.size();
		 //cout << FC_mult << endl;
			int mysize = LG_Tof->size();
			for(int j=0; j<mysize; j++){
					// LG //
					int index_LG = LG_ID->at(j) * LG_Anode_ID->at(j);
					if(LG_ID->at(j)>0 && LG_Anode_ID->at(j)>0){ 
							hLG[index_LG-1]->Fill(LG_Tof->at(j));		
					}
			}
			
			mysize = HG_Tof->size();
			for(int j=0; j<mysize; j++){
					// HG //
					int index_HG = HG_ID->at(j) * HG_Anode_ID->at(j);
					if(HG_ID->at(j)>0 && HG_Anode_ID->at(j)>0){
							hHG[index_HG-1]->Fill(HG_Tof->at(j));			
					}
			}
	}

	/* for(int i=0; i<NumberOfDetectors; i++){ */
	/* 		for(int j=0; j<NumberOfAnodes; j++){ */
	/* 				int anode = 6; */
	/* 				//int index = (i+1) * (j+1); */
	/* 				int index = (i+1) * anode; */
						
	/* 				hLG[index]->Write(); */
	/* 				hHG[index]->Write(); */
	/* 		} */
	/* } */

	//hLG[6]->Draw();
	ofile->Write();
	ofile->Close();

}
