TChain* chain;

int NumberOfDetectors= 72;
int NumberOfAnodes= 11;
int nentries=1e6;

/////////////////////////////////////
void LoadRootFile(){
		chain = new TChain("PhysicsTree");
		chain->Add("/home/faster/nptool/Outputs/Analysis/run21.root");
		//chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_1.root");
		//chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_2.root");
		//chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_3.root");
		//chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_4.root");
}

/////////////////////////////////////
void FillTOFHisto(){

		LoadRootFile();
		nentries = chain->GetEntries();
		cout << "Number of entries: " << nentries << endl;

		TFile* ofile = new TFile("histo_tof_file_run21.root","recreate");
		TH1F* hLG[791];
		TH1F* hHG[791];

		vector<double>* FC_Q1 = new vector<double>();

		vector<double>* LG_Tof = new vector<double>();
		vector<int>* LG_ID = new vector<int>();
		vector<double>* LG_DT = new vector<double>();
		vector<double>* LG_Q1 = new vector<double>();
		vector<double>* LG_Q2 = new vector<double>();
		vector<int>* LG_Anode_ID = new vector<int>();
		vector<bool>* LG_FakeFission = new vector<bool>();

		vector<double>* HG_Tof = new vector<double>();
		vector<int>* HG_ID = new vector<int>();
		vector<double>* HG_DT = new vector<double>();
		vector<double>* HG_Q1 = new vector<double>();
		vector<double>* HG_Q2 = new vector<double>();
		vector<int>* HG_Anode_ID = new vector<int>();
		vector<bool>* HG_FakeFission = new vector<bool>();

		TFissionChamberPhysics* FC = new TFissionChamberPhysics();
		chain->SetBranchAddress("FissionChamber",&FC);

		chain->SetBranchAddress("FC_Q1",&FC_Q1);
		chain->SetBranchAddress("LG_Tof",&LG_Tof);
		chain->SetBranchAddress("LG_ID",&LG_ID);
		chain->SetBranchAddress("LG_DT",&LG_DT);
		chain->SetBranchAddress("LG_Q1",&LG_Q1);
		chain->SetBranchAddress("LG_Q2",&LG_Q2);
		chain->SetBranchAddress("LG_FakeFission",&LG_FakeFission);
		chain->SetBranchAddress("LG_Anode_ID",&LG_Anode_ID);
		chain->SetBranchAddress("HG_Tof",&HG_Tof);
		chain->SetBranchAddress("HG_ID",&HG_ID);
		chain->SetBranchAddress("HG_DT",&HG_DT);
		chain->SetBranchAddress("HG_Q1",&HG_Q1);
		chain->SetBranchAddress("HG_Q2",&HG_Q2);
		chain->SetBranchAddress("HG_Anode_ID",&HG_Anode_ID);
		chain->SetBranchAddress("HG_FakeFission",&HG_FakeFission);

		for(int i=0; i<NumberOfDetectors; i++){
				for(int j=0; j<NumberOfAnodes; j++){
						int index = j + i*NumberOfAnodes;
						TString histo_name = Form("hLG_Det%i_Anode%i",i+1,j+1);
						hLG[index] = new TH1F(histo_name,histo_name,2000,-100,300);

						histo_name = Form("hHG_Det%i_Anode%i",i+1,j+1);
						hHG[index] = new TH1F(histo_name,histo_name,2500,0,500);
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
						int index_LG = (LG_Anode_ID->at(j)-1) + (LG_ID->at(j)-1)*NumberOfAnodes;
						double PSD = LG_Q2->at(j)/LG_Q1->at(j);
						//if(LG_ID->at(j)>0 && LG_Anode_ID->at(j)>0 && LG_FakeFission->at(j)==0 && LG_DT->at(j)>1e7 && PSD>0.7){ 
						if(LG_ID->at(j)>0 && LG_Anode_ID->at(j)>0 && LG_FakeFission->at(j)==0){ 
								hLG[index_LG]->Fill(LG_Tof->at(j));		
						}
				}

				mysize = HG_Tof->size();
				if(mysize==1){
						for(int j=0; j<mysize; j++){
								// HG //
								int index_HG = (HG_Anode_ID->at(j)-1) + (HG_ID->at(j)-1)*NumberOfAnodes;
								double PSD = HG_Q2->at(j)/HG_Q1->at(j);
								//if(HG_ID->at(j)>0 && HG_Anode_ID->at(j)>0 && HG_FakeFission->at(j)==0 && HG_DT->at(j)>1e7 && PSD>0.7){
								if(HG_ID->at(j)>0 && HG_Anode_ID->at(j)>0 && HG_FakeFission->at(j)==0){
										hHG[index_HG]->Fill(HG_Tof->at(j));			
								}
						}
						}
				}

				ofile->Write();
				ofile->Close();

				}
