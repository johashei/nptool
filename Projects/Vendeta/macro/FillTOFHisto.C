TChain* chain;

int NumberOfDetectors= 72;
int NumberOfAnodes= 11;
int nentries=1e6;
int run=59;
/////////////////////////////////////
void LoadRootFile(){
		chain = new TChain("PhysicsTree");
		chain->Add(Form("/home/faster/nptool/Outputs/Analysis/run%i.root",run));
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

		TFile* ofile = new TFile(Form("histo_tof_file_run%i.root",run),"recreate");
		TH1F* hLG[791];
		TH1F* hHG[791];

		vector<double>* FC_Q1 = new vector<double>();
		vector<int>* FC_Anode_ID = new vector<int>();
		vector<bool>* FC_FakeFission = new vector<bool>();

		vector<double>* LG_Tof = new vector<double>();
		vector<int>* LG_ID = new vector<int>();
		vector<double>* LG_Q1 = new vector<double>();
		vector<double>* LG_Q2 = new vector<double>();

		vector<double>* HG_Tof = new vector<double>();
		vector<int>* HG_ID = new vector<int>();
		vector<double>* HG_Q1 = new vector<double>();
		vector<double>* HG_Q2 = new vector<double>();

		TFissionChamberPhysics* FC = new TFissionChamberPhysics();
		chain->SetBranchAddress("FissionChamber",&FC);

		chain->SetBranchAddress("FC_Q1",&FC_Q1);
		chain->SetBranchAddress("FC_Anode_ID",&FC_Anode_ID);
		chain->SetBranchAddress("FC_FakeFission",&FC_FakeFission);
		chain->SetBranchAddress("LG_Tof",&LG_Tof);
		chain->SetBranchAddress("LG_ID",&LG_ID);
		chain->SetBranchAddress("LG_Q1",&LG_Q1);
		chain->SetBranchAddress("LG_Q2",&LG_Q2);
		chain->SetBranchAddress("HG_Tof",&HG_Tof);
		chain->SetBranchAddress("HG_ID",&HG_ID);
		chain->SetBranchAddress("HG_Q1",&HG_Q1);
		chain->SetBranchAddress("HG_Q2",&HG_Q2);

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
				
				if(FC_Anode_ID->size()>0){
						bool Fake = FC_FakeFission->at(0);
						int anode = FC_Anode_ID->at(0);


						int mysize = LG_Tof->size();
						for(int j=0; j<mysize; j++){
								// LG //
								int index_LG = (anode-1) + (LG_ID->at(j)-1)*NumberOfAnodes;
								double PSD = LG_Q2->at(j)/LG_Q1->at(j);
								if(LG_ID->at(j)>0 && anode>0 && Fake==0){
										hLG[index_LG]->Fill(LG_Tof->at(j));		
								}
						}

						mysize = HG_Tof->size();
						for(int j=0; j<mysize; j++){
								// HG //
								int index_HG = (anode-1) + (HG_ID->at(j)-1)*NumberOfAnodes;
								double PSD = HG_Q2->at(j)/HG_Q1->at(j);
								if(HG_ID->at(j)>0 && anode>0 && Fake==0 && HG_ID->size()==1){
										hHG[index_HG]->Fill(HG_Tof->at(j));			
								}
						}
				}
		}

		ofile->Write();
		ofile->Close();

}
