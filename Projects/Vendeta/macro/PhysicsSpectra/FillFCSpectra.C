TChain* chain;

int NumberOfAnodes= 11;
int nentries=1e6;

TH1F* hQ1[11];
TH1F* hQ2[11];
TH1F* hQmax[11];
TH2F* hQ2vsQ1[11];
TH2F* hQmaxvsQ1[11];
TH1F* hinE;
TH1F* hinTof;

/////////////////////////////////////
void LoadRootFile(){
		chain = new TChain("PhysicsTree");
		chain->Add("/home/faster/nptool/Outputs/Analysis/run31.root");
}

/////////////////////////////////////
void FillFCSpectra(){

		LoadRootFile();
		nentries = chain->GetEntries();
		cout << "Number of entries: " << nentries << endl;

		TFile* ofile = new TFile("histo_FC_run31.root","recreate");
		TH1F* hLG[791];
		TH1F* hHG[791];


		vector<double>* inToF = new vector<double>();
		vector<double>* inEnergy = new vector<double>();
		vector<double>* FC_Q1 = new vector<double>();
		vector<double>* FC_Q2 = new vector<double>();
		vector<double>* FC_Qmax = new vector<double>();
		vector<double>* FC_DT = new vector<double>();
		vector<int>* FC_Anode_ID = new vector<int>();
		vector<bool>* FC_FakeFission = new vector<bool>();

		TFissionChamberPhysics* FC = new TFissionChamberPhysics();
		chain->SetBranchAddress("FissionChamber",&FC);

		chain->SetBranchAddress("inToF",&inToF);
		chain->SetBranchAddress("inEnergy",&inEnergy);
		chain->SetBranchAddress("FC_Q1",&FC_Q1);
		chain->SetBranchAddress("FC_Q2",&FC_Q2);
		chain->SetBranchAddress("FC_Qmax",&FC_Qmax);
		chain->SetBranchAddress("FC_DT",&FC_DT);
		chain->SetBranchAddress("FC_Anode_ID",&FC_Anode_ID);
		chain->SetBranchAddress("FC_FakeFission",&FC_FakeFission);

		hinE = new TH1F("Incoming Energy","Incoming Energy",1000,0.5,10);
		hinTof = new TH1F("Incoming Tof","Incoming Tof",7200,0,1800);
		for(int i=0; i<NumberOfAnodes; i++){
				TString histo_name = Form("Q1_Anode%i",i+1);
				hQ1[i] = new TH1F(histo_name,histo_name,500,0,100e3);

				histo_name = Form("Q2_Anode%i",i+1);
				hQ2[i] = new TH1F(histo_name,histo_name,500,0,20e3);

				histo_name = Form("Qmax_Anode%i",i+1);
				hQmax[i] = new TH1F(histo_name,histo_name,500,0,10e3);

				histo_name = Form("Q2vsQ1_Anode%i",i+1);
				hQ2vsQ1[i] = new TH2F(histo_name,histo_name,500,0,100e3,500,0,20e3);

				histo_name = Form("QmaxvsQ1_Anode%i",i+1);
				hQmaxvsQ1[i] = new TH2F(histo_name,histo_name,500,0,100e3,500,0,10e3);

		}

		for(int i=0; i<nentries; i++){
				chain->GetEntry(i);

				if(i%100000==0){
						cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
				}

				if(inToF->size()==1){
						double inT = inToF->at(0);
						double inE = inEnergy->at(0);
	
						int anode = FC_Anode_ID->at(0);
						if(anode>0){
								double q1 = FC_Q1->at(0);
								double q2 = FC_Q2->at(0);
								double qmax = FC_Qmax->at(0);
								double DT = FC_DT->at(0);
								bool FakeFission = FC_FakeFission->at(0);
								if(inT>10 && inT<1790 && FakeFission==0 && DT>7e6){
										hQ1[anode-1]->Fill(q1);
										hQ2[anode-1]->Fill(q2);
										hQmax[anode-1]->Fill(qmax);
										hQ2vsQ1[anode-1]->Fill(q1,q2);
										hQmaxvsQ1[anode-1]->Fill(q1,qmax);
										hinE->Fill(inE);
										hinTof->Fill(inT);
								}
						}
				}
		}

		ofile->Write();		
		ofile->Close();

}
