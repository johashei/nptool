TChain* chain;

int NumberOfDetectors= 72;
int nentries=1e6;

//////////////////////////////////////////////////
void OpenRootFile(){
		chain = new TChain("RawTree");
 	chain->Add("/home/faster/fastertonptool/data/rootfiles/run59_*.root");
}

//////////////////////////////////////////////////
void FillRawPSD(){
		OpenRootFile();
		nentries = chain->GetEntries();
		cout << "Number of entries= " << nentries << endl;

		TFile* ofile = new TFile("PSD_histo_run59.root","recreate");
		TH2F* hLG[72];
		TH2F* hHG[72];

		TString histo_name;
		for(int i=0; i<NumberOfDetectors;i++){
				histo_name = Form("hLG_det%i",i+1);
				hLG[i] = new TH2F(histo_name,histo_name,500,0,500e3,500,0,1.1);
				histo_name = Form("hHG_det%i",i+1);
				hHG[i] = new TH2F(histo_name,histo_name,500,0,900e3,500,0,1.1);
		}

		TVendetaData* Vendeta = new TVendetaData();
		chain->SetBranchAddress("Vendeta",&Vendeta);

		for(int i=0; i<nentries; i++){
			chain->GetEntry(i);
			if(i%100000==0){
				cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
			}

			// LG //
			int mult_LG = Vendeta->GetLGMultEnergy();
			if(mult_LG>0){
					for(int j=0; j<mult_LG; j++){
							double Q1 = Vendeta->GetLGQ1(j);
							double Q2 = Vendeta->GetLGQ2(j);
							double PSD_LG = Q2/Q1;
							int det = Vendeta->GetLGDetectorNbr(j);

							hLG[det-1]->Fill(Q1,PSD_LG);
					}
			}
				
			// HG //
			int mult_HG = Vendeta->GetHGMultEnergy();
			if(mult_HG>0){
					for(int j=0; j<mult_HG; j++){
							double Q1 = Vendeta->GetHGQ1(j);
							double Q2 = Vendeta->GetHGQ2(j);
							double PSD_HG = Q2/Q1;
							int det = Vendeta->GetHGDetectorNbr(j);

							hHG[det-1]->Fill(Q1,PSD_HG);
					}
			}
		

		}
	
		ofile->Write();
		ofile->Close();
}
