TChain* chain;

int NumberOfAnodes= 11;
int nentries= 1e6;

//////////////////////////////////////////////////
void OpenRootFile(){
		chain = new TChain("RawTree");
		chain->Add("/home/faster/fastertonptool/data/rootfiles/run31_*.root");
}

//////////////////////////////////////////////////
void FillRawFCSpectra(){
		OpenRootFile();
		nentries = chain->GetEntries();

		TFile* ofile = new TFile("FC_Raw_spectra_run31.root","recreate");
		TH1F* Q1[11];
		TH1F* Q2[11];
		TH1F* Qmax[11];
		TH2F* Q2vsQ1[11];
		TH2F* QmaxvsQ1[11];

		TFissionChamberData* FC = new TFissionChamberData();
		chain->SetBranchAddress("FissionChamber",&FC);

		TString histo_name;
		for(int i=0; i< NumberOfAnodes; i++){
				histo_name = Form("Q1_Anode%i",i+1);
				Q1[i] = new TH1F(histo_name, histo_name, 500,0,100e3);

				histo_name = Form("Q2_Anode%i",i+1);
				Q2[i] = new TH1F(histo_name, histo_name, 500,0,20e3);

				histo_name = Form("Qmax_Anode%i",i+1);
				Qmax[i] = new TH1F(histo_name, histo_name, 500,0,10e3);

				histo_name = Form("Q2vsQ1_Anode%i",i+1);
				Q2vsQ1[i] = new TH2F(histo_name, histo_name, 500,0,100e3,500,0,20e3);

				histo_name = Form("QmaxvsQ1_Anode%i",i+1);
				QmaxvsQ1[i] = new TH2F(histo_name, histo_name, 500,0,100e3,500,0,10e3);
		}

		for(int i=0; i<nentries; i++){
				chain->GetEntry(i);

				if(i%100000==0){
						cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
				}
				int mult = FC->GetMultiplicity();
				for(int j=0; j<mult; j++){

						int anode = FC->GetAnodeNbr(j);
						if(anode>0){
								double Q1val = FC->GetQ1(j);
								double Q2val = FC->GetQ2(j);
								double Qmaxval = FC->GetQmax(j);

								Q1[anode-1]->Fill(Q1val);
								Q2[anode-1]->Fill(Q2val);
								Qmax[anode-1]->Fill(Qmaxval);
								Q2vsQ1[anode-1]->Fill(Q1val,Q2val);
								QmaxvsQ1[anode-1]->Fill(Q1val,Q2val);
						}

				}
		}

		ofile->Write();
		ofile->Close();
}
