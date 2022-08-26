TChain* chain;

int NumberOfDetectors= 2;
int NumberOfAnodes= 11;

/////////////////////////////////////
void LoadRootFile(){
	chain = new TChain("PhysicsTree");
	chain->Add("/home/faster/nptool/Outputs/Analysis/test_calib.root");
}

/////////////////////////////////////
void GenerateTOFHisto(){
	LoadRootFile();

	TFile* ofile = new TFile("histo_tof_file.root","recreate");

	for(int i=0; i<NumberOfDetectors; i++){
		for(int j=0; j<NumberOfAnodes; j++){
			// LG // 		
			TString to_draw = Form("LG_Tof>>hLG_Det%i_Anode%i(1000,0,100)",i+1,j+1);
			TString condition = Form("LG_ID==%i && LG_Anode_ID==%i",i+1,j+1);
			TString histo_name_LG = Form("hLG_Det%i_Anode%i",i+1,j+1);
			chain->Draw(to_draw,condition,"",1e6);
			TH1F* hLG = (TH1F*) gDirectory->FindObjectAny(histo_name_LG);

			// HG // 		
			to_draw = Form("HG_Tof>>hHG_Det%i_Anode%i(1000,0,100)",i+1,j+1);
			
			cout << to_draw << endl;
			condition = Form("HG_ID==%i && HG_Anode_ID==%i",i+1,j+1);
			TString histo_name_HG = Form("hHG_Det%i_Anode%i",i+1,j+1);
			chain->Draw(to_draw,condition,"",1e6);
			TH1F* hHG = (TH1F*) gDirectory->FindObjectAny(histo_name_HG);
		}
	}

	ofile->Write();
	ofile->Close();
}
