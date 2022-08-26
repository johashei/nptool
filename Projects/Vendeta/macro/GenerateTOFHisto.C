TChain* chain;

int NumberOfDetectors=2;

/////////////////////////////////////
void LoadRootFile(){
	chain = new TChain("PhysicsTree");
	chain->Add("/home/faster/nptool/Outputs/Analysis/test.root");
}

/////////////////////////////////////
void GenerateTOFHisto(){
	LoadRootFile();

	TFile* ofile = new TFile("histo_tof_file.root","recreate");

	for(int i=0; i<NumberOfDetectors; i++){
		// LG // 		
		TString to_draw = Form("LG_Tof>>hLG_%i(1000,0,100)",i+1);
		TString condition = Form("LG_ID==%i",i+1);
		TString histo_name_LG = Form("hLG_%i",i+1);
		chain->Draw(to_draw,condition,"",1e6);
		TH1F* hLG = (TH1F*) gDirectory->FindObjectAny(histo_name_LG);

		// HG // 		
		to_draw = Form("HG_Tof>>hHG_%i(1000,0,100)",i+1);
		condition = Form("HG_ID==%i",i+1);
		TString histo_name_HG = Form("hHG_%i",i+1);
		chain->Draw(to_draw,condition,"",1e6);
		TH1F* hHG = (TH1F*) gDirectory->FindObjectAny(histo_name_HG);
	}

	ofile->Write();
	ofile->Close();
}
