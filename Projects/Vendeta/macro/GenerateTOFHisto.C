TChain* chain;

int NumberOfDetectors= 72;
int NumberOfAnodes= 1;
int nentries=1e6;
/////////////////////////////////////
void LoadRootFile(){
	chain = new TChain("PhysicsTree");
	//chain->Add("/home/faster/nptool/Outputs/Analysis/test_cf4_q1_160ns.root");
	chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_1.root");
	chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_2.root");
	chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_3.root");
	chain->Add("/home/faster/nptool/Outputs/Analysis/test_sampler_qdc_cf_4.root");
}

/////////////////////////////////////
void GenerateTOFHisto(){
	LoadRootFile();

	nentries = chain->GetEntries();
	TFile* ofile = new TFile("histo_tof_file.root","recreate");

	
	for(int i=0; i<NumberOfDetectors; i++){
		for(int j=0; j<NumberOfAnodes; j++){
				j=5;
			// LG // 		
			TString to_draw = Form("LG_Tof>>hLG_Det%i_Anode%i(5000,-200,300)",i+1,j+1);
			TString condition = Form("LG_ID==%i && LG_Anode_ID==%i && FC_Q1>5500",i+1,j+1);
			TString histo_name_LG = Form("hLG_Det%i_Anode%i",i+1,j+1);
			chain->Draw(to_draw,condition,"",1e6);
			TH1F* hLG = (TH1F*) gDirectory->FindObjectAny(histo_name_LG);

			// HG // 		
			to_draw = Form("HG_Tof>>hHG_Det%i_Anode%i(3000,0,300)",i+1,j+1);
			
			cout << to_draw << endl;
			condition = Form("HG_ID==%i && HG_Anode_ID==%i && FC_Q1>5500",i+1,j+1);
			TString histo_name_HG = Form("hHG_Det%i_Anode%i",i+1,j+1);
			chain->Draw(to_draw,condition,"",1e6);
			TH1F* hHG = (TH1F*) gDirectory->FindObjectAny(histo_name_HG);
		}
	}

	ofile->Write();
	ofile->Close();
}
