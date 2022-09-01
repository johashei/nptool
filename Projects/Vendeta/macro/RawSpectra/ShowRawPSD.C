TChain* chain;

int NumberOfDetectors= 72;
int NumberOfAnodes= 1;
int NumberOfEvents= 1e7;

//////////////////////////////////////////////////
void OpenRootFile(){
		chain = new TChain("RawTree");
		//chain->Add("/home/faster/fastertonptool/data/rootfiles/V4B_SAMPLING_6_0001.root");
 	chain->Add("/home/faster/fastertonptool/data/rootfiles/test_cf4_q1_160ns.root");
}

//////////////////////////////////////////////////
void ShowRawPSD(){
		OpenRootFile();

		TFile* ofile = new TFile("PSD_Q1_160ns.root","recreate");

		// Canvas Definition for Low Gain //
		TCanvas* cLG_RI   = new TCanvas("Det-LG RI","Det-LG RI",1800,1800);
		TCanvas* cLG_RII  = new TCanvas("Det-LG RII","Det-LG RII",1800,1800);
		TCanvas* cLG_RIII = new TCanvas("Det-LG RIII","Det-LG RIII",1800,1800);
		TCanvas* cLG_LI   = new TCanvas("Det-LG LI","Det-LG LI",1800,1800);
		TCanvas* cLG_LII  = new TCanvas("Det-LG LII","Det-LG LII",1800,1800);
		TCanvas* cLG_LIII = new TCanvas("Det-LG LIII","Det-LG LIII",1800,1800);

		cLG_RI->Divide(3,4);
		cLG_RII->Divide(3,4);
		cLG_RIII->Divide(3,4);
		cLG_LI->Divide(3,4);
		cLG_LII->Divide(3,4);
		cLG_LIII->Divide(3,4);

		// Canvas Definition for High Gain //
		TCanvas* cHG_RI   = new TCanvas("Det-HG RI","Det-HG RI",1800,1800);
		TCanvas* cHG_RII  = new TCanvas("Det-HG RII","Det-HG RII",1800,1800);
		TCanvas* cHG_RIII = new TCanvas("Det-HG RIII","Det-HG RIII",1800,1800);
		TCanvas* cHG_LI   = new TCanvas("Det-HG LI","Det-HG LI",1800,1800);
		TCanvas* cHG_LII  = new TCanvas("Det-HG LII","Det-HG LII",1800,1800);
		TCanvas* cHG_LIII = new TCanvas("Det-HG LIII","Det-HG LIII",1800,1800);

		cHG_RI->Divide(3,4);
		cHG_RII->Divide(3,4);
		cHG_RIII->Divide(3,4);
		cHG_LI->Divide(3,4);
		cHG_LII->Divide(3,4);
		cHG_LIII->Divide(3,4);

		for(int i=0; i<NumberOfDetectors; i++){
				// LG //	
				TString to_draw_LG = Form("fVendeta_LG_Q2/fVendeta_LG_Q1:fVendeta_LG_Q1>>hLG_det%i(500,0,500e3,500,0,1)",i+1);
				TString condition_LG = Form("fVendeta_LG_DetectorNbr==%i",i+1);

				// HG //
				TString to_draw_HG = Form("fVendeta_HG_Q2/fVendeta_HG_Q1:fVendeta_HG_Q1>>hHG_det%i(500,0,800e3,500,0,1)",i+1);
				TString condition_HG = Form("fVendeta_HG_DetectorNbr==%i",i+1);

				cout << to_draw_LG << endl;
				
				if(i+1<12){
						cLG_LI->cd(i+1);
						chain->Draw(to_draw_LG,condition_LG,"colz",NumberOfEvents);

						cHG_LI->cd(i+1);
						chain->Draw(to_draw_HG,condition_HG,"colz",NumberOfEvents);
				}
				if(i+1>12 && i+1<24){
						cLG_LII->cd(i+1-12);
						chain->Draw(to_draw_LG,condition_LG,"colz",NumberOfEvents);
						cHG_LII->cd(i+1-12);
						chain->Draw(to_draw_HG,condition_HG,"colz",NumberOfEvents);
				}
				if(i+1>24 && i+1<36){
						cLG_LIII->cd(i+1-24);
						chain->Draw(to_draw_LG,condition_LG,"colz",NumberOfEvents);

						cHG_LIII->cd(i+1-24);
						chain->Draw(to_draw_HG,condition_HG,"colz",NumberOfEvents);
				}
				if(i+1>36 && i+1<48){
						cLG_RI->cd(i+1-36);
						chain->Draw(to_draw_LG,condition_LG,"colz",NumberOfEvents);

						cHG_RI->cd(i+1-36);
						chain->Draw(to_draw_HG,condition_HG,"colz",NumberOfEvents);
				}
				if(i+1>48 && i+1<60){
						cLG_RII->cd(i+1-48);
						chain->Draw(to_draw_LG,condition_LG,"colz",NumberOfEvents);

						cHG_RII->cd(i+1-48);
						chain->Draw(to_draw_HG,condition_HG,"colz",NumberOfEvents);
				}
				if(i+1>60){
						cLG_RIII->cd(i+1-60);
						chain->Draw(to_draw_LG,condition_LG,"colz",NumberOfEvents);
				
						cHG_RIII->cd(i+1-60);
						chain->Draw(to_draw_HG,condition_HG,"colz",NumberOfEvents);
				}
		}

		ofile->Write();
		ofile->Close();
}
