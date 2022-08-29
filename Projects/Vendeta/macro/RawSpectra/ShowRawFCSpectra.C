TChain* chain;

int NumberOfAnodes= 1;
int NumberOfEvents= 1e6;

//////////////////////////////////////////////////
void OpenRootFile(){
		chain = new TChain("RawTree");
		chain->Add("/home/faster/fastertonptool/data/rootfiles/V4B_SAMPLING_6_0001.root");
}

//////////////////////////////////////////////////
void ShowRawFCSpectra(string Nucleus="Cf"){
		OpenRootFile();

		if(Nucleus=="Cf" || Nucleus=="252Cf")
				NumberOfAnodes=1;
		if(Nucleus=="U" || Nucleus=="238U")
				NumberOfAnodes=11;

		// Canvas Definition for Low Gain //
		TCanvas* c1   = new TCanvas("Charge Q1","Charge Q1",1800,1800);
		TCanvas* c2   = new TCanvas("Charge Q2","Charge Q2",1800,1800);
		TCanvas* c3   = new TCanvas("Charge Qmax","Charge Qmax",1800,1800);
		TCanvas* c4   = new TCanvas("Q2 vs Q1","Q2 vs Q1",1800,1800);
		TCanvas* c5   = new TCanvas("Qmax vs Q1","Qmax vs Q1",1800,1800);

		if(Nucleus=="U" || Nucleus=="238U"){
				c1->Divide(3,4);
				c2->Divide(3,4);
				c3->Divide(3,4);
				c4->Divide(3,4);
				c5->Divide(3,4);
		}

		for(int i=0; i<NumberOfAnodes; i++){
				// Draw //
				TString draw_Q1 = Form("fFC_Q1>>hQ1_Anode%i(500,0,100e3)",i+1);
				TString draw_Q2 = Form("fFC_Q2>>hQ2_Anode%i(500,0,20e3)",i+1);
				TString draw_Qmax = Form("fFC_Qmax>>hQmax_Anode%i(500,0,10e3)",i+1);
				TString draw_Q2vsQ1 = Form("fFC_Q2:fFC_Q1>>hQ2vsQ1_Anode%i(500,0,100e3,500,0,20e3)",i+1);
				TString draw_QmaxvsQ1 = Form("fFC_Qmax:fFC_Q1>>hmax2vsQ1_Anode%i(500,0,100e3,500,0,10e3)",i+1);

				TString condition;
				if(Nucleus=="Cf" || Nucleus=="252Cf")
						condition = Form("fFC_AnodeNbr==%i",6);

				if(Nucleus=="U" || Nucleus=="238U")
						condition = Form("fFC_AnodeNbr==%i",i+1);

				if(Nucleus=="U")c1->cd(i+1);
				else c1->cd();
				chain->Draw(draw_Q1,condition,"",NumberOfEvents);
				c1->Update();

				if(Nucleus=="U")c2->cd(i+1);
				else c2->cd();
				chain->Draw(draw_Q2,condition,"",NumberOfEvents);

				if(Nucleus=="U")c3->cd(i+1);
				else c3->cd();
				chain->Draw(draw_Qmax,condition,"",NumberOfEvents);

				if(Nucleus=="U")c4->cd(i+1)->SetLogz();
				else c4->cd()->SetLogz();
				chain->Draw(draw_Q2vsQ1,condition,"colz",NumberOfEvents);

				if(Nucleus=="U") c5->cd(i+1)->SetLogz();
				else c5->cd()->SetLogz();
				chain->Draw(draw_QmaxvsQ1,condition,"colz",NumberOfEvents);
		}
}
