TChain* chain;

int NumberOfAnodes=11;
//////////////////////////////////////////////
void OpenRootFile(){
		chain = new TChain("PhysicsTree");
		chain->Add("/home/faster/nptool/Outputs/Analysis/test.root");
}

//////////////////////////////////////////////
void ShowIncomingNeutron(){
		OpenRootFile();

		TCanvas* cToF = new TCanvas("cToF","cToF",1800,1800);
		TCanvas* cE = new TCanvas("cE","cE",1800,1800);

		cToF->Divide(3,4);
		cE->Divide(3,4);

		for(int i=0; i<NumberOfAnodes; i++){
				TString condition = Form("LG_Anode_ID==%i",i+1);
				TString to_draw;

				cToF->cd(i+1);
				to_draw = Form("inToF>>hToF_Anode%i(5400,0,1800)",i+1);
				chain->Draw(to_draw,condition);

				cE->cd(i+1);
				to_draw = Form("inEnergy>>hE_Anode%i(1000,0,700)",i+1);
				chain->Draw(to_draw,condition);
		}
}
