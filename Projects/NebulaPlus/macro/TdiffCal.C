TFile* file=NULL ;

////////////////////////////////////////////////////////////////////////////////
void LoadFile(){
    file = new TFile("../root/analysis/hTdiff_ID_run18.root");
}

////////////////////////////////////////////////////////////////////////////////
void TdiffCal(bool savehist=true){
    LoadFile();
    TH2F* hQu_ID = (TH2F*) gDirectory->FindObjectAny("hTdiff_ID");

    TCanvas* c1 = new TCanvas("Raw","Raw",0,0,1200,600);
    c1->cd();
    //c1->SetLogy();
 
    TH1F* hQu[90];
    Double_t mean=0;
    string histname;
    bool found = false;

    //RANGE FOR THE MUON PEAK SEARCH
    double range_inf= -30;
    double range_sup= 30;
    
    ofstream peaks_file;
    peaks_file.open("tmp.txt");
    for(int i=1; i<93; i++){

	    histname = "hQu_"+to_string(i);
	    cout<<i<<"histname\t"<<histname<<endl;
	    hQu_ID->ProjectionY(histname.c_str(),i+1,i+1);   
	    hQu[i-1] = (TH1F*) gDirectory->FindObjectAny(histname.c_str());
	    hQu[i-1]->Draw("");
	    hQu[i-1]->SetTitle(histname.c_str());
	    hQu[i-1]->GetXaxis()->SetRangeUser(range_inf,range_sup);
	    mean = hQu[i-1]->GetMean(1);
	    //if(!found) cout<<"BarID:"<<i<<"UP, problem in finding peaks"<<endl; 
	    //gPad->WaitPrimitive();
	    peaks_file<<"NebulaPlus_Tdiff_Offset_"<<i<<"\t"<<mean<<endl;
	    cout<<i<<"mean:\t"<<mean<<endl;

    }
    peaks_file.close();
}
