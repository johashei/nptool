TChain* chain=NULL ;

////////////////////////////////////////////////////////////////////////////////
void LoadChain(){
    chain = new TChain("PhysicsTree");
    chain->Add("./root/analysis/testmult16.root");
}

////////////////////////////////////////////////////////////////////////////////
void ShowTrack(){
    LoadChain();
    TH1 *h[20];
    TGraph *g[20];
    string command;
    string histname;

    TCanvas* c1 = new TCanvas("YX","YX",0,0,1200,600);
    TCanvas* c2 = new TCanvas("ZX","ZX",0,0,1200,600);
    c1->Divide(4,4);
    c2->Divide(4,4);
    for(int i=0; i<16; i++){
        c1->cd(i+1);
        command ="PosY:PosX>>hYX_"+to_string(i)+"(100,0,300,400,-200,200)";
        chain->Draw(command.c_str(),"","",1,i); 
        //histname="hYX_"+to_string(i);
        //h[i]= (TH1F*) gDirectory->FindObjectAny(histname.c_str());
        //g[i] = new TGraph(h[i]);
        //g[i]->SetMarkerStyle(3) ;
        //g[i]->Draw("AP") ;

        c2->cd(i+1);
        command ="PosZ/10.:PosX>>hZX_"+to_string(i)+"(100,0,300,120,-10,110)";
        chain->Draw(command.c_str(),"","*",1,i); 
        //histname="hWall2_"+to_string(i);
        //h[i-1]= (TH1F*) gDirectory->FindObjectAny(histname.c_str());
        //h[i-1]->SetMarkerColor(kRed) ;

        //gPad->WaitPrimitive();
    }
}
