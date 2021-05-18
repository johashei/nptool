void CheckPPAC(int id=0){
    auto tree = new TChain("PhysicsTree");
    tree->Add("../../Outputs/Analysis/PhysicsTree.root");

    if(id==0) {  
        TString Cond="";
        auto c = new TCanvas();
        c->Divide(3,2);
        c->cd(1);
        tree->Draw("PPAC.TSumX:PPAC.ID>>hTSumX(40,0,40,400,0,200)",Cond,"colz");
        c->cd(2);
        tree->Draw("PPAC.TSumY:PPAC.ID>>hTSumY(40,0,40,400,0,200)",Cond,"colz");
        c->cd(3);
        tree->Draw("PPAC.TDiffX:PPAC.ID>>hTDiffX(40,0,40,400,-100,100)",Cond,"colz");
        c->cd(4);
        tree->Draw("PPAC.TDiffY:PPAC.ID>>hTDiffY(40,0,40,400,-100,100)",Cond,"colz");
        c->cd(5);
        tree->Draw("PPAC.X:PPAC.ID>>hX(40,0,40,600,-120,120)",Cond,"col");
        c->cd(6);
        tree->Draw("PPAC.Y:PPAC.ID>>hY(40,0,40,600,-120,120)",Cond,"col");
    } else { 
        TString Cond="ID==19";
        auto c = new TCanvas();
        c->Divide(3,2);
        c->cd(1);
        tree->Draw("PPAC.TSumX>>hTSumX(120,165,195)",Cond,"");
        c->cd(2);
        tree->Draw("PPAC.TDiffX>>hTDiffX(400,-40,40)",Cond,"");
        c->cd(3);
        tree->Draw("PPAC.X>>hX(240,-30,30)",Cond,"");
        c->cd(4);
        tree->Draw("PPAC.TSumY>>hTSumY(120,90,120)",Cond,"");
        c->cd(5);
        tree->Draw("PPAC.TDiffY>>hTDiffY(400,-40,40)",Cond,"");
        c->cd(6);
        tree->Draw("PPAC.Y>>hY(240,-30,30)",Cond,"");

    }

}
