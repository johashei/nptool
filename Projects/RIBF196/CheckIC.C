void CheckIC(int id=0){
    auto tree = new TChain("PhysicsTree");
    tree->Add("../../Outputs/Analysis/PhysicsTree.root");

    if(id==0) {  
        //TString Cond="";
        auto c = new TCanvas();
        c->Divide(3,3);
        c->cd(1);
        tree->Draw("IC.E:IC.E_ID>>hEvsID(9,0,9,1000,0,500)","","colz");
        c->cd(2);
        tree->Draw("IC.E:IC.E_Layer>>hEvsLayer_F7(9,0,9,1000,0,500)","IC.E_ID==1","colz");
        c->cd(3);
        tree->Draw("IC.E:IC.E_Layer>>hEvsLayer_F11(9,0,9,1000,0,500)","IC.E_ID==2","colz");
        c->cd(4);
        tree->Draw("IC.NLayerFired>>hFired_F7(9,0,9)","IC.ID==1","");
        c->cd(5);
        tree->Draw("IC.NLayerFired>>hFired_F11(9,0,9)","IC.ID==2","");
        c->cd(6);
        tree->Draw("IC.CalAvSum>>hCalAvSum_F7(3000,0,1000)","IC.ID==1","");
        c->cd(7);
        tree->Draw("IC.CalAvSum>>hCalAvSum_F11(3000,0,1000)","IC.ID==2","");
    }

}
