void CheckPlastic(int id=2){
    auto tree = new TChain("PhysicsTree");
    tree->Add("../../Outputs/Analysis/PhysicsTree.root");

    if(id==0) {  
        TString Cond="";
        auto c = new TCanvas();
        c->Divide(3,2);
        c->cd(1);
        tree->Draw("Plastic.TL:Plastic.ID>>hTL(20,0,20,1400,0,140)",Cond,"colz");
        c->cd(2);
        tree->Draw("Plastic.TR:Plastic.ID>>hTR(20,0,20,1400,0,140)",Cond,"colz");
        c->cd(3);
        tree->Draw("Plastic.QL:Plastic.ID>>hQL(20,0,20,500,0,2000)",Cond,"colz");
        c->cd(4);
        tree->Draw("Plastic.QR:Plastic.ID>>hQR(20,0,20,500,0,2000)",Cond,"colz");
        c->cd(5);
        tree->Draw("Plastic.TLSlew:Plastic.ID>>hTLSlew(20,0,20,1400,0,140)",Cond,"colz");
        c->cd(6);
        tree->Draw("Plastic.TRSlew:Plastic.ID>>hTRSlew(20,0,20,1400,0,140)",Cond,"colz");
    } else { 
        //TString Cond="ID==3";
        auto c = new TCanvas();
        c->Divide(2,2);
        c->cd(1);
        tree->Draw("TMath::Log(Plastic.QL/Plastic.QR):Plastic.TR-Plastic.TL>>hCorrF3(360,-6,6,1000,-2,2)","Plastic.ID==1 && Plastic.fired","colz");
        c->cd(2);
        tree->Draw("TMath::Log(Plastic.QL/Plastic.QR):Plastic.TR-Plastic.TL>>hCorrF7(360,-6,6,1000,-2,2)","Plastic.ID==3 && Plastic.fired","colz");
        c->cd(3);
        tree->Draw("TMath::Log(Plastic.QL/Plastic.QR):Plastic.TR-Plastic.TL>>hCorrF8(360,-6,6,1000,-2,2)","Plastic.ID==4 && Plastic.fired","colz");
        c->cd(4);
        tree->Draw("TMath::Log(Plastic.QL/Plastic.QR):Plastic.TR-Plastic.TL>>hCorrF11(360,-6,6,1000,-2,2)","Plastic.ID==5 && Plastic.fired","colz");
    }

}
