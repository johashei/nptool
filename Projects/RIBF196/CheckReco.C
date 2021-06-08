void CheckReco(){
    auto tree = new TChain("PhysicsTree");
    //tree->Add("../../Outputs/Analysis/PhysicsTree.root");
    tree->Add("../../Outputs/Analysis/ribf196_0030_phy.root");

    //TString Cond="";
    auto c = new TCanvas();
    c->Divide(3,2);
    c->cd(1);
    tree->Draw("Rec35.delta>>hdeltaBrho35(100,-10,10)","","");
    c->cd(2);
    tree->Draw("Rec57.delta>>hdeltaBrho57(100,-10,10)","","");
    c->cd(3);
    tree->Draw("Rec37.z:Rec37.aoq>>haoq37(1000,2.5,2.8,1000,25,35)","","colz");
    c->cd(4);
    tree->Draw("Rec35.z:Rec35.aoq>>haoq35(1000,2.5,2.8,1000,25,35)","","colz");
    c->cd(5);
    tree->Draw("Rec57.z:Rec57.aoq>>haoq57(1000,2.5,2.8,1000,25,35)","","colz");
    c->cd(6);
    tree->Draw("Rec37.z:aoqc37>>haoqc37(1000,2.5,2.8,1000,25,35)","","colz");

    auto c2 = new TCanvas();
    c2->Divide(3,2);
    c2->cd(1);
    tree->Draw("Rec89.delta>>hdeltaBrho89(100,-10,10)","","");
    c2->cd(2);
    tree->Draw("Rec911.delta>>hdeltaBrho911(100,-10,10)","","");
    c2->cd(3);
    tree->Draw("Rec811.z:Rec811.aoq>>haoq811(1000,2.5,2.8,1000,25,35)","","colz");
    c2->cd(4);
    tree->Draw("Rec89.z:Rec89.aoq>>haoq89(1000,2.5,2.8,1000,25,35)","","colz");
    c2->cd(5);
    tree->Draw("Rec911.z:Rec911.aoq>>haoq911(1000,2.5,2.8,1000,25,35)","","colz");
    c2->cd(6);
    tree->Draw("Rec811.z:aoqc811>>haoqc811(1000,2.5,2.8,1000,25,35)","","colz");
}
