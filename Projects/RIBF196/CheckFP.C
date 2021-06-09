void CheckFP(int id=3){
    auto tree = new TChain("PhysicsTree");
    //tree->Add("../../Outputs/Analysis/PhysicsTree.root");
    tree->Add("../../Outputs/Analysis/ribf196_0030_phy.root");

    TString Cond="";
    auto c = new TCanvas();
    c->Divide(3,2);
    string buffer, command;
    buffer = Form("RecFP[%i]",id);
    c->cd(1);
    command = buffer + "[0]>>hfX(300,-150,150)";
    cout<< command<<endl;
    tree->Draw(command.c_str(),Cond,"colz");
    c->cd(2);
    command = buffer + "[2]>>hfY(300,-150,150)";
    cout<< command<<endl;
    tree->Draw(command.c_str(),Cond,"colz");
    c->cd(3);
    command = buffer + "[2]:"+ buffer +"[0]>>hfXY(300,-150,150,300,-150,150)";
    cout<< command<<endl;
    tree->Draw(command.c_str(),Cond,"colz");
    c->cd(4);
    command = buffer + "[1]>>hfA(200,-50,50)";
    cout<< command<<endl;
    tree->Draw(command.c_str(),Cond,"colz");
    c->cd(5);
    command = buffer + "[3]>>hfB(200,-50,50)";
    cout<< command<<endl;
    tree->Draw(command.c_str(),Cond,"colz");
    c->cd(5);
    command = buffer + "[3]:"+ buffer +"[1]>>hfAB(200,-50,50,200,-50,50)";
    c->cd(6);

}
