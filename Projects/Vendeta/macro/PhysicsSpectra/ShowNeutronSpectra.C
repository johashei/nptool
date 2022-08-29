TCutG* cutn_LG=NULL;
TCutG* cutn_HG=NULL;
TCutG* cut_light_LG=NULL;
TCutG* cut_light_HG=NULL;

TChain* chain=NULL ;

///////////////////////////////////////////
void LoadCuts(){
    TFile* File_cutn_LG = new TFile("cuts/cutn_LG.root","READ");
    cutn_LG = (TCutG*) File_cutn_LG->FindObjectAny("cutn_LG");
}

///////////////////////////////////////////
void LoadChain(){
    chain = new TChain("PhysicsTree");
    chain->Add("/home/faster/nptool/Outputs/Analysis/test.root");
}

///////////////////////////////////////////
void ShowNeutronSpectra(){
    LoadChain();
    //LoadCuts();


}
