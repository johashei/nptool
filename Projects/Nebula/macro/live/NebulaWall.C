#include "/local/lemair/nptool/NPLib/Detectors/Nebula/TNebulaPhysics.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::find
#include <vector>   

void BuildBox(double Z, double X, double Zsize, double Xsize, Color_t color){
  TBox* box = new TBox(Z,X,Zsize,Xsize); //(Position of bottom-left corner, Position of top-right corner)
  box->SetFillStyle(1001);
  box->SetFillColor(color);
  box->SetLineColor(kBlack);
  box->SetLineWidth(1);

  box->Draw("l");
}

void DrawWall(std::vector<int> v_ID){
  vector<double> RefZ = {11060, 11181, 11910, 12031};
  double RefX = 0;
  double Z,X;
  int ID;
  Color_t color;

  for(int sublayer = 0; sublayer < 4; sublayer ++){
    for(int XOffset = -15; XOffset<15; XOffset++){
      //CHECK IF BAR WAS HIT
      ID = 30*sublayer + (XOffset+16);
      if(std::find(v_ID.begin(), v_ID.end(), ID)!=v_ID.end()){
        color = kOrange-8;
      }else{
        color = kCyan-8;
      }

      //DRAW BAR
      Z = RefZ[sublayer];
      X = RefX + XOffset*121;
      BuildBox(Z,X,Z+120,X+120, color);
    }
  }
}

void NebulaWall(){
  //SETTING THE VISUAL INTERFACE
  TCanvas* c = new TCanvas("c", "Nebula Walls", 1920, 1200);
  gPad->SetPad("pad", "New Pad", 0, 0, 1, 1, /*backgrd color*/ 0, /*line thickness*/ 0, /*stamping style*/ 0);
  gPad->DrawFrame(10000,-2000,13500,2000);
  c->SetGrid();
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.1);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.05);
  TFrame* frame = new TFrame(0,0,100,100);


  //CONNECTING TO TREE AND DATA ANALYSIS
  TFile *f = TFile::Open("root/analysis/PhysicsTree.root");
  TTree *t = (TTree*)f->Get("PhysicsTree");

  TNebulaPhysics *Neb_Phys = new TNebulaPhysics();
  std::vector<int> v_ID;
  int vec_size;
  TBranch* branch;
  t->SetBranchAddress("Nebula", &Neb_Phys, &branch);
  int nentries = branch->GetEntries();

  for(int i = 0; i < nentries; i++){
    t->GetEntry(i);
    v_ID = Neb_Phys->DetectorNumber;
    vec_size = v_ID.size(); 

    if(vec_size>=4){
      DrawWall(v_ID);
      gPad->Modified(); gPad->Update(); // make sure it is (re)drawn
      gSystem->ProcessEvents(); gSystem->ProcessEvents();
      gSystem->Sleep(1000);
      
      cout << "size " << vec_size << endl;
      cout << "event " << i << endl << endl;

      //for(int j = 0 ; j < vec_size ; j++)
    }
  }
}
