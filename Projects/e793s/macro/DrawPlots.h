#include "NPReaction.h"
#include <string>
#include <sstream>

using namespace std;

TChain* chain=NULL ;
char cond[1000];

NPL::Reaction Kdp("47K(d,p)48K@362");
NPL::Reaction Kdt("47K(d,t)46K@362");
NPL::Reaction Kdd("47K(d,d)47K@362");
NPL::Reaction Kpp("47K(p,p)47K@362");
NPL::Reaction K12C12C("47K(12C,12C)47K@362");
NPL::Reaction Tidp("47Ti(d,p)48Ti@362");
NPL::Reaction Tidt("47Ti(d,t)46Ti@362");
NPL::Reaction Tidd("47Ti(d,d)47Ti@362");
NPL::Reaction Ti12C12C("47Ti(12C,12C)47Ti@362");

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
TChain* Chain(std::string TreeName, std::vector<std::string>& file, bool EventList){
  TChain*  chain = new TChain(TreeName.c_str());
  unsigned int size =file.size(); 
    for(unsigned int i = 0 ; i < size ; i++){
      cout << "Adding " << file[i] << endl;
      chain->Add(file[i].c_str());
    } 
  return chain;
}
/////////////////////////////////////
void LoadChainNP(){
  vector<string> files;
  //files.push_back("../../../Outputs/Analysis/47K_Full_09June_MG3_Target.root");
  files.push_back("../../../Outputs/Analysis/47K_Full_25Jun_TestingAddBack2.root");
  chain = Chain("PhysicsTree",files,true);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void plot_kine(NPL::Reaction r, double Ex,Color_t c,int w, int s){
  r.SetExcitation4(Ex);
  TGraph* g= r.GetKinematicLine3();
  g->SetLineColor(c) ;
  g->SetLineStyle(s) ;
  g->SetLineWidth(w) ;
  g->Draw("c");
}
/////////////////////////////////////
void plot_state(double Ex,double max,Color_t c,int w, int s){
  TLine* line = new TLine(Ex,0,Ex,max) ; 
  line->SetLineColor(c) ;
  line->SetLineStyle(s) ;
  line->SetLineWidth(w) ;
  line->Draw();
}
/////////////////////////////////////
void AddTiStates(double E){
 NPL::Reaction Ti("47Ti(d,p)48Ti@362");
 Ti.SetExcitationHeavy(E);
 auto g = Ti.GetKinematicLine3();
 g->SetLineWidth(1);
 g->SetLineStyle(2);
 g->Draw("c");
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Draw_1DGamma(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  chain->Draw("AddBack_EDC>>Eg(2500,0,5)","abs(T_MUGAST_VAMOS-2777)<600");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->SetTitle("Egamma (using 25Jun Full run)");
  Eg->GetXaxis()->SetTitle("Eg [MeV]");
}
/////////////////////////////////////
void Draw_1DParticle(){
  TCanvas *cEx = new TCanvas("cEx","cEx",1000,1000);
  chain->Draw("Ex>>Ep(220,-1,10)","abs(T_MUGAST_VAMOS-2777)<600");
  TH1F* Ep = (TH1F*) gDirectory->Get("Ep");
  Ep->SetTitle("Ex (using 25Jun Full run)");
  Ep->GetXaxis()->SetTitle("Ex [MeV]");
}
/////////////////////////////////////
void Draw_2DParticleGamma(){
  TCanvas *cExEg = new TCanvas("cExEg","cExEg",1000,1000);
  chain->Draw("AddBack_EDC:Ex>>ExEg(100,-1,5,2500,0,5)",
		  "abs(T_MUGAST_VAMOS-2777)<600","colz");
  TH1F* ExEg = (TH1F*) gDirectory->Get("ExEg");
  ExEg->SetTitle("Ex-Egamma (using 25Jun Full run)");
  ExEg->GetXaxis()->SetTitle("Ex [MeV]");
  ExEg->GetYaxis()->SetTitle("Eg [MeV]");
  TLine *XeqY = new TLine(0,0,9,9);
  XeqY->SetLineColor(kRed);
  XeqY->Draw();
}
/////////////////////////////////////
void Draw_2DGammaGamma(){
  TCanvas *cEgEg = new TCanvas("cEgEg","cEgEg",1000,1000);
  chain->Draw("AddBack_EDC:AddBack_EDC2>>EgEg(590,0.05,3,590,0.05,3)","","colz");
  TH1F* EgEg = (TH1F*) gDirectory->Get("EgEg");
  EgEg->SetTitle("Egamma-Egamma (25Jun Full run)");
  EgEg->GetXaxis()->SetTitle("Eg [MeV]");
  EgEg->GetYaxis()->SetTitle("Eg [MeV]");
  //TLine *XeqY = new TLine(0,0,3,3);
  //XeqY->SetLineColor(kRed);
  //XeqY->Draw();
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void GateGamma_SeeParticle(double gamma, double width){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  chain->Draw("Ex>>ExGate(220,-1,10)",gating.c_str(),"colz");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  //ExGate->SetTitle("Ex gated on (using 09Jun Full run)");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
}
///////////////////////////////////
void GateGamma_SeeParticle_WithBG(double gamma, double width, double bg){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  chain->Draw("Ex>>ExGate(220,-1,10)",gating.c_str(),"");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ExGate->SetLineColor(kGreen);
  ExGate->SetFillColor(kGreen);
  ExGate->SetFillStyle(3154);

  chain->Draw("Ex>>ExBG(220,-1,10)",bggate.c_str(),"same");
  TH1F* ExBG = (TH1F*) gDirectory->Get("ExBG");
  ExBG->SetLineColor(kRed);
  ExBG->SetFillColor(kRed);
  ExBG->SetFillStyle(3345);
}
///////////////////////////////////
void GateParticle_SeeGamma(double particle, double width){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-" 
      + to_string(particle)
      + ")<"
      + to_string(width);

  TCanvas *cEg_Gate = new TCanvas("cEg_Gate","cEg_Gate",1000,1000);
  chain->Draw("AddBack_EDC>>EgGate(1000,0,10)",gating.c_str(),"colz");
  TH1F* EgGate = (TH1F*) gDirectory->Get("EgGate");
  //ExGate->SetTitle("Ex gated on (using 09Jun Full run)");
  EgGate->GetXaxis()->SetTitle("Eg [MeV]");
  EgGate->GetYaxis()->SetTitle("Counts / 10 keV");
}
///////////////////////////////////
void GateParticle_SeeGamma_WithBG(double particle, double width, double bg){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-" 
      + to_string(particle)
      + ")<"
      + to_string(width);
  string bggate = "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

  TCanvas *cEg_Gate = new TCanvas("cEg_Gate","cEg_Gate",1000,1000);
  chain->Draw("AddBack_EDC>>EgGate(1000,0,10)",gating.c_str(),"");
  TH1F* EgGate = (TH1F*) gDirectory->Get("EgGate");
  EgGate->GetXaxis()->SetTitle("Eg [MeV]");
  EgGate->GetYaxis()->SetTitle("Counts / 10 keV");
  EgGate->SetLineColor(kGreen);
  EgGate->SetFillColor(kGreen);
  EgGate->SetFillStyle(3154);

  chain->Draw("AddBack_EDC>>EgBG(1000,0,10)",bggate.c_str(),"same");
  TH1F* EgBG = (TH1F*) gDirectory->Get("EgBG");
  EgBG->SetLineColor(kRed);
  EgBG->SetFillColor(kRed);
  EgBG->SetFillStyle(3345);
}
///////////////////////////////////
void GateGamma_SeeGamma(double gamma, double width){
  string gating = "abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  TCanvas *cEx_Gate = new TCanvas("cggGate","cggGate",1000,1000);
  chain->Draw("AddBack_EDC>>ggGate(590,0.05,3)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  //ExGate->SetTitle("Ex gated on (using 09Jun Full run)");
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
}
///////////////////////////////////
void GateGamma_SeeGamma_WithBG(double gamma, double width, double bg){
  string gating = "abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = "abs(AddBack_EDC2-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

  TCanvas *cggGate = new TCanvas("cggGate","cggGate",1000,1000);
  chain->Draw("AddBack_EDC>>ggGate(590,0.05,3)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ggGate->SetLineColor(kGreen);
  ggGate->SetFillColor(kGreen);
  ggGate->SetFillStyle(3154);

  chain->Draw("AddBack_EDC>>ggBG(590,0.05,3)",bggate.c_str(),"same");
  TH1F* ggBG = (TH1F*) gDirectory->Get("ggBG");
  ggBG->SetLineColor(kRed);
  ggBG->SetFillColor(kRed);
  ggBG->SetFillStyle(3345);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CompareExsAt4MeV(){
  TCanvas *cExCompare = new TCanvas("cExCompare","cExCompare",1000,1000);
  chain->Draw("AddBack_EDC>>gate3p0(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-3.0)<0.1","same");
  chain->Draw("AddBack_EDC>>gate3p5(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-3.5)<0.1","same");
  chain->Draw("AddBack_EDC>>gate3p9(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-3.9)<0.1","same");
  chain->Draw("AddBack_EDC>>gate4p3(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-4.3)<0.1","same");
 
  TH1F* gate3p0 = (TH1F*) gDirectory->Get("gate3p0");
    gate3p0->GetXaxis()->SetTitle("Egamma [MeV]");
    gate3p0->GetYaxis()->SetTitle("Counts");
    gate3p0->SetLineColor(kRed);
  
  TH1F* gate3p5 = (TH1F*) gDirectory->Get("gate3p5");
    gate3p5->SetLineColor(kBlue);

  TH1F* gate3p9 = (TH1F*) gDirectory->Get("gate3p9");
    gate3p9->SetLineColor(kGreen);

  TH1F* gate4p3 = (TH1F*) gDirectory->Get("gate4p3");
    gate4p3->SetLineColor(kViolet);
    gate4p3->GetXaxis()->SetRangeUser(0.,5.);

  cout << " 3.0 - Red\n"
       << " 3.5 - Blue\n"
       << " 3.9 - Green\n"
       << " 4.3 - Violet"
  << endl;
}
///////////////////////////////////
void CompareSimExp(){
  TCanvas *cSimExp = new TCanvas("cSimExp","cSimExp",1000,1000);
  gStyle->SetOptStat(0);
  
  chain->Draw("Ex>>hexp(70,-1,6)","abs(T_MUGAST_VAMOS-2777)<600","");
  TH1F* hexp = (TH1F*) gDirectory->Get("hexp");
  hexp->SetTitle("Comparing Simulation to Experiment");
  hexp->GetXaxis()->SetTitle("Ex [MeV]");
  hexp->GetYaxis()->SetTitle("Counts / 0.1 MeV");
  hexp->SetLineColor(kRed);

  TFile* simfile = new TFile("../../../Outputs/Analysis/SimTest_Jun22_TWOFNR_p32.root", 
		  "READ");
  TTree* simtree = (TTree*) simfile->FindObjectAny("PhysicsTree");
  simtree->Draw("Ex>>hsimMGp32(70,-1,6)",
		  "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8","same");
  TH1F* hsimMGp32 = (TH1F*) gDirectory->Get("hsimMGp32");
  hsimMGp32->SetLineColor(kBlue);
  //simtree->Draw("Ex>>hsimALLp32(70,-1,6)","","same");
  //  TH1F* hsimALLp32 = (TH1F*) gDirectory->Get("hsimALLp32");
  //  hsimALLp32->SetLineColor(kBlue);
  //  hsimALLp32->SetLineStyle(kDashed);

  TFile* simfile3 = new TFile("../../../Outputs/Analysis/SimTest_Jun22_TWOFNR_p12.root",
		  "READ");
  TTree* simtree3 = (TTree*) simfile3->FindObjectAny("PhysicsTree");
  simtree3->Draw("Ex>>hsimMGp12(70,-1,6)",
		  "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8","same");
  TH1F* hsimMGp12 = (TH1F*) gDirectory->Get("hsimMGp12");
  hsimMGp12->SetLineColor(kOrange);
  //simtree3->Draw("Ex>>hsimALLp12(70,-1,6)","","same");
  //  TH1F* hsimALLp12 = (TH1F*) gDirectory->Get("hsimALLp12");
  //  hsimALLp12->SetLineColor(kOrange);
  //  hsimALLp12->SetLineStyle(kDashed);

  TFile* simfile2 = new TFile("../../../Outputs/Analysis/SimTest_Jun22_TWOFNR_f72.root",
		  "READ");
  TTree* simtree2 = (TTree*) simfile2->FindObjectAny("PhysicsTree");
  simtree2->Draw("Ex>>hsimMGf72(70,-1,6)",
		  "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8","same");
  TH1F* hsimMGf72 = (TH1F*) gDirectory->Get("hsimMGf72");
  hsimMGf72->SetLineColor(kGreen);
  //simtree2->Draw("Ex>>hsimALLf72(70,-1,6)","","same");
  //  TH1F* hsimALLf72 = (TH1F*) gDirectory->Get("hsimALLf72");
  //  hsimALLf72->SetLineColor(kGreen);
  //  hsimALLf72->SetLineStyle(kDashed);

  TFile* simfile4 = new TFile("../../../Outputs/Analysis/SimTest_Jun22_TWOFNR_f52.root", 
		  "READ");
  TTree* simtree4 = (TTree*) simfile4->FindObjectAny("PhysicsTree");
  simtree4->Draw("Ex>>hsimMGf52(70,-1,6)",
		  "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8","same");
  TH1F* hsimMGf52 = (TH1F*) gDirectory->Get("hsimMGf52");
  hsimMGf52->SetLineColor(kViolet);
  //simtree4->Draw("Ex>>hsimALLf72(70,-1,6)","","same");
  //  TH1F* hsimALLf52 = (TH1F*) gDirectory->Get("hsimALLf52");
  //  hsimALLf52->SetLineColor(kViolet);
  //  hsimALLf52->SetLineStyle(kDashed);


  auto legend = new TLegend(0.7,0.8,0.9,0.9);
  legend->AddEntry(hexp,       "Experiment","l");
  legend->AddEntry(hsimMGp32,  "Simulation p3/2, MG","l");
  legend->AddEntry(hsimMGp12,  "Simulation p1/2, MG","l");
  //legend->AddEntry(hsimALLp32, "Simulation p3/2, MG+MM","l");
  //legend->AddEntry(hsimALLp12, "Simulation p1/2, MG+MM","l");
  legend->AddEntry(hsimMGf72,  "Simulation f7/2, MG","l");
  legend->AddEntry(hsimMGf52,  "Simulation f5/2, MG","l");
  //legend->AddEntry(hsimALLf72, "Simulation f7/2, MG+MM","l");
  //legend->AddEntry(hsimALLf52, "Simulation f5/2, MG+MM","l");
  legend->Draw();
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
