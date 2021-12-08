#include "NPReaction.h"
#include <string>
#include <sstream>

using namespace std;

TChain* chain=NULL ;
char cond[1000];

NPL::Reaction Kdp("47K(d,p)48K@355");
NPL::Reaction Kdt("47K(d,t)46K@355");
NPL::Reaction Kdd("47K(d,d)47K@355");
NPL::Reaction Kpp("47K(p,p)47K@355");
NPL::Reaction K12C12C("47K(12C,12C)47K@355");
NPL::Reaction Tidp("47Ti(d,p)48Ti@355");
NPL::Reaction Tidt("47Ti(d,t)46Ti@355");
NPL::Reaction Tidd("47Ti(d,d)47Ti@355");
NPL::Reaction Ti12C12C("47Ti(12C,12C)47Ti@355");

void KnownLines_Ex(bool isVertical, double rangemin, double rangemax, Style_t lType, Color_t lColour);

/* BASE FUNCTIONS */

TChain* Chain(std::string TreeName, std::vector<std::string>& file, bool EventList){
  TChain*  chain = new TChain(TreeName.c_str());
  unsigned int size =file.size(); 
    for(unsigned int i = 0 ; i < size ; i++){
      cout << "Adding " << file[i] << endl;
      chain->Add(file[i].c_str());
    } 
  return chain;
}

void LoadChainNP(){
  vector<string> files;
  
  //files.push_back("../../../Outputs/Analysis/OriginalValues_ptI.root");
  //files.push_back("../../../Outputs/Analysis/OriginalValues_ptII.root");
  //files.push_back("../../../Outputs/Analysis/OriginalValues_ptIII.root");

  //files.push_back("../../../Outputs/Analysis/47K_Full_07Sep_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47K_Full_07Sep_PartII.root");

  //files.push_back("../../../Outputs/Analysis/47Kdp_11Oct_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_11Oct_PartII.root");

  //files.push_back("../../../Outputs/Analysis/47Kdp_13Oct_wildTry_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_13Oct_wildTry_PartII.root");

  files.push_back("../../../Outputs/Analysis/47Kdp_08Nov_PartI.root");
  files.push_back("../../../Outputs/Analysis/47Kdp_08Nov_PartII.root");

  chain = Chain("PhysicsTree",files,true);
}

void DrawParticleStates(TCanvas* canvas){
  canvas->cd();
  
  canvas->Update();
  double max = canvas->GetUymax();
  TLine *Sn = new TLine(4.644, 0.0, 4.644, max);
    Sn->SetLineColor(kRed);
    Sn->SetLineStyle(7);
    Sn->Draw();
  TLine *gs = new TLine(0.000, 0.0, 0.000, max);
    gs->SetLineColor(kGreen);
    gs->SetLineStyle(7);
    gs->Draw();
  TLine *l0143 = new TLine(0.143, 0.0, 0.143, max);
    l0143->SetLineStyle(kDashed);
    l0143->Draw();
  TLine *l0728 = new TLine(0.728, 0.0, 0.728, max);
    l0728->SetLineStyle(kDotted);
    l0728->Draw();
  TLine *l0968 = new TLine(0.968, 0.0, 0.968, max);
    l0968->SetLineStyle(kDotted);
    l0968->Draw();
  TLine *l1410 = new TLine(1.410, 0.0, 1.410, max);
    l1410->SetLineStyle(kDotted);
    l1410->Draw();
  TLine *l1981 = new TLine(1.981, 0.0, 1.981, max);
    l1981->SetLineStyle(kDotted);
    l1981->Draw("same");
  TLine *l2410 = new TLine(2.410, 0.0, 2.410, max);
    l2410->SetLineStyle(kDotted);
    l2410->Draw("same");
  TLine *l2907 = new TLine(2.907, 0.0, 2.907, max);
    l2907->SetLineStyle(kDotted);
    l2907->Draw("same");
  TLine *l3600 = new TLine(3.600, 0.0, 3.600, max);
    l3600->SetLineStyle(kDotted);
    l3600->Draw("same");
  TLine *l3800 = new TLine(3.792, 0.0, 3.792, max);
    l3800->SetLineStyle(kDotted);
    l3800->Draw("same");
}

void plot_kine(NPL::Reaction r, double Ex,Color_t c,int w, int s){
  r.SetExcitation4(Ex);
  TGraph* g= r.GetKinematicLine3();
  g->SetLineColor(c) ;
  g->SetLineStyle(s) ;
  g->SetLineWidth(w) ;
  g->Draw("c");
}

void plot_state(double Ex,double max,Color_t c,int w, int s){
  TLine* line = new TLine(Ex,0,Ex,max) ; 
  line->SetLineColor(c) ;
  line->SetLineStyle(s) ;
  line->SetLineWidth(w) ;
  line->Draw();
}

void AddTiStates(double E){
 NPL::Reaction Ti("47Ti(d,p)48Ti@355");
 Ti.SetExcitationHeavy(E);
 auto g = Ti.GetKinematicLine3();
 g->SetLineWidth(1);
 g->SetLineStyle(2);
 g->Draw("c");
}

/* DRAWING FUNCTIONS */

void Draw_1DGamma(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)","abs(T_MUGAST_VAMOS-2777)<600");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
}

void Draw_1DParticle(){
  TCanvas *cEx = new TCanvas("cEx","cEx",1000,1000);
  chain->Draw("Ex>>Ep(120,-1,5)",
	"abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
	"");
  TH1F* Ep = (TH1F*) gDirectory->Get("Ep");
//  Ep->SetTitle("Ex");
  Ep->GetXaxis()->SetTitle("Ex [MeV]");
  Ep->GetYaxis()->SetTitle("Counts / 0.05 MeV");

  DrawParticleStates(cEx);
}

void Draw_2DParticleGamma(){
  TCanvas *cExEg = new TCanvas("cExEg","cExEg",1000,1000);
  chain->Draw("AddBack_EDC:Ex>>ExEg(140,-1,7,2500,0,5)",
        "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
	"colz");
  TH1F* ExEg = (TH1F*) gDirectory->Get("ExEg");
  ExEg->SetTitle("Ex-Egamma");
  ExEg->GetXaxis()->SetTitle("Ex [MeV]");
  ExEg->GetYaxis()->SetTitle("Eg [MeV]");
  TLine *XeqY = new TLine(0,0,9,9);
  XeqY->SetLineColor(kRed);
  XeqY->Draw();
}

void Draw_2DGammaGamma(){
  TCanvas *cEgEg = new TCanvas("cEgEg","cEgEg",1000,1000);
  chain->Draw("AddBack_EDC:AddBack_EDC2>>EgEg(999,0.005,5,999,0.005,5)","","colz");
  TH1F* EgEg = (TH1F*) gDirectory->Get("EgEg");
  chain->Draw("AddBack_EDC2:AddBack_EDC>>EgEg2(999,0.005,5,999,0.005,5)","","colz");
  TH1F* EgEg2 = (TH1F*) gDirectory->Get("EgEg2");
  EgEg->Add(EgEg2,1);
  EgEg->SetTitle("Egamma-Egamma");
  EgEg->GetXaxis()->SetTitle("Eg [MeV]");
  EgEg->GetYaxis()->SetTitle("Eg [MeV]");
  EgEg->Draw("colz");
  TLine *XeqY = new TLine(0,0,5,5);
  XeqY->SetLineColor(kRed);
  XeqY->SetLineStyle(kDashed);
  XeqY->Draw("same");
}

void GateGamma_SeeParticle(double gamma, double width, double binsize){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" 
  //string gating = "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber!=3 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
/*  
cout << " NO MG3 IN THIS ONE!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
cout << " NO MG3 IN THIS ONE!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
cout << " NO MG3 IN THIS ONE!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
cout << " NO MG3 IN THIS ONE!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
cout << " NO MG3 IN THIS ONE!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
cout << " NO MG3 IN THIS ONE!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
*/

  string title = to_string(gamma-width)+" < Eg < "+to_string(gamma+width);
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>ExGate(" + to_string(8.0/binsize) + ",-1,7)";

  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  //chain->Draw("Ex>>ExGate(60,-1,5)",gating.c_str(),"colz");
  chain->Draw(draw.c_str(),gating.c_str(),"colz");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle(ytitle.c_str());
  ExGate->SetTitle(title.c_str());

  DrawParticleStates(cEx_Gate);
}

void GateGamma_SeeParticle_WithBG(double gamma, double width, double bg){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

  string title = "Gate: "+to_string(gamma-width)+" to "+to_string(gamma+width)+"."
	  + "  BG: "+to_string(bg-width)+" to "+to_string(bg+width)+".";
  
  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  chain->Draw("Ex>>ExGate(80,-1,7)",gating.c_str(),"");
  //chain->Draw("Ex>>ExGate(120,-1,5)",gating.c_str(),"");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle("Counts / 0.10 MeV");
  //ExGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ExGate->SetLineColor(kGreen);
  ExGate->SetFillColor(kGreen);
  ExGate->SetFillStyle(3154);
  ExGate->SetTitle(title.c_str());

  chain->Draw("Ex>>ExBG(80,-1,7)",bggate.c_str(),"same");
  //chain->Draw("Ex>>ExBG(120,-1,5)",bggate.c_str(),"same");
  TH1F* ExBG = (TH1F*) gDirectory->Get("ExBG");
  ExBG->SetLineColor(kRed);
  ExBG->SetFillColor(kRed);
  ExBG->SetFillStyle(3345);

  DrawParticleStates(cEx_Gate);
}

void GateParticle_SeeGamma(double particle, double width){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber<8 && Mugast.TelescopeNumber>0 && abs(Ex-" 
      + to_string(particle)
      + ")<"
      + to_string(width);

  string title = to_string(particle-width)+" < Ex < "+to_string(particle+width);
  
  TCanvas *cEg_Gate = new TCanvas("cEg_Gate","cEg_Gate",1000,1000);
  chain->Draw("AddBack_EDC>>EgGate(1000,0,10)",gating.c_str(),"colz");
  TH1F* EgGate = (TH1F*) gDirectory->Get("EgGate");
  //ExGate->SetTitle("Ex gated on (using 09Jun Full run)");
  EgGate->GetXaxis()->SetTitle("Eg [MeV]");
  EgGate->GetYaxis()->SetTitle("Counts / 10 keV");
  EgGate->SetTitle(title.c_str());
  EgGate->GetXaxis()->SetRangeUser(0.0, particle);
  
  cEg_Gate->Update();
  TLine *limit = new TLine(particle, 0.0, particle, cEg_Gate->GetUymax());
  limit->SetLineColor(kRed);

  EgGate->Draw();
  limit->Draw();
}

void GateParticle_SeeGamma_WithBG(double particle, double width, double bg){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-" 
      + to_string(particle)
      + ")<"
      + to_string(width);
  string bggate = "abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

  string title = "Gate: "+to_string(particle-width)+" to "+to_string(particle+width)+"."
	  + "  BG: "+to_string(bg-width)+" to "+to_string(bg+width)+".";

  TCanvas *cEg_Gate = new TCanvas("cEg_Gate","cEg_Gate",1000,1000);
  chain->Draw("AddBack_EDC>>EgGate(1000,0,10)",gating.c_str(),"");
  TH1F* EgGate = (TH1F*) gDirectory->Get("EgGate");
  EgGate->GetXaxis()->SetTitle("Eg [MeV]");
  EgGate->GetYaxis()->SetTitle("Counts / 10 keV");
  EgGate->SetTitle(title.c_str());
  EgGate->SetLineColor(kGreen);
  EgGate->SetFillColor(kGreen);
  EgGate->SetFillStyle(3154);
  EgGate->GetXaxis()->SetRangeUser(0.0, particle);
  EgGate->Draw();

  chain->Draw("AddBack_EDC>>EgBG(1000,0,10)",bggate.c_str(),"same");
  TH1F* EgBG = (TH1F*) gDirectory->Get("EgBG");
  EgBG->SetTitle(title.c_str());
  EgBG->SetLineColor(kRed);
  EgBG->SetFillColor(kRed);
  EgBG->SetFillStyle(3345);
}
                                   
void GateGamma_SeeGamma(double gamma, double width){
  string gating = "abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string gating2 = "abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  string title = to_string(gamma-width) + " < Eg < " + to_string(gamma+width);
  
  TCanvas *cEx_Gate = new TCanvas("cggGate","cggGate",1000,1000);
  //chain->Draw("AddBack_EDC>>ggGate(990,0.05,5)",gating.c_str(),"");
  chain->Draw("AddBack_EDC>>ggGate(999,0.005,5)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.005 MeV");
  ggGate->SetTitle(title.c_str());

  //chain->Draw("AddBack_EDC2>>ggGate2(990,0.05,5)",gating2.c_str(),"");
  chain->Draw("AddBack_EDC2>>ggGate2(999,0.005,5)",gating2.c_str(),"");
  TH1F* ggGate2 = (TH1F*) gDirectory->Get("ggGate2");
  ggGate->Add(ggGate2,1);
  ggGate->Draw();
}
                                   
void GateGamma_SeeGamma_WithBG(double gamma, double width, double bg){
/**/
  string gating = "abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = "abs(AddBack_EDC2-" 
      + to_string(bg)
      + ")<"
      + to_string(width);
  string gating2 = "abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate2 = "abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

TCanvas *cggGate = new TCanvas("cggGate","cggGate",1000,1000);
  chain->Draw("AddBack_EDC>>ggGate(999,0.005,5)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  chain->Draw("AddBack_EDC2>>ggGate2(999,0.005,5)",gating2.c_str(),"");
  TH1F* ggGate2 = (TH1F*) gDirectory->Get("ggGate2");
  ggGate->Add(ggGate2,1);
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ggGate->SetLineColor(kGreen);
  ggGate->SetFillColor(kGreen);
  ggGate->SetFillStyle(3154);

  chain->Draw("AddBack_EDC>>ggBG(999,0.005,5)",bggate.c_str(),"");
  TH1F* ggBG = (TH1F*) gDirectory->Get("ggBG");
  chain->Draw("AddBack_EDC2>>ggBG2(999,0.005,5)",bggate2.c_str(),"");
  TH1F* ggBG2 = (TH1F*) gDirectory->Get("ggBG2");
  ggBG->Add(ggBG2,1);
  ggBG->SetLineColor(kRed);
  ggBG->SetFillColor(kRed);
  ggBG->SetFillStyle(3345);

  ggGate->Draw();
  ggBG->Draw("same");
/**/

/*
  string gating = "abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string gating2 = "abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  string title = to_string(gamma-width)+" < Eg < "+to_string(gamma+width);
  
  TCanvas *cEx_Gate = new TCanvas("cggGate","cggGate",1000,1000);
  //chain->Draw("AddBack_EDC>>ggGate(990,0.05,5)",gating.c_str(),"");
  chain->Draw("AddBack_EDC>>ggGate(999,0.005,5)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.005 MeV");
  ggGate->SetTitle(title.c_str());

  //chain->Draw("AddBack_EDC2>>ggGate2(990,0.05,5)",gating2.c_str(),"");
  chain->Draw("AddBack_EDC2>>ggGate2(999,0.005,5)",gating2.c_str(),"");
  TH1F* ggGate2 = (TH1F*) gDirectory->Get("ggGate2");
  ggGate->Add(ggGate2,1);
  ggGate->SetLineColor(kGreen);
  ggGate->SetFillColor(kGreen);
  ggGate->SetFillStyle(3154);
//  ggGate->Draw();

  string bggating = "abs(AddBack_EDC2-" 
      + to_string(bg)
      + ")<"
      + to_string(width);
  string bggating2 = "abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

  //string title = to_string(bg-width)+" < Eg < "+to_string(gamma+width);
  
  //TCanvas *cEx_Gate = new TCanvas("cggGate","cggGate",1000,1000);
  //chain->Draw("AddBack_EDC>>ggGate(990,0.05,5)",gating.c_str(),"");
  chain->Draw("AddBack_EDC>>ggbgGate(999,0.005,5)",bggating.c_str(),"");
  TH1F* ggbgGate = (TH1F*) gDirectory->Get("ggbgGate");
  //ggbgGate->GetXaxis()->SetTitle("Eg [MeV]");
  //ggbgGate->GetYaxis()->SetTitle("Counts / 0.005 MeV");
  //ggbgGate->SetTitle(title.c_str());

  //chain->Draw("AddBack_EDC2>>ggGate2(990,0.05,5)",gating2.c_str(),"");
  chain->Draw("AddBack_EDC2>>ggbgGate2(999,0.005,5)",bggating2.c_str(),"");
  TH1F* ggbgGate2 = (TH1F*) gDirectory->Get("ggbgGate2");
  ggbgGate->Add(ggbgGate2,1);
  ggbgGate->SetLineColor(kRed);
  ggbgGate->SetFillColor(kRed);
  ggbgGate->SetFillStyle(3345);

  ggGate->Draw();
  ggbgGate->Draw("same");

*/

}

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

void CompareSimExp(){
  TCanvas *cSimExp = new TCanvas("cSimExp","cSimExp",1000,1000);
  gStyle->SetOptStat(0);
  
  chain->Draw("Ex>>hexp(70,-1,6)","abs(T_MUGAST_VAMOS-2777)<600","");
  TH1F* hexp = (TH1F*) gDirectory->Get("hexp");
  hexp->SetTitle("Comparing Simulation to Experiment");
  hexp->GetXaxis()->SetTitle("Ex [MeV]");
  hexp->GetYaxis()->SetTitle("Counts / 0.1 MeV");
  hexp->SetLineColor(kRed);

  TFile* simfile = new TFile("../../../Outputs/Analysis/SimTest_Jun28_p32_ExHeavy0p143.root", 
  //TFile* simfile = new TFile("../../../Outputs/Analysis/SimTest_Jun22_TWOFNR_p32.root", 
		  "READ");
  TTree* simtree = (TTree*) simfile->FindObjectAny("PhysicsTree");
  simtree->Draw("Ex>>hsimMGp32(70,-1,6)",
		  "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8","same");
  TH1F* hsimMGp32 = (TH1F*) gDirectory->Get("hsimMGp32");
  hsimMGp32->SetLineColor(kBlue);
  hsimMGp32->SetFillColor(kBlue);
  hsimMGp32->SetFillStyle(3345);
  //simtree->Draw("Ex>>hsimALLp32(70,-1,6)","","same");
  //  TH1F* hsimALLp32 = (TH1F*) gDirectory->Get("hsimALLp32");
  //  hsimALLp32->SetLineColor(kBlue);
  //  hsimALLp32->SetLineStyle(kDashed);
/*
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
*/
}

void MugastMisalignment(){
  TCanvas *cMisaligned = new TCanvas("cMisaligned","cMisaligned",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("Ex>>hcG_T1(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==1","");
  chain->Draw("Ex>>hcG_T2(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==2","same");
  chain->Draw("Ex>>hcG_T3(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==3","same");
  chain->Draw("Ex>>hcG_T4(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==4","same");
  chain->Draw("Ex>>hcG_T5(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","same");
  chain->Draw("Ex>>hcG_T7(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==7","same");
  TH1F* hcG_T1 = (TH1F*) gDirectory->Get("hcG_T1");
  hcG_T1->SetTitle("Misalignment of MUGAST telescopes");
  hcG_T1->GetXaxis()->SetTitle("E_{x} [MeV]");
  hcG_T1->GetYaxis()->SetTitle("Counts / (50 keV)");
  hcG_T1->SetLineColor(kRed);
  hcG_T1->SetFillStyle(3244);
  hcG_T1->SetFillColor(kRed);
  hcG_T1->GetYaxis()->SetRangeUser(0.,500.);
  TH1F* hcG_T2 = (TH1F*) gDirectory->Get("hcG_T2");
  hcG_T2->SetLineColor(kOrange);
  hcG_T2->SetFillStyle(3344);
  hcG_T2->SetFillColor(kOrange);
  TH1F* hcG_T3 = (TH1F*) gDirectory->Get("hcG_T3");
  hcG_T3->SetLineColor(kGreen);
  hcG_T3->SetFillStyle(3444);
  hcG_T3->SetFillColor(kGreen);
  TH1F* hcG_T4 = (TH1F*) gDirectory->Get("hcG_T4");
  hcG_T4->SetLineColor(kTeal);
  hcG_T4->SetFillStyle(3544);
  hcG_T4->SetFillColor(kTeal);
  TH1F* hcG_T5 = (TH1F*) gDirectory->Get("hcG_T5");
  hcG_T5->SetLineColor(kBlue);
  hcG_T5->SetFillStyle(3644);
  hcG_T5->SetFillColor(kBlue);
  TH1F* hcG_T7 = (TH1F*) gDirectory->Get("hcG_T7");
  hcG_T7->SetLineColor(kViolet);
  hcG_T7->SetFillStyle(3644);
  hcG_T7->SetFillColor(kViolet);
//  hcG_T7->GetYaxis()->SetRangeUser(0.,750.);

  TLine *line = new TLine(0,0,0,500);
  line->SetLineColor(kBlack);
  line->Draw();

  auto legend = new TLegend(0.1,0.7,0.35,0.9);
  legend->AddEntry(hcG_T1,"MUGAST 1","f");
  legend->AddEntry(hcG_T2,"MUGAST 2","f");
  legend->AddEntry(hcG_T3,"MUGAST 3","f");
  legend->AddEntry(hcG_T4,"MUGAST 4","f");
  legend->AddEntry(hcG_T5,"MUGAST 5","f");
  legend->AddEntry(hcG_T7,"MUGAST 7","f");
  legend->Draw();

  DrawParticleStates(cMisaligned);
}

void MugastMisalignment(double gamma, double width){
  TCanvas *cMisaligned = new TCanvas("cMisaligned","cMisaligned",1000,1000);
  gStyle->SetOptStat(0);

  string base = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" + to_string(gamma) 
	      + ")<" + to_string(width) + " && Mugast.TelescopeNumber==";
  string str1 = base + "1";
  string str2 = base + "2";
  string str3 = base + "3";
  string str4 = base + "4";
  string str5 = base + "5";
  string str7 = base + "7";

  chain->Draw("Ex>>hcG_T1(120,-1,5)", str1.c_str(), "");
  chain->Draw("Ex>>hcG_T2(120,-1,5)", str2.c_str(), "same");
  chain->Draw("Ex>>hcG_T3(120,-1,5)", str3.c_str(), "same");
  chain->Draw("Ex>>hcG_T4(120,-1,5)", str4.c_str(), "same");
  chain->Draw("Ex>>hcG_T5(120,-1,5)", str5.c_str(), "same");
  chain->Draw("Ex>>hcG_T7(120,-1,5)", str7.c_str(), "same");
  TH1F* hcG_T1 = (TH1F*) gDirectory->Get("hcG_T1");
  hcG_T1->SetTitle("Misalignment of MUGAST telescopes, gated on gamma");
  hcG_T1->GetXaxis()->SetTitle("E_{x} [MeV]");
  hcG_T1->GetYaxis()->SetTitle("Counts / (50 keV)");
  hcG_T1->SetLineColor(kRed);
  //hcG_T1->SetFillStyle(3244);
  //hcG_T1->SetFillColor(kRed);
  hcG_T1->GetYaxis()->SetRangeUser(0.,500.);
  TH1F* hcG_T2 = (TH1F*) gDirectory->Get("hcG_T2");
  hcG_T2->SetLineColor(kOrange);
  //hcG_T2->SetFillStyle(3344);
  //hcG_T2->SetFillColor(kOrange);
  TH1F* hcG_T3 = (TH1F*) gDirectory->Get("hcG_T3");
  hcG_T3->SetLineColor(kGreen);
  //hcG_T3->SetFillStyle(3444);
  //hcG_T3->SetFillColor(kGreen);
  TH1F* hcG_T4 = (TH1F*) gDirectory->Get("hcG_T4");
  hcG_T4->SetLineColor(kTeal);
  //hcG_T4->SetFillStyle(3544);
  //hcG_T4->SetFillColor(kTeal);
  TH1F* hcG_T5 = (TH1F*) gDirectory->Get("hcG_T5");
  hcG_T5->SetLineColor(kBlue);
  //hcG_T5->SetFillStyle(3644);
  //hcG_T5->SetFillColor(kBlue);
  TH1F* hcG_T7 = (TH1F*) gDirectory->Get("hcG_T7");
  hcG_T7->SetLineColor(kViolet);
  //hcG_T7->SetFillStyle(3644);
  //hcG_T7->SetFillColor(kViolet);
//  hcG_T7->GetYaxis()->SetRangeUser(0.,750.);

  TLine *line = new TLine(0,0,0,500);
  line->SetLineColor(kBlack);
  line->Draw();

  auto legend = new TLegend(0.1,0.7,0.35,0.9);
  legend->AddEntry(hcG_T1,"MUGAST 1","f");
  legend->AddEntry(hcG_T2,"MUGAST 2","f");
  legend->AddEntry(hcG_T3,"MUGAST 3","f");
  legend->AddEntry(hcG_T4,"MUGAST 4","f");
  legend->AddEntry(hcG_T5,"MUGAST 5","f");
  legend->AddEntry(hcG_T7,"MUGAST 7","f");
  legend->Draw();

  DrawParticleStates(cMisaligned);
}

void ExPhiLab(){
  TCanvas *diagnosePhi = new TCanvas("diagnosePhi","diagnosePhi",1000,1000);
  chain->Draw(
    "Ex:PhiLab>>phiHist(180,-180,180,120,-1,5)", 
    "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    "colz");
  TH1F* phiHist = (TH1F*) gDirectory->Get("phiHist");  
  phiHist->GetXaxis()->SetTitle("Phi (degrees)");
  phiHist->GetYaxis()->SetTitle("Ex [MeV]");
  phiHist->SetTitle("Phi dependance testing");

  diagnosePhi->Update();
  TLine *l0143 = new TLine(-180., 0.143, 180., 0.143);
    l0143->SetLineStyle(kDashed);
    l0143->SetLineColor(kRed);
    l0143->Draw();
  TLine *l0968 = new TLine(-180., 0.968, 180., 0.968);
    l0968->SetLineStyle(kDotted);
    l0968->SetLineColor(kRed);
    l0968->Draw();
  TLine *l1410 = new TLine(-180., 1.410, 180., 1.410);
    l1410->SetLineStyle(kDotted);
    l1410->SetLineColor(kRed);
    l1410->Draw();
  TLine *l1981 = new TLine(-180., 1.981, 180., 1.981);
    l1981->SetLineStyle(kDotted);
    l1981->SetLineColor(kRed);
    l1981->Draw("same");
  TLine *l2410 = new TLine(-180., 2.410, 180., 2.410);
    l2410->SetLineStyle(kDotted);
    l2410->SetLineColor(kRed);
    l2410->Draw("same");
  TLine *l2907 = new TLine(-180., 2.907, 180., 2.907);
    l2907->SetLineStyle(kDotted);
    l2907->SetLineColor(kRed);
    l2907->Draw("same");
  TLine *l3600 = new TLine(-180., 3.600, 180., 3.600);
    l3600->SetLineStyle(kDotted);
    l3600->SetLineColor(kRed);
    l3600->Draw("same");
}

void ExPhiLab_ForPoster(){
  TCanvas *diagnosePhi = new TCanvas("diagnosePhi","diagnosePhi",1000,1000);
  diagnosePhi->Divide(1,2);
  diagnosePhi->cd(1);
  chain->Draw(
    "Ex:PhiLab>>phiHist(180,-180,180,40,-1,+1)", 
    "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    "colz");
  TH1F* phiHist = (TH1F*) gDirectory->Get("phiHist");  
  phiHist->GetXaxis()->SetTitle("#Phi_{Lab} [deg]");
  phiHist->GetYaxis()->SetTitle("Ex [MeV]");

  diagnosePhi->Update();
  TLine *l0143 = new TLine(-180., 0.143, 180., 0.143);
    l0143->SetLineColor(kRed);
    l0143->Draw();

  TFile* file = new TFile("../../../Outputs/Analysis/47K_Full_22July.root", "READ");
  TTree* tree = (TTree*) file->FindObjectAny("PhysicsTree");
  diagnosePhi->cd(2);
  tree->Draw("Ex:PhiLab>>phiHist2(180,-180,180,40,-1,+1)", 
    "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    "colz");
  TH1F* phiHist2 = (TH1F*) gDirectory->Get("phiHist2");  
  
  diagnosePhi->Update();
  TLine *l0143_2 = new TLine(-180., 0.143, 180., 0.143);
    l0143_2->SetLineColor(kRed);
    l0143_2->Draw();
}

void ExThetaLab(){
  TCanvas *diagnoseTheta = new TCanvas("diagnoseTheta","diagnoseTheta",1000,1000);
  chain->Draw(
    "Ex:ThetaLab>>thetaHist(60,100,160,120,-1,5)", 
    "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    "colz");
  TH1F* thetaHist = (TH1F*) gDirectory->Get("thetaHist");  
  thetaHist->GetXaxis()->SetTitle("Theta (degrees)");
  thetaHist->GetYaxis()->SetTitle("Ex [MeV]");
  thetaHist->SetTitle("Theta dependance testing");
  
  diagnoseTheta->Update();
  TLine *l0000 = new TLine(100., 0.000, 160., 0.000);
    l0000->SetLineStyle(kDashed);
    l0000->SetLineColor(kRed);
    l0000->Draw();
  TLine *l0143 = new TLine(100., 0.143, 160., 0.143);
    l0143->SetLineStyle(kDashed);
    l0143->SetLineColor(kRed);
    l0143->Draw();
  TLine *l0968 = new TLine(100., 0.968, 160., 0.968);
    l0968->SetLineStyle(kDotted);
    l0968->SetLineColor(kRed);
    l0968->Draw();
  TLine *l1410 = new TLine(100., 1.410, 160., 1.410);
    l1410->SetLineStyle(kDotted);
    l1410->SetLineColor(kRed);
    l1410->Draw();
  TLine *l1981 = new TLine(100., 1.981, 160., 1.981);
    l1981->SetLineStyle(kDotted);
    l1981->SetLineColor(kRed);
    l1981->Draw("same");
  TLine *l2410 = new TLine(100., 2.410, 160., 2.410);
    l2410->SetLineStyle(kDotted);
    l2410->SetLineColor(kRed);
    l2410->Draw("same");
  TLine *l2907 = new TLine(100., 2.907, 160., 2.907);
    l2907->SetLineStyle(kDotted);
    l2907->SetLineColor(kRed);
    l2907->Draw("same");
  TLine *l3600 = new TLine(100., 3.600, 160., 3.600);
    l3600->SetLineStyle(kDotted);
    l3600->SetLineColor(kRed);
    l3600->Draw("same");
}

void ExThetaLab(double gamma, double width){
  TCanvas *diagnoseTheta = new TCanvas("diagnoseTheta","diagnoseTheta",1000,1000);

  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && abs(AddBack_EDC-"
	        + to_string(gamma) + ") < " + to_string(width); 

  chain->Draw("Ex:ThetaLab>>thetaHist(60,100,160,60,-1,5)", gating.c_str(), "colz");
  TH1F* thetaHist = (TH1F*) gDirectory->Get("thetaHist");  
  thetaHist->GetXaxis()->SetTitle("Theta (degrees)");
  thetaHist->GetYaxis()->SetTitle("Ex [MeV]");
  thetaHist->SetTitle("Theta dependance testing w/ gamma gating");
  
  diagnoseTheta->Update();
  //TLine *line = new TLine(100., gamma, 160., gamma);
  //  line->SetLineStyle(kDashed);
  //  line->SetLineColor(kGreen);
  //  line->Draw();
  TLine *l0000 = new TLine(100., 0.000, 160., 0.000);
    l0000->SetLineStyle(kDashed);
    l0000->SetLineColor(kRed);
    l0000->Draw();
  TLine *l0143 = new TLine(100., 0.143, 160., 0.143);
    l0143->SetLineStyle(kDashed);
    l0143->SetLineColor(kRed);
    l0143->Draw();
  TLine *l0968 = new TLine(100., 0.968, 160., 0.968);
    l0968->SetLineStyle(kDotted);
    l0968->SetLineColor(kRed);
    l0968->Draw();
  TLine *l1410 = new TLine(100., 1.410, 160., 1.410);
    l1410->SetLineStyle(kDotted);
    l1410->SetLineColor(kRed);
    l1410->Draw();
  TLine *l1981 = new TLine(100., 1.981, 160., 1.981);
    l1981->SetLineStyle(kDotted);
    l1981->SetLineColor(kRed);
    l1981->Draw("same");
  TLine *l2410 = new TLine(100., 2.410, 160., 2.410);
    l2410->SetLineStyle(kDotted);
    l2410->SetLineColor(kRed);
    l2410->Draw("same");
  TLine *l2907 = new TLine(100., 2.907, 160., 2.907);
    l2907->SetLineStyle(kDotted);
    l2907->SetLineColor(kRed);
    l2907->Draw("same");
  TLine *l3600 = new TLine(100., 3.600, 160., 3.600);
    l3600->SetLineStyle(kDotted);
    l3600->SetLineColor(kRed);
    l3600->Draw("same");

}

void ELabThetaLab(){
  TCanvas *cELabTLaab = new TCanvas("cELabTLab","cELabTLab",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("ELab:ThetaLab>>hKine(360,0,180,450,0,7)","abs(T_MUGAST_VAMOS-2777)<600","col");
  TH2F* hKine = (TH2F*) gDirectory->Get("hKine");
  hKine->SetTitle("");
  hKine->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  hKine->GetYaxis()->SetTitle("E_{lab} [MeV]");
  plot_kine(Kdp, 0.000, kBlack, 2, 1);
  plot_kine(Kdp, 4.644, kBlack, 2, 1);

  plot_kine(Kdp, 0.143, kRed, 1, 2);
  plot_kine(Kdp, 0.968, kRed, 1, 2);
  plot_kine(Kdp, 1.410, kRed, 1, 2);
  plot_kine(Kdp, 1.981, kRed, 1, 2);
  plot_kine(Kdp, 2.410, kRed, 1, 2);
  plot_kine(Kdp, 2.907, kRed, 1, 2);
  plot_kine(Kdp, 3.600, kRed, 1, 2);
  plot_kine(Kdp, 3.792, kRed, 1, 2);

  plot_kine(Kdd, 0.000, kGreen+2, 2, 9);
  plot_kine(Kpp, 0.000, kYellow, 2, 9);

  plot_kine(Tidp, 0.000, kBlack, 2, 1);
  plot_kine(Tidp, 5.652, kBlack, 2, 6); //strongest populated state according to PDBarnes(1965)
}

void XYMust2(){
  TCanvas *cXYMust2 = new TCanvas("cXYMM","cXYMM",1000,1000);
  chain->Draw("Y:X>>hXYMust2(300,-150,+150,300,-150,+150)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && MUST2.TelescopeNumber>0 && MUST2.TelescopeNumber<5",
		  "colz");
}

void XYMugast(){
  TCanvas *cXYMugast = new TCanvas("cXYMG","cXYMG",1000,1000);
  chain->Draw("Y:X>>hXYMugast(150,-150,+150, 150,-150,+150)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
		  "colz");
}

void MM5_ELabThetaLab(){
  chain->Draw("ELab:ThetaLab>>hMM5_el(180,50,95,700,0,7)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && MUST2.TelescopeNumber==5",
		  "colz");
  TH2F* hMM5_el = (TH2F*) gDirectory->Get("hMM5_el");
  hMM5_el->GetXaxis()->SetTitle("Theta Lab");
  hMM5_el->GetYaxis()->SetTitle("Ex");

  plot_kine(Kdd, 0, kGreen+2, 2, 9);
  plot_kine(Kpp, 0, kYellow, 2, 9);
}

void MM5_RawEThetaLab(){
  chain->Draw("RawEnergy:ThetaLab>>hMM5_el(90,50,95,700,0,7)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && MUST2.TelescopeNumber==5",
		  "colz");

  //plot_kine(Kdd, 0, kGreen+2, 2, 9);
  //plot_kine(Kpp, 0, kYellow, 2, 9);
}

void MM5_ExThetaLab(){
  chain->Draw("Ex:ThetaLab>>hMM5_ex(180,0,180,400,0,20)",
		  "abs(T_MUGAST_VAMOS-2777)<600 && Must2.TelescopeNumber==5",
		  "colz");
}

void ExMugast_ForPoster(){

  TCanvas *forPoster = new TCanvas("forPoster","forPoster",1000,1000);
  gStyle->SetOptStat(0);
  forPoster->Divide(1,2);
  forPoster->cd(1);

  chain->Draw("Ex>>hcG_T1(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==1","same");
  chain->Draw("Ex>>hcG_T2(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==2","same");
  chain->Draw("Ex>>hcG_T5(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","same");
  chain->Draw("Ex>>hcG_T7(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==7","same");
  TH1F* hcG_T1 = (TH1F*) gDirectory->Get("hcG_T1");
  hcG_T1->GetXaxis()->SetRangeUser(-1.0,+1.0);
  TH1F* hcG_T2 = (TH1F*) gDirectory->Get("hcG_T2");
  hcG_T2->SetLineColor(kRed);
    hcG_T2->SetTitle("Before Correction");
    hcG_T2->GetXaxis()->SetTitle("E_{x} [MeV]");
    hcG_T2->GetYaxis()->SetTitle("Counts / (50 keV)");
    hcG_T2->GetYaxis()->SetRangeUser(0.,500.);
    hcG_T2->GetXaxis()->SetRangeUser(-1.0,+1.0);
    hcG_T2->GetXaxis()->SetLabelSize(0.05);
    hcG_T2->GetYaxis()->SetLabelSize(0.05);
    hcG_T2->GetXaxis()->SetTitleSize(0.05);
    hcG_T2->GetYaxis()->SetTitleSize(0.05);
  hcG_T2->SetFillStyle(3444);
  hcG_T2->SetFillColor(kRed);
  TH1F* hcG_T5 = (TH1F*) gDirectory->Get("hcG_T5");
  hcG_T5->SetLineColor(kBlue);
  hcG_T5->SetFillStyle(3444);
  hcG_T5->SetFillColor(kBlue);
  TH1F* hcG_T7 = (TH1F*) gDirectory->Get("hcG_T7");
  hcG_T7->SetLineColor(kViolet);
  hcG_T7->SetFillStyle(3444);
  hcG_T7->SetFillColor(kViolet);

  hcG_T1->Draw();
  hcG_T2->Draw();
  hcG_T5->Draw("same");

  TLine *line = new TLine(0.143,0,0.143,500);
  line->SetLineColor(kBlack);
  line->SetLineStyle(kDashed);
  line->Draw();

  auto legend = new TLegend(0.1,0.7,0.35,0.9);
  legend->AddEntry(hcG_T2,"Detector #2","f");
  legend->AddEntry(hcG_T5,"Detector #5","f");
  legend->Draw();


  TFile* file = new TFile("../../../Outputs/Analysis/47K_Full_22July.root", "READ");
  TTree* tree = (TTree*) file->FindObjectAny("PhysicsTree");
  forPoster->cd(2);


  tree->Draw("Ex>>h_T1(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==1","same");
  tree->Draw("Ex>>h_T2(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==2","same");
  tree->Draw("Ex>>h_T5(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","same");
  tree->Draw("Ex>>h_T7(120,-1,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==7","same");
  TH1F* h_T1 = (TH1F*) gDirectory->Get("h_T1");
  h_T1->GetXaxis()->SetRangeUser(-1.0,+1.0);
  TH1F* h_T2 = (TH1F*) gDirectory->Get("h_T2");
  h_T2->SetLineColor(kRed);
    h_T2->SetTitle("After Correction");
    h_T2->GetXaxis()->SetTitle("E_{x} [MeV]");
    h_T2->GetYaxis()->SetTitle("Counts / (50 keV)");
    h_T2->GetYaxis()->SetRangeUser(0.,500.);
    h_T2->GetXaxis()->SetRangeUser(-1.0,+1.0);
    h_T2->GetXaxis()->SetLabelSize(0.05);
    h_T2->GetYaxis()->SetLabelSize(0.05);
    h_T2->GetXaxis()->SetTitleSize(0.05);
    h_T2->GetYaxis()->SetTitleSize(0.05);
  h_T2->SetFillStyle(3444);
  h_T2->SetFillColor(kRed);
  TH1F* h_T5 = (TH1F*) gDirectory->Get("h_T5");
  h_T5->SetLineColor(kBlue);
  h_T5->SetFillStyle(3444);
  h_T5->SetFillColor(kBlue);
  TH1F* h_T7 = (TH1F*) gDirectory->Get("h_T7");
  h_T7->SetLineColor(kViolet);
  h_T7->SetFillStyle(3444);
  h_T7->SetFillColor(kViolet);

  h_T1->Draw();
  h_T2->Draw();
  h_T5->Draw("same");

  line->Draw();
}

void AGATA_efficiency(){
  TF1 *fit_1 = new TF1("fit_1","TMath::Exp([0]+[1]*TMath::Log(x)+[2]*pow(TMath::Log(x),2.0)+[3]*pow(TMath::Log(x),3.0)+[4]*pow(TMath::Log(x),4.0))",10,5000);

  fit_1->SetParameters(-6.34543e+01,
		       +4.24746e+01,
		       -1.00304e+01,
		       +1.03468e+00,
		       -3.97076e-02);
  fit_1->Draw();
}

void AGATA_efficiency(double Energy_keV){
  TF1 *fit_1 = new TF1("fit_1","TMath::Exp([0]+[1]*TMath::Log(x)+[2]*pow(TMath::Log(x),2.0)+[3]*pow(TMath::Log(x),3.0)+[4]*pow(TMath::Log(x),4.0))",10,5000);
  fit_1->SetParameters(-6.34543e+01,
		       +4.24746e+01,
		       -1.00304e+01,
		       +1.03468e+00,
		       -3.97076e-02);
  fit_1->Draw();
  cout << "At E = " 
       << Energy_keV 
       << " keV, AGATA efficiency = " 
       << fit_1->Eval(Energy_keV)
       << " %" << endl;
}

void ElasticsGate(double EMin, double EMax){
  string gates = "abs(T_MUGAST_VAMOS-2777)<600 && MUST2.TelescopeNumber==5 && ELab > " 
	       + to_string(EMin) 
	       + " && ELab < " 
	       + to_string(EMax);

  chain->Draw("ThetaLab>>hist(80,50,90)", gates.c_str(), "");
}

void GateThetaCM(double minTheta, double maxTheta, double binsize){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && ThetaCM > " 
      + to_string(minTheta)
      + " && ThetaCM < "
      + to_string(maxTheta);

  string title = to_string(minTheta)+" < ThetaCM < "+to_string(maxTheta);
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>Ex_ThetaCMGate(" + to_string(8.0/binsize) + ",-1,7)";

  TCanvas *cEx_ThetaCMGate = new TCanvas("cEx_ThetaCMGate","cEx_ThetaCMGate",1000,1000);
  chain->Draw(draw.c_str(),gating.c_str(),"colz");
  TH1F* Ex_ThetaCMGate = (TH1F*) gDirectory->Get("Ex_ThetaCMGate");
  Ex_ThetaCMGate->GetXaxis()->SetTitle("Ex [MeV]");
  Ex_ThetaCMGate->GetYaxis()->SetTitle(ytitle.c_str());
  Ex_ThetaCMGate->SetTitle(title.c_str());

  DrawParticleStates(cEx_ThetaCMGate);
}

void GateThetaLab(double minTheta, double maxTheta, double binsize){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && ThetaLab > " 
      + to_string(minTheta)
      + " && ThetaLab < "
      + to_string(maxTheta);

  string title = to_string(minTheta)+" < ThetaLab < "+to_string(maxTheta);
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>Ex_ThetaLabGate(" + to_string(8.0/binsize) + ",-1,7)";

  TCanvas *cEx_ThetaLabGate = new TCanvas("cEx_ThetaLabGate","cEx_ThetaLabGate",1000,1000);
  chain->Draw(draw.c_str(),gating.c_str(),"colz");
  TH1F* Ex_ThetaLabGate = (TH1F*) gDirectory->Get("Ex_ThetaLabGate");
  Ex_ThetaLabGate->GetXaxis()->SetTitle("Ex [MeV]");
  Ex_ThetaLabGate->GetYaxis()->SetTitle(ytitle.c_str());
  Ex_ThetaLabGate->SetTitle(title.c_str());

  DrawParticleStates(cEx_ThetaLabGate);
}


