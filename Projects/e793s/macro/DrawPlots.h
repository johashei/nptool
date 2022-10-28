#include <math.h>
#include "NPReaction.h"
#include <string>
#include <sstream>
using namespace std;

TChain* chain=NULL ;
char cond[1000];

NPL::Reaction Cadp("47Ca(d,p)48Ca@355");
NPL::Reaction Scdp("47Sc(d,p)48Sc@355");
NPL::Reaction K46dp("46K(d,p)47K@355");

NPL::Reaction Kdp("47K(d,p)48K@355");
NPL::Reaction Kdt("47K(d,t)46K@355");//@355");
NPL::Reaction Kdd("47K(d,d)47K@355");
NPL::Reaction Kpp("47K(p,p)47K@355");
NPL::Reaction K12C12C("47K(12C,12C)47K@355");
NPL::Reaction Tidp("47Ti(d,p)48Ti@355");
NPL::Reaction Tidt("47Ti(d,t)46Ti@355");
NPL::Reaction Tidd("47Ti(d,d)47Ti@355");
NPL::Reaction Ti12C12C("47Ti(12C,12C)47Ti@355");

void KnownLines_Ex(bool isVertical, double rangemin, double rangemax, Style_t lType, Color_t lColour);
void AddGammaLines(TH1F* hist, double particle, double ymax);
void AddPlacedGammas(TH1F* hist, double ymax);

double tCentre;
double tRange;

/* BASE FUNCTIONS */
TF1* f_efficAGATA(){
  TF1 *f_E = new TF1("fit_1","TMath::Exp([0]+[1]*TMath::Log(x)+[2]*pow(TMath::Log(x),2.0)+[3]*pow(TMath::Log(x),3.0)+[4]*pow(TMath::Log(x),4.0))",10,6000);
  f_E->SetParameters(-6.34543e+01, +4.24746e+01, -1.00304e+01, +1.03468e+00, -3.97076e-02);
  return f_E; 
}

TChain* Chain(std::string TreeName, std::vector<std::string>& file, bool EventList){
  TChain*  chain = new TChain(TreeName.c_str());
  unsigned int size =file.size(); 
    for(unsigned int i = 0 ; i < size ; i++){
      cout << "Adding " << file[i] << endl;
      chain->Add(file[i].c_str());
    } 
  return chain;
}

void LoadChain47Kdp(){
  vector<string> files;
  
  //files.push_back("../../../Outputs/Analysis/OriginalValues_ptI.root");
  //files.push_back("../../../Outputs/Analysis/OriginalValues_ptII.root");
  //files.push_back("../../../Outputs/Analysis/OriginalValues_ptIII.root");

  //files.push_back("../../../Outputs/Analysis/47Kdp_08Nov_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_08Nov_PartII.root");
  
  /* With thresholds, strip matching, and bad strips out */
  //files.push_back("../../../Outputs/Analysis/47Kdp_11Apr22_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_11Apr22_PartII.root");

  /* New target thickness analysis */
  //files.push_back("../../../Outputs/Analysis/47Kdp_11Jul22_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_11Jul22_PartII.root");

  //files.push_back("../../../Outputs/Analysis/47Kdp_10Aug22_TrueStripRemoval_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_10Aug22_TrueStripRemoval_PartII.root");
  
  //files.push_back("../../../Outputs/Analysis/47Kdp_22Sep22_RmvMM5_NoRun51-52_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_22Sep22_RmvMM5_NoRun51-52_PartII.root");
  
  /*************************/

  files.push_back("../../../Outputs/Analysis/47Kdp_18Oct22_PartI.root");
  files.push_back("../../../Outputs/Analysis/47Kdp_18Oct22_PartII.root");

  chain = Chain("PhysicsTree",files,true);
}

void LoadChain47Kdt(){
  vector<string> files;
  
  /* Offset MM1 timing by -5 */
  /* Push MM1-4 +200mm in Z */

  files.push_back("../../../Outputs/Analysis/47Kdt_18Oct22_PartI.root");
  files.push_back("../../../Outputs/Analysis/47Kdt_18Oct22_PartII.root");

  chain = Chain("PhysicsTree",files,true);
}

void LoadChain47Kdd(){
  vector<string> files;
  //files.push_back("../../../Outputs/Analysis/47Kdd_08Nov_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdd_08Nov_PartII.root");

  //files.push_back("../../../Outputs/Analysis/47Kdd_11Jul22_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdd_11Jul22_PartII.root");
 
  /* Removing MM5x103+ */
  //files.push_back("../../../Outputs/Analysis/47Kdd_01Sep22_RemoveMoreMM5_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdd_01Sep22_RemoveMoreMM5_PartII.root");

  /* Now without runs 51&52, due to CATS not being in the trigger for these runs */
  files.push_back("../../../Outputs/Analysis/47Kdd_21Sep22_RmvMM5_NoRun51-52_PartI.root");
  files.push_back("../../../Outputs/Analysis/47Kdd_21Sep22_RmvMM5_NoRun51-52_PartII.root");
  //cout << RED << " ONLY USING ONE PART OF THE SORT" << RESET << endl;


  chain = Chain("PhysicsTree",files,true);
}

void LoadChain47Kpp(){
  vector<string> files;
  //files.push_back("../../../Outputs/Analysis/47Kpp_08Nov_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kpp_08Nov_PartII.root");

  //files.push_back("../../../Outputs/Analysis/47Kpp_11Jul22_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kpp_11Jul22_PartII.root");

  //files.push_back("../../../Outputs/Analysis/24Oct22_47Kpp_PartI.root");
  //files.push_back("../../../Outputs/Analysis/24Oct22_47Kpp_PartII.root");

  files.push_back("../../../Outputs/Analysis/25Oct22_47Kpp_PartI.root");
  files.push_back("../../../Outputs/Analysis/25Oct22_47Kpp_PartII.root");
  

  chain = Chain("PhysicsTree",files,true);
}

void CorrectForAGATAEffic(TH1F* hist){
  int nbins = hist->GetNbinsX();
  TF1* effic_keV = f_efficAGATA();

  for(int i=1;i<nbins+1;i++){
    double x = hist->GetBinCenter(i);
    double val = effic_keV->Eval(x*1000.);
//    cout << x << "  " << val << endl;

    double bincontent = hist->GetBinContent(i);
    bincontent = bincontent/(val/100.);
    hist->SetBinContent(i,bincontent);
  }
  hist->SetFillStyle(3244);
  hist->SetFillColor(kBlue);
  hist->Draw();
}

void DrawParticleStates(TCanvas* canvas){
  canvas->cd();
  
  canvas->Update();
  double max = canvas->GetUymax();
  TLine *gs = new TLine(0.000, 0.0, 0.000, max);
    gs->SetLineColor(kGreen);
    gs->SetLineStyle(7);
    gs->Draw();
  TLine *Sn = new TLine(4.644, 0.0, 4.644, max);
    Sn->SetLineColor(kRed);
    Sn->SetLineStyle(7);
    Sn->Draw();

  
/*  if(reactionName=="47K(d,p)"){
    TLine** lines = new TLine*[numPeaks];
    for(int j = 0; j < numPeaks; j++){
      lines[j] = new TLine(means[j],0.0,means[j],max);
      lines[j]->SetLineStyle(8);
      lines[j]->SetLineColor(kGray);
      lines[j]->Draw();
    }
  }
*/
//  if(reactionName=="47K(d,t)"){
    TLine *l1945 = new TLine(1.945, 0.0, 1.945, max);
    l1945->SetLineStyle(kDashed);
    l1945->Draw();
 
    TLine *l2233 = new TLine(2.233, 0.0, 2.233, max);
    l2233->SetLineStyle(kDashed);
    l2233->Draw();

    TLine *l3344 = new TLine(3.344, 0.0, 3.344, max);
    l3344->SetLineStyle(kDashed);
    l3344->Draw();
//  }
}

void plot_kine(NPL::Reaction r, double Ex,Color_t c,int w, int s){
  r.SetExcitation4(Ex);
  TGraph* g= r.GetKinematicLine3();
  g->SetLineColor(c) ;
  g->SetLineStyle(s) ;
  g->SetLineWidth(w) ;
  g->Draw("c same");
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

string timegate; /* defined by choice of dp or dt */
string det_gate; /* defined by choice of dp or dt */
string reactionName; /* defined by choice of dp or dt */
string exclBmDcy = "abs(AGATA_GammaE-2.013)>0.004 && abs(AGATA_GammaE-0.511)>0.003 && abs(AGATA_GammaE-0.564)>0.004 && abs(AGATA_GammaE-0.586)>0.003";

void Draw_1DGamma(){
  string gate = timegate;
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)",gate.c_str(),"");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
}

void Load_1DGamma(){
  TH1F *hEg = new TH1F("hEg","Loaded 1D Gamma Spectrum",600,-15,15);
  TFile *file = new TFile("LoadHistograms/Load_1DGamma.root","READ");
  hEg = (TH1F*)file->Get("Eg");
  hEg->Draw();
}

void Draw_1DGamma_DetGate(){
  string gate = timegate 
	      + " && " + det_gate
	      + " && " + exclBmDcy;
  if(reactionName=="47K(d,t)"){
    gate = gate + " && cutTritons && cutTime";
  }

  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)",gate.c_str(),"");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
}

void Draw_1DGamma_DetGate_x10x100(){
  string gate = timegate 
	      + " && " + det_gate
	      + " && " + exclBmDcy;

  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("AddBack_EDC>>Eg(1250,0,5)",gate.c_str(),"");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.004 MeV");
  
  for(int b = 126; b<501; b++){
    Eg->SetBinContent(b,(Eg->GetBinContent(b) * 10.));
  }
  for(int b = 501; b<1250; b++){
    Eg->SetBinContent(b,(Eg->GetBinContent(b) * 100.));
  }

  Eg->Draw();

  TLine *x10 = new TLine(0.5, 0.0, 0.5, 3000.);
    x10->SetLineColor(kBlack); x10->SetLineStyle(7);
    x10->Draw("same");
  TLine *x100 = new TLine(2.0, 0.0, 2.0, 3000.);
    x100->SetLineColor(kBlack); x100->SetLineStyle(7);
    x100->Draw("same");

}

void Draw_1DGamma_DetGate_SplitCanv(){
  string gate = timegate 
	      + " && " + det_gate
	      + " && " + exclBmDcy;

  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  gStyle->SetOptStat(0);
  cEg->Divide(2,1,0.005,0.005,0);
  cEg->cd(1); 
    gStyle->SetPadLeftMargin(0.0);
    gStyle->SetPadRightMargin(0.0);
    gPad->SetTickx();
    gPad->SetTicky();
  cEg->cd(2); 
    gStyle->SetPadLeftMargin(0.0);
    gStyle->SetPadRightMargin(0.0);
    gPad->SetTickx();
    gPad->SetTicky();

  cEg->cd(1);
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)",gate.c_str(),"");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
  Eg->GetXaxis()->SetRangeUser(0.,2.);

  cEg->cd(2);
  chain->Draw("AddBack_EDC>>Eg2(1250,0,5)",gate.c_str(),"");
  TH1F* Eg2 = (TH1F*) gDirectory->Get("Eg2");
  Eg2->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg2->GetYaxis()->SetTitle("Counts / 0.004 MeV");
  Eg2->GetXaxis()->SetRangeUser(2.,4.5);
}

void Load_1DGamma_MG(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  TH1F *hEgMG = new TH1F("hEg","Loaded 1D Gamma Spectrum, MG gated",600,-15,15);
  TFile *file = new TFile("LoadHistograms/Load_1DGamma_MG.root","READ");
  hEgMG = (TH1F*)file->Get("Eg");
  hEgMG->Draw();
}

void Load_1DParticle(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  TH1F *hEx = new TH1F("hEx","Loaded 1D Particle Spectrum",600,-15,15);
  TFile *file = new TFile("LoadHistograms/Load_1DParticle.root","READ");
  hEx = (TH1F*)file->Get("Ep");
  hEx->Draw();
}

void Load_1DParticle_SubPhaseSpace(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  TH1F *hSub = new TH1F("hSubtracted",
		  "Loaded 1D Particle Spectrum, Phase Space & Flat BG Subtracted (10May)",600,-15,15);
  TFile *file = new TFile("LoadHistograms/Load_1DParticle_SubPhaseSpace.root","READ");
  hSub = (TH1F*)file->Get("ExSubPSpace");
  hSub->Draw();
}

void Draw_1DParticle(){
  string gate2;
  string gate3;
  string gate = timegate 
	      + " && " + det_gate
              + " && Ex@.size()==1";
  if(reactionName=="47K(d,t)"){
    gate = gate + " && cutTritons && cutTime";
  }

  TCanvas *cEx = new TCanvas("cEx","cEx",1000,1000);
  chain->Draw("Ex>>Ep(400,-2,8)", gate.c_str(),"");
  TH1F* Ep = (TH1F*) gDirectory->Get("Ep");
  Ep->GetXaxis()->SetTitle("Ex [MeV]");
  Ep->GetYaxis()->SetTitle("Counts / 0.025 MeV");

  if(reactionName=="47K(d,t)"){
    Ep->Rebin(2); 
    Ep->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  }

  if(reactionName=="47K(d,p)"){DrawParticleStates(cEx);}
}

void Load_2DParticleGamma(){
  TCanvas *cExEg = new TCanvas("cExEg","cExEg",1000,1000);
  TH2F *hExEg = new TH2F("hExEg","Loaded 2D Particle-Gamma",600,-15,15,2500,0,5);
  TFile *file = new TFile("LoadHistograms/Load_2DParticleGamma.root","READ");
  hExEg = (TH2F*)file->Get("ExEg");
  hExEg->Draw("colz");
}

void Draw_2DParticleGamma(){
  string gate = timegate 
	      + " && " + det_gate
              + " && Ex@.size()==1";
  if(reactionName=="47K(d,t)"){
    gate = gate + " && cutTritons && cutTime";
  }

  TCanvas *cExEg = new TCanvas("cExEg","cExEg",1000,1000);
  gStyle->SetOptStat(0);
  //chain->Draw("AddBack_EDC:Ex>>ExEg(600,-15,15,2500,0,5)", gate.c_str(), "colz");
  chain->Draw("Ex:AddBack_EDC>>ExEg(5000,0,5,3000,-15,15)", gate.c_str(), "");
  TH1F* ExEg = (TH1F*) gDirectory->Get("ExEg");
  ExEg->SetTitle("");
  ExEg->GetYaxis()->SetTitle("Ex [MeV]");
  ExEg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  ExEg->GetYaxis()->SetRangeUser(-1.0,7.0);
  ExEg->Draw();
  if(reactionName=="47K(d,p)"){
    TLine *Sn = new TLine(0,4.644,4.644,4.644);
      Sn->SetLineColor(kRed);Sn->SetLineStyle(2);
      Sn->Draw("same");
    TLatex *TSn = new TLatex(.5,.5,"S_{n}");
      TSn->SetTextColor(kRed);
      TSn->SetTextSize(0.05);
      TSn->SetX(2.50);
      TSn->SetY(4.90);
      TSn->Draw("same");
  }
  TLine *XeqY = new TLine(0,0,5,5);
    XeqY->SetLineColor(kRed);
    XeqY->Draw("same");
  TLatex *Texeg = new TLatex(.5,.5,"Ex = E_{#gamma}");
    Texeg->SetTextColor(kRed);
    Texeg->SetTextSize(0.05);
    Texeg->SetX(2.35);
    Texeg->SetY(1.50);
    Texeg->Draw("same");
}

void Load_2DGammaGamma(){
  TCanvas *cEgEg = new TCanvas("cEgEg","cEgEg",1000,1000);
  TH2F *hEgEg = new TH2F("hEgEg","Loaded 2D Gamma-Gamma",600,-15,15,2500,0,5);
  TFile *file = new TFile("LoadHistograms/Load_2DGammaGamma.root","READ");
  hEgEg = (TH2F*)file->Get("gg");
  hEgEg->SetName("hEgEg");
  hEgEg->Draw("colz");
}

void Draw_2DGammaGamma_ExcludeBeam(){
  string gate = timegate 
              + " && " + exclBmDcy;

  TCanvas *cEgEg = new TCanvas("cEgEg","cEgEg",1000,1000);
  chain->Draw("AddBack_EDC:AddBack_EDC2>>EgEg(2500,0,5,2500,0,5)",gate.c_str(),"colz");
  TH1F* EgEg = (TH1F*) gDirectory->Get("EgEg");
  chain->Draw("AddBack_EDC2:AddBack_EDC>>EgEg2(2500,0,5,2500,0,5)",gate.c_str(),"colz");
  TH1F* EgEg2 = (TH1F*) gDirectory->Get("EgEg2");
  EgEg->Add(EgEg2,1);
  EgEg->SetTitle("Egamma-Egamma");
  EgEg->GetXaxis()->SetTitle("Eg [Counts / 0.002 MeV]");
  EgEg->GetYaxis()->SetTitle("Eg [Counts / 0.002 MeV]");
  EgEg->GetXaxis()->SetRangeUser(0.005,5.0);
  EgEg->GetYaxis()->SetRangeUser(0.005,5.0);
  EgEg->Draw("colz");
  TLine *XeqY = new TLine(0,0,5,5);
  XeqY->SetLineColor(kRed);
  XeqY->SetLineStyle(kDashed);
  XeqY->Draw("same");
}

/*
void gg(){
  // Example filling of GG matrix 
  unsigned int size = AddBack_EDC.size();
  double e1,e2;
  auto h=new TH2("gg","gg",1000,0,10,1000,0,10);
  for(unsigned int i = 0 ; i < size-1 ; i++){
    e1=AddBack_EDC[i] ; e2 = AddBack_EDC[i+1];
    // Folding of the matrix, always fill big first
    if (e1>e2) {h->Fill(e1,e2);}
    else {h->Fill(e2,e1);}
  }
}
*/

void Draw_2DGammaGamma_TimeGated(){
  string gate = timegate;
  TCanvas *cEgEg = new TCanvas("cEgEg","cEgEg",1000,1000);
  chain->Draw("AddBack_EDC:AddBack_EDC2>>EgEg(2500,0,5,2500,0,5)",gate.c_str(),"colz");
  TH1F* EgEg = (TH1F*) gDirectory->Get("EgEg");
  chain->Draw("AddBack_EDC2:AddBack_EDC>>EgEg2(2500,0,5,2500,0,5)",gate.c_str(),"colz");
  TH1F* EgEg2 = (TH1F*) gDirectory->Get("EgEg2");
  EgEg->Add(EgEg2,1);
  EgEg->SetTitle("Egamma-Egamma");
  EgEg->GetXaxis()->SetTitle("Eg [Counts / 0.002 MeV]");
  EgEg->GetYaxis()->SetTitle("Eg [Counts / 0.002 MeV]");
  EgEg->GetXaxis()->SetRangeUser(0.005,5.0);
  EgEg->GetYaxis()->SetRangeUser(0.005,5.0);
  EgEg->Draw("colz");
  TLine *XeqY = new TLine(0,0,5,5);
  XeqY->SetLineColor(kRed);
  XeqY->SetLineStyle(kDashed);
  XeqY->Draw("same");
}

void GateGamma_SeeParticle(double gamma, double width, double binsize){
  string gating = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  if(reactionName=="47K(d,t)"){
    gating = gating + " && cutTritons && cutTime";
  }

  string title = to_string(gamma-width)+" < Eg < "+to_string(gamma+width);
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>ExGate(" + to_string(10.0/binsize) + ",-2,8)";

  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  //chain->Draw("Ex>>ExGate(60,-1,5)",gating.c_str(),"colz");
  chain->Draw(draw.c_str(),gating.c_str(),"colz");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle(ytitle.c_str());
  ExGate->SetTitle(title.c_str());
  
  if(reactionName=="47K(d,p)"){DrawParticleStates(cEx_Gate);}
}

void GateGamma_SeeParticle_WithBG(double gamma, double width, double bg){
  string gating = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(width);
  if(reactionName=="47K(d,t)"){
    gating = gating + " && cutTritons && cutTime";
    bggate = bggate + " && cutTritons && cutTime";
  }

  string title = "Gate: "+to_string(gamma-width)+" to "+to_string(gamma+width)+"."
	  + "  BG: "+to_string(bg-width)+" to "+to_string(bg+width)+".";
  
  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  chain->Draw("Ex>>ExGate(200,-2,8)",gating.c_str(),"");
  //chain->Draw("Ex>>ExGate(120,-1,5)",gating.c_str(),"");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle("Counts / 0.10 MeV");
  //ExGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ExGate->SetLineColor(kGreen);
  ExGate->SetFillColor(kGreen);
  ExGate->SetFillStyle(3154);
  ExGate->SetTitle(title.c_str());

  chain->Draw("Ex>>ExBG(600,-15,15)",bggate.c_str(),"same");
  //chain->Draw("Ex>>ExBG(120,-1,5)",bggate.c_str(),"same");
  TH1F* ExBG = (TH1F*) gDirectory->Get("ExBG");
  ExBG->SetLineColor(kRed);
  ExBG->SetFillColor(kRed);
  ExBG->SetFillStyle(3345);

  DrawParticleStates(cEx_Gate);
}

void GateGamma_SeeParticle_WithBG(double gamma, double width, double bg, double widthbg){
  string gating = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(widthbg);
  if(reactionName=="47K(d,t)"){
    gating = gating + " && cutTritons && cutTime";
    bggate = bggate + " && cutTritons && cutTime";
  }

  double ratio = width/widthbg;

  string title = "Gate: "+to_string(gamma-width)+" to "+to_string(gamma+width)+"."
	  + "  BG: "+to_string(bg-width)+" to "+to_string(bg+width)+".";
  
  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  chain->Draw("Ex>>ExGate(200,-2,8)",gating.c_str(),"");
  //chain->Draw("Ex>>ExGate(120,-1,5)",gating.c_str(),"");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle("Counts / 0.10 MeV");
  //ExGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ExGate->SetLineColor(kGreen);
  ExGate->SetFillColor(kGreen);
  ExGate->SetFillStyle(3154);
  ExGate->SetTitle(title.c_str());

  chain->Draw("Ex>>ExBG(200,-2,8)",bggate.c_str(),"same");
  //chain->Draw("Ex>>ExBG(120,-1,5)",bggate.c_str(),"same");
  TH1F* ExBG = (TH1F*) gDirectory->Get("ExBG");
  ExBG->Scale(ratio);
  ExBG->SetLineColor(kRed);
  ExBG->SetFillColor(kRed);
  ExBG->SetFillStyle(3345);
  //ExBG->Draw("BSAME");

  DrawParticleStates(cEx_Gate);
}

void GateParticle_SeeGamma(double particle, double width){ 
  gStyle->SetOptStat("nemMrRi");

  string gating = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(Ex-" 
      + to_string(particle)
      + ")<"
      + to_string(width);
  if(reactionName=="47K(d,t)"){
    gating = gating + " && cutTritons && cutTime";
  }

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

  AddGammaLines(EgGate, particle, cEg_Gate->GetUymax());

}

void GateParticle_SeeGamma_WithBG(double particle, double width, double bg, double width2){
  string gating = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(Ex-" 
      + to_string(particle)
      + ")<"
      + to_string(width);
  string bggate = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(Ex-" 
      + to_string(bg)
      + ")<"
      + to_string(width2);
  if(reactionName=="47K(d,t)"){
    gating = gating + " && cutTritons && cutTime";
    bggate = bggate + " && cutTritons && cutTime";
  }

  double ratio = width/width2;

  string title = "Gate: "+to_string(particle-width)+" to "+to_string(particle+width)+"."
	  + "  BG: "+to_string(bg-width2)+" to "+to_string(bg+width2)+".";

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
  EgBG->Scale(ratio);
  EgBG->SetTitle(title.c_str());
  EgBG->SetLineColor(kRed);
  EgBG->SetFillColor(kRed);
  EgBG->SetFillStyle(3345);
}
 
/*
void GateGamma_SeeGamma_ExcludeBeamDecay(double gamma, double width){
  string gating 
      = "abs(AGATA_GammaE-2.013)>0.004 && abs(AGATA_GammaE-0.511)>0.003 && abs(AGATA_GammaE-0.564)>0.004 && abs(AGATA_GammaE-0.586)>0.003 && abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string gating2
      = "abs(AGATA_GammaE-2.013)>0.004 && abs(AGATA_GammaE-0.511)>0.003 && abs(AGATA_GammaE-0.564)>0.004 && abs(AGATA_GammaE-0.586)>0.003 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  string title = to_string(gamma-width) + " < Eg < " + to_string(gamma+width); 
  TCanvas *cEx_Gate = new TCanvas("cggGate","cggGate",1000,1000);

  chain->Draw("AddBack_EDC>>ggGate(999,0.005,5)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.005 MeV");
  ggGate->SetTitle(title.c_str());

  chain->Draw("AddBack_EDC2>>ggGate2(999,0.005,5)",gating2.c_str(),"");
  TH1F* ggGate2 = (TH1F*) gDirectory->Get("ggGate2");
  ggGate->Add(ggGate2,1);
  ggGate->Draw();
}
*/

void GateGamma_SeeGamma(double gamma, double width){
  string gating = "abs(AGATA_GammaE-2.013)>0.004 && abs(AGATA_GammaE-0.511)>0.003 && abs(AGATA_GammaE-0.564)>0.004 && abs(AGATA_GammaE-0.586)>0.00 && abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string gating2 = "abs(AGATA_GammaE-2.013)>0.004 && abs(AGATA_GammaE-0.511)>0.003 && abs(AGATA_GammaE-0.564)>0.004 && abs(AGATA_GammaE-0.586)>0.00 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  string title = to_string(gamma-width) + " < Eg < " + to_string(gamma+width);
  TCanvas *cEx_Gate = new TCanvas("cggGate","cggGate",1000,1000);

  chain->Draw("AddBack_EDC>>ggGate(999,0.005,5)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.005 MeV");
  ggGate->SetTitle(title.c_str());

  chain->Draw("AddBack_EDC2>>ggGate2(999,0.005,5)",gating2.c_str(),"");
  TH1F* ggGate2 = (TH1F*) gDirectory->Get("ggGate2");
  ggGate->Add(ggGate2,1);
  ggGate->Draw();
}

void GateGamma_SeeGamma_TimeGate(double gamma, double width){
  string gating = timegate + " && abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string gating2 = timegate + " && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  string title = to_string(gamma-width) + " < Eg < " + to_string(gamma+width);
  
  TCanvas *cEx_Gate = new TCanvas("cggGate","cggGate",1000,1000);
  chain->Draw("AddBack_EDC>>ggGate(999,0.005,5)",gating.c_str(),"");
  TH1F* ggGate = (TH1F*) gDirectory->Get("ggGate");
  ggGate->GetXaxis()->SetTitle("Eg [MeV]");
  ggGate->GetYaxis()->SetTitle("Counts / 0.005 MeV");
  ggGate->SetTitle(title.c_str());

  chain->Draw("AddBack_EDC2>>ggGate2(999,0.005,5)",gating2.c_str(),"");
  TH1F* ggGate2 = (TH1F*) gDirectory->Get("ggGate2");
  ggGate->Add(ggGate2,1);
  ggGate->Draw();
}

void GateGamma_SeeGamma_WithBG(double gamma, double width, double bg, double width2){
  string gating = "abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = "abs(AddBack_EDC2-" 
      + to_string(bg)
      + ")<"
      + to_string(width2);
  string gating2 = "abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate2 = "abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(width2);

  double ratio = width/width2;
  
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
  ggBG->Scale(ratio);
  ggBG->SetLineColor(kRed);
  ggBG->SetFillColor(kRed);
  ggBG->SetFillStyle(3345);

  ggGate->Draw();
  ggBG->Draw("same");
}

void CompareExsAt4MeV(){
  TCanvas *cExCompare = new TCanvas("cExCompare","cExCompare",1000,1000);
  chain->Draw("AddBack_EDC>>gate3p0(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && abs(Ex-3.0)<0.1","same");
  chain->Draw("AddBack_EDC>>gate3p5(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && abs(Ex-3.5)<0.1","same");
  chain->Draw("AddBack_EDC>>gate3p9(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && abs(Ex-3.9)<0.1","same");
  chain->Draw("AddBack_EDC>>gate4p3(1000,0,10)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && abs(Ex-4.3)<0.1","same");
 
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
  
  chain->Draw("Ex>>hexp(600,-15,15)","abs(T_MUGAST_VAMOS-2700)<400","");
  TH1F* hexp = (TH1F*) gDirectory->Get("hexp");
  hexp->SetTitle("Comparing Simulation to Experiment");
  hexp->GetXaxis()->SetTitle("Ex [MeV]");
  hexp->GetYaxis()->SetTitle("Counts / 0.1 MeV");
  hexp->SetLineColor(kRed);

  TFile* simfile = new TFile("../../../Outputs/Analysis/SimTest_Jun28_p32_ExHeavy0p143.root", 
  //TFile* simfile = new TFile("../../../Outputs/Analysis/SimTest_Jun22_TWOFNR_p32.root", 
		  "READ");
  TTree* simtree = (TTree*) simfile->FindObjectAny("PhysicsTree");
  simtree->Draw("Ex>>hsimMGp32(600,-15,15)",
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

void ForPoster_DiffCrossSec(){
  ifstream infile("DiffCrossSecInputfile.txt");
  vector<double> theta, p32, p32Pos, p32Neg, f72, f72Pos, f72Neg, zero;

  double a, b, c, d, e, f, g;
  while(infile){
    infile >> a >> b >> c >> d >> e >> f >> g;

    theta.push_back(a);
    zero.push_back(0.0);
    
    p32.push_back(b);
    p32Pos.push_back(c);
    p32Neg.push_back(d);
    
    f72.push_back(e);
    f72Pos.push_back(f);
    f72Neg.push_back(g);
  }

  TGraph *graphp32 = new TGraph(theta.size(), &theta[0], &p32[0]);
  TGraph *graphf72 = new TGraph(theta.size(), &theta[0], &f72[0]);

  int num = theta.size();

  TGraphAsymmErrors *graphp32error 
  = new TGraphAsymmErrors(num, &theta[0],  &p32[0], &zero[0], &zero[0], &p32Neg[0], &p32Pos[0]);

  TGraphAsymmErrors *graphf72error 
  = new TGraphAsymmErrors(num, &theta[0],  &f72[0], &zero[0], &zero[0], &f72Neg[0], &f72Pos[0]);


  graphp32->SetFillColor(kRed);
  graphp32->SetFillStyle(3005);
  graphf72->SetFillColor(kBlue);
  graphf72->SetFillStyle(3005);
  
  graphp32->Draw();
  graphf72->Draw("same");
  graphp32error->Draw("same");
  graphf72error->Draw("same");

/*
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(graphp32);
  mg->Add(graphf72);
  mg->Add(graphp32error);
  mg->Add(graphf72error);
  mg->Draw();
*/


}

void MugastMisalignment(){
  TCanvas *cMisaligned = new TCanvas("cMisaligned","cMisaligned",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("Ex>>hcG_T1(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==1","");
  chain->Draw("Ex>>hcG_T2(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==2","same");
  chain->Draw("Ex>>hcG_T3(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==3","same");
  chain->Draw("Ex>>hcG_T4(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==4","same");
  chain->Draw("Ex>>hcG_T5(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==5","same");
  chain->Draw("Ex>>hcG_T7(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==7","same");
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

  string base = "abs(T_MUGAST_VAMOS-2700)<400 && abs(AddBack_EDC-" + to_string(gamma) 
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
    "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
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
    "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
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
    "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    "colz");
  TH1F* phiHist2 = (TH1F*) gDirectory->Get("phiHist2");  
  
  diagnosePhi->Update();
  TLine *l0143_2 = new TLine(-180., 0.143, 180., 0.143);
    l0143_2->SetLineColor(kRed);
    l0143_2->Draw();
}

void ExThetaLab(){
  string gate = timegate 
	      + " && " + det_gate
              + " && Ex@.size()==1";
  if(reactionName=="47K(d,t)"){
    gate = gate + " && cutTritons && cutTime";
  }

  TCanvas *diagnoseTheta = new TCanvas("diagnoseTheta","diagnoseTheta",1000,1000);
  chain->Draw(
    "Ex:ThetaLab>>thetaHist(360,0,180,180,-1,8)", 
    gate.c_str(), "colz");
  TH1F* thetaHist = (TH1F*) gDirectory->Get("thetaHist");  
  thetaHist->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  thetaHist->GetYaxis()->SetTitle("Ex [MeV]");
  thetaHist->SetTitle("Theta dependance testing");
  
  diagnoseTheta->Update();
  
  TLine *l0000 = new TLine(100., 0.000, 160., 0.000);
    l0000->SetLineStyle(kDashed);
    l0000->SetLineColor(kRed);
    l0000->Draw();
  /*
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
  */
  TLine *Sn = new TLine(100., 4.644, 160., 4.644);
    Sn->SetLineStyle(kDashed);
    Sn->SetLineColor(kRed);
    Sn->Draw();
}

void ExThetaLab(double gamma, double width){
  TCanvas *diagnoseTheta = new TCanvas("diagnoseTheta","diagnoseTheta",1000,1000);
  string gate = timegate 
	      + " && " + det_gate
              + " && Ex@.size()==1";
  if(reactionName=="47K(d,t)"){
    gate = gate + " && cutTritons && cutTime";
  }
  gate = gate + "&& abs(AddBack_EDC-"
	      + to_string(gamma) + ") < " + to_string(width); 

  chain->Draw("Ex:ThetaLab>>thetaHist(360,0,180,100,-1,9)", gate.c_str(), "colz");
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
  TCanvas *cELabTLab = new TCanvas("cELabTLab","cELabTLab",1000,1000);
  gStyle->SetOptStat(0);

  string gate = timegate 
	      + " && " + det_gate;
  if(reactionName=="47K(d,t)"){
    gate = gate + " && cutTritons && cutTime";
  }


  chain->Draw("ELab:ThetaLab>>hKine(360,0,180,500,0,10)",gate.c_str(),"col");
  TH2F* hKine = (TH2F*) gDirectory->Get("hKine");
  hKine->SetTitle("");
  hKine->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  hKine->GetYaxis()->SetTitle("E_{lab} [MeV]");

  plot_kine(K12C12C, 0.000, kRed, 2, 1);

  if(reactionName=="47K(d,t)"){
    cout << "  Trying to draw lines for " << reactionName << endl;
    plot_kine(Kdt, 0.000, kBlack, 2, 1);
    plot_kine(Kdt, 1.944, kBlack, 2, 1);
    plot_kine(Kdt, 3.340, kBlack, 2, 1);
    plot_kine(Kdt, 4.3  , kBlack, 2, 1);
    plot_kine(Kdt, 5.8  , kBlack, 2, 1);
  }
  
  if(reactionName=="47K(d,d)"){
    cout << "  Trying to draw lines for " << reactionName << endl;
    plot_kine(Kdd, 0.000, kBlack, 2, 1);
    plot_kine(Kpp, 0.000, kBlack, 2, 6);
  }
  
  if(reactionName=="47K(p,p)"){
    cout << "  Trying to draw lines for " << reactionName << endl;
    plot_kine(Kdd, 0.000, kBlack, 2, 6);
    plot_kine(Kpp, 0.000, kBlack, 2, 1);
  }

  if(reactionName=="47K(d,p)"){
    cout << "  Trying to draw lines for " << reactionName << endl;
    plot_kine(Kdp, 0.000, kBlack, 2, 1);
    plot_kine(Kdp, 4.644, kBlack, 2, 1);
    plot_kine(Cadp, 0.000, kRed, 2, 1);
    plot_kine(Tidp, 0.000, kBlue, 2, 1);
    plot_kine(Scdp, 0.000, kGreen, 2, 1);
    plot_kine(K46dp, 0.000, kViolet, 2, 1);
  }
}

void XYMust2(){
  TCanvas *cXYMust2 = new TCanvas("cXYMM","cXYMM",1000,1000);
  chain->Draw("Y:X>>hXYMust2(300,-150,+150,300,-150,+150)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && MUST2.TelescopeNumber>0 && MUST2.TelescopeNumber<5",
		  "colz");
}

void XYMugast(){
  TCanvas *cXYMugast = new TCanvas("cXYMG","cXYMG",1000,1000);
  chain->Draw("Y:X>>hXYMugast(150,-150,+150, 150,-150,+150)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
		  "colz");
}

void MM5_ELabThetaLab(){
  chain->Draw("ELab:ThetaLab>>hMM5_el(180,50,95,700,0,7)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && MUST2.TelescopeNumber==5",
		  "colz");
  TH2F* hMM5_el = (TH2F*) gDirectory->Get("hMM5_el");
  hMM5_el->GetXaxis()->SetTitle("Theta Lab");
  hMM5_el->GetYaxis()->SetTitle("Ex");

  plot_kine(Kdd, 0, kGreen+2, 2, 9);
  plot_kine(Kpp, 0, kYellow, 2, 9);
}

void MM5_RawEThetaLab(){
  chain->Draw("RawEnergy:ThetaLab>>hMM5_el(90,50,95,700,0,7)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && MUST2.TelescopeNumber==5",
		  "colz");

  //plot_kine(Kdd, 0, kGreen+2, 2, 9);
  //plot_kine(Kpp, 0, kYellow, 2, 9);
}

void MM5_ExThetaLab(){
  chain->Draw("Ex:ThetaLab>>hMM5_ex(180,0,180,400,0,20)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && Must2.TelescopeNumber==5",
		  "colz");
}

void thickness(){

  std::ifstream infile("thicknessTheory4.txt");
  
  // THEORY -------------------------------------
  double x, pCH, dHSS, pBG, pPer, dKD, d79D, dPer, dLH, dBel;
  vector<double> vx, vpCH, vdHSS, vpBG, vpPer, vdKD, vd79D, vdPer, vdLH, vdBel;
  int count = 0;

  while (infile >> x >> pCH >> dHSS >> pBG >> pPer >> dKD >> d79D >> dPer >> dLH >> dBel)
  {
    vx.push_back(x);
    vpBG.push_back(pBG);
    vdHSS.push_back(dHSS);
    vpCH.push_back(pCH);
    vpPer.push_back(pPer);
    vdKD.push_back(dKD);
    vd79D.push_back(d79D);
    vdPer.push_back(dPer);
    vdLH.push_back(dLH);
    vdBel.push_back(dBel);

    count++;
  }
  // --------------------------------------------

  // EXPERIMENT ---------------------------------

double expDx[20]= {22.5	,
23.5	,
24.5	,
25.5	,
26.5	,
27.5	,
28.5	,
29.5	,
30.5	,
31.5	,
32.5	,
33.5	,
34.5	,
35.5	,
36.5	,
37.5	,
38.5	,
39.5	,
40.5	,
41.5	};

double expDy[20]={
0.595726876922256	,
0.543213473232366	,
0.494283286886849	,
0.421774881771817	,
0.341156378590349	,
0.307099360773692	,
0.262569839264486	,
0.232682384534885	,
0.216888228577843	,
0.23470124429325	,
0.240399882846714	,
0.228134082861555	,
0.32920533438895	,
0.296543336595945	,
0.367649491313769	,
0.413065289661424	,
0.437192259660495	,
0.440938347607344	,
0.378811580496408	,
0.404394227640296	};



double expDyErr[20]={
0.028324914039265	,
0.0218144191759532	,
0.0142504076860264	,
0.0146454365657758	,
0.0088989334482194	,
0.00824524006184739	,
0.0108064634609478	,
0.0111948065344291	,
0.00904509721171868	,
0.0119548919107791	,
0.0180227167275249	,
0.0201388245171666	,
0.0213346009360056	,
0.0193056550679402	,
0.0155440903197556	,
0.0136119544115286	,
0.0204840651640718	,
0.0248743996273666	,
0.023151640426499	,
0.0183223089693594	};

double expDyErr2[20]={0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5};

double expPx[15]={31.5	,
32.5	,
33.5	,
34.5	,
35.5	,
36.5	,
37.5	,
38.5	,
39.5	,
40.5	,
41.5	,
42.5	,
43.5	,
44.5	,
45.5	};

double expPy[15]={
0.927940963711002	,
0.886373091269099	,
0.997232654140559	,
0.850828825193078	,
0.975515571786442	,
0.799588217203809	,
0.816607420458171	,
0.982225026964146	,
0.888146409552698	,
0.901367241759583	,
1.07321471789101	,
1.0011298892341	,
0.966230221622026	,
1.1926071484492	,
0.865315617070875	};





double expPyErr[15]={
0.0990008372934202	,
0.14013781860333	,
0.0725936416502272	,
0.066957912885687	,
0.0518874166956934	,
0.0314030147725784	,
0.0396294675991484	,
0.0474887746673471	,
0.0504518388646933	,
0.0432543357898785	,
0.0509856961378907	,
0.0376160621523655	,
0.043968910197424	,
0.0649913895880728	,
0.0576066735490418	};







double expPyErr2[15]={0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5};




  /*  double expDx[16] =   {22.8659168020674,
			24.4932387071222,
			25.4884605201683,
			26.4478963605459,
			27.8768790582849,
			30.0521503123866,
			31.6308395766488,
			32.4198704973062,

			34.0083546157897,
			35.4165962236503,
			36.4835938795446,
			37.5238982330089,
			38.5397371697557,
			39.5330583015003,
			40.5055760670275,
			41.3};

  double expDy[16] =   {9.52350198819289,
			5.64694896708287,
			3.57621524475264,
			2.22306813461932,
			1.84242236351724,
			0.839564977947224,
			0.533641048642503,
			0.811298146594814,
			0.452491201918824,
			0.628258484221101,
			0.700018256426215,
			0.598634895567586,
			0.510176367468261,
			0.457345826343762,
			0.341874900915227,
			0.41491476578731};

  double expPx[9] =    {30.6251132222237,
			32.9493281027741,
			34.6378972639558,
			37.6049166106808,
			39.1281033807373,
			41.0541710577087,
			43.6797461320023,
			43.5133105051628,
			45.6589118858263};

  double expPy[9] =    {15.4053624168101,
			10.0326697910848,
			11.2010890871452,
			6.28484827873840,
			5.93048525466625,
			5.18590383783819,
			4.93601112583176,
			3.77917031974004,
			6.32306635135636};
  */
  // --------------------------------------------

  //cout << vx.front() << " to " << vx.back() << "   count " << count << endl;

  TCanvas *canvThick = new TCanvas("canvThick","canvThick",1000,1000);
  canvThick->SetLogy();

  int HSS = 10, Perey = 4, BGreen = 1;


  TGraph* gp1 = new TGraph(vx.size(), &vx[0], &vpBG[0]);
    gp1->GetXaxis()->SetLimits(10.,80.);
    gp1->SetLineColor(kRed);
    gp1->SetLineStyle(BGreen);
    gp1->SetLineWidth(2);
    gp1->SetTitle("(p,p) Bechetti-Greenlees");

  TGraph* gp2 = new TGraph(vx.size(), &vx[0], &vpCH[0]);
    gp2->GetXaxis()->SetLimits(10.,80.);
    gp2->SetLineColor(kRed);
    gp2->SetLineStyle(Perey);
    gp2->SetLineWidth(2);
    gp2->SetTitle("(p,p) Chapel-Hill");

  TGraph* gp3 = new TGraph(vx.size(), &vx[0], &vpPer[0]);
    gp3->GetXaxis()->SetLimits(10.,80.);
    gp3->SetLineColor(kRed);
    gp3->SetLineStyle(HSS);
    gp3->SetLineWidth(2);
    gp3->SetTitle("(p,p) Perey");

  TGraph* gd1 = new TGraph(vx.size(), &vx[0], &vdHSS[0]);
    gd1->GetXaxis()->SetLimits(10.,80.);
    gd1->SetLineColor(kBlue);
    gd1->SetLineStyle(HSS);
    gd1->SetLineWidth(2);
    gd1->SetTitle("(d,d) HSS");

  TGraph* gd2 = new TGraph(vx.size(), &vx[0], &vdKD[0]);
    gd2->GetXaxis()->SetLimits(10.,80.);
    gd2->SetLineColor(kBlue);
    gd2->SetLineStyle(1);
    gd2->SetLineWidth(2);
    gd2->SetTitle("(d,d) Koning-Delaroche");

  TGraph* gd3 = new TGraph(vx.size(), &vx[0], &vd79D[0]);
    gd3->GetXaxis()->SetLimits(10.,80.);
    gd3->SetLineColor(kBlue);
    gd3->SetLineStyle(6);
    gd3->SetLineWidth(2);
    gd3->SetTitle("(d,d) 79DCV");

  TGraph* gd4 = new TGraph(vx.size(), &vx[0], &vdPer[0]);
    gd4->GetXaxis()->SetLimits(10.,80.);
    gd4->SetLineColor(kBlue);
    gd4->SetLineStyle(Perey);
    gd4->SetLineWidth(2);
    gd4->SetTitle("(d,d) Perey");

  TGraph* gd5 = new TGraph(vx.size(), &vx[0], &vdLH[0]);
    gd5->GetXaxis()->SetLimits(10.,80.);
    gd5->SetLineColor(kBlue);
    gd5->SetLineStyle(2);
    gd5->SetLineWidth(2);
    gd5->SetTitle("(d,d) Lohr-Haeberli");

  TGraph* gd6 = new TGraph(vx.size(), &vx[0], &vdBel[0]);
    gd6->GetXaxis()->SetLimits(10.,80.);
    gd6->SetLineColor(kBlue);
    gd6->SetLineStyle(9);
    gd6->SetLineWidth(2);
    gd6->SetTitle("(d,d) Belote 48Ca");

  TGraphErrors* expP = new TGraphErrors( 15, expPx, expPy, expPyErr2, expPyErr);
  expP->SetTitle("(p,p) Experiment");
    //expP->SetMarkerStyle(22);
  TGraphErrors* expD = new TGraphErrors(20, expDx, expDy, expDyErr2, expDyErr);
  expD->SetTitle("(d,d) Experiment");
    //expD->SetMarkerStyle(20);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(gp1);
  mg->Add(gp2);
  mg->Add(gp3);
  mg->Add(gd1);
  mg->Add(gd2);
  mg->Add(gd3);
  mg->Add(gd4);
  mg->Add(gd5);
  mg->Add(gd6);
  //mg->GetXaxis()->SetTitle("{#theta}_{CM}");
  //mg->GetYaxis()->SetTitle("Elastic counts [mb/sr]");
  mg->Draw("AC");
  mg->GetXaxis()->SetTitle("#theta_{CM} [deg]");
  mg->GetYaxis()->SetTitle("Ratio #sigma/#sigma_{Rutherford}");
  expP->Draw("same*");
  expD->Draw("same*");

  canvThick->BuildLegend();

}

/*
void ExThetaAnalysis(double gamma, double width, int version){

  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && abs(AddBack_EDC-"
	        + to_string(gamma) + ") < " + to_string(width); 

  TCanvas *diagnoseTheta2 = new TCanvas("diagnoseTheta2","diagnoseTheta2",1000,1000);
  chain->Draw(
    "Ex:ThetaLab>>thetaHist(60,100,160,120,-1,5)", 
    //"abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    gating.c_str(),
    "colz");
  TH2F* thetaHist = (TH2F*) gDirectory->Get("thetaHist");  


//  while(running){
//    cout << "Overlay projections (1) or candle plots (2)?" << endl;
//    cin >> version;
    cout << "Processing..." << endl;
    if(version==1){
      thetaHist->ProjectionY("tpy1",06.,15.); 
      TH1F* tpy1 = (TH1F*) gDirectory->Get("tpy1");  
      thetaHist->ProjectionY("tpy2",16.,25.); 
      TH1F* tpy2 = (TH1F*) gDirectory->Get("tpy2");  
      thetaHist->ProjectionY("tpy3",26.,35.); 
      TH1F* tpy3 = (TH1F*) gDirectory->Get("tpy3");  
      thetaHist->ProjectionY("tpy4",36.,45.); 
      TH1F* tpy4 = (TH1F*) gDirectory->Get("tpy4");  
      thetaHist->ProjectionY("tpy5",46.,55.); 
      TH1F* tpy5 = (TH1F*) gDirectory->Get("tpy5");  

      tpy1->SetLineColor(kRed);
      tpy2->SetLineColor(kOrange);
      tpy3->SetLineColor(kGreen);
      tpy4->SetLineColor(kBlue);
      tpy5->SetLineColor(kViolet);

      tpy1->Rebin(2);
      tpy2->Rebin(2);
      tpy3->Rebin(2);
      tpy4->Rebin(2);
      tpy5->Rebin(2);

      tpy1->Draw();
      tpy2->Draw("same");
      tpy3->Draw("same");
      tpy4->Draw("same");
      tpy5->Draw("same");
    }else if (version==2){
      thetaHist->GetXaxis()->SetRangeUser(105.,155.);
      thetaHist->RebinX(5);
      thetaHist->GetXaxis()->SetRangeUser(105.,155.);
      thetaHist->Draw("candlex6");

    }//else{running=0;}
//  }
}
*/

void ExTheta_Analysis(double gamma, double width){
  string gate = timegate 
	      + " && " + det_gate
              + " && Ex@.size()==1"
              + " && abs(AddBack_EDC-"
	      + to_string(gamma) + ")<"
	      + to_string(width);
  if(reactionName=="47K(d,t)"){
    gate = gate + " && cutTritons && cutTime";
  }

  string gateLow = gate + " && abs(ThetaLab-117.5)<12.5";
  string gateHigh = gate + " && abs(ThetaLab-142.5)<12.5";

  TCanvas *diagnosis_ExTheta = new TCanvas("diagnosis_ExTheta","diagnosis_ExTheta",1000,1000);
  chain->Draw("Ex>>thetaHistLow(200,-1,9)", gateLow.c_str(),"");
  TH2F* thetaHistLow = (TH2F*) gDirectory->Get("thetaHistLow");  
  thetaHistLow->SetLineColor(kBlue);

  chain->Draw("Ex>>thetaHistHigh(200,-1,9)", gateHigh.c_str(),"");
  TH2F* thetaHistHigh = (TH2F*) gDirectory->Get("thetaHistHigh");  
  thetaHistHigh->SetLineColor(kRed);

  thetaHistLow->Draw();
  thetaHistHigh->Draw("SAME");
}

void ExMugast_ForPoster(){

  TCanvas *forPoster = new TCanvas("forPoster","forPoster",1000,1000);
  gStyle->SetOptStat(0);
  forPoster->Divide(1,2);
  forPoster->cd(1);

  chain->Draw("Ex>>hcG_T1(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==1","same");
  chain->Draw("Ex>>hcG_T2(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==2","same");
  chain->Draw("Ex>>hcG_T5(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==5","same");
  chain->Draw("Ex>>hcG_T7(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==7","same");
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


  tree->Draw("Ex>>h_T1(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==1","same");
  tree->Draw("Ex>>h_T2(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==2","same");
  tree->Draw("Ex>>h_T5(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==5","same");
  tree->Draw("Ex>>h_T7(120,-1,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber==7","same");
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
  TF1* func = f_efficAGATA();
  func->Draw();
  cout << "At E = " << Energy_keV 
       << " keV, AGATA efficiency = " << func->Eval(Energy_keV)
       << " %" << endl;
}

void ElasticsGate(double EMin, double EMax){
  double width = EMax-EMin;
  double centre = EMin+(0.5*width);
//  string gates = "abs(T_MUGAST_VAMOS-2700)<400 && MUST2.TelescopeNumber==5 && abs(ELab - " 
  string gate = timegate + "&&" + det_gate 
               + "&& abs(ELab - "
	       + to_string(centre) 
	       + ")< " 
	       + to_string(width);

  chain->Draw("ThetaLab>>hist(80,50,90)", gate.c_str(), "");
}

void GateThetaCM(double minTheta, double maxTheta, double binsize){
  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ThetaCM > " 
      + to_string(minTheta)
      + " && ThetaCM < "
      + to_string(maxTheta);

  string title = to_string(minTheta)+" < ThetaCM < "+to_string(maxTheta);
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>Ex_ThetaCMGate(" + to_string(10.0/binsize) + ",-1,9)";

  TCanvas *cEx_ThetaCMGate = new TCanvas("cEx_ThetaCMGate","cEx_ThetaCMGate",1000,1000);
  chain->Draw(draw.c_str(),gating.c_str(),"colz");
  TH1F* Ex_ThetaCMGate = (TH1F*) gDirectory->Get("Ex_ThetaCMGate");
  Ex_ThetaCMGate->GetXaxis()->SetTitle("Ex [MeV]");
  Ex_ThetaCMGate->GetYaxis()->SetTitle(ytitle.c_str());
  Ex_ThetaCMGate->SetTitle(title.c_str());

  DrawParticleStates(cEx_ThetaCMGate);
}

void GateThetaLab(double minTheta, double maxTheta, double binsize){
  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ThetaLab > " 
      + to_string(minTheta)
      + " && ThetaLab < "
      + to_string(maxTheta);

  string title = to_string(minTheta)+" < ThetaLab < "+to_string(maxTheta);
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>Ex_ThetaLabGate(" + to_string(10.0/binsize) + ",-1,9)";

  TCanvas *cEx_ThetaLabGate = new TCanvas("cEx_ThetaLabGate","cEx_ThetaLabGate",1000,1000);
  chain->Draw(draw.c_str(),gating.c_str(),"colz");
  TH1F* Ex_ThetaLabGate = (TH1F*) gDirectory->Get("Ex_ThetaLabGate");
  Ex_ThetaLabGate->GetXaxis()->SetTitle("Ex [MeV]");
  Ex_ThetaLabGate->GetYaxis()->SetTitle(ytitle.c_str());
  Ex_ThetaLabGate->SetTitle(title.c_str());

  DrawParticleStates(cEx_ThetaLabGate);
}

void GateThetaLab_AllOverlaid(){
  double binsize = 0.1;
  string basegate = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ";
      //+ to_string(minTheta)
      //+ " && ThetaLab < "
      //+ to_string(maxTheta);

  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>Ex_ThetaLabGate(" + to_string(10.0/binsize) + ",-1,9)";

  TCanvas *cThetaLabGates = new TCanvas("cThetaLabGates","cThetaLabGates",1000,1000);
  
  /* 105 to 110 */
  
  for(int i=0; i<9;i++){
    int min = 105+(i*5);
    int max = 110+(i*5);


    string gate = basegate + "ThetaLab > " + to_string(min) + " && ThetaLab < " + to_string(max);
    string histname = "Gate" + to_string(min) + "to" + to_string(max); 
    string draw = "Ex>>" + histname + "(600,-15,15)";

    
    chain->Draw(draw.c_str(),gate.c_str(),"same");
    TH1F* hist = (TH1F*) gDirectory->Get(histname.c_str());
    hist->SetLineColor(i+1);

  }

}

void GateThetaLab_MultiWrite(double startTheta, double finishTheta, int numGates, double binsize){
  string core = timegate 
	      + " && " + det_gate;
  if(reactionName=="47K(d,t)"){
    core = core + " && cutTritons && cutTime";
  }
  core = core + " && Ex@.size()==1 && ThetaLab > ";

  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  double gatesize = (finishTheta-startTheta)/numGates;
  TList* list = new TList();

  for (int i=0; i<numGates; i++){
    cout << GREEN << "Writing gate " << i+1 << "/" << numGates << RESET << endl;
    double minTheta = startTheta + (i * gatesize);
    string title = to_string((int) minTheta)+" < ThetaLab < "+to_string((int) (minTheta+gatesize));
    string gating = core
        + to_string(minTheta)
        + " && ThetaLab < "
        + to_string(minTheta+gatesize);
    string histname = "cThetaLabGate_" + to_string((int) minTheta) + "-" + to_string((int) (minTheta+gatesize));
    string draw = "Ex>>" + histname + "(" + to_string(30.0/binsize) + ",-15,15)";

    TCanvas *cEx_ThetaLabGate = new TCanvas(histname.c_str(),histname.c_str(),1000,1000);
    chain->Draw(draw.c_str(),gating.c_str(),"colz");
    TH1F* Ex_ThetaLabGate = (TH1F*) gDirectory->Get(histname.c_str());
    Ex_ThetaLabGate->GetXaxis()->SetTitle("Ex [MeV]");
    Ex_ThetaLabGate->GetYaxis()->SetTitle(ytitle.c_str());
    Ex_ThetaLabGate->Sumw2();
    Ex_ThetaLabGate->SetTitle(title.c_str());
    list->Add(Ex_ThetaLabGate);
    delete cEx_ThetaLabGate;
  }

  TFile* file = new TFile("GateThetaLabHistograms.root","RECREATE");
  list->Write("GateThetaLabHistograms",TObject::kSingleKey);
  file->ls();
}

void GateThetaCM_MultiWrite(double startTheta, double finishTheta, int numGates, double binsize){
  string core = timegate 
	      + " && " + det_gate;
  if(reactionName=="47K(d,t)"){
    core = core + " && cutTritons && cutTime";
  }
  core = core + " && Ex@.size()==1 && ThetaCM > ";

  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  double gatesize = (finishTheta-startTheta)/numGates;
  TList* list = new TList();

  for (int i=0; i<numGates; i++){
    cout << GREEN << "Writing gate " << i+1 << "/" << numGates << RESET << endl;
    double minTheta = startTheta + (i * gatesize);
    string title = to_string((int) minTheta)+" < ThetaCM < "+to_string((int) (minTheta+gatesize));
    string gating = core
        + to_string(minTheta)
        + " && ThetaCM < "
        + to_string(minTheta+gatesize);
    string histname = "cThetaCMGate_" + to_string((int) minTheta) + "-" + to_string((int) (minTheta+gatesize));
    string draw = "Ex>>" + histname + "(" + to_string(30.0/binsize) + ",-15,15)";

    TCanvas *cEx_ThetaCMGate = new TCanvas(histname.c_str(),histname.c_str(),1000,1000);
    chain->Draw(draw.c_str(),gating.c_str(),"colz");
    TH1F* Ex_ThetaCMGate = (TH1F*) gDirectory->Get(histname.c_str());
    Ex_ThetaCMGate->GetXaxis()->SetTitle("Ex [MeV]");
    Ex_ThetaCMGate->GetYaxis()->SetTitle(ytitle.c_str());
    Ex_ThetaCMGate->Sumw2();
    Ex_ThetaCMGate->SetTitle(title.c_str());
    list->Add(Ex_ThetaCMGate);
    delete cEx_ThetaCMGate;
  }

  TFile* file = new TFile("GateThetaCMHistograms.root","RECREATE");
  list->Write("GateThetaCMHistograms",TObject::kSingleKey);
  file->ls();
}



void GateThetaLab_MultiWrite(double startTheta, double finishTheta, int numGates, double binsize, int MGX){
   string core = timegate 
	      + " && " + det_gate
	      + " && Mugast.TelescopeNumber==" + to_string(MGX)
              + " && Ex@.size()==1 && ThetaLab > ";
// string core = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ThetaLab > ";
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  double gatesize = (finishTheta-startTheta)/numGates;
  TList* list = new TList();

  for (int i=0; i<numGates; i++){
    double minTheta = startTheta + (i * gatesize);
    string title = to_string((int) minTheta)+" < ThetaLab < "+to_string((int) (minTheta+gatesize));
    string gating = core
        + to_string(minTheta)
        + " && ThetaLab < "
        + to_string(minTheta+gatesize);
    string histname = "cThetaLabGate_" + to_string((int) minTheta) + "-" + to_string((int) (minTheta+gatesize));
    string draw = "Ex>>" + histname + "(" + to_string(30.0/binsize) + ",-15,15)";

    TCanvas *cEx_ThetaLabGate = new TCanvas(histname.c_str(),histname.c_str(),1000,1000);
    chain->Draw(draw.c_str(),gating.c_str(),"colz");
    TH1F* Ex_ThetaLabGate = (TH1F*) gDirectory->Get(histname.c_str());
    Ex_ThetaLabGate->GetXaxis()->SetTitle("Ex [MeV]");
    Ex_ThetaLabGate->GetYaxis()->SetTitle(ytitle.c_str());
    Ex_ThetaLabGate->Sumw2();
    Ex_ThetaLabGate->SetTitle(title.c_str());
    list->Add(Ex_ThetaLabGate);
    delete cEx_ThetaLabGate;
  }

  TFile* file = new TFile("GateThetaLabHistograms_MGX.root","RECREATE");
  list->Write("GateThetaLabHistograms_MGX",TObject::kSingleKey);
  file->ls();
}

void GateThetaLab_MultiWrite(double startTheta, double finishTheta, int numGates, double binsize, double gateGammaE, double gateGammaWdth){
    string core = timegate 
	      + " && " + det_gate
	      + " && abs(AddBack_EDC-" + to_string(gateGammaE) + ")<" + to_string(gateGammaWdth)
              + " && Ex@.size()==1 && ThetaLab > ";
// string core = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && abs(AddBack_EDC-" + to_string(gateGammaE) + ")<" + to_string(gateGammaWdth) + " && Mugast.TelescopeNumber<8 && ThetaLab > ";
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  double gatesize = (finishTheta-startTheta)/numGates;
  TList* list = new TList();

  for (int i=0; i<numGates; i++){
    double minTheta = startTheta + (i * gatesize);
    string title = to_string((int) minTheta)+" < ThetaLab < "+to_string((int) (minTheta+gatesize));
    string gating = core
        + to_string(minTheta)
        + " && ThetaLab < "
        + to_string(minTheta+gatesize);
    string histname = "cThetaLabGate_" + to_string((int) minTheta) + "-" + to_string((int) (minTheta+gatesize));
    string draw = "Ex>>" + histname + "(" + to_string(30.0/binsize) + ",-15,15)";

    TCanvas *cEx_ThetaLabGate = new TCanvas(histname.c_str(),histname.c_str(),1000,1000);
    chain->Draw(draw.c_str(),gating.c_str(),"colz");
    TH1F* Ex_ThetaLabGate = (TH1F*) gDirectory->Get(histname.c_str());
    Ex_ThetaLabGate->GetXaxis()->SetTitle("Ex [MeV]");
    Ex_ThetaLabGate->GetYaxis()->SetTitle(ytitle.c_str());
    Ex_ThetaLabGate->Sumw2();
    Ex_ThetaLabGate->SetTitle(title.c_str());
    list->Add(Ex_ThetaLabGate);
    delete cEx_ThetaLabGate;
  }

  TFile* file = new TFile("GateThetaLabHistograms_GammaGate.root","RECREATE");
  list->Write("GateThetaLabHistograms",TObject::kSingleKey);
  file->ls();
}

void CompareThetaLabGatesFromFile(){
  //TFile* oldF = new TFile("GateThetaLabHistograms_11Apr22_20angles.root","READ");
  //TFile* oldF = new TFile("GateThetaLabHistograms_11Jul22.root","READ");
  //TFile* oldF = new TFile("GateThetaLabHistograms_29Aug22_TrueStripRemoval_0p05.root","READ");
  TFile* oldF = new TFile("GateThetaLabHistograms_14Oct22_bin0p05.root","READ");
  TList* oldL = (TList*) oldF->FindObjectAny("GateThetaLabHistograms");
  
  //TFile* newF = new TFile("GateThetaLabHistograms_11Jul22.root","READ");
  //TFile* newF = new TFile("GateThetaLabHistograms_Test2um.root","READ");
  //TFile* newF = new TFile("GateThetaLabHistograms_30Aug22_2p06um.root","READ");
  TFile* newF = new TFile("GateThetaLabHistograms_14Oct22_2_bin0p05.root","READ");
  TList* newL = (TList*) newF->FindObjectAny("GateThetaLabHistograms");

  double minTheta=105.;
  double gatesize=2.5;
  int numGates=20;

  TCanvas* cOldNew = new TCanvas("cOldNew","cOldNew",2000,1000);
  gStyle->SetOptStat(0);
  cOldNew->Divide(3,2,0.01,0.01,0);
  /*SetPad(xlow, ylow, xupp, yupp) */
  cOldNew->cd(1)->SetPad(0., 0.5, 0.5, 1.0);
  cOldNew->cd(2)->SetPad(0.5, 0.5, 0.75, 1.0 );
  cOldNew->cd(3)->SetPad(0.75, 0.5, 1, 1.0 );
  cOldNew->cd(4)->SetPad(0., 0., 0.5, 0.5);
  cOldNew->cd(5)->SetPad(0.5, 0., 0.75, 0.5 );
  cOldNew->cd(6)->SetPad(0.75, 0., 1, 0.5 );

  for (int i=0; i<numGates; i++){ 
    string hname = "cThetaLabGate_" 
	         + to_string((int) (minTheta+(i*gatesize))) 
	         + "-" 
	         + to_string((int) (minTheta+((i+1.)*gatesize)));

    cout << " Looking for " << hname << endl;
    
    TH1F* oldH = (TH1F*) oldL->FindObject(hname.c_str());
    TH1F* newH = (TH1F*) newL->FindObject(hname.c_str());

    //would be faster to set some of these on just the last one but speed not issue here
    cOldNew->cd(1);
    oldH->SetLineColor(30+i);
    //oldH->Rebin(2);
    oldH->GetYaxis()->SetTitle("Counts / 0.100000 MeV");
    oldH->GetXaxis()->SetRangeUser(-1,7);
    oldH->GetYaxis()->SetRangeUser(0,400);
    //oldH->SetTitle("08Nov22: 1.3008 um CD2");
    //oldH->SetTitle("11Jul22: 2.8798 um CD2");
    oldH->SetTitle("14Oct22: 2.8798 um CD2");
    oldH->Draw("HIST SAME");

    cOldNew->cd(2);
    TH1F* oldH1 = (TH1F*) oldH->Clone();
    oldH1->SetName("oldH1");
    oldH1->GetXaxis()->SetRangeUser(-1,1);
    oldH1->SetTitle("");
    oldH1->Draw("HIST SAME");
 
    cOldNew->cd(3);
    TH1F* oldH2 = (TH1F*) oldH->Clone();
    oldH2->GetYaxis()->SetRangeUser(0,250);
    oldH2->SetName("oldH2");
    oldH2->GetXaxis()->SetRangeUser(1.5,2.5);
    oldH2->SetTitle("");
    oldH2->Draw("HIST SAME");
   
    cOldNew->cd(4);
    newH->GetXaxis()->SetRangeUser(-1,7);
    newH->GetYaxis()->SetRangeUser(0,400);
    newH->SetLineColor(30+i);
    //newH->SetTitle("11Jul22: 2.8798 um CD2");
    //newH->SetTitle("AugTest: 2.000 um CD2");
    newH->SetTitle("14Oct22_2: 2.9992 um CD2");
    newH->Draw("HIST SAME");

    cOldNew->cd(5);
    TH1F* newH1 = (TH1F*) newH->Clone();
    newH1->SetName("newH1");
    newH1->GetXaxis()->SetRangeUser(-1,1);
    newH1->SetTitle("");
    newH1->Draw("HIST SAME");
 
    cOldNew->cd(6);
    TH1F* newH2 = (TH1F*) newH->Clone();
    newH2->GetYaxis()->SetRangeUser(0,250);
    newH2->SetName("newH2");
    newH2->GetXaxis()->SetRangeUser(1.5,2.5);
    newH2->SetTitle("");
    newH2->Draw("HIST SAME");

  }


}

void CompareTritonsFomFile(){
  string gate = timegate 
	      + " && " + det_gate
              + " && cutTritons && cutTime";
  string gate2 = gate + " && Ex@.size()==1 && abs(AddBack_EDC-1.94)<0.05";

  /* Push MM1-4 +190mm in Z */
  vector<string> f190;
  f190.push_back("../../../Outputs/Analysis/47Kdt_17Aug22_MM+190mm_PartI.root");
  f190.push_back("../../../Outputs/Analysis/47Kdt_17Aug22_MM+190mm_PartII.root");
  TChain* c190 = Chain("PhysicsTree",f190,true);

  /* Push MM1-4 +200mm in Z */
  vector<string> f200;
  f200.push_back("../../../Outputs/Analysis/47Kdt_17Aug22_MM+200mm_PartI.root");
  f200.push_back("../../../Outputs/Analysis/47Kdt_17Aug22_MM+200mm_PartII.root");
  TChain* c200 = Chain("PhysicsTree",f200,true);

  /* Push MM1-4 +210mm in Z */
  vector<string> f210;
  f210.push_back("../../../Outputs/Analysis/47Kdt_17Aug22_MM+210mm_PartI.root");
  f210.push_back("../../../Outputs/Analysis/47Kdt_17Aug22_MM+210mm_PartII.root");
  TChain* c210 = Chain("PhysicsTree",f210,true);

  auto cTritons = new TCanvas("cTritons","cTritons",1000,1000);
  
  c190->Draw("Ex>>t190(200,-2,8)", gate2.c_str(),"");
  TH1F* t190 = (TH1F*) gDirectory->Get("t190");
  t190->GetXaxis()->SetTitle("Ex [MeV]");
  t190->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  t190->SetLineColor(kRed); 

  c200->Draw("Ex>>t200(200,-2,8)", gate2.c_str(),"");
  TH1F* t200 = (TH1F*) gDirectory->Get("t200");
  t200->SetLineColor(kBlue); 

  c210->Draw("Ex>>t210(200,-2,8)", gate2.c_str(),"");
  TH1F* t210 = (TH1F*) gDirectory->Get("t210");
  t210->SetLineColor(kGreen);

  t190->Draw();
  t200->Draw("same");
  t210->Draw("same");

  auto cTest = new TCanvas("cTest","cTest",1000,1000);
  cTest->Divide(3);

  cTest->cd(1);
  c190->Draw("ELab:ThetaLab>>h190(90,0,45,500,0,10)",gate.c_str(),"colz");
  TH2F* h190 = (TH2F*) gDirectory->Get("h190");
  plot_kine(Kdt, 0.000, kBlack, 2, 1);
  plot_kine(Kdt, 1.944, kBlack, 2, 5);
  plot_kine(Kdt, 3.290, kBlack, 2, 5);

  cTest->cd(2);
  c200->Draw("ELab:ThetaLab>>h200(90,0,45,500,0,10)",gate.c_str(),"colz");
  TH2F* h200 = (TH2F*) gDirectory->Get("h200");
  plot_kine(Kdt, 0.000, kBlack, 2, 1);
  plot_kine(Kdt, 1.944, kBlack, 2, 5);
  plot_kine(Kdt, 3.290, kBlack, 2, 5);

  cTest->cd(3);
  c210->Draw("ELab:ThetaLab>>h210(90,0,45,500,0,10)",gate.c_str(),"colz");
  TH2F* h210 = (TH2F*) gDirectory->Get("h210");
  plot_kine(Kdt, 0.000, kBlack, 2, 1);
  plot_kine(Kdt, 1.944, kBlack, 2, 5);
  plot_kine(Kdt, 3.290, kBlack, 2, 5);


}

void GatePhaseSpaceByThetaLab_MultiWrite(double startTheta, double finishTheta, int numGates, double binsize){
  string core = "EventWeight*(Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ThetaLab > ";
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  double gatesize = (finishTheta-startTheta)/numGates;
  TList* list = new TList();

  //TFile* psfile = new TFile("../../../Outputs/Analysis/Sim_02Mar_47Kdp_PhaseSpace.root","READ");
  TFile* psfile = new TFile("../../../Outputs/Analysis/Sim_PhaseSpace_11Jul22.root","READ");
  TTree* PSTree = (TTree*) psfile->FindObjectAny("PhysicsTree");

  for (int i=0; i<numGates; i++){
    double minTheta = startTheta + (i * gatesize);
    string title = to_string((int) minTheta)+" < ThetaLab < "+to_string((int) (minTheta+gatesize));
    string gating = core
        + to_string(minTheta)
        + " && ThetaLab < "
        + to_string(minTheta+gatesize)
	+ ")";
    string histname = "cPSpaceThetaLabGate_" + to_string((int) minTheta) + "-" + to_string((int) (minTheta+gatesize));
    string draw = "Ex>>" + histname + "(" + to_string(30.0/binsize) + ",-15,15)";

    TCanvas *cPSpace_ThetaLabGate = new TCanvas(histname.c_str(),histname.c_str(),1000,1000);
    PSTree->Draw(draw.c_str(),gating.c_str(),"colz");
    TH1F* PSpace_ThetaLabGate = (TH1F*) gDirectory->Get(histname.c_str());
    PSpace_ThetaLabGate->GetXaxis()->SetTitle("Ex [MeV]");
    PSpace_ThetaLabGate->GetYaxis()->SetTitle(ytitle.c_str());
    PSpace_ThetaLabGate->Sumw2();
    PSpace_ThetaLabGate->SetTitle(title.c_str());
    list->Add(PSpace_ThetaLabGate);
    delete cPSpace_ThetaLabGate;
  }

  TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms.root","RECREATE");
  list->Write("GatePhaseSpaceThetaLabHistograms",TObject::kSingleKey);
  file->ls();
}

void GammaSub_NoDoppler(){
  TCanvas *cGammaSubNoDopp = new TCanvas("cGammaSubNoDopp","cGammaSubNoDopp",1000,1000);
  chain->Draw("AGATA_GammaE>>hNoDoppler(5000,0,5)","abs(T_MUGAST_VAMOS-7000)<3000","");
  TH1F* hNoDoppler = (TH1F*) gDirectory->Get("hNoDoppler");
  hNoDoppler->SetTitle("AGATA, Timing gate 4k-10k, no doppler");
}

void GammaSub_WithDoppler(){
  TCanvas *cGammaSubWithDopp = new TCanvas("cGammaSubWithDopp","cGammaSubWithDopp",1000,1000);
  chain->Draw("AddBack_EDC>>hWithDoppler(5000,0,5)","abs(T_MUGAST_VAMOS-7000)<3000","");
  TH1F* hWithDoppler = (TH1F*) gDirectory->Get("hWithDoppler");
  hWithDoppler->SetTitle("AGATA, Timing gate 4k-10k, with doppler");
}

/*
void GammaSub_Subbed(){
  TCanvas *cGammaSub = new TCanvas("cGammaSub","cGammaSub",1000,1000);
  chain->Draw("AddBack_EDC>>Eg_Sub(5000,0,5)","abs(T_MUGAST_VAMOS-2700)<400","");
  TH1F* Eg_Sub = (TH1F*) gDirectory->Get("Eg_Sub");
  chain->Draw("AddBack_EDC>>Eg_True(5000,0,5)","abs(T_MUGAST_VAMOS-2700)<400","");
  TH1F* Eg_True = (TH1F*) gDirectory->Get("Eg_True");
  //chain->Draw("AddBack_EDC>>Eg_False(5000,0,5)","abs(T_MUGAST_VAMOS-1250)<850 || abs(T_MUGAST_VAMOS-4650)<850");
  chain->Draw("AddBack_EDC>>Eg_False(5000,0,5)",
		  "abs(T_MUGAST_VAMOS-4050)<250 || abs(T_MUGAST_VAMOS-4900)<200 || abs(T_MUGAST_VAMOS-650)<250 || abs(T_MUGAST_VAMOS-1550)<250");
  TH1F* Eg_False = (TH1F*) gDirectory->Get("Eg_False");

  double scale = (400.*2.)/(500.+400.+500.+500.);
  Eg_False->Scale(scale);
  Eg_Sub->Add(Eg_False,-1);
  Eg_Sub->Draw();
}
*/

void GammaSub_Actual_ExcludeBeamDecay(){
  TCanvas *cGammaSub = new TCanvas("cGammaSub","cGammaSub",1000,1000);
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)",
		  "abs(T_MUGAST_VAMOS-2700)<400 && abs(AGATA_GammaE-2.013)>0.004 && abs(AGATA_GammaE-0.511)>0.003 && abs(AGATA_GammaE-0.564)>0.004 && abs(AGATA_GammaE-0.586)>0.003");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->Draw();
}

/*
void gg(){
  static vector<double> *AddBack_EDC;
  auto Gamma_Branch = chain->GetBranch("AddBack_EDC");
  Gamma_Branch->SetAddress(&AddBack_EDC);

  auto h=new TH2F("gg","gg",1000,0,10,1000,0,10);

//cout << "here1" << endl;

  unsigned int numEntries = chain->GetEntries();
//cout << numEntries << endl;
  for(unsigned int i=0; i<numEntries; i++){
//cout << "here2" << endl;
    chain->GetEntry(i);
//cout << "here3" << endl;
    // Example filling of GG matrix 
    unsigned int size = AddBack_EDC->size();
//cout << size << endl;
    if(size>0){
	    double e1,e2;
      for(unsigned int s = 0 ; s < size-1 ; s++){
//cout << "here4, s = " << s << " & size = " << size << endl;
        e1=AddBack_EDC->at(s) ; e2 = AddBack_EDC->at(s+1);
//cout << "here5" << endl;
       // Folding of the matrix, always fill big first
       if(e1>e2)
         h->Fill(e1,e2);
       else
         h->Fill(e2,e1);
      }
    }
//cout << "here6" << endl;
  }
  h->Draw();
}
*/

void ggLoad(TTree* chain, TH2F* h){

   // Initilise the Mugast branch
   auto Mugast = new TMugastPhysics();
 
   // Initilise access variables for chain
   static double T_MUGAST_VAMOS;
   static vector<double> *AddBack_EDC, *AGATA_GammaE;

   // Pull chain branches
   auto Gamma_Branch = chain->GetBranch("AddBack_EDC");
   auto RawGamma_Branch = chain->GetBranch("AGATA_GammaE");
   auto MugVam_Branch = chain->GetBranch("T_MUGAST_VAMOS");

   // Set chain variable addresses
   Gamma_Branch->SetAddress(&AddBack_EDC);
   RawGamma_Branch->SetAddress(&AGATA_GammaE);
   MugVam_Branch->SetAddress(&T_MUGAST_VAMOS);
 
   // Build loop variables
   unsigned int numEntries = chain->GetEntries();
   unsigned int multiplicity = 0;
 
   // Loop on entries
   for(unsigned int i=0; i<numEntries; i++){
     chain->GetEntry(i);
       // Gate on Timing
       if(abs(T_MUGAST_VAMOS-2700)<400){
         int gammaMultip = AddBack_EDC->size();
	 //no muliplicity 1
         if(gammaMultip>=1){
	   double e1,e2;
	   //loop through events
           for(unsigned int s=0 ; s<gammaMultip-1 ; s++){
	     //remove beam decay gammas
             if(abs(AGATA_GammaE->at(s)-2.013)>0.004 && abs(AGATA_GammaE->at(s)-0.511)>0.003 
	     && abs(AGATA_GammaE->at(s)-0.564)>0.004 && abs(AGATA_GammaE->at(s)-0.586)>0.003){
               e1=AddBack_EDC->at(s); e2 = AddBack_EDC->at(s+1);
               // Folding of the matrix, always fill big first
               if(e1>e2){
                 h->Fill(e1,e2);
	       }
	       else{
                 h->Fill(e2,e1);
	       }
	     }
           }
         }//if gamma
      }//timing
   }//for i
}

void gggLoad(TTree* chain, THnSparseF* h){

cout << "THIS IS OLD!!!! UPDATE WITH THE BEAM EXCLUSION!!!" << endl;
	// Initilise the Mugast branch
   auto Mugast = new TMugastPhysics();
 
   // Initilise access variables for chain
//   static double T_MUGAST_VAMOS;
   static vector<double> //*X, *Y, *Z, *RawEnergy, 
	   *AddBack_EDC;
 
   // Pull chain branches
//   auto Energy_Branch = chain->GetBranch("RawEnergy");
   auto Gamma_Branch = chain->GetBranch("AddBack_EDC");
//   auto X_Branch = chain->GetBranch("X");
//   auto Y_Branch = chain->GetBranch("Y");
//   auto Z_Branch = chain->GetBranch("Z");
//   auto MugVam_Branch = chain->GetBranch("T_MUGAST_VAMOS");
 
   // Set Mugast branch address
//   chain->SetBranchAddress("Mugast",&Mugast);
 
   // Set chain variable addresses
//   Energy_Branch->SetAddress(&RawEnergy);
   Gamma_Branch->SetAddress(&AddBack_EDC);
//   X_Branch->SetAddress(&X);
//   Y_Branch->SetAddress(&Y);
//   Z_Branch->SetAddress(&Z);
//   MugVam_Branch->SetAddress(&T_MUGAST_VAMOS);
  
   // Build loop variables
   unsigned int numEntries = chain->GetEntries();
   unsigned int multiplicity = 0;
 
   // Loop on entries
   for(unsigned int i=0; i<numEntries; i++){
     chain->GetEntry(i);
//     multiplicity = Mugast->TelescopeNumber.size(); 

     // Loop on MUGAST multiplicity
//     for(int m=0; m<multiplicity; m++){
 
       // Gate on Timing
//       if(abs(T_MUGAST_VAMOS-2700)<400){
         int gammaMultip = AddBack_EDC->size();
         if(gammaMultip>=2){
	   double e1,e2,e3;
           for(unsigned int s=0; s<gammaMultip-2; s++){
             e1 = AddBack_EDC->at(s);
	     e2 = AddBack_EDC->at(s+1);
	     e3 = AddBack_EDC->at(s+2);

	     double arr[] = {e1, e2, e3};
             int n = sizeof(arr)/sizeof(arr[0]);
	     sort(arr, arr+n, greater<double>());

             //h->Fill(arr[0], arr[1], arr[2]);
             h->Fill(arr, 1.);

/*
             // Folding of the matrix, always fill big first
             if(e1>e2){
	       if(e2>e3){
		 //e1 > e2 > e3
                 h->Fill(e1,e2,e3);
	       }
	       else{
		 //e1 > e2 > e3
	       }
	     }
	     else{
               h->Fill(e2,e1,e3);
	     }
*/
           }
         }//if gamma
//       }//if timing
//     }//for m
   }//for i
}

void gg(){
 
   //cout << "LOADING FILES: 47Kdp_11Apr22_PartI & II" << endl;
   cout << "LOADING FILES: 47Kdp_11Jul22_PartI & II" << endl;

   auto h=new TH2F("gg","gg",1000,0,10,1000,0,10);
   //auto DataFile = new TFile("../../../Outputs/Analysis/47Kdp_11Apr22_PartI.root", "READ");
   auto DataFile = new TFile("../../../Outputs/Analysis/47Kdp_11Jul22_PartI.root", "READ");
   auto chain = (TTree*) DataFile->FindObjectAny("PhysicsTree");
   ggLoad(chain, h);

   auto h2=new TH2F("gg","gg",1000,0,10,1000,0,10);
   //auto DataFile2 = new TFile("../../../Outputs/Analysis/47Kdp_11Apr22_PartII.root", "READ");
   auto DataFile2 = new TFile("../../../Outputs/Analysis/47Kdp_11Jul22_PartII.root", "READ");
   auto chain2 = (TTree*) DataFile->FindObjectAny("PhysicsTree");
   ggLoad(chain2, h2);

   h->Add(h2,1);

   TFile* file = new TFile("GGMatrix.root","RECREATE");
   h->Write();
   file->Close();

   //h->Draw("colz");
}

void ggg(){
   int bins[3] = {1000,1000,1000};
   double min[3] = {0.,0.,0.};
   double max[3] = {10.,10.,10.};

   //auto h3d=new TH3F("ggg","ggg",1000,0,10,1000,0,10,1000,0,10);
   auto h3d=new THnSparseF("hggg","hggg",3,bins,min,max);
   auto DataFile = new TFile("../../../Outputs/Analysis/47Kdp_08Nov_PartI.root", "READ");
   auto chain = (TTree*) DataFile->FindObjectAny("PhysicsTree");
 
   gggLoad(chain, h3d);

   //auto h3d2=new TH3F("gg","gg",1000,0,10,1000,0,10,1000,0,10);
   auto h3d2=new THnSparseF("hggg","hggg",3,bins,min,max);
   auto DataFile2 = new TFile("../../../Outputs/Analysis/47Kdp_08Nov_PartII.root", "READ");
   auto chain2 = (TTree*) DataFile->FindObjectAny("PhysicsTree");

   gggLoad(chain2, h3d2);

   h3d->Add(h3d2,1);
   TFile* file = new TFile("GGGMatrix.root","RECREATE");
   h3d->Write();
   file->Close();

//   h3d->Draw();
   double gate = 0.10;
   int xmin = (int)((0.66-gate)*1000.);
   int xmax = (int)((0.66+gate)*1000.);
   int ymin = (int)((2.28-gate)*1000.);
   int ymax = (int)((2.28+gate)*1000.);
   cout << "GATING X AXIS FROM " << 0.66-gate << " - " << 0.66+gate << " -> bins " << xmin << " to " << xmax << endl;
   cout << "GATING Y AXIS FROM " << 2.28-gate << " - " << 2.28+gate << " -> bins " << ymin << " to " << ymax << endl;
   
   h3d->GetAxis(0)->SetRange(xmin,xmax);
   h3d->GetAxis(1)->SetRange(ymin,ymax);
   TH1D* projZ = h3d->Projection(2);
   projZ->Draw();
   //h->SaveAs("Save3Dgammas.root");
}

void gggGater(THnSparseF* h3d, double xE, double xgate, double yE, double ygate){
  int xmin = (int)((xE-xgate)*1000.);
  int xmax = (int)((xE+xgate)*1000.);
  int ymin = (int)((yE-ygate)*1000.);
  int ymax = (int)((yE+ygate)*1000.);

  cout << "GATING X AXIS FROM " << xE-xgate << " - " << xE+xgate << " -> bins " << xmin << " to " << xmax << endl;
  cout << "GATING Y AXIS FROM " << yE-ygate << " - " << yE+ygate << " -> bins " << ymin << " to " << ymax << endl;

  h3d->GetAxis(0)->SetRange(xmin,xmax);
  h3d->GetAxis(1)->SetRange(ymin,ymax);
  TH1D* projZ = h3d->Projection(2);
  projZ->Draw();

}

void ggGater(TH2F* h, double E, double gate){
  int binmin = (int)((E-gate)*100.)+1;//h->GetXaxis()->GetBin(E-gate);
  int binmax = (int)((E+gate)*100.)+1;//h->GetXaxis()->GetBin(E+gate);

  TH1D* h1 = h->ProjectionX("_px",binmin,binmax);
  TH1D* h2 = h->ProjectionY("_py",binmin,binmax);

  h1->SetTitle("gg Gate, ASSUMING 1000 bins from 0 to 10 MeV");
  h1->Add(h2,1);
  h1->Draw();


}

void Figure_Eg_MG(){
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
  Eg->SetTitle("#gamma-ray spectrum, requiring MUGAST upstream coincidence");

  TH1F* Eg2 = (TH1F*) Eg->Clone();
  Eg2->SetTitle("");

  TCanvas* canv = new TCanvas("canv","canv",1000,1000);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetOptStat(0);


  TPad *pad1 = new TPad("pad1","pad1",0,0.5,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.5);
  pad1->SetTopMargin(0.08);
  pad1->SetBottomMargin(0.08);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.08);
  pad2->SetBottomMargin(0.08);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  Eg->GetXaxis()->SetRangeUser(0.,2.);
  Eg->Draw();


  pad2->cd();
 

  Eg2->Rebin(2);
  Eg2->GetYaxis()->SetTitle("Counts / 0.001 MeV");
  Eg2->GetXaxis()->SetRangeUser(2.,4.5);
  Eg2->Draw();


}

void Figure_GateGamma_SeeParticle(double gamma, double width, double bg, double widthbg, double gammaBinWidth, double particleBinWidth){
  string gating = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = timegate + "&&" + det_gate + " && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(widthbg);

  string gammagate = timegate 
	      + " && " + det_gate
	      + " && " + exclBmDcy;

  double ratio = width/widthbg;

  ////////////////////////////////////////////////

  TCanvas *cFig_GGSP = new TCanvas("cFig_GGSP","cFig_GGSP",1000,1000);
  gStyle->SetOptStat(0);
  cFig_GGSP->Divide(1,2,0.005,0.005,0);
  cFig_GGSP->cd(1); 
    gStyle->SetPadLeftMargin(0.10);
    gStyle->SetPadRightMargin(0.01);
    gPad->SetTickx();
    gPad->SetTicky();
  cFig_GGSP->cd(2); 
    gStyle->SetPadLeftMargin(0.10);
    gStyle->SetPadRightMargin(0.01);
    gPad->SetTickx();
    gPad->SetTicky();
 
  cFig_GGSP->cd(2);
  string draw1 = "Ex>>hEx(" + to_string(30./particleBinWidth) + ",-15,15)";
  chain->Draw(draw1.c_str(),gating.c_str(),"");
  TH1F* ExGate = (TH1F*) gDirectory->Get("hEx");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  string nameExYaxis = "Counts / " + to_string(particleBinWidth) + " MeV";
  ExGate->GetYaxis()->SetTitle(nameExYaxis.c_str());
  ExGate->SetLineColor(kGreen);
  ExGate->SetFillColor(kGreen);
  ExGate->SetFillStyle(3154);
  ExGate->SetTitle("");
  ExGate->GetXaxis()->SetTitleSize(0.05);
  ExGate->GetXaxis()->SetLabelSize(0.05);
  ExGate->GetYaxis()->SetTitleSize(0.05);
  ExGate->GetYaxis()->SetLabelSize(0.05);
  ExGate->GetYaxis()->SetTitleOffset(0.5);
  ExGate->GetYaxis()->SetNdivisions(520);
  ExGate->GetXaxis()->SetRangeUser(-1.,7.);

  string draw2 = "Ex>>hExBG(" + to_string(30./particleBinWidth) + ",-15,15)";
  chain->Draw(draw2.c_str(),bggate.c_str(),"same hist");
  TH1F* ExBG = (TH1F*) gDirectory->Get("hExBG");
  ExBG->Scale(ratio);
  ExBG->SetLineColor(kRed);
  ExBG->SetFillColor(kRed);
  ExBG->SetFillStyle(3345);
  ExBG->SetTitle("");
  ExBG->GetXaxis()->SetRangeUser(-1.,7.);
  ExBG->Draw("same hist");

  cFig_GGSP->Update();
  double maxEx = cFig_GGSP->cd(2)->GetUymax();
  TLine *Sn = new TLine(4.644, 0.0, 4.644, maxEx);
    Sn->SetLineColor(kBlack); Sn->SetLineStyle(7);
    Sn->Draw();
  TLine *gs = new TLine(0.000, 0.0, 0.000, maxEx);
    gs->SetLineColor(kBlack); gs->SetLineStyle(1);
    gs->Draw();

  cFig_GGSP->cd(1);
  string drawg = "AddBack_EDC>>hEg(" + to_string(5./gammaBinWidth) + ",0,5)";
  chain->Draw(drawg.c_str(),gammagate.c_str(),"");
  TH1F* hEg = (TH1F*) gDirectory->Get("hEg");
  hEg->SetLineColor(kBlack);
  string nameEgYaxis = "Counts / " + to_string(gammaBinWidth) + " MeV";
  hEg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg->GetYaxis()->SetTitle(nameEgYaxis.c_str());
  double zoomMin = min((double) floor((gamma-width)*10)/10., (double) floor((bg-widthbg)*10)/10.);
  double zoomMax = max((double) ceil((gamma+width)*10)/10.,  (double) ceil((bg+widthbg)*10)/10.);
  hEg->GetXaxis()->SetRangeUser(zoomMin,zoomMax);
  hEg->SetTitle("");
  hEg->GetXaxis()->SetTitleSize(0.05);
  hEg->GetXaxis()->SetLabelSize(0.05);
  hEg->GetYaxis()->SetTitleSize(0.05);
  hEg->GetYaxis()->SetLabelSize(0.05);
  hEg->GetYaxis()->SetTitleOffset(0.5);
  hEg->Draw();

  cFig_GGSP->Update();
  double max = cFig_GGSP->cd(1)->GetUymax();
  TLine* gateL = new TLine(gamma-width,0,gamma-width,max);
    gateL->SetLineColor(kGreen); gateL->SetLineStyle(10);
  TLine* gateH = new TLine(gamma+width,0,gamma+width,max);
    gateH->SetLineColor(kGreen); gateH->SetLineStyle(10);
  TLine* bgL = new TLine(bg-widthbg,0,bg-widthbg,max);
    bgL->SetLineColor(kRed); bgL->SetLineStyle(10);
  TLine* bgH = new TLine(bg+widthbg,0,bg+widthbg,max);
    bgH->SetLineColor(kRed); bgH->SetLineStyle(10);
  gateL->Draw("SAME");
  gateH->Draw("SAME");
  bgL->Draw("SAME");
  bgH->Draw("SAME");

}

void Figure_TopGamma_BottomParticle(double gammaBinWidth, double particleBinWidth){
  string gating = timegate + "&&" + det_gate + " && Ex@.size()==1";
  string gammagate = timegate 
	      + " && " + det_gate
	      + " && " + exclBmDcy;

  ////////////////////////////////////////////////

  TCanvas *cFig_GGSP = new TCanvas("cFig_GGSP","cFig_GGSP",1000,1000);
  gStyle->SetOptStat(0);
  cFig_GGSP->Divide(1,2,0.005,0.005,0);
  cFig_GGSP->cd(1); 
    gStyle->SetPadLeftMargin(0.10);
    gStyle->SetPadRightMargin(0.01);
    gPad->SetTickx();
    gPad->SetTicky();
  cFig_GGSP->cd(2); 
    gStyle->SetPadLeftMargin(0.10);
    gStyle->SetPadRightMargin(0.01);
    gPad->SetTickx();
    gPad->SetTicky();
 
  cFig_GGSP->cd(2);
  string draw1 = "Ex>>hEx(" + to_string(30./particleBinWidth) + ",-15,15)";
  chain->Draw(draw1.c_str(),gating.c_str(),"");
  TH1F* hEx = (TH1F*) gDirectory->Get("hEx");
  hEx->GetXaxis()->SetTitle("Ex [MeV]");
  string nameExYaxis = "Counts / " + to_string(particleBinWidth) + " MeV";
  hEx->GetYaxis()->SetTitle(nameExYaxis.c_str());
  //ExGate->SetLineColor(kGreen);
  //ExGate->SetFillColor(kGreen);
  //ExGate->SetFillStyle(3154);
  hEx->SetTitle("");
  hEx->GetXaxis()->SetTitleSize(0.05);
  hEx->GetXaxis()->SetLabelSize(0.05);
  hEx->GetYaxis()->SetTitleSize(0.05);
  hEx->GetYaxis()->SetLabelSize(0.05);
  hEx->GetYaxis()->SetTitleOffset(0.5);
  hEx->GetYaxis()->SetNdivisions(520);
  hEx->GetXaxis()->SetRangeUser(-1.,7.);
  hEx->Draw();
  FitKnownPeaks(hEx);

  cFig_GGSP->Update();
  double maxEx = cFig_GGSP->cd(2)->GetUymax();
  TLine *Sn = new TLine(4.644, 0.0, 4.644, maxEx);
    Sn->SetLineColor(kBlack); Sn->SetLineStyle(7);
    Sn->Draw();
  TLine *gs = new TLine(0.000, 0.0, 0.000, maxEx);
    gs->SetLineColor(kBlack); gs->SetLineStyle(1);
    gs->Draw();

  cFig_GGSP->cd(1);
  string drawg = "AddBack_EDC>>hEg(" + to_string(5./gammaBinWidth) + ",0,5)";
  chain->Draw(drawg.c_str(),gammagate.c_str(),"");
  TH1F* hEg = (TH1F*) gDirectory->Get("hEg");
  hEg->SetLineColor(kBlack);
  string nameEgYaxis = "Counts / " + to_string(gammaBinWidth) + " MeV";
  hEg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg->GetYaxis()->SetTitle(nameEgYaxis.c_str());
  hEg->SetTitle("");
  hEg->GetXaxis()->SetTitleSize(0.05);
  hEg->GetXaxis()->SetLabelSize(0.05);
  hEg->GetYaxis()->SetTitleSize(0.05);
  hEg->GetYaxis()->SetLabelSize(0.05);
  hEg->GetYaxis()->SetTitleOffset(0.5);
  hEg->Draw();

  cFig_GGSP->Update();
}

