#include <math.h>
#include "NPReaction.h"
#include <string>
#include <sstream>

using namespace std;

TChain* chain=NULL ;
char cond[1000];

NPL::Reaction Cadp("47Ca(d,p)48Ca@355");
NPL::Reaction Scdp("47Sc(d,p)48Sc@355");

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

static double tCentre = 2700.;
static double tRange =   400.;

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

void LoadChainNP(){
  vector<string> files;
  
  //files.push_back("../../../Outputs/Analysis/47Kdp_08Nov_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_08Nov_PartII.root");
  
  //files.push_back("../../../Outputs/Analysis/47Kdp_08Apr_PartI.root");
  //files.push_back("../../../Outputs/Analysis/47Kdp_08Apr_PartII.root");

  /* With thresholds, strip mathcing, and bad strips out */
  files.push_back("../../../Outputs/Analysis/47Kdp_11Apr22_PartI.root");
  files.push_back("../../../Outputs/Analysis/47Kdp_11Apr22_PartII.root");

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
  TLine *l3200 = new TLine(3.2, 0.0, 3.2, max);
    l3200->SetLineStyle(kDotted);
    l3200->Draw("same");
  TLine *l3605 = new TLine(3.605, 0.0, 3.605, max);
    l3605->SetLineStyle(kDotted);
    l3605->Draw("same");
  TLine *l3800 = new TLine(3.792, 0.0, 3.792, max);
    l3800->SetLineStyle(kDotted);
    l3800->Draw("same");
  TLine *l3870 = new TLine(3.876, 0.0, 3.876, max);
    l3870->SetLineStyle(kDotted);
    l3870->Draw("same");
  TLine *l4100 = new TLine(4.1, 0.0, 4.1, max);
    l4100->SetLineStyle(kDotted);
    l4100->Draw("same");
 TLine *l4510 = new TLine(4.51, 0.0, 4.51, max);
    l4510->SetLineStyle(kDotted);
    l4510->Draw("same");
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
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)","abs(T_MUGAST_VAMOS-2700)<400");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
}

void Load_1DGamma(){
  TH1F *hEg = new TH1F("hEg","Loaded 1D Gamma Spectrum",200,-1,9);
  TFile *file = new TFile("LoadHistograms/Load_1DGamma.root","READ");
  hEg = (TH1F*)file->Get("Eg");
  hEg->Draw();
}

void Draw_1DGamma_MG(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)",
    "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && abs(AGATA_GammaE-2.013)>0.004 && abs(AGATA_GammaE-0.511)>0.003 && abs(AGATA_GammaE-0.564)>0.004 && abs(AGATA_GammaE-0.586)>0.003");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
}

void Load_1DGamma_MG(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  TH1F *hEgMG = new TH1F("hEg","Loaded 1D Gamma Spectrum, MG gated",200,-1,9);
  TFile *file = new TFile("LoadHistograms/Load_1DGamma_MG.root","READ");
  hEgMG = (TH1F*)file->Get("Eg");
  hEgMG->Draw();
}

void Draw_1DGamma_MM(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  gStyle->SetOptStat(0);
  chain->Draw("AddBack_EDC>>Eg(5000,0,5)","abs(T_MUGAST_VAMOS-2700)<400 && MUST2.TelescopeNumber>0 && MUST2.TelescopeNumber<5");
  TH1F* Eg = (TH1F*) gDirectory->Get("Eg");
  Eg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  Eg->GetYaxis()->SetTitle("Counts / 0.001 MeV");
}

void Load_1DParticle(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  TH1F *hEx = new TH1F("hEx","Loaded 1D Particle Spectrum",200,-1,9);
  TFile *file = new TFile("LoadHistograms/Load_1DParticle.root","READ");
  hEx = (TH1F*)file->Get("Ep");
  hEx->Draw();
}

void Load_1DParticle_SubPhaseSpace(){
  TCanvas *cEg = new TCanvas("cEg","cEg",1000,1000);
  TH1F *hSubtracted = new TH1F("hSubtracted",
		  "Loaded 1D Particle Spectrum, Phase Spacae Subtracted (17Mar)",200,-1,9);
  TFile *file = new TFile("LoadHistograms/Load_1DParticle_SubPhaseSpace.root","READ");
  hSubtracted = (TH1F*)file->Get("subtracted");
  hSubtracted->Draw();
}

void Draw_1DParticle(){
  TCanvas *cEx = new TCanvas("cEx","cEx",1000,1000);
  chain->Draw("Ex>>Ep(200,-1,9)",
	"abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && Ex@.size()==1",
	"");
  TH1F* Ep = (TH1F*) gDirectory->Get("Ep");
//  Ep->SetTitle("Ex");
  Ep->GetXaxis()->SetTitle("Ex [MeV]");
  Ep->GetYaxis()->SetTitle("Counts / 0.05 MeV");

  DrawParticleStates(cEx);
}

void Draw_1DParticleMUST2(){
  TCanvas *cEx = new TCanvas("cEx","cEx",1000,1000);
  chain->Draw("Ex>>Ep(200,-1,9)",
	"abs(T_MUGAST_VAMOS-2700)<400 && MUST2.TelescopeNumber>0 && Ex@.size()==1 && MUST2.TelescopeNumber<4",
	"");
  TH1F* Ep = (TH1F*) gDirectory->Get("Ep");
//  Ep->SetTitle("Ex");
  Ep->GetXaxis()->SetTitle("Ex [MeV]");
  Ep->GetYaxis()->SetTitle("Counts / 0.05 MeV");
}

void Load_2DParticleGamma(){
  TCanvas *cExEg = new TCanvas("cExEg","cExEg",1000,1000);
  TH2F *hExEg = new TH2F("hExEg","Loaded 2D Particle-Gamma",200,-1,9,2500,0,5);
  TFile *file = new TFile("LoadHistograms/Load_2DParticleGamma.root","READ");
  hExEg = (TH2F*)file->Get("ExEg");
  hExEg->Draw("colz");
}

void Draw_2DParticleGamma(){
  TCanvas *cExEg = new TCanvas("cExEg","cExEg",1000,1000);
  chain->Draw("AddBack_EDC:Ex>>ExEg(200,-1,9,2500,0,5)",
        "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && Ex@.size()==1 ",
	"colz");
  TH1F* ExEg = (TH1F*) gDirectory->Get("ExEg");
  ExEg->SetTitle("Ex-Egamma");
  ExEg->GetXaxis()->SetTitle("Ex [MeV]");
  ExEg->GetYaxis()->SetTitle("Eg [MeV]");
  TLine *XeqY = new TLine(0,0,9,9);
  XeqY->SetLineColor(kRed);
  XeqY->Draw();
}

void Load_2DGammaGamma(){
  TCanvas *cEgEg = new TCanvas("cEgEg","cEgEg",1000,1000);
  TH2F *hEgEg = new TH2F("hEgEg","Loaded 2D Gamma-Gamma",200,-1,9,2500,0,5);
  TFile *file = new TFile("LoadHistograms/Load_2DGammaGamma.root","READ");
  hEgEg = (TH2F*)file->Get("gg");
  hEgEg->SetName("hEgEg");
  hEgEg->Draw("colz");
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
  TCanvas *cEgEg = new TCanvas("cEgEg","cEgEg",1000,1000);
  chain->Draw("AddBack_EDC:AddBack_EDC2>>EgEg(999,0.005,5,999,0.005,5)","abs(T_MUGAST_VAMOS-2700)<400","colz");
  TH1F* EgEg = (TH1F*) gDirectory->Get("EgEg");
  chain->Draw("AddBack_EDC2:AddBack_EDC>>EgEg2(999,0.005,5,999,0.005,5)","abs(T_MUGAST_VAMOS-2700)<400","colz");
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
  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && Ex@.size()==1 && abs(AddBack_EDC-" 
  //string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber!=3 && abs(AddBack_EDC-" 
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
  string draw = "Ex>>ExGate(" + to_string(10.0/binsize) + ",-1,9)";

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
  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(width);

  string title = "Gate: "+to_string(gamma-width)+" to "+to_string(gamma+width)+"."
	  + "  BG: "+to_string(bg-width)+" to "+to_string(bg+width)+".";
  
  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  chain->Draw("Ex>>ExGate(100,-1,9)",gating.c_str(),"");
  //chain->Draw("Ex>>ExGate(120,-1,5)",gating.c_str(),"");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle("Counts / 0.10 MeV");
  //ExGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ExGate->SetLineColor(kGreen);
  ExGate->SetFillColor(kGreen);
  ExGate->SetFillStyle(3154);
  ExGate->SetTitle(title.c_str());

  chain->Draw("Ex>>ExBG(100,-1,9)",bggate.c_str(),"same");
  //chain->Draw("Ex>>ExBG(120,-1,5)",bggate.c_str(),"same");
  TH1F* ExBG = (TH1F*) gDirectory->Get("ExBG");
  ExBG->SetLineColor(kRed);
  ExBG->SetFillColor(kRed);
  ExBG->SetFillStyle(3345);

  DrawParticleStates(cEx_Gate);
}

void GateGamma_SeeParticle_WithBG(double gamma, double width, double bg, double widthbg){
  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string bggate = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && Ex@.size()==1 && abs(AddBack_EDC-" 
      + to_string(bg)
      + ")<"
      + to_string(widthbg);

  double ratio = width/widthbg;

  string title = "Gate: "+to_string(gamma-width)+" to "+to_string(gamma+width)+"."
	  + "  BG: "+to_string(bg-width)+" to "+to_string(bg+width)+".";
  
  TCanvas *cEx_Gate = new TCanvas("cEx_Gate","cEx_Gate",1000,1000);
  chain->Draw("Ex>>ExGate(100,-1,9)",gating.c_str(),"");
  //chain->Draw("Ex>>ExGate(120,-1,5)",gating.c_str(),"");
  TH1F* ExGate = (TH1F*) gDirectory->Get("ExGate");
  ExGate->GetXaxis()->SetTitle("Ex [MeV]");
  ExGate->GetYaxis()->SetTitle("Counts / 0.10 MeV");
  //ExGate->GetYaxis()->SetTitle("Counts / 0.05 MeV");
  ExGate->SetLineColor(kGreen);
  ExGate->SetFillColor(kGreen);
  ExGate->SetFillStyle(3154);
  ExGate->SetTitle(title.c_str());

  chain->Draw("Ex>>ExBG(100,-1,9)",bggate.c_str(),"same");
  //chain->Draw("Ex>>ExBG(120,-1,5)",bggate.c_str(),"same");
  TH1F* ExBG = (TH1F*) gDirectory->Get("ExBG");
  ExBG->Scale(ratio);
  ExBG->SetLineColor(kRed);
  ExBG->SetFillColor(kRed);
  ExBG->SetFillStyle(3345);
  //ExBG->Draw("BSAME");

  DrawParticleStates(cEx_Gate);
}

void AddGammaLines(TH1F* hist, double particle, double ymax){
  string base = "sub ";

  for(int i=1; i<means.size();i++){
    string name = base + to_string(means.at(i));
    TLine *line = new TLine(particle-means.at(i), 0.0, particle-means.at(i), ymax);
    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
    line->Draw();
    TText *text = new TText((1.-(means.at(i)/particle))*particle,0.8*ymax,name.c_str());
    text->SetTextAngle(90);
    //text->SetTextSize(40);
    text->Draw();
  }

}

void AddPlacedGammas(TH1F* hist, double ymax){
  hist->Draw();
  for(int i=0; i<knowngammas.size();i++){
    TLine *line = new TLine(knowngammas.at(i), 0.0, knowngammas.at(i), ymax);
    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
    line->Draw();
  }
}

void GateParticle_SeeGamma(double particle, double width){ 
  gStyle->SetOptStat("nemMrRi");

  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber<8 && Mugast.TelescopeNumber>0 && Ex@.size()==1 && abs(Ex-" 
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
  AddGammaLines(EgGate, particle, cEg_Gate->GetUymax());
/*
  TLine *sub0143 = new TLine(particle-0.143, 0.0, particle-0.143, cEg_Gate->GetUymax());
  sub0143->SetLineColor(kBlack); sub0143->SetLineStyle(kDotted);
  TLine *sub0279 = new TLine(particle-0.279, 0.0, particle-0.279, cEg_Gate->GetUymax());
  sub0279->SetLineColor(kBlack); sub0279->SetLineStyle(kDotted);
  TLine *sub0728 = new TLine(particle-0.728, 0.0, particle-0.728, cEg_Gate->GetUymax());
  sub0728->SetLineColor(kBlack); sub0728->SetLineStyle(kDotted);

  TLine *g0143 = new TLine(0.143, 0.0, 0.143, cEg_Gate->GetUymax());
  g0143->SetLineColor(kBlack); g0143->SetLineStyle(kDotted);
  TLine *g0279 = new TLine(0.279, 0.0, 0.279, cEg_Gate->GetUymax());
  g0279->SetLineColor(kBlack); g0279->SetLineStyle(kDotted);
  TLine *g0449 = new TLine(0.449, 0.0, 0.449, cEg_Gate->GetUymax());
  g0449->SetLineColor(kBlack); g0449->SetLineStyle(kDotted);

  EgGate->Draw();
  limit->Draw();
  sub0143->Draw();
  sub0279->Draw();
  sub0728->Draw();
*/

}

void GateParticle_SeeGamma_WithBG(double particle, double width, double bg, double width2){
  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Ex@.size()==1 && abs(Ex-" 
      + to_string(particle)
      + ")<"
      + to_string(width);
  string bggate = "abs(T_MUGAST_VAMOS-2700)<400 && Ex@.size()==1 && abs(Ex-" 
      + to_string(bg)
      + ")<"
      + to_string(width2);

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

void GateGamma_SeeGamma_TimeGate(double gamma, double width){
  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && abs(AddBack_EDC2-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);
  string gating2 = "abs(T_MUGAST_VAMOS-2700)<400 && abs(AddBack_EDC-" 
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





void GateGamma_SeeGamma_WithBG(double gamma, double width, double bg, double width2){
/**/
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
  
  chain->Draw("Ex>>hexp(100,-1,9)","abs(T_MUGAST_VAMOS-2700)<400","");
  TH1F* hexp = (TH1F*) gDirectory->Get("hexp");
  hexp->SetTitle("Comparing Simulation to Experiment");
  hexp->GetXaxis()->SetTitle("Ex [MeV]");
  hexp->GetYaxis()->SetTitle("Counts / 0.1 MeV");
  hexp->SetLineColor(kRed);

  TFile* simfile = new TFile("../../../Outputs/Analysis/SimTest_Jun28_p32_ExHeavy0p143.root", 
  //TFile* simfile = new TFile("../../../Outputs/Analysis/SimTest_Jun22_TWOFNR_p32.root", 
		  "READ");
  TTree* simtree = (TTree*) simfile->FindObjectAny("PhysicsTree");
  simtree->Draw("Ex>>hsimMGp32(100,-1,9)",
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
  TCanvas *diagnoseTheta = new TCanvas("diagnoseTheta","diagnoseTheta",1000,1000);
  chain->Draw(
    "Ex:ThetaLab>>thetaHist(120,100,160,180,-1,8)", 
    "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    "colz");
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

  string gating = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && abs(AddBack_EDC-"
	        + to_string(gamma) + ") < " + to_string(width); 

  chain->Draw("Ex:ThetaLab>>thetaHist(60,100,160,100,-1,9)", gating.c_str(), "colz");
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
  chain->Draw("ELab:ThetaLab>>hKine(360,0,180,500,0,10)","abs(T_MUGAST_VAMOS-2700)<400","col");
  TH2F* hKine = (TH2F*) gDirectory->Get("hKine");
  hKine->SetTitle("");
  hKine->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  hKine->GetYaxis()->SetTitle("E_{lab} [MeV]");
  plot_kine(Kdt, 0.000, kBlack, 2, 1);
  
  plot_kine(Kdp, 0.000, kBlack, 2, 1);
  plot_kine(Kdp, 4.644, kBlack, 2, 1);

  /**
  plot_kine(Kdp, 0.143, kRed, 1, 2);
  plot_kine(Kdp, 0.968, kRed, 1, 2);
  plot_kine(Kdp, 1.410, kRed, 1, 2);
  plot_kine(Kdp, 1.981, kRed, 1, 2);
  plot_kine(Kdp, 2.410, kRed, 1, 2);
  plot_kine(Kdp, 2.907, kRed, 1, 2);
  plot_kine(Kdp, 3.600, kRed, 1, 2);
  plot_kine(Kdp, 3.8  , kRed, 1, 2);
  plot_kine(Kdp, 4.3  , kRed, 1, 2);
  plot_kine(Kdp, 4.507, kRed, 1, 2);
  **/

  plot_kine(Kdd, 0.000, kBlack, 2, 9);
  plot_kine(Kpp, 0.000, kBlack, 2, 9);

  plot_kine(Cadp, 0.000, kRed, 2, 1);
  plot_kine(Tidp, 0.000, kBlue, 2, 1);
  plot_kine(Scdp, 0.000, kGreen, 2, 1);
  //plot_kine(Tidp, 5.652, kBlack, 2, 6); //strongest populated state according to PDBarnes(1965)
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
  string gates = "abs(T_MUGAST_VAMOS-2700)<400 && MUST2.TelescopeNumber==5 && ELab > " 
	       + to_string(EMin) 
	       + " && ELab < " 
	       + to_string(EMax);

  chain->Draw("ThetaLab>>hist(80,50,90)", gates.c_str(), "");
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
    string draw = "Ex>>" + histname + "(100,-1,9)";

    
    chain->Draw(draw.c_str(),gate.c_str(),"same");
    TH1F* hist = (TH1F*) gDirectory->Get(histname.c_str());
    hist->SetLineColor(i+1);

  }

}

void GateThetaLab_MultiWrite(double startTheta, double finishTheta, int numGates, double binsize){
  string core = "abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ThetaLab > ";
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
    string draw = "Ex>>" + histname + "(" + to_string(10.0/binsize) + ",-1,9)";

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

void GatePhaseSpaceByThetaLab_MultiWrite(double startTheta, double finishTheta, int numGates, double binsize){
  string core = "EventWeight*(Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ThetaLab > ";
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  double gatesize = (finishTheta-startTheta)/numGates;
  TList* list = new TList();

  TFile* psfile = new TFile("../../../Outputs/Analysis/Sim_02Mar_47Kdp_PhaseSpace.root","READ");
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
    string draw = "Ex>>" + histname + "(" + to_string(10.0/binsize) + ",-1,9)";

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

//void gggLoad(TTree* chain, TH3F* h){
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
 
   cout << "LOADING FILES: 47Kdp_11Apr22_PartI & II" << endl;

   auto h=new TH2F("gg","gg",1000,0,10,1000,0,10);
   auto DataFile = new TFile("../../../Outputs/Analysis/47Kdp_11Apr22_PartI.root", "READ");
   auto chain = (TTree*) DataFile->FindObjectAny("PhysicsTree");
   ggLoad(chain, h);

   auto h2=new TH2F("gg","gg",1000,0,10,1000,0,10);
   auto DataFile2 = new TFile("../../../Outputs/Analysis/47Kdp_11Apr22_PartII.root", "READ");
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

