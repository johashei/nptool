#include "DefineColours.h"
#include "GausFit.h"
#include "KnownPeakFitter.h"
#include "DrawPlots.h"

#include "CS2_dt.h"
//#include "ThreeBodyBreakup.h"
//#include "ThreeBodyBreakup_FitPhaseSpace.h"


void AddGammaLines(TH1F* hist, double particle, double ymax){
//  string base = "sub ";
//
//  for(int i=1; i<means.size();i++){
//    string name = base + to_string(means.at(i));
//    TLine *line = new TLine(particle-means.at(i), 0.0, particle-means.at(i), ymax);
//    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
//    line->Draw();
//    TText *text = new TText((1.-(means.at(i)/particle))*particle,0.8*ymax,name.c_str());
//    text->SetTextAngle(90);
//    //text->SetTextSize(40);
//    text->Draw();
//  }
}

void AddPlacedGammas(TH1F* hist, double ymax){
//  hist->Draw();
//  for(int i=0; i<knowngammas.size();i++){
//    TLine *line = new TLine(knowngammas.at(i), 0.0, knowngammas.at(i), ymax);
//    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
//    line->Draw();
//  }
}

/* MAIN FUNCTION */

void MM_Timing_Comparison(){
  TCanvas* c = new TCanvas("cMMT","cMMT",1000,1000);


  chain->Draw("MUST2.Si_T>>h1(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==1 && MUST2.CsI_E>0",
	      "same");
               TH1F* h1 = (TH1F*) gDirectory->Get("h1");
  chain->Draw("MUST2.Si_T>>h2(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==2 && MUST2.CsI_E>0",
	      "same");
               TH1F* h2 = (TH1F*) gDirectory->Get("h2");
  chain->Draw("MUST2.Si_T>>h3(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==3 && MUST2.CsI_E>0",
	      "same");
               TH1F* h3 = (TH1F*) gDirectory->Get("h3");
  chain->Draw("MUST2.Si_T>>h4(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==4 && MUST2.CsI_E>0",
	      "same");
               TH1F* h4 = (TH1F*) gDirectory->Get("h4");

  h1->SetLineColor(kRed);
  h2->SetLineColor(kGreen);
  h3->SetLineColor(kBlue);
  h4->SetLineColor(kViolet);

//  h1->Draw();
//  h2->Draw("same");
//  h3->Draw("same");
//  h4->Draw("same");
}

void Plots_47Kdt(){

  /* Load graphical cut */
  //TFile gcIn("GraphicalCut_22Jun22.root");
  //TCutG* cutTritons = (TCutG*) gcIn.FindObjectAny("cutTritons");
  
  //TFile gcIn("cutTritonsWide.root");
  TFile gcIn("cutTritons_26Aug22Long.root");
  //TCutG* cutTritons = (TCutG*) gcIn.FindObjectAny("cutTritonsWide");
  TCutG* cutTritons = (TCutG*) gcIn.FindObjectAny("cutTritons");
  cutTritons->SetName("cutTritons");

  TFile gcIn2("cutTime.root");
  TCutG* cutTime = (TCutG*) gcIn2.FindObjectAny("cutTime");
  cutTime->SetName("cutTime");

  TFile gcIn3("cutDoublePeakGarbage.root");
  TCutG* cutGarbage = (TCutG*) gcIn3.FindObjectAny("cutDoublePeakGarbage");
  cutGarbage->SetName("cutGarbage");


  /**************/
  TFile gcInA("cutTritons_HighTLowE.root");
  TCutG* cutHighTLowE = (TCutG*) gcInA.FindObjectAny("cutTritons");
  cutHighTLowE->SetName("cutHighTLowE");

  TFile gcInB("cutTritons_SlimGate.root");
  TCutG* cutSlim = (TCutG*) gcInB.FindObjectAny("cutTritons");
  cutSlim->SetName("cutSlim");

  TFile gcInC("cutTritons_HighELowT.root");
  TCutG* cutHighELowT = (TCutG*) gcInC.FindObjectAny("cutTritons");
  cutHighELowT->SetName("cutHighELowT");


  /**************/

  LoadChain47Kdt();
  gStyle->SetOptStat("nemMrRi");

  tCentre = 2750;  tRange = 350; //Wide is fine because I use the 2D time gate
  //tCentre = 2550;  tRange = 150;
  timegate = "abs(T_MUGAST_VAMOS-" + to_string(tCentre) + ")<" + to_string(tRange);
  det_gate = "MUST2.TelescopeNumber>0 && MUST2.TelescopeNumber<5";
  reactionName = "47K(d,t)";
  
  cout << "==============================================" << endl;
  cout << "=============== (d,t) reaction ===============" << endl;
  cout << "==============================================" << endl;
  cout << " - (d,t) selection: 'cutTritons'  " << endl;


}
