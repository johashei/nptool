#include "DefineColours.h"
#include "GausFit.h"
#include "KnownPeakFitter.h"
#include "DrawPlots.h"
#include "ElasticsFitELabGates.h"

//#include "CS2.h"
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

void Plots_47Kdd(){

  LoadChain47Kdd();
  gStyle->SetOptStat("nemMrRi");

  tCentre = 2750;  tRange = 200;
  timegate = "abs(T_MUGAST_VAMOS-" + to_string(tCentre) + ")<" + to_string(tRange);
  det_gate = "MUST2.TelescopeNumber==5";
  reactionName = "47K(d,d)";
  
  cout << "==============================================" << endl;
  cout << "=============== (d,d) reaction ===============" << endl;
  cout << "==============================================" << endl;

}
