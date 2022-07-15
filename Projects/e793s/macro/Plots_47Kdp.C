#include "GausFit.h"
#include "KnownPeakFitter.h"
#include "DrawPlots.h"

#include "CS2.h"
#include "ThreeBodyBreakup.h"
#include "ThreeBodyBreakup_FitPhaseSpace.h"


void AddGammaLinesMG(TH1F* hist, double particle, double ymax){
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

void AddPlacedGammasMG(TH1F* hist, double ymax){
  hist->Draw();
  for(int i=0; i<knowngammas.size();i++){
    TLine *line = new TLine(knowngammas.at(i), 0.0, knowngammas.at(i), ymax);
    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
    line->Draw();
  }
}

/* MAIN FUNCTION */

void Plots_47Kdp(){

  LoadChain47Kdp();
  gStyle->SetOptStat("nemMrRi");

  tCentre = 2700;  tRange = 400;
  timegate = "abs(T_MUGAST_VAMOS-" + to_string(tCentre) + ")<" + to_string(tRange);
  det_gate = "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8";
  
  cout << "==============================================" << endl;
  cout << "=============== (d,p) reaction ===============" << endl;
  cout << "==============================================" << endl;
  cout << ""<< endl;
  cout << "- CS(stateE, stateSp, orb_l, orb_j, nodes) "<< endl;
  cout << "---- 0.143, p3/2 = CS(0.143, 2, 1, 1.5) "<< endl;
  cout << "---- 0.279, p3/2 = CS(0.279, 2, 1, 1.5) "<< endl;
  cout << "---- 0.728, f7/2 = CS(0.728, 3, 3, 3.5) "<< endl;
  cout << "---- 0.968, p1/2 = CS(0.968, 0, 1, 0.5) "<< endl;
  cout << "---- 1.410, p3/2 = CS(1.410, 2, 1, 1.5) "<< endl;
  cout << "---- 1.981, p3/2 = CS(1.981, 2, 1, 1.5) "<< endl;
  cout << "---- 2.410, p3/2 = CS(2.410, 0, 1, 0.5) "<< endl;
  cout << "---- 3.2  , f7/2 = CS(3.2  , 3, 3, 3.5) "<< endl;
  cout << "---- 3.6  , f5/2 = CS(3.6  , 3, 3, 2.5) "<< endl;
  cout << "---- 3.8  , f5/2 = CS(3.8  , 3, 3, 2.5) "<< endl;
  cout << "---- 4.1  , f5/2 = CS(4.1  , 3, 3, 2.5) "<< endl;
  cout << "---- 4.4  , f5/2 = CS(4.4  , 3, 3, 2.5) "<< endl;
  cout << ""<< endl;
  cout << "==============================================" << endl;

}
