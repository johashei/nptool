#include "DefineColours.h"
#include "GausFit.h"
#include "KnownPeakFitter.h"
#include "DrawPlots.h"

#include "CS2.h"
//#include "CS2_MGX.h"
#include "ThreeBodyBreakup.h"
#include "ThreeBodyBreakup_FitPhaseSpace.h"
#include "20Oct22_CompareYield.h"


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

void CompareCountsInThetaLab(){
  auto canv = new TCanvas("cCompareExpSim","cCompareExpSim",1500,1500);
  auto file = new TFile("../../../Outputs/Analysis/Sim_47Kdp_10Aug22_TrueStripRemoval.root");  
  auto tree = (TTree*) file->FindObjectAny("PhysicsTree");

  canv->Divide(2,3);
  canv->cd(1);
    chain->Draw("ThetaLab>>exp1(120,100,160)","Mugast.TelescopeNumber==1","");
    tree->Draw("ThetaLab>>sim1(120,100,160)","Mugast.TelescopeNumber==1","same hist");
    auto exp1 = (TH1F*) gDirectory->Get("exp1");
    auto sim1 = (TH1F*) gDirectory->Get("sim1");
    exp1->SetLineColor(kRed);
    exp1->GetYaxis()->SetRangeUser(0.,700.);
    sim1->SetLineColor(kBlue);
    sim1->Scale(0.05);
  canv->cd(2);
    chain->Draw("ThetaLab>>exp2(120,100,160)","Mugast.TelescopeNumber==2","");
    tree->Draw("ThetaLab>>sim2(120,100,160)","Mugast.TelescopeNumber==2","same hist");
    auto exp2 = (TH1F*) gDirectory->Get("exp2");
    auto sim2 = (TH1F*) gDirectory->Get("sim2");
    exp2->SetLineColor(kRed);
    exp2->GetYaxis()->SetRangeUser(0.,700.);
    sim2->SetLineColor(kBlue);
    sim2->Scale(0.05);
  canv->cd(3);
    chain->Draw("ThetaLab>>exp3(120,100,160)","Mugast.TelescopeNumber==3","");
    tree->Draw("ThetaLab>>sim3(120,100,160)","Mugast.TelescopeNumber==3","same hist");
    auto exp3 = (TH1F*) gDirectory->Get("exp3");
    auto sim3 = (TH1F*) gDirectory->Get("sim3");
    exp3->SetLineColor(kRed);
    exp3->GetYaxis()->SetRangeUser(0.,700.);
    sim3->SetLineColor(kBlue);
    sim3->Scale(0.05);
  canv->cd(4);
    chain->Draw("ThetaLab>>exp4(120,100,160)","Mugast.TelescopeNumber==4","");
    tree->Draw("ThetaLab>>sim4(120,100,160)","Mugast.TelescopeNumber==4","same hist");
    auto exp4 = (TH1F*) gDirectory->Get("exp4");
    auto sim4 = (TH1F*) gDirectory->Get("sim4");
    exp4->SetLineColor(kRed);
    exp4->GetYaxis()->SetRangeUser(0.,700.);
    sim4->SetLineColor(kBlue);
    sim4->Scale(0.05);
  canv->cd(5);
    chain->Draw("ThetaLab>>exp5(120,100,160)","Mugast.TelescopeNumber==5","");
    tree->Draw("ThetaLab>>sim5(120,100,160)","Mugast.TelescopeNumber==5","same hist");
    auto exp5 = (TH1F*) gDirectory->Get("exp5");
    auto sim5 = (TH1F*) gDirectory->Get("sim5");
    exp5->SetLineColor(kRed);
    exp5->GetYaxis()->SetRangeUser(0.,700.);
    sim5->SetLineColor(kBlue);
    sim5->Scale(0.05);
  canv->cd(6);
    chain->Draw("ThetaLab>>exp7(120,100,160)","Mugast.TelescopeNumber==7","");
    tree->Draw("ThetaLab>>sim7(120,100,160)","Mugast.TelescopeNumber==7","same hist");
    auto exp7 = (TH1F*) gDirectory->Get("exp7");
    auto sim7 = (TH1F*) gDirectory->Get("sim7");
    exp7->SetLineColor(kRed);
    exp7->GetYaxis()->SetRangeUser(0.,700.);
    sim7->SetLineColor(kBlue);
    sim7->Scale(0.05);

}

/* MAIN FUNCTION */

void Plots_47Kdp(){

  LoadChain47Kdp();
  gStyle->SetOptStat("nemMrRi");

  tCentre = 2700;  tRange = 400;
  timegate = "abs(T_MUGAST_VAMOS-" + to_string(tCentre) + ")<" + to_string(tRange);
  det_gate = "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8";
  reactionName = "47K(d,p)";

  cout << "==============================================" << endl;
  cout << "=============== (d,p) reaction ===============" << endl;
  cout << "==============================================" << endl;
  cout << ""<< endl;
  CS();
  cout << ""<< endl;
  cout << "==============================================" << endl;

}
