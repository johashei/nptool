#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

double refE = 0.143; // the energy of the selected states
vector<TVector3*> pos;
vector<double> energy;
vector<int> detnum;
NPL::EnergyLoss CD2("proton_CD2.G4table","G4Table",100);
NPL::EnergyLoss Al("proton_Al.G4table","G4Table",100);
using namespace std;

bool flagDraw = 0;
static auto h = new TH1D("h","h", 80,-1.,1.);
static auto h1 = new TH1D("h1","h1", 40,-1.,1.);
static auto h2 = new TH1D("h2","h2", 40,-1.,1.);
static auto h3 = new TH1D("h3","h3", 40,-1.,1.);
static auto h4 = new TH1D("h4","h4", 40,-1.,1.);
static auto h5 = new TH1D("h5","h5", 40,-1.,1.);
static auto h7 = new TH1D("h7","h7", 40,-1.,1.);

////////////////////////////////////////////////////////////////////////////////
void LoadFile(){
  // Open XYZE gamma-gated file
  ifstream file("XYZE_Full_02June.txt");
  if(!file.is_open()){
    cout << "fail to load file" << endl;
    exit(1);
  }
  else {
    cout <<  "Success opening file" << endl;
  }

  // Read in
  int mg;
  double x,y,z,e;
  while(file >> mg >> x >> y >> z >> e ){
    auto p = new TVector3(x,y,z);
    detnum.push_back(mg);
    pos.push_back(p);
    energy.push_back(e);
  }
  file.close();
}
////////////////////////////////////////////////////////////////////////////////
void InitiliseCanvas(){
  TCanvas *canv = new TCanvas("canv","Ex Histograms",20,20,1600,800);
  gStyle->SetOptStat(0);
  canv->Divide(2,1,0.005,0.005,0);
  canv->cd(1)->SetLeftMargin(0.15);
  canv->cd(1)->SetBottomMargin(0.15);
  gPad->SetTickx();
  gPad->SetTicky();
  canv->cd(2)->SetLeftMargin(0.15);
  canv->cd(2)->SetBottomMargin(0.15);
  gPad->SetTickx();
  gPad->SetTicky();
      
  canv->cd(1);
  h1->SetMaximum(75.);
  h1->GetXaxis()->SetTitle("Ex [MeV]");
  h1->GetYaxis()->SetTitle("Counts");
      
  // ----- MG1 -----
  h1->SetStats(0);
  h1->SetLineColor(kRed);
  h1->SetFillStyle(3244);
  h1->SetFillColor(kRed);
  h1->Draw();
  h1->Fit("gaus","WQ"); //add N to stop it drawing
  
  // ----- MG2 -----
  h2->SetStats(0);
  h2->SetLineColor(kOrange);
  h2->SetFillStyle(3244);
  h2->SetFillColor(kOrange);
  h2->Draw("same");
  h2->Fit("gaus","WQ"); //add N to stop it drawing
          
  // ----- MG3 -----
  h3->SetStats(0);
  h3->SetLineColor(kGreen);
  h3->SetFillStyle(3344);
  h3->SetFillColor(kGreen);
  h3->Draw("same");
  h3->Fit("gaus","WQ"); //add N to stop it drawing
          
  // ----- MG4 -----
  h4->SetStats(0);
  h4->SetLineColor(kTeal);
  h4->SetFillStyle(3444);
  h4->SetFillColor(kTeal);
  h4->Draw("same");
  h4->Fit("gaus","WQ"); //add N to stop it drawing
          
  // ----- MG5 -----
  h5->SetStats(0);
  h5->SetLineColor(kBlue);
  h5->SetFillStyle(3544);
  h5->SetFillColor(kBlue);
  h5->Draw("same");
  h5->Fit("gaus","WQ"); //add N to stop it drawing
	  
  // ----- MG7 -----
  h7->SetStats(0);
  h7->SetLineColor(kViolet);
  h7->SetFillStyle(3644);
  h7->SetFillColor(kViolet);
  h7->Draw("same");
  h7->Fit("gaus","WQ"); //add N to stop it drawing
          
  // Format legend
  auto legend = new TLegend(0.15,0.7,0.35,0.9);
  legend->AddEntry(h1,"MUGAST 1","f");
  legend->AddEntry(h2,"MUGAST 2","f");
  legend->AddEntry(h3,"MUGAST 3","f");
  legend->AddEntry(h4,"MUGAST 4","f");
  legend->AddEntry(h5,"MUGAST 5","f");
  legend->AddEntry(h7,"MUGAST 7","f");
  legend->Draw();
  
  // ----- ALL -----
  canv->cd(2);
  h->SetStats(0);
  h->GetXaxis()->SetTitle("Ex [MeV]");
  h->GetYaxis()->SetTitle("Counts");
  h->Draw();
  h->Fit("gaus", "WQ");
  gPad->Update();
}
////////////////////////////////////////////////////////////////////////////////

