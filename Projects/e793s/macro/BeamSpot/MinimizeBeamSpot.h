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
unsigned int mgSelect = 10;
NPL::EnergyLoss CD2("proton_CD2.G4table","G4Table",100);
NPL::EnergyLoss Al("proton_Al.G4table","G4Table",100);
using namespace std;

bool flagDraw = 0;

//double FitResultMatrix[7][5];
// 7 => Sum in 0 and them MG's in 1-6
// 5 => Mean, MeanErr, StdDev, StdDevErr, Chi2/NDF


//TCanvas *canv = new TCanvas("canv","Ex Histograms",20,20,1600,800);

static auto h = new TH1D("h","All MG#'s", 80,-1.,1.);
static auto h1 = new TH1D("h1","Individual MG#'s", 40,-1.,1.);
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
void FillMatrix(double* matrix, TFitResultPtr fit){
  matrix[0] = fit->Parameter(1);    //Mean
  matrix[1] = fit->ParError(1);
  matrix[2] = fit->Parameter(2);    //StdDev
  matrix[3] = fit->ParError(2);
  matrix[4] = fit->Chi2()/fit->Ndf(); //Chi2/NDF

  if(flagDraw){
    cout << "\n        Mean = " << matrix[0] << " +/- " << matrix[1] << endl;
    cout << "      StdDev = " << matrix[2] << " +/- " << matrix[3] << endl;
    cout << "    Chi2/NDF = " << matrix[4] << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//[[[[ UPDATE WITH NEW MG POSITIONS FROM SURVEY ]]]]
//Overloaded function definiton; this is for MUGAST Normal vectors
void DetectorSwitch(int MG, TVector3& Normal ){
  switch(MG){
      case 1:
        Normal.SetXYZ(-0.453915, +0.455463, -0.765842);
        break;
      case 2:
	Normal.SetXYZ(-0.642828, +0.000000, -0.766010);
        break;
      case 3:
	Normal.SetXYZ(-0.454594, -0.450670, -0.768271);
        break;
      case 4:
	Normal.SetXYZ(-0.002437, -0.638751, -0.769409);
        break;
      case 5:
	Normal.SetXYZ(+0.452429, -0.454575, -0.767248);
        break;
      case 7:
	Normal.SetXYZ(+0.443072, +0.443265, -0.779232);
        break;
      default:
	cout << "ERROR:: Invalid DetNum " << MG << endl;
        return 1; // Exit code
    }
}
////////////////////////////////////////////////////////////////////////////////
//Overloaded function definiton; this is for filling individual Ex histograms
void DetectorSwitch(int MG, double Ex){
  switch(MG){
      case 1:
        h1->Fill(Ex); 
        break;
      case 2:
        h2->Fill(Ex); 
        break;
      case 3:
        h3->Fill(Ex); 
        break;
      case 4:
        h4->Fill(Ex); 
        break;
      case 5:
        h5->Fill(Ex); 
        break;
      case 7:
        h7->Fill(Ex); 
        break;
      default:
        cout << "ERROR:: Invalid DetNum " << MG << endl;
        return 1; // Exit code
    }
}
////////////////////////////////////////////////////////////////////////////////
void DrawOneHistogram(TH1D* hist, int mg, int colour, int fill, double *FitResultMatrixMG){
  //Hist settings
  hist->SetStats(0);
  hist->SetLineColor(colour);
  hist->SetFillStyle(fill);
  hist->SetFillColor(colour);
  hist->Draw("same");

  if (flagDraw){
    //Header
    cout << noshowpos;
    cout << "\n==================================================" << endl;
    if (mg==6){
      cout << "=---------------------- MG7 ---------------------=" << endl;
    } else if (mg==0) {
      cout << "=---------------------- SUM ---------------------=" << endl;
    } else {
      cout << "=---------------------- MG" << mg << " ---------------------=" << endl;
    }
    cout << showpos;
  }

  TFitResultPtr fit = hist->Fit("gaus","WQS"); //N=stop drawing, Q=stop writing
  FillMatrix(FitResultMatrixMG,fit);
} 
////////////////////////////////////////////////////////////////////////////////
void InitiliseCanvas(double FitResultMatrix[7][5]){

  //Canvas setup
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

  //Histogram setup - Individual
  canv->cd(1);
  h1->SetMaximum(75.);
  h1->GetXaxis()->SetTitle("Ex [MeV]");
  h1->GetYaxis()->SetTitle("Counts");
  
  //Histogram draw - Individual
  DrawOneHistogram(h1, 1, 632, 3244, FitResultMatrix[1]);
  DrawOneHistogram(h2, 2, 800, 3244, FitResultMatrix[2]);
  DrawOneHistogram(h3, 3, 416, 3344, FitResultMatrix[3]);
  DrawOneHistogram(h4, 4, 840, 3444, FitResultMatrix[4]);
  DrawOneHistogram(h5, 5, 600, 3544, FitResultMatrix[5]);
  DrawOneHistogram(h7, 6, 880, 3644, FitResultMatrix[6]);
   
  //Format legend
  auto legend = new TLegend(0.15,0.7,0.35,0.9);
  legend->AddEntry(h1,"MUGAST 1","f");
  legend->AddEntry(h2,"MUGAST 2","f");
  legend->AddEntry(h3,"MUGAST 3","f");
  legend->AddEntry(h4,"MUGAST 4","f");
  legend->AddEntry(h5,"MUGAST 5","f");
  legend->AddEntry(h7,"MUGAST 7","f");
  legend->Draw();
  
  //Histogram setup - Sum
  canv->cd(2);
  h->GetXaxis()->SetTitle("Ex [MeV]");
  h->GetYaxis()->SetTitle("Counts");
  
  //Histogram draw - Sum
  DrawOneHistogram(h, 0, 1, 0, FitResultMatrix[0]);

  //Refresh
  gPad->Update();

  if(!flagDraw){delete canv;}
}
////////////////////////////////////////////////////////////////////////////////

