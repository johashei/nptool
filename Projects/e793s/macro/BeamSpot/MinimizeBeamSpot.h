////////////////////////////////////////////////////////////////////////////////
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

////////////////////////////////////////////////////////////////////////////////
/*   Global   */
//Various numbers and objects

vector<TVector3*> pos;
vector<double> energy;
vector<int> detnum;
unsigned int mgSelect = 10;
NPL::EnergyLoss CD2("proton_CD2.G4table","G4Table",100);
NPL::EnergyLoss Al("proton_Al.G4table","G4Table",100);
NPL::EnergyLoss BeamTarget("K47_CD2.G4table","G4Table",100);
bool flagDraw = 0;
bool allButMG3 = 0;

//Output files
TFile *histfile = new TFile("./gridHistograms.root", "RECREATE");
ofstream file;
int writeCount = 0;

//Histograms
string filename;
double refE; // the energy of the selected states
static auto h = new TH1D("h","All MG#'s", 400,-1.0,3.0);
static auto h1 = new TH1D("h1","Individual MG#'s", 80,-1.0,3.0);
static auto h2 = new TH1D("h2","h2", 80,-1.0,3.0);
static auto h3 = new TH1D("h3","h3", 80,-1.0,3.0);
static auto h4 = new TH1D("h4","h4", 80,-1.0,3.0);
static auto h5 = new TH1D("h5","h5", 80,-1.0,3.0);
static auto h7 = new TH1D("h7","h7", 80,-1.0,3.0);
static auto hT = new TH2F("hT","hT", 60,100.,160.,80,-1.0,3.0);

Double_t f_full(Double_t *x, Double_t *par) {
  float xx = x[0];
  double result, norm;
  // Flat background
  //result = par[0];
  result = 0;
  // Add N peaks
  for(int pk=0; pk<3; pk++){
    result += (par[3+(pk*3)]/(par[1+(pk*3)]*sqrt(2*pi)))
	      * exp(-0.5*pow((xx-par[2+(pk*3)])/par[1+(pk*3)],2));
  }
  return result;
}


//static auto hT = new TH2F("hT","hT", 20,100.,160.,20,-1.0,1.0);
////////////////////////////////////////////////////////////////////////////////
void LoadFile(){
  // Open XYZE gamma-gated file
  ifstream file(filename.c_str());
  if(!file.is_open()){
    cout << "fail to load file " << filename << endl;
  }
  else {
    cout <<  "Success opening file " << filename << endl;
  }

  // Read in
  int mg;
  double x,y,z,e;
  while(file >> mg >> x >> y >> z >> e ){
    auto p = new TVector3(x,y,z);
//    if(mg==3){
//      p->SetZ(z+1.0); //Edit MG3 position directly
//      cout << "EDITING MG3 POSITION" << endl;
//    }
    detnum.push_back(mg);
    pos.push_back(p);
    energy.push_back(e);
  }
  file.close();
}
////////////////////////////////////////////////////////////////////////////////
void FillMatrix(double* matrix, TFitResultPtr fit){
  matrix[0] = fit->Parameter(2);    //Mean
  matrix[1] = fit->ParError(2);
  matrix[2] = fit->Parameter(1);    //StdDev
  matrix[3] = fit->ParError(1);
  matrix[4] = fit->Chi2()/fit->Ndf(); //Chi2/NDF
  matrix[5] = fit->Parameter(5);//(8);    //Mean2
  matrix[6] = fit->Parameter(4);//(7);    //StdDev2
  matrix[7] = fit->Parameter(8);//(8);    //Mean2
  matrix[8] = fit->Parameter(7);//(7);    //StdDev2

  if(flagDraw){
    cout << "\n        Mean = " << matrix[0] << " +/- " << matrix[1] << endl;
    cout << "      StdDev = " << matrix[2] << " +/- " << matrix[3] << endl;
    cout << "    Chi2/NDF = " << matrix[4] << endl;
    cout << "    Mean2 = " << matrix[5] << " StdDev2 = " << matrix[6] << endl;
    cout << "    Mean3 = " << matrix[7] << " StdDev2 = " << matrix[8] << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//Overloaded function definiton; this is for MUGAST Normal vectors
void DetectorSwitch(int MG, TVector3& Normal ){
  switch(MG){
      case 1:
        Normal.SetXYZ(-0.454552, +0.454996, -0.765742);
        break;
      case 2:
	Normal.SetXYZ(-0.641920, +0.002239, -0.766769);
        break;
      case 3:
	Normal.SetXYZ(-0.455476, -0.452406, -0.766727);
        break;
      case 4:
	Normal.SetXYZ(-0.003212, -0.641831, -0.766839);
        break;
      case 5:
	Normal.SetXYZ(+0.452522, -0.455595, -0.766588);
        break;
      case 7:
	Normal.SetXYZ(+0.454034, +0.458527, -0.763941);
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
void WriteToCout(double* result, double thick, double metric){
  cout 
    << "Mean: " 
    << result[0] 
    << " +/- " 
    << result[1] 
    << "    StdDev: " 
    << result[2] 
    << " +/- " 
    << result[3] 
    << "    Thick: " 
    << thick 
    << " um" 
    << "    Chi2/NDF = " 
    << result[4]
    << "    Metric: " 
    << metric
    << "    Mean2: " 
    << result[5] 
    << "    StdDev2: " 
    << result[6] 
    << "    Mean3: " 
    << result[7] 
    << "    StdDev3: " 
    << result[8] 
    << endl;
}
////////////////////////////////////////////////////////////////////////////////
void WriteToFile(double* result, const double* parameter, double metric){
  file 
    << parameter[0] << "\t"
    << parameter[1] << "\t"
    << parameter[2] << "\t"
    << parameter[3] << "\t"
    << metric << "\t"
    << result[0] << "\t" 
    << result[1] << "\t"
    << result[2] << "\t"
    << result[3] << "\t"
    << result[4] 
    << endl;
}
////////////////////////////////////////////////////////////////////////////////
void DrawOneHistogram(TH1D* hist, int mg, int colour, int fill, double *FitResultMatrixMG, bool drawFit){
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

  const char* settings;
  if (drawFit){
    settings = "RBWQS";
  } else {
    settings = "RBWQSN";
  }  


  TF1 *full = new TF1("fitThreePeaks", f_full, -1.0, +2.6, (int) 1+(3*3));
  for(int i=0; i<3; i++) {
    full->SetParameter((i*3)+1,0.14);
    full->SetParameter((i*3)+3,1e3);
  }
  full->SetParameter((0*3)+2,0.143);
  full->SetParameter((1*3)+2,1.410);
  full->SetParameter((2*3)+2,1.981);

  TFitResultPtr fit = hist->Fit(full, settings, "", -1.0, +2.6);

  

  //TFitResultPtr fit = hist->Fit("gaus",settings); //N=stop drawing, Q=stop writing
  FillMatrix(FitResultMatrixMG,fit);
} 
////////////////////////////////////////////////////////////////////////////////
void InitiliseCanvas(double FitResultMatrix[7][9]){

  //Canvas setup
  TCanvas *canv = new TCanvas("canv","Ex Histograms",20,20,1600,800);
  gStyle->SetOptStat(0);
  //canv->Divide(2,1,0.005,0.005,0);
  canv->Divide(3,1,0.005,0.005,0);
  canv->cd(1)->SetLeftMargin(0.15);
  canv->cd(1)->SetBottomMargin(0.15);
  gPad->SetTickx();
  gPad->SetTicky();
  canv->cd(2)->SetLeftMargin(0.15);
  canv->cd(2)->SetBottomMargin(0.15);
  gPad->SetTickx();
  gPad->SetTicky();

  canv->cd(3)->SetLeftMargin(0.15);
  canv->cd(3)->SetBottomMargin(0.15);
  gPad->SetTickx();
  gPad->SetTicky();

  //Histogram setup - Individual
  canv->cd(1);
  h1->SetMaximum(150.);
  h1->GetXaxis()->SetTitle("Ex [MeV]");
  h1->GetYaxis()->SetTitle("Counts");
  
  //Histogram draw - Individual
  DrawOneHistogram(h1, 1, 632, 0,    FitResultMatrix[1], 0);
  DrawOneHistogram(h2, 2, 800, 3001,    FitResultMatrix[2], 0);
  if(allButMG3){
    DrawOneHistogram(h3, 3, 416, 0,    FitResultMatrix[3], 0); //3344
  } else {
    DrawOneHistogram(h3, 3, 416, 0,    FitResultMatrix[3], 0); //3344
  }
  DrawOneHistogram(h4, 4, 840, 0,    FitResultMatrix[4], 0);
  DrawOneHistogram(h5, 5, 600, 3001,    FitResultMatrix[5], 0);
  DrawOneHistogram(h7, 6, 880, 0,    FitResultMatrix[6], 0);

  canv->Update();
  TLine *line=new TLine(refE,0.0,refE,150.0);
  line->SetLineColor(kBlack);
  line->SetLineStyle(7);
  line->Draw();

  //Format legend
  auto legend = new TLegend(0.15,0.7,0.45,0.9);
  legend->AddEntry(h1,"MG #1","f");
  legend->AddEntry(h2,"MG #2","f");
  legend->AddEntry(h3,"MG #3","f");
  legend->AddEntry(h4,"MG #4","f");
  legend->AddEntry(h5,"MG #5","f");
  legend->AddEntry(h7,"MG #7","f");
  //legend->SetTextSize(20);
  legend->Draw();
  
  //Histogram setup - Sum
  canv->cd(2);
  h->GetXaxis()->SetTitle("Ex [MeV]");
  h->GetYaxis()->SetTitle("Counts");
  
  //Histogram draw - Sum
  DrawOneHistogram(h, 0, 1, 0, FitResultMatrix[0],1);

  canv->Update();
  TLine *line2=new TLine(refE,0.0,refE,h->GetMaximum()+10.);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(7);
  line2->Draw();

  canv->cd(3);
  hT->GetXaxis()->SetTitle("Theta (degrees)");
  hT->GetYaxis()->SetTitle("Ex [MeV]");
  hT->Draw("colz");
  canv->Update();
  TLine *l0143 = new TLine(100., 0.143, 160., 0.143);
  TLine *l1981 = new TLine(100., 1.981, 160., 1.981);
  l0143->SetLineStyle(kDashed);
  l0143->SetLineColor(kRed);
  l0143->Draw("same");

  //Refresh
  gPad->Update();

  if(flagDraw){
    //writeCount++;
    histfile->cd();
    canv->Write("Minimum#");
  }
  else { 
    delete canv;
  }
}
////////////////////////////////////////////////////////////////////////////////

