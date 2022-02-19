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

vector<TVector3*> GammaP;
vector<TVector3*> Beta;
vector<double> GammaE;
unsigned int mgSelect = 10;
bool flagDraw = 0;

//Output files
TFile *histfile = new TFile("./gridHistograms.root", "RECREATE");
ofstream file;
int writeCount = 0;

//Histograms
string filename;
double refE; // the energy of the selected states
static auto h = new TH1D("h","EGamma", 10000,0.0,5.0);
////////////////////////////////////////////////////////////////////////////////
void LoadFile(){
  // Open EventReader output file
  ifstream file(filename.c_str());
  if(!file.is_open()){
    cout << "fail to load file " << filename << endl;
  }
  else {
    cout <<  "Success opening file " << filename << endl;
  }

  // Read in
  //int mg;
  double GPx,GPy,GPz,GE,Bx,By,Bz;
  while(file >> GPx >> GPy >> GPz >> GE >> Bx >> By >> Bz ){
    auto GP = new TVector3(GPx,GPy,GPz);
    auto B = new TVector3(Bx,By,Bz);
    GammaP.push_back(GP);
    Beta.push_back(B);
    GammaE.push_back(GE);
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
//Overloaded function definiton; this is for MUGAST Normal vectors
/*void DetectorSwitch(int MG, TVector3& Normal ){
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
}*/
////////////////////////////////////////////////////////////////////////////////
//Overloaded function definiton; this is for filling individual Ex histograms
/*void DetectorSwitch(int MG, double Ex){
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
}*/
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
/*void DrawOneHistogram(TH1D* hist, int mg, int colour, int fill, double *FitResultMatrixMG, bool drawFit){
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
    settings = "WQS";
  } else {
    settings = "WQSN";
  }  

  TFitResultPtr fit = hist->Fit("gaus",settings); //N=stop drawing, Q=stop writing
  FillMatrix(FitResultMatrixMG,fit);
}*/ 
////////////////////////////////////////////////////////////////////////////////
void InitiliseCanvas(double FitResultMatrix[5]){
  //Canvas setup
  TCanvas *canv = new TCanvas("canv","Egamma",20,20,1600,800);
  gStyle->SetOptStat(0);
  h->Draw();

  h->GetXaxis()->SetRangeUser(0.05,refE+0.5);
  //h->GetXaxis()->SetRangeUser(0.1,2.1);

  TF1 *gauswithbg = new TF1("gauswithbg",
		  "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",
		  refE-0.005,refE+0.005);
  gauswithbg->SetParNames("Const","Mean","Sigma","Flat");
  gauswithbg->SetParameter(0,300.);  gauswithbg->SetParLimits(0,5.,800.);
  gauswithbg->SetParameter(1,refE);  gauswithbg->SetParLimits(1,refE-0.001,refE+0.001);
  gauswithbg->SetParameter(2,0.005); gauswithbg->SetParLimits(2,0.0005,0.010);
  gauswithbg->SetParameter(3,10.);   gauswithbg->SetParLimits(3,0.,100.);

  TFitResultPtr fit;
  //if(refE<1.0){
    fit = h->Fit("gauswithbg","WQS","",refE-0.05,refE+0.05);
  //} else {
  //  fit = h->Fit("gauswithbg","WQS","",refE-0.05,refE+0.05);
  //}
  FillMatrix(FitResultMatrix,fit);


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
