// STL
#include<stdlib.h>
#include<stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

// Root
#include "TString.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"

//NPTool
#include "NPEnergyLoss.h"
using namespace NPL;

/// DEFINING GLOBAL VARIABLE

double mean_extrapolation=100;

/// Parameter used in the macro
// Micrometer
Double_t AlThickness;
Double_t SiThickness;
EnergyLoss EL_Al("./EnergyLossTable/alpha_Al.G4table" , "G4Table", 100) ;
EnergyLoss EL_Si("./EnergyLossTable/alpha_Si.G4table" , "G4Table", 100) ;
// Information about the calibration condition (use Latex marks-up)

const TString xy                  = "FRONT";

const TString Experiment          = "ccam";
const TString Run_Period          = "2020";
const TString Operator            = "Anne MEYER";
const TString Source              = "207Bi";
const TString sourceName          = "207Bi";
const TString Comment             = "Source at 0$^{\\circ}$";
const char* frun_207Bi            = "20200128_10h44_bi207_conv";
const char* frun_pulserFront      = "run37";
const char* frun_pulserBack       = "run40";

int Telescope_Number=0;
const int Strip_Start=1;
const int Strip_End=32;

// choosing a method for the fit
const TString method = "ZeroExtrapolation" ;
const bool RefitWithSatellite = false ;
const bool Pedestals_Aligned = false ;   
const bool Pedestals_pulser = false ;
Int_t CurrentTelescope = 0 ;
Int_t CurrentStrip     = 0 ;
TString folder;
TString main_name;
TCanvas* Tsummary;
TCanvas* Buffer;

map<int,string> BadStrip;

// Defining the array used after (order needs to be diffent for X and Y )
Int_t NumberOfIsotope;
Int_t NumberOfPeaks;

// Source original value
Int_t Source_Number_Peak;
TString*  Source_isotope;
Double_t* Source_branching_ratio;
Double_t* Source_E;
Double_t* Source_Sig;

// Source corrected value
Double_t particle_angle;
Double_t* energyFront;
Double_t* errorsFront;
Double_t* energyBack;
Double_t* errorsBack;

// Calibration Coefficient
Double_t a ;
Double_t b ;
Double_t pedestal;

Double_t* mean       = new Double_t[6]; 
Double_t* error_mean = new Double_t[6]; 
Double_t* sigma      = new Double_t[6];
Double_t* error_par  = new Double_t[6];

TGraph* ZeroDispersion ;
ofstream peaks_file, calib_file, dispersion_file , calibError_file, calib_online_file, latex_file;; 

TH1F* sigma_fit  ;
TH1F* Dispersion ;
TGraph* coeff_a ;
TGraph* coeff_b ;

TFile *inFile;
TFile *f_pf;
TFile *f_pb;

/// Function Header
void AutoCalibration(int,int);
void EnergyCalibrator();
Double_t Pedestals(TH1F *);
void Alpha(TH1F*, TString, Double_t);
bool Finder(TH1F*, TString , Double_t*, Double_t*);
Double_t Calib_ZeroForceMethod(string ,TGraphErrors*,float, Double_t*, Double_t*);
Double_t Calib_ZeroExtrapolationMethod(TH1F* hist ,string ,TGraphErrors*,float, Double_t*, Double_t*, Double_t*, Double_t &a , Double_t &b);
void LatexSummaryHeader(TString xy);
void LatexSummaryEnder();
void LatexSummaryTelescope();
void DefineSource(TString sourceName, Double_t);
//void DefineSource(TString sourceName="3 alphas");


void Find_Satellites(TH1F *h);
Double_t source_Pu(Double_t *x, Double_t *par);
Double_t source_Am(Double_t *x, Double_t *par);
Double_t source_Cm(Double_t *x, Double_t *par);

