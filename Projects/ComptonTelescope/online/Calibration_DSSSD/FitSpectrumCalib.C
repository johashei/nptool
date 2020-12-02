#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TMatrixTSym.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>

using namespace std;

Double_t fpeaks(Double_t*, Double_t*);
Double_t gaussianPeak(Double_t*, Double_t*);
Double_t indivFcn(Double_t*, Double_t*);
Double_t indivFcn_BG(Double_t*, Double_t*);
Double_t fit_pol1(Double_t*, Double_t*);
void     AddPeak(Double_t, Double_t, Double_t);
void     RemovePeak(Double_t, Double_t);
Double_t CalibUncertainty(Double_t, TF1*, TMatrixDSym);

static vector<Double_t>  peakList;
static Int_t            npeaks;
static Int_t            pmin, pmax, rmin, rmax;
Bool_t Fit_withBGLin;
Bool_t reject;

#define NBSTRIPS  32


void FitSpectrumCalib(Bool_t isFitBGLin = 1)
{

  const int dimE = 6;
  double Energy[dimE] = {481.6935e-3, 553.8372e-3, 565.8473e-3, 975.651e-3, 1047.795e-3, 1059.805e-3}; // MeV
  double error_E[dimE] = {0.0021e-3, 0.0021e-3, 0.0021e-3, 0.003e-3, 0.003e-3, 0.003e-3};
  //double Energy[dimE] = {481.6935, 553.8372, 975.651, 1047.795};
  //double error_E[dimE] = {0.0021, 0.0021, 0.003, 0.003};


  // fit range
  pmin = 40;
  pmax = 980;
  // range to remove from the fit
  rmin = 70;
  rmax = 300;


  // misc.
  Fit_withBGLin = isFitBGLin;
  peakList.clear();

  // define canvas
  TCanvas *c1 = new TCanvas("c1", "resultats", 1000, 700);
  c1->Draw();    

  // open file and get histo
  TFile *inFile = new TFile("./Histograms/20200128_10h44_bi207_conv_RawDSSSDHistos.root");

  // Front strips
  for (int k = 0; k < NBSTRIPS; k++)
  {

    if (k != 1) // remove strip front 2 with no data
    {
      // open histo
      TH1F *hist = (TH1F*) inFile->Get(Form("h_D1_FRONT_E%d", k+1));
      cout << "Fitting histo " << Form("h_D1_FRONT_E%d", k+1) << endl;
      hist->GetXaxis()->SetRangeUser(rmax, pmax);

      // peak search
      TSpectrum *s = new TSpectrum();
      npeaks = s->Search(hist, 8, "", 0.03); // hist, sigma, option, threshold
      cout << "Found " << npeaks << " peaks to fit" << endl;

      // get list of peaks
      Double_t *xpeaks = s->GetPositionX();
      peakList.clear();

      // print in terminal list of found peaks
      // loop on sorted xpeaks
      for (Int_t i = 0; i < npeaks; ++i) { 
        peakList.push_back(xpeaks[i]);
        cout << i << "\t" << xpeaks[i] << endl;
      }

      // order peaks and get number of peaks
      sort(peakList.begin(), peakList.end());

      // remove peaks between range 1 and 2 and recalculate number of peaks
      //       reject = kTRUE;
      //       if (reject) RemovePeak(rmin, rmax);
      //       npeaks = peakList.size();
      //       cout << "Remove range from " << rmin << " to " << rmax << endl;

      // Add peaks when doublet
      AddPeak(430, rmax, 500); 
      AddPeak(885, 700, pmax); 

      // 7 peaks found for strip 31
      if (k == 30) RemovePeak(700,770);

      // order peaks and get number of peaks
      sort(peakList.begin(), peakList.end());
      npeaks = peakList.size();
      cout << endl;
      cout << "Found " << npeaks << " peaks to fit" << endl;

      // Loop on all found peaks
      Double_t par[3000];
      par[0] = 8;  // common width for all states
      par[1] = 0;  // common exponential decay for all states
      for (Int_t p = 0; p < npeaks; p++) {
        Double_t xp = peakList[p];
        Int_t bin  = hist->GetXaxis()->FindBin(xp);
        Double_t yp = hist->GetBinContent(bin);
        par[2*p+2] = xp;
        par[2*p+3] = yp;
      }

      if(Fit_withBGLin) {
        par[2*npeaks+2] = 150.; // y-intercept
        par[2*npeaks+3] = -0.15; // slope
      }

      // define fit function and parameters
      Int_t nparams = 2*npeaks + 2; // gaussian width and exponential tail parameters
      if (Fit_withBGLin) nparams += 2;
      TF1 *fit = new TF1("fit", fpeaks, rmax, pmax, nparams);
      TVirtualFitter::Fitter(hist, nparams);
      fit->SetParameters(par);

      // constrain parameters
      fit->SetParLimits(0, 1, 20); // width
      fit->FixParameter(1, 0); // to fit only with gaussian
      for (Int_t p = 0; p < npeaks; p++) {
        cout << "niveau " << par[2*p+2] << endl;
        /*      if (par[2*p+2] < (rmin+rmax)/2) {
                fit->SetParLimits(2*p+2, pmin, pmax);    // position(mean parameter)
                }
                else {
                fit->SetParLimits(2*p+2, pmin2, pmax2);    // position(mean parameter)
                }*/

        if (p==2 || p==5) { // doublet
          fit->SetParLimits(2*p+2, par[2*(p-1)+2]+5, par[2*(p-1)+2]+12); // position(mean parameter) 
          fit->SetParLimits(2*p+3, par[2*(p-1)+3]/6., par[2*(p-1)+3]/4.);        // Amplitude
        }
        else {
          //fit->FixParameter(2*p+2, par[2*p+2]); // position(mean parameter) 
          fit->SetParLimits(2*p+2, par[2*p+2]-2, par[2*p+2]+2); // position(mean parameter) 
          fit->SetParLimits(2*p+3, 0, 1e5);        // Amplitude
        }

      }

      // constrain for background
      if (Fit_withBGLin) {
        fit->SetParLimits(2*npeaks+2, 100, 500);
        fit->SetParLimits(2*npeaks+3, -0.4, -0.1);
      }

      // fit spectrum
      cout << endl;
      cout << "Now fitting: be patient.... " << endl;
      cout << endl;

      fit->SetNpx(1000);
      TFitResultPtr r = hist->Fit("fit", "RS");
      Int_t fitStatus = r;
      cout << "status " << fitStatus << endl;
      cout << endl;
      TMatrixDSym cov = r->GetCovarianceMatrix();
      //   r->Print("V");

      // individual skewed gaussian function
      TF1 *signalFcn = new TF1("signalFcn", indivFcn, rmax, pmax, 4); // "4": nb of parameters
      signalFcn->SetLineColor(kRed);
      signalFcn->SetLineStyle(2);
      signalFcn->SetLineWidth(2);
      signalFcn->SetNpx(500);

      // individual background
      TF1 *signalFcn_BG = new TF1("signalFcn_BG", indivFcn_BG, rmax, pmax, 2);
      signalFcn_BG->SetLineColor(kMagenta);
      signalFcn_BG->SetLineStyle(2);
      signalFcn_BG->SetLineWidth(2);
      signalFcn_BG->SetNpx(500);

      // open output files
      ofstream myfile("peak_param.txt");

      // declare TGraph for calibration
      TGraphErrors* gr_calib= new TGraphErrors(4);

      // draw individual contributions
      Double_t param[8];
      param[2] = fit->GetParameter(0); // gaussian width 
      param[6] = fit->GetParError(0);  // gaussian width uncertainty
      param[3] = fit->GetParameter(1); // expo tail
      param[7] = fit->GetParError(1);  // expo tail uncertainty
      for (Int_t i = 0; i < npeaks; ++i) {   // loop on number of peaks
        // parameters from drawing individual contribution
        param[1] = fit->GetParameter(2 + 2*i);  // position
        param[5] = fit->GetParError(2 + 2*i);   // Erreur position
        param[0] = fit->GetParameter(3 + 2*i);  // amplitude
        param[4] = fit->GetParError(3 + 2*i);   // Erreur amplitude

        // draw
        signalFcn->SetParameters(param);
        signalFcn->DrawCopy("same");

        // write peak parameters in myfile
        myfile << param[1] << "\t"<< param[5] << endl;

        // Fill TGraph
        gr_calib->SetPoint(i, param[1], Energy[i]);
        gr_calib->SetPointError(i, param[5], error_E[i]);

      } // end loop on number of peaks

      // remove small doublet peak
      gr_calib->RemovePoint(2);
      gr_calib->RemovePoint(4);

      // fit 
      TF1 *fit1 = new TF1("fit1", fit_pol1, 300, 1000, 2);
      fit1->SetParNames("p0","p1");
      fit1->SetParameters(0,1.1);
      TFitResultPtr graph=gr_calib->Fit(fit1,"RS");
      graph->Print("V");

      // write fit param in txt file
      ofstream fitParam_file("DSSSD_calibration.txt", ios::app);
      fitParam_file << Form("COMPTONTELESCOPE_D1_STRIP_FRONT%d_E",k) << " " << fit1->GetParameter(0) << " " << fit1->GetParameter(1) << endl;
      //fitParam_file << Form("COMPTONTELESCOPE_D1_STRIP_FRONT%d_E",k) << " " << fit1->GetParameter(1) << " " << fit1->GetParameter(0) << endl;
      fitParam_file.close();


      // write graph in root file
      TFile *calibFile = new TFile("graph_calib_207Bi_spectrum.root","update");
      gr_calib->SetNameTitle(Form("grCalib_D1_Front%d", k), Form("D1_Front%d", k));
      gr_calib->GetXaxis()->SetTitle("position (channel)");
      gr_calib->GetYaxis()->SetTitle("Energy (MeV)");
      gr_calib->Write();
      calibFile->Close();

      Double_t paramBG[4];
      paramBG[0] = fit->GetParameter(2*npeaks+2); // y-intercept 
      paramBG[2] = fit->GetParError(2*npeaks+2);  // y-intercept uncertainty
      paramBG[1] = fit->GetParameter(2*npeaks+3); // slope
      paramBG[3] = fit->GetParError(2*npeaks+3);  // slope uncertainty

      signalFcn_BG->SetParameters(paramBG);
      signalFcn_BG->DrawCopy("same");

      // close files
      myfile.close();

      //write in a root file
      TFile *f = new TFile("./Histograms/Fit_207Bi_spectrum.root","UPDATE");
      c1->SetName(Form("h_D1_FRONT_E%d", k));
      c1->Write();
      //hist->Write();
      //signalFcn->Write();
      f->Close();
    }
}


  // Back strips
  for (int k = 0; k < NBSTRIPS; k++)
  {

    // change range
    rmax = 370;

    // open histo
    TH1F *hist = (TH1F*) inFile->Get(Form("h_D1_BACK_E%d", k+1));
    cout << "Fitting histo " << Form("h_D1_BACK_E%d", k+1) << endl;
    hist->GetXaxis()->SetRangeUser(rmax, pmax);

    // peak search
    TSpectrum *s = new TSpectrum();
    npeaks = s->Search(hist, 8, "", 0.03); // hist, sigma, option, threshold
    cout << "Found " << npeaks << " peaks to fit" << endl;

    // get list of peaks
    Double_t *xpeaks = s->GetPositionX();
    peakList.clear();

    // print in terminal list of found peaks
    // loop on sorted xpeaks
    for (Int_t i = 0; i < npeaks; ++i) { 
      peakList.push_back(xpeaks[i]);
      cout << i << "\t" << xpeaks[i] << endl;
    }

    // order peaks and get number of peaks
    sort(peakList.begin(), peakList.end());

    // remove peaks between range 1 and 2 and recalculate number of peaks
    //       reject = kTRUE;
    //       if (reject) RemovePeak(rmin, rmax);
    //       npeaks = peakList.size();
    //       cout << "Remove range from " << rmin << " to " << rmax << endl;

    // Add peaks when doublet
    AddPeak(480, rmax, 500); 
    AddPeak(930, 700, pmax); 

    // 7 peaks found for strip 26
    if (k == 24 || k == 25) RemovePeak(700,830);

    // order peaks and get number of peaks
    sort(peakList.begin(), peakList.end());
    npeaks = peakList.size();
    cout << endl;
    cout << "Found " << npeaks << " peaks to fit" << endl;

    // Loop on all found peaks
    Double_t par[3000];
    par[0] = 8;  // common width for all states
    par[1] = 0;  // common exponential decay for all states
    for (Int_t p = 0; p < npeaks; p++) {
      Double_t xp = peakList[p];
      Int_t bin  = hist->GetXaxis()->FindBin(xp);
      Double_t yp = hist->GetBinContent(bin);
      par[2*p+2] = xp;
      par[2*p+3] = yp;
    }

    if(Fit_withBGLin) {
      par[2*npeaks+2] = 150.; // y-intercept
      par[2*npeaks+3] = -0.15; // slope
    }

    // define fit function and parameters
    Int_t nparams = 2*npeaks + 2; // gaussian width and exponential tail parameters
    if (Fit_withBGLin) nparams += 2;
    TF1 *fit = new TF1("fit", fpeaks, rmax, pmax, nparams);
    TVirtualFitter::Fitter(hist, nparams);
    fit->SetParameters(par);

    // constrain parameters
    fit->SetParLimits(0, 1, 20); // width
    fit->FixParameter(1, 0); // to fit only with gaussian
    for (Int_t p = 0; p < npeaks; p++) {
      cout << "niveau " << par[2*p+2] << endl;
      /*      if (par[2*p+2] < (rmin+rmax)/2) {
              fit->SetParLimits(2*p+2, pmin, pmax);    // position(mean parameter)
              }
              else {
              fit->SetParLimits(2*p+2, pmin2, pmax2);    // position(mean parameter)
              }*/

      if (p==2 || p==5) { // doublet
        fit->SetParLimits(2*p+2, par[2*(p-1)+2]+5, par[2*(p-1)+2]+12); // position(mean parameter) 
        fit->SetParLimits(2*p+3, par[2*(p-1)+3]/6., par[2*(p-1)+3]/4.);        // Amplitude
      }
      else {
        //fit->FixParameter(2*p+2, par[2*p+2]); // position(mean parameter) 
        fit->SetParLimits(2*p+2, par[2*p+2]-2, par[2*p+2]+2); // position(mean parameter) 
        fit->SetParLimits(2*p+3, 0, 1e5);        // Amplitude
      }

    }

    // constrain for background
    if (Fit_withBGLin) {
      fit->SetParLimits(2*npeaks+2, 100, 500);
      fit->SetParLimits(2*npeaks+3, -0.4, -0.1);
    }

    // fit spectrum
    cout << endl;
    cout << "Now fitting: be patient.... " << endl;
    cout << endl;

    fit->SetNpx(1000);
    TFitResultPtr r = hist->Fit("fit", "RS");
    Int_t fitStatus = r;
    cout << "status " << fitStatus << endl;
    cout << endl;
    TMatrixDSym cov = r->GetCovarianceMatrix();
    //   r->Print("V");

    // individual skewed gaussian function
    TF1 *signalFcn = new TF1("signalFcn", indivFcn, rmax, pmax, 4); // "4": nb of parameters
    signalFcn->SetLineColor(kRed);
    signalFcn->SetLineStyle(2);
    signalFcn->SetLineWidth(2);
    signalFcn->SetNpx(500);

    // individual background
    TF1 *signalFcn_BG = new TF1("signalFcn_BG", indivFcn_BG, rmax, pmax, 2);
    signalFcn_BG->SetLineColor(kMagenta);
    signalFcn_BG->SetLineStyle(2);
    signalFcn_BG->SetLineWidth(2);
    signalFcn_BG->SetNpx(500);

    // open output files
    ofstream myfile("peak_param.txt");

    // declare TGraph for calibration
    TGraphErrors* gr_calib= new TGraphErrors(4);

    // draw individual contributions
    Double_t param[8];
    param[2] = fit->GetParameter(0); // gaussian width 
    param[6] = fit->GetParError(0);  // gaussian width uncertainty
    param[3] = fit->GetParameter(1); // expo tail
    param[7] = fit->GetParError(1);  // expo tail uncertainty
    for (Int_t i = 0; i < npeaks; ++i) {   // loop on number of peaks
      // parameters from drawing individual contribution
      param[1] = fit->GetParameter(2 + 2*i);  // position
      param[5] = fit->GetParError(2 + 2*i);   // Erreur position
      param[0] = fit->GetParameter(3 + 2*i);  // amplitude
      param[4] = fit->GetParError(3 + 2*i);   // Erreur amplitude

      // draw
      signalFcn->SetParameters(param);
      signalFcn->DrawCopy("same");

      // write peak parameters in myfile
      myfile << param[1] << "\t"<< param[5] << endl;

      // Fill TGraph
      gr_calib->SetPoint(i, param[1], Energy[i]);
      gr_calib->SetPointError(i, param[5], error_E[i]);

    } // end loop on number of peaks

    // remove small doublet peak
    gr_calib->RemovePoint(2);
    gr_calib->RemovePoint(4);

    // fit 
    TF1 *fit1 = new TF1("fit1", fit_pol1, 300, 1000, 2);
    fit1->SetParNames("p0","p1");
    fit1->SetParameters(0.025,0.0011);
    TFitResultPtr graph=gr_calib->Fit(fit1,"RS");
    graph->Print("V");

    // write fit param in txt file
    ofstream fitParam_file("DSSSD_calibration.txt", ios::app);
    fitParam_file << Form("COMPTONTELESCOPE_D1_STRIP_BACK%d_E",k) << " " << fit1->GetParameter(0) << " " << fit1->GetParameter(1) << endl;
    //fitParam_file << Form("COMPTONTELESCOPE_D1_STRIP_BACK%d_E",k) << " " << fit1->GetParameter(1) << " " << fit1->GetParameter(0) << endl;
    fitParam_file.close();

 
    // write graph in root file
    TFile *calibFile = new TFile("graph_calib_207Bi_spectrum.root","update");
    gr_calib->SetNameTitle(Form("grCalib_D1_Back%d", k), Form("D1_Back%d", k));
    gr_calib->GetXaxis()->SetTitle("position (channel)");
    gr_calib->GetYaxis()->SetTitle("Energy (MeV)");
    gr_calib->Write();
    calibFile->Close();
    
    Double_t paramBG[4];
    paramBG[0] = fit->GetParameter(2*npeaks+2); // y-intercept 
    paramBG[2] = fit->GetParError(2*npeaks+2);  // y-intercept uncertainty
    paramBG[1] = fit->GetParameter(2*npeaks+3); // slope
    paramBG[3] = fit->GetParError(2*npeaks+3);  // slope uncertainty

    signalFcn_BG->SetParameters(paramBG);
    signalFcn_BG->DrawCopy("same");

    // close files
    myfile.close();

    //write in a root file
    TFile *f = new TFile("./Histograms/Fit_207Bi_spectrum.root","UPDATE");
    c1->SetName(Form("h_D1_BACK_E%d", k));
    c1->Write();
    //hist->Write();
    //signalFcn->Write();
    f->Close();

  }

}


// add peak at pos if not detected for the fit
void AddPeak(Double_t pos, Double_t pmin, Double_t pmax)
{
  if (pos>pmin && pos<pmax) {
    peakList.push_back(pos);
    cout << "Added 1 peak to fit: \t" << pos << endl;
  }
}


// remove peaks between brmin and brmax
void RemovePeak(Double_t brmin, Double_t brmax)
{
   Int_t last = 0;
   Int_t number = 0;
   for (UInt_t i = 0; i < peakList.size(); ++i) {   // loop on found peaks
      if (peakList[i]>brmin && peakList[i]<brmax) {
         cout << "Removed 1 peak to fit: \t" << peakList[i] << endl;
         last = i;
         number++;
      }
   } // end loop on found peaks

   peakList.erase(peakList.begin()+last-number+1, peakList.begin()+last+1);
}


// fit fonction without background
Double_t fpeaks(Double_t *x, Double_t *par) 
{
   if (reject && x[0] > (rmin) && x[0] < (rmax)) {
      TF1::RejectPoint();
      return 0;
   }

   Double_t result = 0;

   // skewed gaussian
   Double_t width = par[0];
   Double_t tail = par[1];
   for (Int_t p = 0; p < npeaks; p++) {
      Double_t mean  = par[2*p+2];
      Double_t norm  = par[2*p+3];
      Double_t arg = (x[0] - mean) / width;
      // if "small tail" : use a simple gaussian fonction
      if (tail < 1e-1) {
      //if (tail < 1e-6) {
         result += norm * TMath::Gaus(x[0], mean, width);
      }
      else {
         result += norm * TMath::Exp((x[0] - mean)/tail) * TMath::Erfc(arg + width/2/tail);
      }
   }

   // Linear background
   if (Fit_withBGLin) {
     result += par[2*npeaks+2] + par[2*npeaks+3] * x[0];
   }

   return result;
}


// gaussian fit function
Double_t gaussianPeak(Double_t *x, Double_t *par)
{ 
   Double_t arg = (x[0] - par[1]) / par[2];
   return par[0] * TMath::Exp(-0.5 * arg*arg);
}

// individual fit
Double_t indivFcn(Double_t *x, Double_t *par)
{
   Double_t arg = (x[0] - par[1]) / par[2];

   Double_t result = 0;
   if (par[3] < 1e-1) {    // gaussian
     result = par[0] * TMath::Gaus(x[0], par[1], par[2]);
   }
   else {   // skewed gaussian
      result = par[0] * TMath::Exp((x[0] - par[1])/par[3]) * TMath::Erfc(arg + par[2]/2/par[3]);
   }

   return result;
}

Double_t indivFcn_BG(Double_t *x, Double_t *par)
{
   Double_t result = 0;
   result = par[0] + par[1]*x[0];

   return result;
}

//// fit polynomial degree 1
Double_t fit_pol1(Double_t *x, Double_t *par)
{
  Double_t fit_pol1 = par[0] + par[1]*x[0];
  return fit_pol1;
}

Double_t CalibUncertainty(Double_t x, TF1* f, TMatrixDSym cov)
{
   // get number of TF1 parameters
   Int_t npar = f->GetNpar();

   // define array with derivative wrt parameter number
   Double_t grad[npar];

   // estimate derivative at x value
   f->GradientPar(&x, grad);

   // calculate error
   Double_t result = 0;
   for (Int_t i = 0; i < npar; ++i) {   // loop on parameters
      //      cout << grad[i] << endl;
      for (Int_t j = 0; j < npar; ++j) {   // loop on parameters
         //         cout << i << "\t" << j << "\t" << cov[i][j] << "\t" << grad[i]*grad[j] << "\t" << grad[i]*grad[j] * cov[i][j] << endl;
         result += grad[i]*grad[j] * cov[i][j];
         //         if (i == j) result += grad[i]*grad[j] * cov[i][j];
      } // end loop on parameters
   } // end loop on parameters

   return sqrt(result);
}



