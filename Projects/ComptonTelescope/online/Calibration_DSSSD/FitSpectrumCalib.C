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
Double_t fpeaks_BG(Double_t*, Double_t*);
Double_t gaussianPeak(Double_t*, Double_t*);
Double_t indivFcn(Double_t*, Double_t*);
Double_t indivFcn_BG(Double_t*, Double_t*);
Double_t fit_pol1(Double_t*, Double_t*);
Double_t fit_pol2(Double_t*, Double_t*);
void     AddPeak(Double_t, Double_t, Double_t);
void     RemovePeak(Double_t, Double_t);
Double_t CalibUncertainty(Double_t, TF1*, TMatrixDSym);

static vector<Double_t>  peakList, peakListP;
static Int_t            npeaks;
static Int_t            pmin, pmax, rmin, rmax;
Bool_t reject;

#define NBSTRIPS  32


void FitSpectrumCalib()
{

  const int dimE = 6;
  double Energy[dimE] = {481.6935e-3, 553.8372e-3, 565.8473e-3, 975.651e-3, 1047.795e-3, 1059.805e-3}; // MeV
  double error_E[dimE] = {0.0021e-3, 0.0021e-3, 0.0021e-3, 0.003e-3, 0.003e-3, 0.003e-3};

  // range for pedestal
  pmin = 40;
  rmin = 70;
  // range for spectrum
  rmax = 300;
  pmax = 980;

  // misc.
  peakList.clear();
  peakListP.clear();

  // define canvas
  TCanvas *c1 = new TCanvas("c1", "spectrum", 1000, 700);
  TCanvas *c2 = new TCanvas("c2", "pedestal", 1000, 700);
  // declare TGraph for calibration
  TGraphErrors* gr_calib = new TGraphErrors();

  // open file and get histo
  TFile *inFile = new TFile("20200128_10h44_bi207_conv_RawDSSSDHistos.root");

  // Front strips
  for (int k = 0; k < NBSTRIPS; k++)
  {

    cout << "\n Fitting histo " << Form("h_D1_FRONT_E%d", k+1) << endl;

    ///// Pedestal fit /////
    cout << "\n" << "/////////// Pedestal fit  //////////////" << endl;
    c2->cd();

    // peak search
    TH1F *histP = (TH1F*) inFile->Get(Form("h_D1_FRONT_E%d", k+1));
    histP->GetXaxis()->SetRangeUser(pmin, rmin);
    TSpectrum *sP = new TSpectrum();
    npeaks = sP->Search(histP, 3, "nobackground", 0.5); // hist, sigma, option, threshold
    cout << "Found " << npeaks << " peaks to fit" << endl;
    Double_t *xpeaksP = sP->GetPositionX();
    peakListP.clear();
    for (Int_t i = 0; i < npeaks; ++i) {
      peakListP.push_back(xpeaksP[i]);
      cout << i << "\t" << xpeaksP[i] << endl;
    }
    sort(peakListP.begin(), peakListP.end());

    // fit first peak with gaussian function, no linear background
    Double_t parP[100];      
    parP[0] = 2.;  // width        
    parP[1] = 0;  // exponential decay
    //parP[1] = 0.4;  // exponential decay
    for (Int_t p = 0; p < 1; p++) {
      Double_t xp = peakListP[p];
      cout << "xp " << xp << endl;
      Int_t bin  = histP->GetXaxis()->FindBin(xp);
      Double_t yp = histP->GetBinContent(bin);           
      parP[2] = xp; // mean
      parP[3] = yp; // amplitude
    }
    TF1 *fitP = new TF1("fitP", fpeaks, parP[2]-3, parP[2]+3, 4);
    TVirtualFitter::Fitter(histP, 4);
    fitP->SetParameters(parP);
    fitP->SetParLimits(0, 1, 20); // width
    fitP->FixParameter(1, 0); // no exponential decay
    //fitP->SetParLimits(1, 0.01, 10); // exponential decay
    fitP->SetParLimits(2, parP[2]-10, parP[2]+10); // mean
    fitP->SetParLimits(3, 1e2, 1e10); // amplitude
    fitP->SetNpx(1000);
    TFitResultPtr rP = histP->Fit("fitP", "RS");
    Int_t fitStatusP = rP;
    cout << "status " << fitStatusP << endl;
    cout << endl;
    TMatrixDSym covP = rP->GetCovarianceMatrix();

    // write pedestal fit in root file
    TFile *fP = new TFile("./fit/Fit_pedestal_207Bi_spectrum.root","UPDATE");
    c2->SetName(Form("h_D1_FRONT_E%d_ped", k));
    c2->Write();
    fP->Close();

    // Add in TGraph
    Double_t paramP[2];
    paramP[0] = fitP->GetParameter(2); // mean
    paramP[1] = fitP->GetParError(2); // error mean
    gr_calib->SetPoint(0, paramP[0], 0);
    gr_calib->SetPointError(0, paramP[1], 0);
    gr_calib->Draw();


    if (k != 1) // remove strip front 2 with no data
    {

      //////// Spectrum fit //////
      cout << "\n" << "/////////// Spectrum //////////////" << endl;
      c1->cd();

      // open histo
      TH1F *hist = (TH1F*) inFile->Get(Form("h_D1_FRONT_E%d", k+1));
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

      par[2*npeaks+2] = 150.; // y-intercept
      par[2*npeaks+3] = -0.15; // slope

      // define fit function and parameters
      Int_t nparams = 2*npeaks + 4; // gaussian width and exponential tail parameters + background
      TF1 *fit = new TF1("fit", fpeaks_BG, rmax, pmax, nparams);
      TVirtualFitter::Fitter(hist, nparams);
      fit->SetParameters(par);

      // constrain parameters
      fit->SetParLimits(0, 1, 20); // width
      fit->FixParameter(1, 0); // to fit only with gaussian
      for (Int_t p = 0; p < npeaks; p++) {
        cout << "niveau " << par[2*p+2] << endl;
        //if (par[2*p+2] < (rmin+rmax)/2) {
        //fit->SetParLimits(2*p+2, pmin, pmax);    // position(mean parameter)
        //}
        //else {
        //fit->SetParLimits(2*p+2, pmin2, pmax2);    // position(mean parameter)
        //}

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
      fit->SetParLimits(2*npeaks+2, 100, 500);
      fit->SetParLimits(2*npeaks+3, -0.4, -0.1);

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

        // Fill TGraph
        gr_calib->SetPoint(i+1, param[1], Energy[i]);
        gr_calib->SetPointError(i+1, param[5], error_E[i]);

      } // end loop on number of peaks

      // remove small doublet peak
      gr_calib->RemovePoint(3);
      gr_calib->RemovePoint(5);

      // Background parameters
      Double_t paramBG[4];
      paramBG[0] = fit->GetParameter(2*npeaks+2); // y-intercept 
      paramBG[2] = fit->GetParError(2*npeaks+2);  // y-intercept uncertainty
      paramBG[1] = fit->GetParameter(2*npeaks+3); // slope
      paramBG[3] = fit->GetParError(2*npeaks+3);  // slope uncertainty

      signalFcn_BG->SetParameters(paramBG);
      signalFcn_BG->DrawCopy("same");

      //write in a root file
      TFile *f = new TFile("./fit/Fit_noPedestal_207Bi_spectrum.root","UPDATE");
      c1->SetName(Form("h_D1_FRONT_E%d", k));
      c1->Write();
      f->Close();


      ////////// calibration //////////
      cout << "\n" << "/////////// calibration //////////////" << endl;

      // fit TGraph pol2
      TF1 *fitP2 = new TF1("fitP2", fit_pol2, 0, 1000, 3);
      fitP2->SetParNames("p0","p1","p2");
      //      fitP2->SetParLimits(0,
      fitP2->SetParameters(-0.09,0.001,-5e-11);
      TFitResultPtr graphP2 = gr_calib->Fit(fitP2, "RS");
      graphP2->Print("V");

      // write fit param in txt file
      ofstream fitParam_fileP2("./calib/DSSSD_calibration_withPed_pol2.txt", ios::app);
      fitParam_fileP2 << Form("COMPTONTELESCOPE_D1_STRIP_FRONT%d_E",k) << " " << fitP2->GetParameter(0) << " " << fitP2->GetParameter(1) << " " << fitP2->GetParameter(2) << endl;
      fitParam_fileP2.close();

      // write graph in root file
      TFile *calibFileP2 = new TFile("./calib/graph_calib_207Bi_spectrum_withPed_pol2.root", "update");
      gr_calib->SetNameTitle(Form("grCalib_D1_Front%d", k), Form("D1_Front%d", k));
      gr_calib->GetXaxis()->SetTitle("position (channel)");
      gr_calib->GetYaxis()->SetTitle("Energy (MeV)");
      gr_calib->Write();
      calibFileP2->Close();

    }

  }


  /////// Back strips
  for (int k = 0; k < NBSTRIPS; k++)
  {

    // change range
    pmin = 30;
    rmax = 370; 

    cout << "\n Fitting histo " << Form("h_D1_BACK_E%d", k+1) << endl;

    ///// Pedestal fit /////
    cout << "\n" << "/////////// Pedestal fit  //////////////" << endl;
    c2->cd();

    // peak search
    TH1F *histP = (TH1F*) inFile->Get(Form("h_D1_BACK_E%d", k+1));
    histP->GetXaxis()->SetRangeUser(pmin, rmin);
    TSpectrum *sP = new TSpectrum();
    npeaks = sP->Search(histP, 3, "nobackground", 0.5); // hist, sigma, option, threshold
    cout << "Found " << npeaks << " peaks to fit" << endl;
    Double_t *xpeaksP = sP->GetPositionX();
    peakListP.clear();
    for (Int_t i = 0; i < npeaks; ++i) {
      peakListP.push_back(xpeaksP[i]);
      cout << i << "\t" << xpeaksP[i] << endl;
    }
    sort(peakListP.begin(), peakListP.end());

    // fit first peak with gaussian function, no linear background
    Double_t parP[100];      
    parP[0] = 2.;  // width        
    parP[1] = 0;  // exponential decay
    //parP[1] = 0.4;  // exponential decay
    for (Int_t p = 0; p < 1; p++) {
      Double_t xp = peakListP[p];
      cout << "xp " << xp << endl;
      Int_t bin  = histP->GetXaxis()->FindBin(xp);
      Double_t yp = histP->GetBinContent(bin);           
      parP[2] = xp; // mean
      parP[3] = yp; // amplitude
    }
    TF1 *fitP = new TF1("fitP", fpeaks, parP[2]-6, parP[2]+6, 4);
    TVirtualFitter::Fitter(histP, 4);
    fitP->SetParameters(parP);
    fitP->SetParLimits(0, 1, 20); // width
    fitP->FixParameter(1, 0); // no exponential decay
    //fitP->SetParLimits(1, 0.01, 10); // exponential decay
    fitP->SetParLimits(2, parP[2]-10, parP[2]+10); // mean
    fitP->SetParLimits(3, 1e2, 1e10); // amplitude
    fitP->SetNpx(1000);
    TFitResultPtr rP = histP->Fit("fitP", "RS");
    Int_t fitStatusP = rP;
    cout << "status " << fitStatusP << endl;
    cout << endl;
    TMatrixDSym covP = rP->GetCovarianceMatrix();

    // write pedestal fit in root file
    TFile *fP = new TFile("./fit/Fit_pedestal_207Bi_spectrum.root","UPDATE");
    c2->SetName(Form("h_D1_BACK_E%d_ped", k));
    c2->Write();
    fP->Close();

    // Add in TGraph
    Double_t paramP[2];
    paramP[0] = fitP->GetParameter(2); // mean
    paramP[1] = fitP->GetParError(2); // error mean
    gr_calib->SetPoint(0, paramP[0], 0);
    gr_calib->SetPointError(0, paramP[1], 0);
    gr_calib->Draw();


    //////// Spectrum fit //////
    cout << "\n" << "/////////// Spectrum //////////////" << endl;
    c1->cd();

    // open histo
    TH1F *hist = (TH1F*) inFile->Get(Form("h_D1_BACK_E%d", k+1));
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

    // 7 peaks found  
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

    par[2*npeaks+2] = 150.; // y-intercept
    par[2*npeaks+3] = -0.15; // slope

    // define fit function and parameters
    Int_t nparams = 2*npeaks + 4; // gaussian width and exponential tail parameters + background
    TF1 *fit = new TF1("fit", fpeaks_BG, rmax, pmax, nparams);
    TVirtualFitter::Fitter(hist, nparams);
    fit->SetParameters(par);

    // constrain parameters
    fit->SetParLimits(0, 1, 20); // width
    fit->FixParameter(1, 0); // to fit only with gaussian
    for (Int_t p = 0; p < npeaks; p++) {
      cout << "niveau " << par[2*p+2] << endl;
      //if (par[2*p+2] < (rmin+rmax)/2) {
      //fit->SetParLimits(2*p+2, pmin, pmax);    // position(mean parameter)
      //}
      //else {
      //fit->SetParLimits(2*p+2, pmin2, pmax2);    // position(mean parameter)
      //}

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
    fit->SetParLimits(2*npeaks+2, 100, 500);
    fit->SetParLimits(2*npeaks+3, -0.4, -0.1);

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

      // Fill TGraph
      gr_calib->SetPoint(i+1, param[1], Energy[i]);
      gr_calib->SetPointError(i+1, param[5], error_E[i]);

    } // end loop on number of peaks

    // remove small doublet peak
    gr_calib->RemovePoint(3);
    gr_calib->RemovePoint(5);

    // Background parameters
    Double_t paramBG[4];
    paramBG[0] = fit->GetParameter(2*npeaks+2); // y-intercept 
    paramBG[2] = fit->GetParError(2*npeaks+2);  // y-intercept uncertainty
    paramBG[1] = fit->GetParameter(2*npeaks+3); // slope
    paramBG[3] = fit->GetParError(2*npeaks+3);  // slope uncertainty

    signalFcn_BG->SetParameters(paramBG);
    signalFcn_BG->DrawCopy("same");

    //write in a root file
    TFile *f = new TFile("./fit/Fit_noPedestal_207Bi_spectrum.root","UPDATE");
    c1->SetName(Form("h_D1_BACK_E%d", k));
    c1->Write();
    f->Close();


    ////////// calibration //////////
    cout << "\n" << "/////////// calibration //////////////" << endl;

    // fit TGraph pol2
    TF1 *fitP2 = new TF1("fitP2", fit_pol2, 0, 1000, 3);
    fitP2->SetParNames("p0","p1","p2");
    //      fitP2->SetParLimits(0,
    fitP2->SetParameters(-0.09,0.001,-5e-11);
    TFitResultPtr graphP2 = gr_calib->Fit(fitP2, "RS");
    graphP2->Print("V");

    // write fit param in txt file
    ofstream fitParam_fileP2("./calib/DSSSD_calibration_withPed_pol2.txt", ios::app);
    fitParam_fileP2 << Form("COMPTONTELESCOPE_D1_STRIP_BACK%d_E",k) << " " << fitP2->GetParameter(0) << " " << fitP2->GetParameter(1) << " " << fitP2->GetParameter(2) << endl;
    fitParam_fileP2.close();

    // write graph in root file
    TFile *calibFileP2 = new TFile("./calib/graph_calib_207Bi_spectrum_withPed_pol2.root", "update");
    gr_calib->SetNameTitle(Form("grCalib_D1_Back%d", k), Form("D1_Back%d", k));
    gr_calib->GetXaxis()->SetTitle("position (channel)");
    gr_calib->GetYaxis()->SetTitle("Energy (MeV)");
    gr_calib->Write();
    calibFileP2->Close();


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

   return result;
}


// fit fonction with background
Double_t fpeaks_BG(Double_t *x, Double_t *par) 
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
   result += par[2*npeaks+2] + par[2*npeaks+3] * x[0];

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

//// fit polynomial degree 2
Double_t fit_pol2(Double_t *x, Double_t *par)
{
  Double_t fit_pol2 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  return fit_pol2;
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



