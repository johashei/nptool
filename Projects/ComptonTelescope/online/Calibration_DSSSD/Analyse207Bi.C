///////////////////////////////////////////////////////////////////////////////
// This macro calibrates DSSSDs with a 207Bi source.                         //
//                                                                           //
// It treats either all ASIC channels (channel = -1) or a single one         //
// (specify channel number). A boolean (isPside) should indicate whether     //
// the input file corresponds to p-side (1) or n-side (0). A pdf file is     //
// generated with the relevant information.                                  //
//                                                                           //
// Use: .L Analyse207Bi.C+                                                   //
//      Analyse207Bi(name, isPside, channel)                                 //
//                                                                           //
// N. de Sereville (September 2019)                                          //
// A. Meyer (January 2021)                                                   //
///////////////////////////////////////////////////////////////////////////////

// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLine.h"
#include "TError.h"
#include "TROOT.h"
#include "TMath.h"

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

// STL headers
#include <vector>
using namespace std;

#define DETECTOR_ID       1
#define NCHANNELS        32
#define BACKGROUND_MIN  300
#define BACKGROUND_MAX 1000

// variables
static std::vector<Double_t> peakListForFit;
static const double branchingLM = 4; //L/M branching ratios

// functions
void AddPeak(Double_t, Double_t, Double_t);
Double_t fpeaks(Double_t*, Double_t*);
Double_t fpeaks2(Double_t*, Double_t*);
Double_t gaussianPeak(Double_t*, Double_t*);


void Analyse207Bi(const char* name = "bb7_3309-7_bi207_20210126_13h09_run5_conv_RawDSSSDHistos.root",
                  Bool_t isPside = 1, Int_t channel = 20)
{
   // no statistical box, no histogram name, no fit info
   if (gStyle->GetOptStat())  gStyle->SetOptStat(0);
   if (gStyle->GetOptTitle()) gStyle->SetOptTitle(0);
   if (gStyle->GetOptFit())   gStyle->SetOptFit(0);
   // remove verbosity when creating pdf file
   gErrorIgnoreLevel = kWarning;

   ///////////////////////////////////////////////////////////////////////////
	// 207Bi information
   std::vector<Double_t> mainLineEnergy = {481.69,  975.65};
   std::vector<Double_t> LLineEnergy    = {553.84, 1047.80};
   std::vector<Double_t> MLineEnergy    = {565.85, 1059.81};
   
   ///////////////////////////////////////////////////////////////////////////
	// open source file
   string sname(name);
	auto f = new TFile(Form("./%s", sname.c_str()));

   ///////////////////////////////////////////////////////////////////////////
   // define main canvas
   auto can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("can");
   if (can) can->Clear();                                                               
   else can = new TCanvas("can", "can", 737, 847);
   Float_t small = 1e-5;
   can->Divide(1, 4, small, small); 

   ///////////////////////////////////////////////////////////////////////////
   // open output pdf file
   string label = sname.substr(4, 6);
   label += sname.substr(16, 21);
   label += (isPside) ? "_pside" : "_nside";
   TString pdfname = Form("Analyse207Bi_%s.pdf", label.c_str());
   can->Print(Form("%s[", pdfname.Data()));

   ///////////////////////////////////////////////////////////////////////////
   // open output calib file
   string fname = (isPside) ? Form("DSSSD_D%d_Calibration_Front_E.txt",DETECTOR_ID) 
                            : Form("DSSSD_D%d_Calibration_Back_E.txt",DETECTOR_ID);
   ofstream calibFile(fname.c_str(), std::ios::out);
   
   ///////////////////////////////////////////////////////////////////////////
   // define some functions
   auto fcalibFull  = new TF1("fcalibFull",  "pol1", 1024);
   auto fcalibLocal = new TF1("fcalibLocal", "pol1", 1024);

   ///////////////////////////////////////////////////////////////////////////
   // define vectors for calibration
   std::vector<Double_t> position, positionError, energy, energyError;

   ///////////////////////////////////////////////////////////////////////////
   // prepare summary information
   // FWHM in keV
   auto grFWHMlow   = new TGraph();
   grFWHMlow->SetMarkerStyle(kFullSquare);
   grFWHMlow->SetMarkerColor(kRed);
   auto grFWHMhigh = new TGraph();
   grFWHMhigh->SetMarkerStyle(kFullCircle);
   grFWHMhigh->SetMarkerColor(kBlue);
   // energy calibration parameters
   auto grOffset = new TGraph();
   grOffset->SetMarkerStyle(kFullSquare);
   grOffset->SetMarkerColor(kRed);
   auto grGain   = new TGraph();
   grGain->SetMarkerStyle(kFullCircle);
   grGain->SetMarkerColor(kBlue);

   ///////////////////////////////////////////////////////////////////////////
   // declare some display histograms
   // sigma histogram
   auto hframe1 = new TH2F("hframe1", "hframe1", 32, 0, 32, 10, 0, 40);
   hframe1->GetXaxis()->SetTitle("strip number");
   hframe1->GetYaxis()->SetTitle("FWHM (keV)");
   hframe1->GetXaxis()->CenterTitle();
   hframe1->GetYaxis()->CenterTitle();
   // calibration histogram
   auto hframe2 = new TH2F("hframe2", "hframe2", 32, 0, 32, 100, 0, 70);
   hframe2->GetXaxis()->SetTitle("strip number");
   hframe2->GetXaxis()->CenterTitle();

   ///////////////////////////////////////////////////////////////////////////
	// treat all channels 
   for (Int_t n = 0; n < NCHANNELS; ++n) {   // loop on channels
      // treat all channels or specified one
      if (channel < 0 || channel == n) {
         cout << "\n--------------- Analysing channel " << n << " ---------------\n";

         ///////////////////////////////////////////////////////////////////////////
         // get histogram
         string hname = (isPside) ? Form("h_D%d_FRONT_E%d",DETECTOR_ID,n+1) 
                                  : Form("h_D%d_BACK_E%d",DETECTOR_ID,n+1);
         auto h = (TH1F*) f->Get(hname.c_str());
         h->Sumw2();
         // draw histogram
         can->cd(1);
         auto pad = (TPad*) can->FindObject("can_1");                                  
         pad->SetLogy();
//         pad->SetBottomMargin(1e-5);
         h->GetYaxis()->CenterTitle();
         h->DrawCopy();

         ///////////////////////////////////////////////////////////////////////////
         // clear variables
         position.clear(); positionError.clear();
         energy.clear();   energyError.clear();

         ///////////////////////////////////////////////////////////////////////////
         // find pedestal position 
         // search
         auto s = new TSpectrum();
         Int_t npeaks = s->Search(h, 6, "", 0.6);
         h->DrawCopy();
         if (npeaks == 1) cout << ". Pedestal found -> OK\n";
         else cout << ". Pedestal PROBLEM\n";
         Double_t *xpeaks = s->GetPositionX();
         // fit
         h->Fit("gaus", "RQ", "", xpeaks[0]-20, xpeaks[0]+20);
         h->DrawCopy();
         // parameters
         auto fitF = h->GetFunction("gaus");
         position.push_back(fitF->GetParameter(1));
         positionError.push_back(fitF->GetParError(1));
         energy.push_back(0 * 1e-3);
         energyError.push_back(1 * 1e-3);
         
         ///////////////////////////////////////////////////////////////////////////
         // substract background
         h->GetXaxis()->SetRangeUser(BACKGROUND_MIN, BACKGROUND_MAX);
         auto background = (TH1F*) h->ShowBackground(50, "");
         background->SetLineColor(kRed);
         background->Draw("same");

         ///////////////////////////////////////////////////////////////////////////
         // draw subtracted spectrum 
         h->Add(background, -1);
         can->cd(2);
         pad = (TPad*) can->FindObject("can_2");                                  
//         pad->SetTopMargin(1e-5);
         h->DrawCopy();
         
         ///////////////////////////////////////////////////////////////////////////
         // search peaks after background subtraction
         npeaks = s->Search(h, 6, "", 0.3);
         h->DrawCopy();
         cout << ". " << npeaks << " main transitions found ";
         if (npeaks == 2) cout << "\t-> OK!\n";
         else cout << "\t-> PROBLEM!\n";
         // get peak position
         xpeaks = s->GetPositionX();

         // order peaks position
         vector<Double_t> peakList;
         for (Int_t i = 0; i < npeaks; ++i) {   // loop on sorted xpeaks 
            peakList.push_back(xpeaks[i]); 
         } // end loop on sorted xpeaks 
         // try to guess which peak is detected and add other one
         // uses calibration parameters from previous strip calibration
         if (npeaks == 1) {
            Double_t energy = fcalibFull->GetParameter(0) + fcalibFull->GetParameter(1)*xpeaks[0];
            Double_t energyGuess = 0;
            if (energy > mainLineEnergy[0]+200) {
               energyGuess = (mainLineEnergy[0]-fcalibFull->GetParameter(0))/fcalibFull->GetParameter(1);;
            }
            else {
               energyGuess = (mainLineEnergy[1]-fcalibFull->GetParameter(0))/fcalibFull->GetParameter(1);;
            }
            AddPeak(energyGuess, BACKGROUND_MIN, BACKGROUND_MAX);
         }
         // order peak position
         sort(peakList.begin(), peakList.end());

         ///////////////////////////////////////////////////////////////////////////
         // rough energy calibration
         auto grcalib = new TGraph(peakList.size(), &peakList[0], &mainLineEnergy[0]);
         grcalib->Fit("fcalibFull", "Q0");
         grOffset->SetPoint(n, n, fcalibFull->GetParameter(0));
         grGain  ->SetPoint(n, n, fcalibFull->GetParameter(1)*10);
         
         // loop on "480"- and "975"-keV features
         // p = 0 -> 480 keV
         // p = 1 -> 975 keV
         Double_t sigma480 = 0;
         for (UInt_t p = 0; p < mainLineEnergy.size(); ++p) {
            can->cd(3+p);
            pad = (TPad*) can->FindObject(Form("can_%d", 3+p));                                  
            // determine range for fitting based on rough calibration
            Double_t cmin = (mainLineEnergy[p]-40 -fcalibFull->GetParameter(0)) / fcalibFull->GetParameter(1);
            Double_t cmax = (mainLineEnergy[p]+130-fcalibFull->GetParameter(0)) / fcalibFull->GetParameter(1);
            if (cmax > 1023) cmax = 1023; // prevents detection of spurious peak at channel 1024
//            std::cout << cmin << "\t" << cmax << "\n";
            h->GetXaxis()->SetRangeUser(cmin, cmax);
            h->DrawCopy();

            // search satellite peaks
            npeaks = s->Search(h, 4, "", 0.1);
            if (p == 0) {
               cout << ". Low  energy transitions\n";
            }
            else {
               cout << ". High energy transitions\n";
            }
            cout << "\t. found " << npeaks << " peaks";
            if (npeaks == 2) cout << "\t-> OK!\n";
            else cout << "\t-> PROBLEM!\n";
            // get peak position
            xpeaks = s->GetPositionX();
            // fill peakListForFit vector with peak position
            peakListForFit.clear();
            for (Int_t i = 0; i < npeaks; ++i) {
               peakListForFit.push_back(xpeaks[i]);
            }
            // order peak list
            sort(peakListForFit.begin(), peakListForFit.end());
            // case where 3 peaks are found corresponding to the typical case 
            // where a "replica" peak close to the main peak is observed, and
            // should be removed
            if (npeaks == 3) peakListForFit.erase(peakListForFit.begin()+1);
            // add CE-M component based on local linear fit
            vector<double> energyLocal = {mainLineEnergy[p], LLineEnergy[p]};
            auto grcalibLocal = new TGraph(peakListForFit.size(), &peakListForFit[0], &energyLocal[0]);
            grcalibLocal->Fit("fcalibLocal", "Q0");
            Double_t channel = (MLineEnergy[p] -fcalibLocal->GetParameter(0)) / fcalibLocal->GetParameter(1);
//            Double_t channel = (MLineEnergy[p] -fcalibFull->GetParameter(0)) / fcalibFull->GetParameter(1);
            AddPeak(channel, cmin, cmax);

            // define parameter list for fit
            std::vector<Double_t> paramList;
            // width (sigma); common for all states
            paramList.push_back(6);
            // position and amplitude
            for (UInt_t i = 0; i < peakListForFit.size(); ++i) {
               paramList.push_back(peakListForFit[i]);
               // force content to be positive (can be negative because of
               // background subtraction)
               double binContent = max(0.0, h->GetBinContent(peakListForFit[i]));
               paramList.push_back(binContent);
            }

            // define fit function 
            auto fit = new TF1("fit", fpeaks2, cmin, cmax, paramList.size());
            fit->SetParameters(&paramList[0]);

            // positive amplitudes and position in range
            int dChannel = 5;
            for (UInt_t i = 0; i < peakListForFit.size(); ++i) {
               // +/- 10 channels prevents L and M inversion
               fit->SetParLimits(1 + 2*i, fit->GetParameter(1+2*i)-dChannel,
                                          fit->GetParameter(1+2*i)+dChannel);
               fit->SetParLimits(2 + 2*i, 0, 1e5);
            }

            // fix arbitrary amplitude for M component 
            // parameter is not used in the fit
            fit->FixParameter(6, 50);

            // width of 975 keV line constrainted wrt 480 keV line
            double factor = 1.2; // 20%
            if (p == 1) {
               fit->SetParLimits(0, sigma480/factor, sigma480*factor);
            }

            // fit spectrum
            h->Fit("fit", "RBQ0");
//            h->Fit("fit", "RB0");
            h->DrawCopy();

            // draw individual components
            // function
            auto signalFcn = new TF1("signalFcn", gaussianPeak, cmin, cmax, 3);
            signalFcn->SetLineColor(kRed);
            signalFcn->SetLineStyle(2);
            signalFcn->SetLineWidth(2);
            signalFcn->SetNpx(500);
            // parameters
            Double_t param[3];
            param[2] = fit->GetParameter(0);
            for (UInt_t i = 0; i < peakListForFit.size(); ++i) {
//               param[0] = fit->GetParameter(2 + 2*i);
               param[1] = fit->GetParameter(1 + 2*i);
               if (i == 2) { // M component
                  param[0] = fit->GetParameter(2*i) / branchingLM;
               }
               else {
                  param[0] = fit->GetParameter(2 + 2*i);
               }
               signalFcn->SetParameters(param);
               signalFcn->DrawCopy("same");
            }

            // get resolution for the 480 keV line
            sigma480 = param[2];

            ///////////////////////////////////////////////////////////////////////////
            // fill calibration arrays
            // K-transition
            position.push_back(fit->GetParameter(1));
            positionError.push_back(fit->GetParError(1));
            energy.push_back(mainLineEnergy[p] * 1e-3);
            energyError.push_back(1 * 1e-3);
            // L-transition
            position.push_back(fit->GetParameter(3));
            positionError.push_back(fit->GetParError(3));
            energy.push_back(LLineEnergy[p] * 1e-3);
            energyError.push_back(1 * 1e-3);

            ///////////////////////////////////////////////////////////////////////////
            // fill fwhm data
            if (p) {
               grFWHMhigh->SetPoint(n, n, fcalibFull->GetParameter(1)*2.35*fit->GetParameter(0));
            }
            else {
               grFWHMlow ->SetPoint(n, n, fcalibFull->GetParameter(1)*2.35*fit->GetParameter(0));
            }
         }

         ///////////////////////////////////////////////////////////////////////////
         // display calibration fit
         // search for existing canvas                                                
         auto canC = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("canC");
         if (canC) canC->Clear();                                                               
         else canC = new TCanvas("canC", "canC");
         // fit
         auto grFullCalib= new TGraphErrors(position.size(), 
                                            &position[0], &energy[0],
                                            &positionError[0], &energyError[0]);
         grFullCalib->Draw("a*");
         grFullCalib->GetXaxis()->SetTitle("ADC channel");
         grFullCalib->GetYaxis()->SetTitle("Energy (MeV)");
         grFullCalib->Draw("same");
         auto fitFullCalib = new TF1("fitFullCalib", "pol3", 0, 1024);
         fitFullCalib->SetParameters(-1.35e-1, 2e-3, -1.3e-6, 5e-10);
         grFullCalib->Fit("fitFullCalib", "QR");
         canC->Update();

         ///////////////////////////////////////////////////////////////////////////
         // write calibration coefficients
         string token = (isPside) ? Form("COMPTONTELESCOPE_D%d_STRIP_FRONT%d_E",DETECTOR_ID,n) 
                                  : Form("COMPTONTELESCOPE_D%d_STRIP_BACK%d_E",DETECTOR_ID,n);
         calibFile << token;
         for (int p = 0; p < fitFullCalib->GetNpar(); ++p) {
            calibFile << "\t" << fitFullCalib->GetParameter(p);
         }
         calibFile << "\n";

         ///////////////////////////////////////////////////////////////////////////
         // channel information
         can->cd(1);
         auto tex = new TLatex(0.7, 0.7, Form("channel N%d", n));
         if (isPside) tex = new TLatex(0.7, 0.7, Form("channel P%d", n));
         tex->SetNDC(); tex->SetTextSize(0.1); tex->SetTextFont(42); tex->Draw();

         can->cd(2);
         tex = new TLatex(0.15, 0.8, "background subtracted");
         tex->SetNDC(); tex->SetTextSize(0.1); tex->SetTextFont(42); tex->Draw();

         can->cd(3);
         tex = new TLatex(0.65, 0.7, Form("%3.0f keV + CE L,M", mainLineEnergy[0]));
         tex->SetNDC(); tex->SetTextSize(0.1); tex->SetTextFont(42); tex->Draw();

         can->cd(4);
         tex = new TLatex(0.65, 0.7, Form("%3.0f keV + CE L,M", mainLineEnergy[1]));
         tex->SetNDC(); tex->SetTextSize(0.1); tex->SetTextFont(42); tex->Draw();

         canC->cd();
         tex = new TLatex(0.2, 0.8, Form("channel N%d", n));
         if (isPside) tex = new TLatex(0.2, 0.8, Form("channel P%d", n));
         tex->SetNDC(); tex->SetTextSize(0.05); tex->SetTextFont(42); tex->Draw();

         // fill pdf file
         can->Update();
         can->Print(pdfname.Data());
         canC->Print(pdfname.Data());
      }
   } // end loop on channels

   ///////////////////////////////////////////////////////////////////////////
   // close calibration file
   calibFile.close();

   ///////////////////////////////////////////////////////////////////////////
   // draw analysis result
   // search for existing canvas                                                
   auto canS = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("canS");
   if (canS) canS->Clear();                                                               
   else canS = new TCanvas("canS", "canS", 737, 847);
   canS->Divide(1, 2);

   // resolution vs strip number
   canS->cd(1);
   hframe1->Draw();
   grFWHMlow->Draw("p");
   grFWHMhigh->Draw("p");
	auto leg = new TLegend(0.22, 0.71, 0.38, 0.86);                        
	leg->AddEntry(grFWHMlow,  Form("%3.0f keV", mainLineEnergy[0]), "P");                                 
	leg->AddEntry(grFWHMhigh, Form("%3.0f keV", mainLineEnergy[1]), "P");                                 
	leg->SetBorderSize(1);                                                 
	leg->Draw();            

   // calibration parameters
   canS->cd(2);
   hframe2->Draw();
   grOffset->Draw("p");
   grGain  ->Draw("p");

   // fill pdf file
   canS->Print(pdfname.Data());

   // close pdf file
   can->Print(Form("%s]", pdfname.Data()));
}



void AddPeak(Double_t channel, Double_t cmin, Double_t cmax)
{
   if (channel>cmin && channel<cmax) {
      peakListForFit.push_back(channel);
      std::cout << "\t. Adding 1 peak at ch " << channel << "\n";
   }
}



Double_t fpeaks(Double_t *x, Double_t *par)         
{                                         
   Double_t result = 0;                     

   // gaussian width                              
   Double_t sigma = par[0];
   for (UInt_t p = 0; p < peakListForFit.size(); p++) {
      // case where brho is allowed to vary in fit
      Double_t mean  = par[2*p+1];
      Double_t norm  = par[2*p+2];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }                          

   return result;
}                                  



Double_t fpeaks2(Double_t *x, Double_t *par)         
{                                         
   Double_t result = 0;                     

   // gaussian width                              
   Double_t sigma = par[0];

   // K component
   result += par[2]*TMath::Gaus(x[0],par[1],sigma);

   // L component
   result += par[4]*TMath::Gaus(x[0],par[3],sigma);

   // M component
   result += par[4]/branchingLM*TMath::Gaus(x[0],par[5],sigma);

   return result;
}                                  



// Gaussian Peak function                                           
Double_t gaussianPeak(Double_t *x, Double_t *par)
{                                                    
   Double_t arg = (x[0] - par[1]) / par[2];                 
   return par[0] * TMath::Exp(-0.5 * arg*arg);          
}
