#include "TMath.h"
#include "math.h"
#include <cmath>
#include "stdlib.h"

const int numPeaks = 17; 
array<double,numPeaks> means = { 0.000,
                           0.143,
                           0.279,
			   0.728,
			   0.968,
			   1.410,
			   1.981,//1.952,//1.981,
			   2.412,
			   2.910,
			   3.253,
			   3.605,
			   3.795,//Split in two?
			   3.870,//Split in two? 
			   4.045,//4.1,
			   4.393,
			   //4.51//,
			   5.15,
			   5.82
                           };

array<double,27> knowngammas = { 0.143,
					0.279,
					0.449,
					0.968,
					1.130,
					1.410,
					1.267,
					//0.575,
					1.013,
					1.838,
					1.981,
					1.000,
					2.412,
					2.767,
					2.518,
					3.325,
					2.878,
					3.605,
					2.839,
					2.734,
					3.522,
					3.076,
					3.875,
					0.834,
					3.325,
					3.77,
					4.037,
					4.364					
					};

/*
Double_t f_bg(Double_t *x, Double_t *par){
  // Flat bg [0] + semicircle [1]*sqrt(6.183^2 - (x-10.829)^2) 
  Float_t xx = x[0];
  Double_t f;
  Double_t a = TMath::Power(6.183,2);
  Double_t b = TMath::Power(xx-10.829,2);
  if(a > b){ f = par[0] + (par[1]*TMath::Sqrt(a-b)); }
  else{ f = par[0]; }
  return f;
}

Double_t f_peak(Double_t *x, Double_t *par){
  float xx = x[0];
  double f = (par[2]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[1])/par[0],2));
  return f;
}
*/

Double_t f_full(Double_t *x, Double_t *par) {
  float xx = x[0];
  double result, norm;
  // Flat background
//  result = par[0];
  result = 0;
  // Add N peaks
  for(int pk=0; pk<numPeaks; pk++){
    result += (par[3+(pk*3)]/(par[1+(pk*3)]*sqrt(2*pi)))
	      //* exp(-0.5*pow((xx-par[2+(pk*3)])/par[1+(pk*3)],2));
	      * exp(-0.5*pow((xx-par[2+(pk*3)]-par[0])/par[1+(pk*3)],2)); //added par 0 as shift in energy
  }
  return result;
}

vector<vector<double>> FitKnownPeaks_RtrnArry(TH1F* hist, double slideshift){
  double minFit=-1.0, maxFit=8.0; 
  double binWidth = hist->GetXaxis()->GetBinWidth(3);
  double sigma = 0.14;

  hist->Sumw2();

  /* Construct flat BG to subtract */ 
  /**
  cout << " REMOVING FLAT BG OF 36 COUNTS!!!!" << endl;
  cout << " REMOVING FLAT BG OF 36 COUNTS!!!!" << endl;
  cout << " REMOVING FLAT BG OF 36 COUNTS!!!!" << endl;
  double ConstBG = 36.0; double ErrBG = 1.0;
  int xbins = hist->GetXaxis()->GetNbins();
  double xmin = hist->GetXaxis()->GetXmin();
  double xmax = hist->GetXaxis()->GetXmax();
  TH1F *FlatBG = new TH1F("FlatBG","FlatBG", xbins, xmin, xmax);
  for(int i=0; i<xbins;i++){
    FlatBG->SetBinContent(i,ConstBG);
    FlatBG->SetBinError(i,ErrBG);
  }
  hist->Add(FlatBG,-1);
  **/

  //Build individual peak fit functions
  string nameBase = "Peak ";
  string function = "([2]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[0],2))";
  TF1 **allPeaks = new TF1*[numPeaks];
  for(int i=0; i<numPeaks; i++) {
    string nameHere = nameBase;
    nameHere +=to_string(i);

    allPeaks[i] = new TF1(nameHere.c_str(), function.c_str(), minFit, maxFit);
    //allPeaks[i] = new TF1(nameHere.c_str(), f_peak, -1, 5);
    allPeaks[i]->SetLineColor(kBlack);  
    allPeaks[i]->SetLineStyle(7);  
    allPeaks[i]->SetParNames("Sigma", "Mean", "Area*BinWidth");
  } 

  //Subtract flat background equal to smallest bin in range
  
  hist->GetXaxis()->SetRange(hist->FindBin(-0.9), hist->FindBin(-0.4));
  double bgmin = hist->GetBinContent(hist->GetMinimumBin());
  hist->GetXaxis()->UnZoom();
  cout << "Subtracting background of " << bgmin << endl;
  for(int b=1; b<hist->GetNbinsX() ; b++){
      hist->SetBinContent(b,hist->GetBinContent(b)-bgmin);
  }
 

  //Build background function
  TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
  bg->SetLineColor(kGreen);
  bg->SetLineStyle(9);
  bg->SetParNames("Background");

  //Build IMPROVED total function
  TF1 *full = new TF1("fitAllPeaks", f_full, minFit, maxFit, (int) 1+(3*numPeaks));
  full->SetLineColor(kRed);
  const int numParams = (numPeaks*3)+1;
  for(int i=0; i<numPeaks; i++) {
    full->FixParameter((i*3)+1,sigma);
    full->FixParameter((i*3)+2,means.at(i));
    full->SetParameter((i*3)+3,1e1);
    full->SetParLimits((i*3)+3,0.0,1e5);
    //full->SetParLimits((i*3)+3,10.0,1e5);
  }
  //full->SetParLimits(0,0.,40.); /* FOR TOTAL SPECTRUM FITTING */
  //full->SetParLimits(0,0.,10.); /* FOR ANGLE GATED FITTING */
  //full->FixParameter(0,0.); /* FOR ANGLE GATED FITTING WITH BG SUBTRACTED */
  full->FixParameter(9,0.); /* FIX 0.279 AREA TO ZERO */
  full->SetParLimits(0,-0.5,+0.5); /* FOR WHEN PAR[0] IS VARIABLE ENERGY CENTRIOD SLIDER */
  //full->FixParameter(0,slideshift); /* FOR WHEN PAR[0] IS FIXED ENERGY CENTRIOD SLIDER */
  full->SetParameter(52,0.39); /* SET 5.8MeV SIGMA */
  full->SetParLimits(52,0.34,0.44); /* SET 5.8MeV SIGMA */
  
  allPeaks[14]->SetLineColor(kOrange);

  //Fit full function to histogram
  hist->Fit(full, "RWQB", "", minFit, maxFit);
  hist->Draw();

  //Extract fitted variables, assign them to individual fits, and draw them
  const Double_t* finalPar = full->GetParameters();
  const Double_t* finalErr = full->GetParErrors();
  for (int i=0; i<numPeaks; i++){
    allPeaks[i]->SetParameters(sigma, means.at(i), finalPar[3+(i*3)]);
  }
  bg->SetParameter(0,finalPar[0]);
  bg->Draw("SAME");
  full->Draw("SAME");

  for (int i=0; i<numPeaks; i++){
    allPeaks[i]->Draw("SAME");
  }

 /* Error propogation:
  * (Abin) +- deltaAbin, B+-0 (no uncertainty)
  * A = Abin/B
  * deltaA/A = deltaAbin/Abin
  * deltaA = A x deltaAbin/Abin
  */

  //Write to screen
  cout << "===========================" << endl;
  cout << "== PEAK =========== AREA ==" << endl;
  
  vector<vector<double>> allpeaks;
  for(int i=0; i<numPeaks; i++){
    double A = finalPar[(i*3)+3]/binWidth;
    double deltaA = A *  (finalErr[(i*3)+3]/finalPar[(i*3)+3]);

    cout << fixed << setprecision(3) 
	 << " #" << i << "  " 
	 << finalPar[(i*3)+2] << "\t" << setprecision(0)
	 << A << "\t+- " 
	 << deltaA << setprecision(3)
	 << endl;

    vector<double> onepeak; //energy, area and error for one peak
    onepeak.push_back(finalPar[(i*3)+2]);
    onepeak.push_back(A);
    onepeak.push_back(deltaA);
    allpeaks.push_back(onepeak);
  }
  cout << " BG  " << full->GetParameter(0) 
       << " +- " << full->GetParError(0) << endl;

  return allpeaks;
}

void FitKnownPeaks(TH1F* hist){
  //Shell function to call Rtrn_Arry without writing vector<vector<double>> to screen
  vector<vector<double>> shell = FitKnownPeaks_RtrnArry(hist,0.0);
}
