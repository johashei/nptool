#include "TMath.h"
#include "math.h"
#include <cmath>
#include "stdlib.h"

const int numPeaks = 15; 
array<double,numPeaks> means = { 0.000,
                           0.143,
                           0.279,
			   0.728,
			   0.968,
			   1.410,
			   1.981,
			   2.410,
			   2.910,
			   3.2,
			   3.605,
			   3.87,//3.792, 
			   4.04,//4.1, 
			   4.4,
			   5.24
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
  result = par[0];
  // Add N peaks
  for(int pk=0; pk<numPeaks; pk++){
    result += (par[3+(pk*3)]/(par[1+(pk*3)]*sqrt(2*pi)))
	      * exp(-0.5*pow((xx-par[2+(pk*3)])/par[1+(pk*3)],2));
  }
  return result;
}

vector<vector<double>> FitKnownPeaks_RtrnArry(TH1F* hist){
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
      // Set max 279 counts to 100
      //full->SetParameter((i*3)+3,0.);//1e1);
      //if(i==2){full->SetParLimits((i*3)+3,0.0,1e2);}
      //else{full->SetParLimits((i*3)+3,0.0,1e5);}
    full->SetParameter((i*3)+3,1e1);
    full->SetParLimits((i*3)+3,0.0,1e5);
  }
  //full->SetParameter(0,30.);
  full->SetParLimits(0,0.,40.); /* FOR TOTAL SPECTRUM FITTING */
  //full->SetParLimits(0,0.,10.); /* FOR ANGLE GATED FITTING */
  //full->FixParameter(0,0.);
  //full->FixParameter(9,0.); //??

  // Specific limits
  // Set max 279 counts to 100
  //full->FixParameter(9,0.0);
  //full->SetParLimits(9,0.0,1e2); // Doesnt work???
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
  vector<vector<double>> shell = FitKnownPeaks_RtrnArry(hist);
}
