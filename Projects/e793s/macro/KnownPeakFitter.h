#include "TMath.h"
#include "math.h"
#include <cmath>
#include "stdlib.h"

void FitKnownPeaks(TH1F* hist){
  double minFit=-1.0, maxFit=5.0; 
  double binWidth = hist->GetXaxis()->GetBinWidth(3);
  double sigma = 0.14;

  string nameBase = "Peak ";
  string function = "([2]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[0],2))";

  /* 11 KNOWN PEAKS as of 12/10/21
   * 0.000
   * 0.143
   * 0.279
   * 0.728
   * 0.968
   * 1.410
   * 1.981
   * 2.410
   * 2.910
   * 3.600
   * 3.792
   * ...
   */

  //Assign peak centroids
  const int numPeaks = 11;
  array<double,11> means = { 0.000,
                             0.143,
                             0.279,
			     0.728,
			     0.968,
			     1.410,
			     1.981,
			     2.410,
			     2.910,
			     3.600,
			     3.792 
                           };

  //Build individual peak fit functions
  TF1 **allPeaks = new TF1*[numPeaks];
  for(int i=0; i<numPeaks; i++) {
    string nameHere = nameBase;
    nameHere +=to_string(i);

    allPeaks[i] = new TF1(nameHere.c_str(), function.c_str(), -1, 5);
    allPeaks[i]->SetLineColor(kRed);  
    allPeaks[i]->SetLineStyle(7);  
    allPeaks[i]->SetParNames("Sigma", "Mean", "Area*BinWidth");
  } 

  //Build background function
  TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
  bg->SetLineColor(kGreen);
  bg->SetLineStyle(9);
  bg->SetParNames("Background");

  //Build total function
  TF1 *full = new TF1("fitAllPeaks", 
    "  ([02]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[01])/[0],2)) " 
    "+ ([04]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[03])/[0],2)) "
    "+ ([06]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[05])/[0],2)) "
    "+ ([08]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[07])/[0],2)) "
    "+ ([10]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[09])/[0],2)) "
    "+ ([12]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[11])/[0],2)) "
    "+ ([14]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[13])/[0],2)) "
    "+ ([16]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[15])/[0],2)) "
    "+ ([18]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[17])/[0],2)) "
    "+ ([20]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[19])/[0],2)) "
    "+ [21]" , minFit, maxFit);
  full->SetLineColor(kBlack);

  //Annoyingly long parameter name assignment 
  //(SetParNames only works for up to 11 variables)
  full->SetParName(0,"Sigma");
  full->SetParName(1,"Mean 01");  full->SetParName(2,"Area*BinWidth 01");
  full->SetParName(3,"Mean 02");  full->SetParName(4,"Area*BinWidth 02");
  full->SetParName(5,"Mean 03");  full->SetParName(6,"Area*BinWidth 03");
  full->SetParName(7,"Mean 04");  full->SetParName(8,"Area*BinWidth 04");
  full->SetParName(9,"Mean 05");  full->SetParName(10,"Area*BinWidth 05");
  full->SetParName(11,"Mean 06"); full->SetParName(12,"Area*BinWidth 06");
  full->SetParName(13,"Mean 07"); full->SetParName(14,"Area*BinWidth 07");
  full->SetParName(15,"Mean 08"); full->SetParName(16,"Area*BinWidth 08");
  full->SetParName(17,"Mean 09"); full->SetParName(18,"Area*BinWidth 09");
  full->SetParName(19,"Mean 10"); full->SetParName(20,"Area*BinWidth 10");
  full->SetParName(21,"Background");

  //Fix sigma & centroid, only allow area to vary  
  full->FixParameter(0, sigma);
  for(int i=0; i<numPeaks; i++) {
    full->FixParameter(2*i+1,means.at(i));
    full->SetParameter(2*i+2,1.0);
    full->SetParLimits(2*i+2,0.0,1e5);
  }

  //Fit full function to histogram
  hist->Fit(full, "WWR", "", minFit, maxFit);
  hist->Draw();
 
  //Extract fitted variables, assign them to individual fits, and draw them
  Double_t finalPar[22];
  Double_t finalErr[22];
  full->GetParameters(&finalPar[0]);
  for (int i=0; i<numPeaks; i++){
    finalErr[2*i+1] = full->GetParError(2*i+1);
    finalErr[2*i+2] = full->GetParError(2*i+2);
      
    allPeaks[i]->SetParameters(sigma, means.at(i), finalPar[2*i+2]);
  }
  bg->SetParameter(0,finalPar[21]);
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
  
  for(int i=0; i<numPeaks; i++){
    cout << fixed << setprecision(3) << finalPar[2*i+1] 
	 << "\t" << setprecision(0)<< finalPar[2*i+2]/binWidth 
	 << "\t+- " << (finalPar[2*i+2]/binWidth)*(finalErr[2*i+2]/finalPar[2*i+2]) 
	 << endl;
  }

}
