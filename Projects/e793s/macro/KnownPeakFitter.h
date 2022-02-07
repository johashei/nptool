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
   * 3.605
   * 3.792
   * ...
   * NEW: 
   * 3.2
   * 4.1
   * 4.4
   *
   */

  //Assign peak centroids
  const int numPeaks = 14;
  array<double,14> means = { 0.000,
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
			     3.792, 
			     4.1, 
			     4.4 
                           };

  //Build individual peak fit functions
  TF1 **allPeaks = new TF1*[numPeaks];
  for(int i=0; i<numPeaks; i++) {
    string nameHere = nameBase;
    nameHere +=to_string(i);

    allPeaks[i] = new TF1(nameHere.c_str(), function.c_str(), -1, 5);
    allPeaks[i]->SetLineColor(kBlack);  
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
    "  ([2]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[0],2)) " 
    "+ ([4]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[3])/[0],2)) "
    "+ ([6]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[5])/[0],2)) "
    "+ ([8]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[7])/[0],2)) "
    "+ ([10]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[9])/[0],2)) "
    "+ ([12]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[11])/[0],2)) "
    "+ ([14]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[13])/[0],2)) "
    "+ ([16]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[15])/[0],2)) "
    "+ ([18]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[17])/[0],2)) "
    "+ ([20]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[19])/[0],2)) "
    "+ ([22]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[21])/[0],2)) "
    "+ ([24]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[23])/[0],2)) "
    "+ ([26]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[25])/[0],2)) "
    "+ ([28]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[27])/[0],2)) "
    "+ [29]" , minFit, maxFit);
  full->SetLineColor(kRed);

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
  full->SetParName(21,"Mean 11"); full->SetParName(22,"Area*BinWidth 11");
  full->SetParName(23,"Mean 12"); full->SetParName(24,"Area*BinWidth 12");
  full->SetParName(25,"Mean 13"); full->SetParName(26,"Area*BinWidth 13");
  full->SetParName(27,"Mean 14"); full->SetParName(28,"Area*BinWidth 14");
  full->SetParName(29,"Background");

  //Fix sigma & centroid, only allow area to vary  
  const int numParams = (numPeaks*2)+2;
  full->FixParameter(0, sigma);
  for(int i=0; i<numPeaks; i++) {
    full->FixParameter(2*i+1,means.at(i));
    full->SetParameter(2*i+2,1.0);
    full->SetParLimits(2*i+2,0.0,1e5);
  }
  //full->FixParameter(numParams-1,0.0);
  full->SetParameter(numParams-1,1.0);
  full->SetParLimits(numParams-1,0.0,1e1);

  //Fit full function tohistogram
  hist->Fit(full, "WWR", "", minFit, maxFit);
  hist->Draw();
 
  //Extract fitted variables, assign them to individual fits, and draw them
  Double_t finalPar[numParams];
  Double_t finalErr[numParams];
  full->GetParameters(&finalPar[0]);
  for (int i=0; i<numPeaks; i++){
    finalErr[2*i+1] = full->GetParError(2*i+1);
    finalErr[2*i+2] = full->GetParError(2*i+2);
      
    allPeaks[i]->SetParameters(sigma, means.at(i), finalPar[2*i+2]);
  }
  bg->SetParameter(0,finalPar[numParams-1]);
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

Double_t f_full(Double_t *x, Double_t *par){
  Float_t xx = x[0];
  Double_t f;
  Double_t peaks = (par[2]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[1])/par[0],2))  
    + (par[4]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[3])/par[0],2)) 
    + (par[6]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[5])/par[0],2)) 
    + (par[8]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[7])/par[0],2)) 
    + (par[10]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[9])/par[0],2)) 
    + (par[12]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[11])/par[0],2))
    + (par[14]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[13])/par[0],2)) 
    + (par[16]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[15])/par[0],2)) 
    + (par[18]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[17])/par[0],2)) 
    + (par[20]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[19])/par[0],2)) 
    + (par[22]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[21])/par[0],2)) 
    + (par[24]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[23])/par[0],2)) 
    + (par[26]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[25])/par[0],2)) 
    + (par[28]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[27])/par[0],2));

  Double_t bg;
  Double_t a = TMath::Power(6.183,2);
  Double_t b = TMath::Power(xx-10.829,2);
  if(a > b){ bg = par[29] + (par[30]*TMath::Sqrt(a-b)); }
  else{ bg = par[29]; }

  f = peaks + bg;
  return f;
}



vector<vector<double>> FitKnownPeaks_RtrnArry(TH1F* hist){
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
   * 3.605
   * 3.792
   * ...
   * NEW: 
   * 3.2
   * 4.1
   * 4.4
   *
   */
  //Assign peak centroids
  const int numPeaks = 14;
  array<double,14> means = { 0.000,
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
			     3.792, 
			     4.1, 
			     4.4 
                           };

  //Build individual peak fit functions
  TF1 **allPeaks = new TF1*[numPeaks];
  for(int i=0; i<numPeaks; i++) {
    string nameHere = nameBase;
    nameHere +=to_string(i);

    allPeaks[i] = new TF1(nameHere.c_str(), function.c_str(), -1, 5);
    allPeaks[i]->SetLineColor(kBlack);  
    allPeaks[i]->SetLineStyle(7);  
    allPeaks[i]->SetParNames("Sigma", "Mean", "Area*BinWidth");
  } 

  //Build background function
  TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
  bg->SetLineColor(kGreen);
  bg->SetLineStyle(9);
  bg->SetParNames("Background");

  //Build IMPROVED background function
  //TF1 *bg = new TF1 ("bg",f_bg, minFit, maxFit);
  //bg->SetLineColor(kGreen);
  //bg->SetLineStyle(9);
  //bg->SetParNames("Background","BreakupScale");

  //Build total function
  TF1 *full = new TF1("fitAllPeaks", 
    "  ([2]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[0],2)) " 
    "+ ([4]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[3])/[0],2)) "
    "+ ([6]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[5])/[0],2)) "
    "+ ([8]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[7])/[0],2)) "
    "+ ([10]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[9])/[0],2)) "
    "+ ([12]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[11])/[0],2)) "
    "+ ([14]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[13])/[0],2)) "
    "+ ([16]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[15])/[0],2)) "
    "+ ([18]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[17])/[0],2)) "
    "+ ([20]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[19])/[0],2)) "
    "+ ([22]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[21])/[0],2)) "
    "+ ([24]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[23])/[0],2)) "
    "+ ([26]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[25])/[0],2)) "
    "+ ([28]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[27])/[0],2)) "
    "+ [29]" , minFit, maxFit);
  full->SetLineColor(kRed);

  //Build IMPROVED total function
  //TF1 *full = new TF1("fitAllPeaks", f_full, minFit, maxFit);
  //full->SetLineColor(kRed);


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
  full->SetParName(21,"Mean 11"); full->SetParName(22,"Area*BinWidth 11");
  full->SetParName(23,"Mean 12"); full->SetParName(24,"Area*BinWidth 12");
  full->SetParName(25,"Mean 13"); full->SetParName(26,"Area*BinWidth 13");
  full->SetParName(27,"Mean 14"); full->SetParName(28,"Area*BinWidth 14");
  full->SetParName(29,"Background");
  full->SetParName(30,"BreakupScale");

  //Fix sigma & centroid, only allow area to vary  
  const int numParams = (numPeaks*2)+2;
  full->FixParameter(0, sigma);
  for(int i=0; i<numPeaks; i++) {
    full->FixParameter(2*i+1,means.at(i));
    full->SetParameter(2*i+2,1.0);
    full->SetParLimits(2*i+2,0.0,1e5);
  }
  //full->FixParameter(numParams-1,0.0);
  full->SetParameter(numParams-1,1.0);
  full->SetParLimits(numParams-1,0.0,1e1);

  //Fit full function to histogram
  hist->Fit(full, "WWRQ", "", minFit, maxFit);
  hist->Draw();
 
  //Extract fitted variables, assign them to individual fits, and draw them
  Double_t finalPar[numParams];
  Double_t finalErr[numParams];
  full->GetParameters(&finalPar[0]);
  for (int i=0; i<numPeaks; i++){
    finalErr[2*i+1] = full->GetParError(2*i+1);
    finalErr[2*i+2] = full->GetParError(2*i+2);
      
    allPeaks[i]->SetParameters(sigma, means.at(i), finalPar[2*i+2]);
  }
  bg->SetParameter(0,finalPar[numParams-1]);
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
    cout << fixed << setprecision(3) << finalPar[2*i+1] 
	 << "\t" << setprecision(0)<< finalPar[2*i+2]/binWidth 
	 << "\t+- " << (finalPar[2*i+2]/binWidth)*(finalErr[2*i+2]/finalPar[2*i+2]) 
	 << setprecision(3) << endl;

    vector<double> onepeak; //energy, area and error for one peak
    onepeak.push_back(finalPar[2*i+1]);
    onepeak.push_back(finalPar[2*i+2]/binWidth);
    onepeak.push_back((finalPar[2*i+2]/binWidth)*(finalErr[2*i+2]/finalPar[2*i+2]));
    allpeaks.push_back(onepeak);
  }
  return allpeaks;
}

 
