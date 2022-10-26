////////////////////////////////////////////////////////
/* FUNCTIONS TO FIT ONE, TWO AND THREE GAUSSIAN PEAKS */
/*      SingleGaus(), DoubleGaus(), TripleGaus()      */
////////////////////////////////////////////////////////

#include "TMath.h"
#include "math.h"
#include <cmath>
#include "stdlib.h"

Double_t pi = 3.14159265358979323846;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/*	
void DoubleGaus(TH1F* hist){
  bool repeat=true, bgbool = true;
  int repeatInt;
  double minFit, maxFit, mean1, mean2; 

  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak 1 = ";
      cin >> mean1;
    cout << " Peak 2 = ";
      cin >> mean2;
    cout << " Background, yes or no?" << endl;
      cin >> bgbool;

  DoubleGausNumbs(hist, minFit, maxFit, mean1, mean2, bgbool);

}
*/

vector<double> DoubleGausNumbs(TH1F* hist, double minFit, double maxFit, double mean1, double mean2, bool bgbool){
  bool repeat=true;
  int repeatInt;
  vector<double> areasOut;

//  while (repeat){ /* Comment out for mass runs */
    double binWidth = hist->GetXaxis()->GetBinWidth(3);

    //TF1 *g1 = new TF1 ("m1", equation, minFit, maxFit);
    TF1 *g1 = new TF1 ("m1", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    //TF1 *g2 = new TF1 ("m1", equation, minFit, maxFit);
    TF1 *g2 = new TF1 ("m2", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);

    g1->SetParNames("Area*BinWidth", "Mean", "Sigma");
    g2->SetParNames("Area*BinWidth", "Mean", "Sigma");
    bg->SetParNames("Background");

    TF1 *f1 = new TF1("double_gaus", 
		      "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2)) + ([3]/([5]*sqrt(2*pi)))*exp(-0.5*pow((x-[4])/[5],2)) + [6]",
		      minFit, maxFit);

    f1->SetParNames("Area*BinWidth 1", "Mean 1", "Sigma 1",
                    "Area*BinWidth 2", "Mean 2", "Sigma 2", 
		    "Background");
    f1->SetLineColor(kBlack);

    /* OPTION */ //bg->FixParameter(0,0.0); f1->FixParameter(6,0.0);

    /* OPTION */ //g1->FixParameter(0,0.0); f1->FixParameter(0,0.0);
    /* OPTION */ g1->FixParameter(2,0.6); f1->FixParameter(2,0.6); 
    g1->SetParameter(0, 100);
      g1->SetParLimits(0, 0.0, 500.0);
    g1->SetParameter(1, mean1);
      g1->SetParLimits(1, mean1-0.5, mean1+0.5);
    //g1->SetParameter(2, 0.13);
    //g1->SetParameter(2, 0.6);//FOR ELASTICS
    //  //g1->SetParLimits(2, 0.05, 0.20);
    //  g1->SetParLimits(2, 0.4, 1.0);//FOR ELASTICS
    g2->SetParameter(0, 100);
      g2->SetParLimits(0, 0.0, 500.0);
    g2->SetParameter(1, mean2);
      g2->SetParLimits(1, mean2-1.0, mean2+1.0);
    //g2->SetParameter(2, 0.13);
    g2->SetParameter(2, 0.6);//FOR ELASTICS
      //g2->SetParLimits(2, 0.05, 0.20);
      g2->SetParLimits(2, 0.4, 1.0);//FOR ELASTICS

    bg->SetParameter(0, 5.0);//FOR ELASTICS
      bg->SetParLimits(0, 0.0, 20.0);//FOR ELASTICS

    hist->Fit(g1, "WWR", "", minFit, mean1+5);//maxFit);
    hist->Fit(g2, "WWR", "", mean2-5, maxFit);//minFit, maxFit);

    Double_t par[7];
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    bg->GetParameters(&par[6]);
    f1->SetParameters(par);

    /* JUST FOR FITTING 0.36-GATED 5.3MeV PEAK! */
    //  f1->SetParLimits(2, 0.137, 0.40);
    //  f1->SetParLimits(5, 0.137, 0.40);

    if(bgbool==false){bg->FixParameter(0,0.); f1->FixParameter(6,0.);}

    hist->Fit(f1, "WWR", "", minFit, maxFit);
    hist->Draw();
 
    Double_t finalPar[7];
    Double_t finalErr[7];
    f1->GetParameters(&finalPar[0]);
    for (int i=0; i<7; i++){finalErr[i] = f1->GetParError(i);}
    g1->SetParameters(finalPar[0], finalPar[1], finalPar[2]);
    g2->SetParameters(finalPar[3], finalPar[4], finalPar[5]);
    bg->SetParameter(0,finalPar[6]);

    g1->Draw("SAME");
    g2->Draw("SAME");
    bg->Draw("SAME");
    f1->Draw("SAME");


/* Error propogation:
 * (Abin) +- deltaAbin, B+-0 (no uncertainty)
 * A = Abin/B
 * deltaA/A = deltaAbin/Abin
 * deltaA = A x deltaAbin/Abin
 */
/*
    cout << "Area Red : " << finalPar[0]/binWidth 
	    << "  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
	 << "\nArea Blue : " << finalPar[3]/binWidth
	    << "  +-  " << (finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])
	 << endl;
*/

    cout << fixed << setprecision(5);

    areasOut.push_back(finalPar[0]/binWidth);
    areasOut.push_back((finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])); 
    areasOut.push_back(finalPar[3]/binWidth);
    areasOut.push_back((finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])); 

    cout << RED;
    cout << " Mean: \t" << finalPar[1]
	    << "\t +- " << finalErr[1]
	    << endl;
    cout << " Sigm: \t" << finalPar[2]
	    << "\t +- " << finalErr[2]
	    << endl;
    cout << " Area: \t" << finalPar[0]/binWidth 
	    << "\t  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
            << endl;

    cout << BLUE;
    cout << " Mean: \t" << finalPar[4]
	    << "\t +- " << finalErr[4]
	    << endl;
    cout << " Sigm: \t" << finalPar[5]
	    << "\t +- " << finalErr[5]
	    << endl;
    cout << " Area: \t" << finalPar[3]/binWidth 
	    << "\t  +-  " << (finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])
            << endl;
    cout << RESET;


    TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(f1,"Total fit","l");
    legend->AddEntry(g1,"Peak 1","l");
    legend->AddEntry(g2,"Peak 2","l");
    legend->AddEntry(bg,"Background","l");
    legend->Draw();

    //cGate->Draw("SAME");
    gPad->Modified();
    gPad->Update();

//    cout << "\33[37m Refit? " << endl;
//    cin >> repeatInt;
//    if(repeatInt!=1){ repeat=false; }
//  } /* Comment out for mass runs */

  return areasOut; 

}


	
void DoubleGaus(TH1F* hist){
  bool repeat=true;
  int repeatInt;
  bool bgbool = true;
  double minFit, maxFit, mean1, mean2; 

  TCanvas* canvGausFit = new TCanvas("canvGausFit","canvGausFit",1000,1000);
  hist->Draw();

  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak 1 = ";
      cin >> mean1;
    cout << " Peak 2 = ";
      cin >> mean2;
    cout << " Background, yes or no?" << endl;
      cin >> bgbool;

    vector<double> areas = DoubleGausNumbs(hist, minFit, maxFit, mean1, mean2, bgbool);

    cout << "\33[37m Refit? " << endl;
    cin >> repeatInt;
    if(repeatInt!=1){ repeat=false; }
  }
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/*
void SingleGaus(TH1F* hist){
  bool repeat=true;
  int repeatInt;
  double minFit, maxFit, mean; 

  double binWidth = hist->GetXaxis()->GetBinWidth(3);

  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak = ";
      cin >> mean;

    TF1 *g1 = new TF1 ("m1", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    
    TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);

    g1->SetParNames("Area*BinWidth", "Mean", "Sigma");
    bg->SetParNames("Background");

    TF1 *f1 = new TF1("single_gaus", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))+[3]", 
		      minFit, maxFit);
    f1->SetParNames("Area*BinWidth 1", "Mean 1", "Sigma 1",
		    "Background");
    f1->SetLineColor(kBlack);


    g1->SetParameter(0, 100);
      g1->SetParLimits(0, 0.0, 500.0);
    g1->SetParameter(1, mean);
      g1->SetParLimits(1, mean-1.0, mean+1.0);
    g1->SetParameter(2, 0.15);
      g1->SetParLimits(2, 0.05, 1.00);
//    bg->FixParameter(0, 0);
    bg->SetParameter(0, 0);
      bg->SetParLimits(0, 0., 50.);

    hist->Fit(g1, "WWR", "", minFit, maxFit);//maxFit);

    Double_t par[4];
    g1->GetParameters(&par[0]);
    bg->GetParameters(&par[3]);
    f1->SetParameters(par);

    hist->Fit(f1, "WWR", "", minFit, maxFit);
    hist->Draw();
 
    Double_t finalPar[4];
    Double_t finalErr[4];
    f1->GetParameters(&finalPar[0]);
    for (int i=0; i<4; i++){finalErr[i] = f1->GetParError(i);}
    g1->SetParameters(finalPar[0], finalPar[1], finalPar[2]);
    bg->SetParameter(0,finalPar[3]);

    g1->Draw("SAME");
    bg->Draw("SAME");
    f1->Draw("SAME");

    cout << fixed << setprecision(5);

    cout << "\033[91m Mean: \t" << finalPar[1]
	    << "\t +- " << finalErr[1]
	    << endl;
    cout << "\033[91m Sigm: \t" << finalPar[2]
	    << "\t +- " << finalErr[2]
	    << endl;
    cout << "\033[91m Area: \t" << finalPar[0]/binWidth 
	    << "\t  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
            << endl;

    TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(f1,"Total fit","l");
    legend->AddEntry(g1,"Peak","l");
    legend->AddEntry(bg,"Background","l");
    legend->Draw();

    //cGate->Draw("SAME");
    gPad->Modified();
    gPad->Update();

    cout << "\033[37m Refit? " << endl;
    cin >> repeatInt;
    if(repeatInt!=1){ repeat=false; }
  }

}
*/
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void TripleGaus(TH1F* hist){
  bool repeat=true;
  int repeatInt;
  double minFit, maxFit, mean1, mean2, mean3; 

  double binWidth = hist->GetXaxis()->GetBinWidth(3);

  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak 1 = ";
      cin >> mean1;
    cout << " Peak 2 = ";
      cin >> mean2;
    cout << " Peak 3 = ";
      cin >> mean3;

    TF1 *g1 = new TF1 ("m1", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    TF1 *g2 = new TF1 ("m2", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    TF1 *g3 = new TF1 ("m3", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g3->SetLineColor(kViolet);
    g3->SetLineStyle(2);

    TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);

    g1->SetParNames("Area*BinWidth", "Mean", "Sigma");
    g2->SetParNames("Area*BinWidth", "Mean", "Sigma");
    g3->SetParNames("Area*BinWidth", "Mean", "Sigma");
    bg->SetParNames("Background");

    TF1 *f1 = new TF1("double_gaus", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2)) + ([3]/([5]*sqrt(2*pi)))*exp(-0.5*pow((x-[4])/[5],2)) + ([6]/([8]*sqrt(2*pi)))*exp(-0.5*pow((x-[7])/[8],2)) + [9]" , minFit, maxFit);

    f1->SetParNames("Area*BinWidth 1", "Mean 1", "Sigma 1",
                    "Area*BinWidth 2", "Mean 2", "Sigma 2", 
                    "Area*BinWidth 3", "Mean 3", "Sigma 3", 
		    "Background");
    f1->SetLineColor(kBlack);

/*
    TF1 *g1 = new TF1 ("m1", "gaus", minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    TF1 *g2 = new TF1 ("m2", "gaus", minFit, maxFit);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    TF1 *bg = new TF1 ("bg", "[0]", minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);
    TF1 *f1 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6]", minFit, maxFit);
    f1->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                    "Constant 2", "Mean 2", "Sigma 2", 
		    "Background");
    f1->SetLineColor(kBlack);
*/

    g1->SetParameter(0, 100);
      g1->SetParLimits(0, 0.0, 1000.0);
    g1->SetParameter(1, mean1);
      g1->SetParLimits(1, mean1-1.0, mean1+1.0);
    g1->SetParameter(2, 0.13);
      g1->SetParLimits(2, 0.05, 0.30);
    g2->SetParameter(0, 100);
      g2->SetParLimits(0, 0.0, 1000.0);
    g2->SetParameter(1, mean2);
      g2->SetParLimits(1, mean2-1.0, mean2+1.0);
    g2->SetParameter(2, 0.13);
      g2->SetParLimits(2, 0.05, 0.30);
    g3->SetParameter(0, 100);
      g3->SetParLimits(0, 0.0, 1000.0);
    g3->SetParameter(1, mean3);
      g3->SetParLimits(1, mean3-1.0, mean3+1.0);
    g3->SetParameter(2, 0.13);
      g3->SetParLimits(2, 0.05, 0.30);


    hist->Fit(g1, "WWR", "", minFit, mean1+5);//maxFit);
    hist->Fit(g2, "WWR", "", mean2-5, maxFit);//minFit, maxFit);

    Double_t par[10];
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    g3->GetParameters(&par[6]);
    bg->GetParameters(&par[9]);
    f1->SetParameters(par);

    hist->Fit(f1, "WWR", "", minFit, maxFit);
    //hist->GetXaxis()->SetTitle("ThetaLab (degrees)");
    //hist->GetYaxis()->SetTitle("Counts (per 0.5 degrees)");
    //hist->SetTitle(titleString.c_str());
    hist->Draw();
 
    Double_t finalPar[10];
    Double_t finalErr[10];
    f1->GetParameters(&finalPar[0]);
    for (int i=0; i<10; i++){finalErr[i] = f1->GetParError(i);}
    g1->SetParameters(finalPar[0], finalPar[1], finalPar[2]);
    g2->SetParameters(finalPar[3], finalPar[4], finalPar[5]);
    g3->SetParameters(finalPar[6], finalPar[7], finalPar[8]);
    bg->SetParameter(0,finalPar[9]);

    g1->Draw("SAME");
    g2->Draw("SAME");
    g3->Draw("SAME");
    bg->Draw("SAME");
    f1->Draw("SAME");

//    cout << "Area Red : " << finalPar[0]/binWidth 
//	    << "  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
//	 << "\nArea Blue : " << finalPar[3]/binWidth 
//	    << "  +-  " << (finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])
//	 << "\nArea Violet : " << finalPar[6]/binWidth 
//	    << "  +-  " << (finalPar[6]/binWidth) * (finalErr[6]/finalPar[6])
//         << endl;

    cout << fixed << setprecision(5);

    cout << "\033[91m Mean: \t" << finalPar[1]
	    << "\t +- " << finalErr[1]
	    << endl;
    cout << "\033[91m Sigm: \t" << finalPar[2]
	    << "\t +- " << finalErr[2]
	    << endl;
    cout << "\033[91m Area: \t" << finalPar[0]/binWidth 
	    << "\t  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
            << endl;

    cout << "\033[36m Mean: \t" << finalPar[4]
	    << "\t +- " << finalErr[4]
	    << endl;
    cout << "\033[36m Sigm: \t" << finalPar[5]
	    << "\t +- " << finalErr[5]
	    << endl;
    cout << "\033[36m Area: \t" << finalPar[3]/binWidth 
	    << "\t  +-  " << (finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])
            << endl;

    cout << "\033[95m Mean: \t" << finalPar[7]
	    << "\t +- " << finalErr[7]
	    << endl;
    cout << "\033[95m Sigm: \t" << finalPar[8]
	    << "\t +- " << finalErr[8]
	    << endl;
    cout << "\033[95m Area: \t" << finalPar[6]/binWidth 
	    << "\t  +-  " << (finalPar[6]/binWidth) * (finalErr[6]/finalPar[6])
            << endl;


    TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(f1,"Total fit","l");
    legend->AddEntry(g1,"Peak 1","l");
    legend->AddEntry(g2,"Peak 2","l");
    legend->AddEntry(g3,"Peak 3","l");
    legend->AddEntry(bg,"Background","l");
    legend->Draw();

    //cGate->Draw("SAME");
    gPad->Modified();
    gPad->Update();

    cout << "\33[37m Refit? " << endl;
    cin >> repeatInt;
    if(repeatInt!=1){ repeat=false; }
  }

}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void SingleGausNoBG(TH1F* hist){
  bool repeat=true;
  int repeatInt;
  double minFit, maxFit, mean; 

  double binWidth = hist->GetXaxis()->GetBinWidth(3);

  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak = ";
      cin >> mean;

    TF1 *g1 = new TF1 ("m1", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    
    g1->SetParNames("Area*BinWidth", "Mean", "Sigma");

    TF1 *f1 = new TF1("single_gaus", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))", 
		      minFit, maxFit);
    f1->SetParNames("Area*BinWidth 1", "Mean 1", "Sigma 1");
    f1->SetLineColor(kBlack);

/*
    TF1 *g1 = new TF1 ("m1", "gaus", minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    TF1 *g2 = new TF1 ("m2", "gaus", minFit, maxFit);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    TF1 *bg = new TF1 ("bg", "[0]", minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);
    TF1 *f1 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6]", minFit, maxFit);
    f1->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                    "Constant 2", "Mean 2", "Sigma 2", 
		    "Background");
    f1->SetLineColor(kBlack);
*/

    g1->SetParameter(0, 100);
      g1->SetParLimits(0, 0.0, 500.0);
    g1->SetParameter(1, mean);
      g1->SetParLimits(1, mean-1.0, mean+1.0);
    g1->SetParameter(2, 0.15);
      g1->SetParLimits(2, 0.05, 1.00);

    hist->Fit(g1, "WWR", "", minFit, maxFit);//maxFit);

    Double_t par[4];
    g1->GetParameters(&par[0]);
    f1->SetParameters(par);

    hist->Fit(f1, "WWR", "", minFit, maxFit);
    hist->Draw();
 
    Double_t finalPar[4];
    Double_t finalErr[4];
    f1->GetParameters(&finalPar[0]);
    for (int i=0; i<4; i++){finalErr[i] = f1->GetParError(i);}
    g1->SetParameters(finalPar[0], finalPar[1], finalPar[2]);

    g1->Draw("SAME");
    f1->Draw("SAME");

    cout << fixed << setprecision(5);

    cout << "\033[91m Mean: \t" << finalPar[1]
	    << "\t +- " << finalErr[1]
	    << endl;
    cout << "\033[91m Sigm: \t" << finalPar[2]
	    << "\t +- " << finalErr[2]
	    << endl;
    cout << "\033[91m Area: \t" << finalPar[0]/binWidth 
	    << "\t  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
            << endl;

    TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(f1,"Total fit","l");
    legend->AddEntry(g1,"Peak","l");
    legend->Draw();

    //cGate->Draw("SAME");
    gPad->Modified();
    gPad->Update();

    cout << "\033[37m Refit? " << endl;
    cin >> repeatInt;
    if(repeatInt!=1){ repeat=false; }
  }

}



////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void SingleGaus(TH1F* hist, bool isGamma){
  bool repeat=true;
  int repeatInt;
  double minFit, maxFit, mean; 

  double binWidth = hist->GetXaxis()->GetBinWidth(3);

  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak = ";
      cin >> mean;

    TF1 *g1 = new TF1 ("m1", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    
    TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);

    g1->SetParNames("Area*BinWidth", "Mean", "Sigma");
    bg->SetParNames("Background");

    TF1 *f1 = new TF1("single_gaus", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))+[3]", 
		      minFit, maxFit);
    f1->SetParNames("Area*BinWidth 1", "Mean 1", "Sigma 1",
		    "Background");
    f1->SetLineColor(kBlack);

    double areaSet, areaMax, meanRange;
    double sigSet, sigMin, sigMax;
    double bgMax;

    if(isGamma){
      areaSet = 1000.; areaMax = 10000.; meanRange = 0.005;
      sigSet = 0.001; sigMin = 0.0005; sigMax = 0.015;
      bgMax = 1000.;
    } else {
      areaSet = 100.; areaMax = 1000.; meanRange = 1.0;
      sigSet = 0.15; sigMin = 0.05; sigMax = 0.3;
      bgMax = 50.;
    }

    g1->SetParameter(0, areaSet);
      g1->SetParLimits(0, 0.0, areaMax);
    g1->SetParameter(1, mean);
      g1->SetParLimits(1, mean-meanRange, mean+meanRange);
    g1->SetParameter(2, sigSet);
      g1->SetParLimits(2, sigMin, sigMax);
//    bg->FixParameter(0, 0);
    bg->SetParameter(0, 0);
      bg->SetParLimits(0, 0., bgMax);

    hist->Fit(g1, "WWR", "", minFit, maxFit);//maxFit);

    Double_t par[4];
    g1->GetParameters(&par[0]);
    bg->GetParameters(&par[3]);
    f1->SetParameters(par);

    hist->Fit(f1, "WWR", "", minFit, maxFit);
    hist->Draw();
 
    Double_t finalPar[4];
    Double_t finalErr[4];
    f1->GetParameters(&finalPar[0]);
    for (int i=0; i<4; i++){finalErr[i] = f1->GetParError(i);}
    g1->SetParameters(finalPar[0], finalPar[1], finalPar[2]);
    bg->SetParameter(0,finalPar[3]);

    g1->Draw("SAME");
    bg->Draw("SAME");
    f1->Draw("SAME");

    cout << fixed << setprecision(5);

    cout << "\033[91m Mean: \t" << finalPar[1]
	    << "\t +- " << finalErr[1]
	    << endl;
    cout << "\033[91m Sigm: \t" << finalPar[2]
	    << "\t +- " << finalErr[2]
	    << endl;
    cout << "\033[91m Area: \t" << finalPar[0]/binWidth 
	    << "\t  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
            << endl;

    TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(f1,"Total fit","l");
    legend->AddEntry(g1,"Peak","l");
    legend->AddEntry(bg,"Background","l");
    legend->Draw();

    //cGate->Draw("SAME");
    gPad->Modified();
    gPad->Update();

    cout << "\033[37m Refit? " << endl;
    cin >> repeatInt;
    if(repeatInt!=1){ repeat=false; }
  }

}



void DoubleGausForElasticsFitting(TH1F* hist){
  bool repeat=true, bgbool = true;
  int repeatInt;
  double minFit, maxFit, mean1, mean2; 

  double binWidth = hist->GetXaxis()->GetBinWidth(3);
  cout << "Bin Width: " << binWidth << endl;


  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak 1 = ";
      cin >> mean1;
    cout << " Peak 2 = ";
      cin >> mean2;
    cout << " Background, yes or no?" << endl;
      cin >> bgbool;

    //TF1 *g1 = new TF1 ("m1", equation, minFit, maxFit);
    TF1 *g1 = new TF1 ("m1", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    //TF1 *g2 = new TF1 ("m1", equation, minFit, maxFit);
    TF1 *g2 = new TF1 ("m2", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);

    g1->SetParNames("Area*BinWidth", "Mean", "Sigma");
    g2->SetParNames("Area*BinWidth", "Mean", "Sigma");
    bg->SetParNames("Background");

    TF1 *f1 = new TF1("double_gaus", 
		      "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2)) + ([3]/([5]*sqrt(2*pi)))*exp(-0.5*pow((x-[4])/[5],2)) + [6]",
		      minFit, maxFit);

    f1->SetParNames("Area*BinWidth 1", "Mean 1", "Sigma 1",
                    "Area*BinWidth 2", "Mean 2", "Sigma 2", 
		    "Background");
    f1->SetLineColor(kBlack);

    g1->SetParameter(0, 100);
      g1->SetParLimits(0, 0.0, 500.0);
    g1->SetParameter(1, mean1);
      g1->SetParLimits(1, mean1-0.5, mean1+0.5);
    g1->SetParameter(2, 0.30);
    //g1->SetParameter(2, 0.14);
      g1->SetParLimits(2, 0.10, 0.80);
      //g1->SetParLimits(2, 0.13, 0.30);
    g2->SetParameter(0, 100);
      g2->SetParLimits(0, 0.0, 500.0);
    g2->SetParameter(1, mean2);
      g2->SetParLimits(1, mean2-1.0, mean2+1.0);
    g2->SetParameter(2, 0.30);
    //g2->SetParameter(2, 0.14);
      g2->SetParLimits(2, 0.10, 0.80);
      //g2->SetParLimits(2, 0.13, 0.30);
    bg->SetParameter(0,25.);
      bg->SetParLimits(0,1.,75.);

    hist->Fit(g1, "WWR", "", minFit, mean1+5);//maxFit);
    hist->Fit(g2, "WWR", "", mean2-5, maxFit);//minFit, maxFit);

    Double_t par[7];
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    bg->GetParameters(&par[6]);
    f1->SetParameters(par);

    /* JUST FOR FITTING 0.36-GATED 5.3MeV PEAK! */
    //  f1->SetParLimits(2, 0.137, 0.40);
    //  f1->SetParLimits(5, 0.137, 0.40);

    if(bgbool==false){bg->FixParameter(0,0.); f1->FixParameter(6,0.);}

    hist->Fit(f1, "WWR", "", minFit, maxFit);
    hist->Draw();
 
    Double_t finalPar[7];
    Double_t finalErr[7];
    f1->GetParameters(&finalPar[0]);
    for (int i=0; i<7; i++){finalErr[i] = f1->GetParError(i);}
    g1->SetParameters(finalPar[0], finalPar[1], finalPar[2]);
    g2->SetParameters(finalPar[3], finalPar[4], finalPar[5]);
    bg->SetParameter(0,finalPar[6]);

    g1->Draw("SAME");
    g2->Draw("SAME");
    bg->Draw("SAME");
    f1->Draw("SAME");


/* Error propogation:
 * (Abin) +- deltaAbin, B+-0 (no uncertainty)
 * A = Abin/B
 * deltaA/A = deltaAbin/Abin
 * deltaA = A x deltaAbin/Abin
 */
/*
    cout << "Area Red : " << finalPar[0]/binWidth 
	    << "  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
	 << "\nArea Blue : " << finalPar[3]/binWidth
	    << "  +-  " << (finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])
	 << endl;
*/

    cout << fixed << setprecision(5);

    cout << "\033[91m Mean: \t" << finalPar[1]
	    << "\t +- " << finalErr[1]
	    << endl;
    cout << "\033[91m Sigm: \t" << finalPar[2]
	    << "\t +- " << finalErr[2]
	    << endl;
    cout << "\033[91m Area: \t" << finalPar[0]/binWidth 
	    << "\t  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
            << endl;

    cout << "\033[36m Mean: \t" << finalPar[4]
	    << "\t +- " << finalErr[4]
	    << endl;
    cout << "\033[36m Sigm: \t" << finalPar[5]
	    << "\t +- " << finalErr[5]
	    << endl;
    cout << "\033[36m Area: \t" << finalPar[3]/binWidth 
	    << "\t  +-  " << (finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])
            << endl;


    TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(f1,"Total fit","l");
    legend->AddEntry(g1,"Peak 1","l");
    legend->AddEntry(g2,"Peak 2","l");
    legend->AddEntry(bg,"Background","l");
    legend->Draw();

    //cGate->Draw("SAME");
    gPad->Modified();
    gPad->Update();

    cout << "\33[37m Refit? " << endl;
    cin >> repeatInt;
    if(repeatInt!=1){ repeat=false; }
  }

}





void SingleGausForElasticsFitting(TH1F* hist){
  bool repeat=true, bgbool = true;
  int repeatInt;
  double minFit, maxFit, mean1; 

  double binWidth = hist->GetXaxis()->GetBinWidth(3);
  cout << "Bin Width: " << binWidth << endl;


  while (repeat){
    cout << "====================================================================" << endl;
    cout << " Input range to fit:" << endl;
    cout << " Min = ";
      cin >> minFit;
    cout << " Max = ";
      cin >> maxFit;
    cout << " Peak 1 = ";
      cin >> mean1;
    cout << " Background, yes or no?" << endl;
      cin >> bgbool;

    //TF1 *g1 = new TF1 ("m1", equation, minFit, maxFit);
    TF1 *g1 = new TF1 ("m1", "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2))",
		       minFit, maxFit);
    g1->SetLineColor(kRed);
    g1->SetLineStyle(2);
    //TF1 *g2 = new TF1 ("m1", equation, minFit, maxFit);
    TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
    bg->SetLineColor(kGreen);
    bg->SetLineStyle(9);

    g1->SetParNames("Area*BinWidth", "Mean", "Sigma");
    bg->SetParNames("Background");

    TF1 *f1 = new TF1("double_gaus", 
		      "([0]/([2]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[2],2)) + [3]",
		      minFit, maxFit);

    f1->SetParNames("Area*BinWidth 1", "Mean 1", "Sigma 1",
		    "Background");
    f1->SetLineColor(kBlack);

    g1->SetParameter(0, 100);
      g1->SetParLimits(0, 0.0, 500.0);
    g1->SetParameter(1, mean1);
      g1->SetParLimits(1, mean1-0.5, mean1+0.5);
    g1->SetParameter(2, 0.13);
    //g1->SetParameter(2, 0.14);
      g1->SetParLimits(2, 0.05, 0.30);
      //g1->SetParLimits(2, 0.13, 0.30);
    bg->SetParameter(0,25.);
      bg->SetParLimits(0,1.,50.);

    hist->Fit(g1, "WWR", "", minFit, mean1+5);//maxFit);

    Double_t par[7];
    g1->GetParameters(&par[0]);
    bg->GetParameters(&par[3]);
    f1->SetParameters(par);

    /* JUST FOR FITTING 0.36-GATED 5.3MeV PEAK! */
    //  f1->SetParLimits(2, 0.137, 0.40);
    //  f1->SetParLimits(5, 0.137, 0.40);

    if(bgbool==false){bg->FixParameter(0,0.); f1->FixParameter(3,0.);}

    hist->Fit(f1, "WWR", "", minFit, maxFit);
    hist->Draw();
 
    Double_t finalPar[4];
    Double_t finalErr[4];
    f1->GetParameters(&finalPar[0]);
    for (int i=0; i<4; i++){finalErr[i] = f1->GetParError(i);}
    g1->SetParameters(finalPar[0], finalPar[1], finalPar[2]);
    bg->SetParameter(0,finalPar[3]);

    g1->Draw("SAME");
    bg->Draw("SAME");
    f1->Draw("SAME");


/* Error propogation:
 * (Abin) +- deltaAbin, B+-0 (no uncertainty)
 * A = Abin/B
 * deltaA/A = deltaAbin/Abin
 * deltaA = A x deltaAbin/Abin
 */
/*
    cout << "Area Red : " << finalPar[0]/binWidth 
	    << "  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
	 << "\nArea Blue : " << finalPar[3]/binWidth
	    << "  +-  " << (finalPar[3]/binWidth) * (finalErr[3]/finalPar[3])
	 << endl;
*/

    cout << fixed << setprecision(5);

    cout << "\033[91m Mean: \t" << finalPar[1]
	    << "\t +- " << finalErr[1]
	    << endl;
    cout << "\033[91m Sigm: \t" << finalPar[2]
	    << "\t +- " << finalErr[2]
	    << endl;
    cout << "\033[91m Area: \t" << finalPar[0]/binWidth 
	    << "\t  +-  " << (finalPar[0]/binWidth) * (finalErr[0]/finalPar[0])
            << endl;

    TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(f1,"Total fit","l");
    legend->AddEntry(g1,"Peak 1","l");
    legend->AddEntry(bg,"Background","l");
    legend->Draw();

    //cGate->Draw("SAME");
    gPad->Modified();
    gPad->Update();

    cout << "\33[37m Refit? " << endl;
    cin >> repeatInt;
    if(repeatInt!=1){ repeat=false; }
  }

}
