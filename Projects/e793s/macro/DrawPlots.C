#include "GausFit.h"
#include "KnownPeakFitter.h"
#include "DrawPlots.h"

#include "CS2.h"
#include "ThreeBodyBreakup.h"
/* USE THIS SPACE TO TEST NEW FEATURES */

void thickness(){

  std::ifstream infile("thicknessTheory4.txt");
  
  // THEORY -------------------------------------
  double x, pCH, dHSS, pBG, pPer, dKD, d79D, dPer, dLH, dBel;
  vector<double> vx, vpCH, vdHSS, vpBG, vpPer, vdKD, vd79D, vdPer, vdLH, vdBel;
  int count = 0;

  while (infile >> x >> pCH >> dHSS >> pBG >> pPer >> dKD >> d79D >> dPer >> dLH >> dBel)
  {
    vx.push_back(x);
    vpBG.push_back(pBG);
    vdHSS.push_back(dHSS);
    vpCH.push_back(pCH);
    vpPer.push_back(pPer);
    vdKD.push_back(dKD);
    vd79D.push_back(d79D);
    vdPer.push_back(dPer);
    vdLH.push_back(dLH);
    vdBel.push_back(dBel);

    count++;
  }
  // --------------------------------------------

  // EXPERIMENT ---------------------------------

double expDx[20]= {22.5	,
23.5	,
24.5	,
25.5	,
26.5	,
27.5	,
28.5	,
29.5	,
30.5	,
31.5	,
32.5	,
33.5	,
34.5	,
35.5	,
36.5	,
37.5	,
38.5	,
39.5	,
40.5	,
41.5	};

double expDy[20]={
0.595726876922256	,
0.543213473232366	,
0.494283286886849	,
0.421774881771817	,
0.341156378590349	,
0.307099360773692	,
0.262569839264486	,
0.232682384534885	,
0.216888228577843	,
0.23470124429325	,
0.240399882846714	,
0.228134082861555	,
0.32920533438895	,
0.296543336595945	,
0.367649491313769	,
0.413065289661424	,
0.437192259660495	,
0.440938347607344	,
0.378811580496408	,
0.404394227640296	};



double expDyErr[20]={
0.028324914039265	,
0.0218144191759532	,
0.0142504076860264	,
0.0146454365657758	,
0.0088989334482194	,
0.00824524006184739	,
0.0108064634609478	,
0.0111948065344291	,
0.00904509721171868	,
0.0119548919107791	,
0.0180227167275249	,
0.0201388245171666	,
0.0213346009360056	,
0.0193056550679402	,
0.0155440903197556	,
0.0136119544115286	,
0.0204840651640718	,
0.0248743996273666	,
0.023151640426499	,
0.0183223089693594	};






double expDyErr2[20]={0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5};

double expPx[15]={31.5	,
32.5	,
33.5	,
34.5	,
35.5	,
36.5	,
37.5	,
38.5	,
39.5	,
40.5	,
41.5	,
42.5	,
43.5	,
44.5	,
45.5	};

double expPy[15]={
0.927940963711002	,
0.886373091269099	,
0.997232654140559	,
0.850828825193078	,
0.975515571786442	,
0.799588217203809	,
0.816607420458171	,
0.982225026964146	,
0.888146409552698	,
0.901367241759583	,
1.07321471789101	,
1.0011298892341	,
0.966230221622026	,
1.1926071484492	,
0.865315617070875	};





double expPyErr[15]={
0.0990008372934202	,
0.14013781860333	,
0.0725936416502272	,
0.066957912885687	,
0.0518874166956934	,
0.0314030147725784	,
0.0396294675991484	,
0.0474887746673471	,
0.0504518388646933	,
0.0432543357898785	,
0.0509856961378907	,
0.0376160621523655	,
0.043968910197424	,
0.0649913895880728	,
0.0576066735490418	};







double expPyErr2[15]={0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5,
	0.5};




  /*  double expDx[16] =   {22.8659168020674,
			24.4932387071222,
			25.4884605201683,
			26.4478963605459,
			27.8768790582849,
			30.0521503123866,
			31.6308395766488,
			32.4198704973062,

			34.0083546157897,
			35.4165962236503,
			36.4835938795446,
			37.5238982330089,
			38.5397371697557,
			39.5330583015003,
			40.5055760670275,
			41.3};

  double expDy[16] =   {9.52350198819289,
			5.64694896708287,
			3.57621524475264,
			2.22306813461932,
			1.84242236351724,
			0.839564977947224,
			0.533641048642503,
			0.811298146594814,
			0.452491201918824,
			0.628258484221101,
			0.700018256426215,
			0.598634895567586,
			0.510176367468261,
			0.457345826343762,
			0.341874900915227,
			0.41491476578731};

  double expPx[9] =    {30.6251132222237,
			32.9493281027741,
			34.6378972639558,
			37.6049166106808,
			39.1281033807373,
			41.0541710577087,
			43.6797461320023,
			43.5133105051628,
			45.6589118858263};

  double expPy[9] =    {15.4053624168101,
			10.0326697910848,
			11.2010890871452,
			6.28484827873840,
			5.93048525466625,
			5.18590383783819,
			4.93601112583176,
			3.77917031974004,
			6.32306635135636};
  */
  // --------------------------------------------

  //cout << vx.front() << " to " << vx.back() << "   count " << count << endl;

  TCanvas *canvThick = new TCanvas("canvThick","canvThick",1000,1000);
  canvThick->SetLogy();

  int HSS = 10, Perey = 4, BGreen = 1;


  TGraph* gp1 = new TGraph(vx.size(), &vx[0], &vpBG[0]);
    gp1->GetXaxis()->SetLimits(10.,80.);
    gp1->SetLineColor(kRed);
    gp1->SetLineStyle(BGreen);
    gp1->SetLineWidth(2);
    gp1->SetTitle("(p,p) Bechetti-Greenlees");

  TGraph* gp2 = new TGraph(vx.size(), &vx[0], &vpCH[0]);
    gp2->GetXaxis()->SetLimits(10.,80.);
    gp2->SetLineColor(kRed);
    gp2->SetLineStyle(Perey);
    gp2->SetLineWidth(2);
    gp2->SetTitle("(p,p) Chapel-Hill");

  TGraph* gp3 = new TGraph(vx.size(), &vx[0], &vpPer[0]);
    gp3->GetXaxis()->SetLimits(10.,80.);
    gp3->SetLineColor(kRed);
    gp3->SetLineStyle(HSS);
    gp3->SetLineWidth(2);
    gp3->SetTitle("(p,p) Perey");

  TGraph* gd1 = new TGraph(vx.size(), &vx[0], &vdHSS[0]);
    gd1->GetXaxis()->SetLimits(10.,80.);
    gd1->SetLineColor(kBlue);
    gd1->SetLineStyle(HSS);
    gd1->SetLineWidth(2);
    gd1->SetTitle("(d,d) HSS");

  TGraph* gd2 = new TGraph(vx.size(), &vx[0], &vdKD[0]);
    gd2->GetXaxis()->SetLimits(10.,80.);
    gd2->SetLineColor(kBlue);
    gd2->SetLineStyle(1);
    gd2->SetLineWidth(2);
    gd2->SetTitle("(d,d) Koning-Delaroche");

  TGraph* gd3 = new TGraph(vx.size(), &vx[0], &vd79D[0]);
    gd3->GetXaxis()->SetLimits(10.,80.);
    gd3->SetLineColor(kBlue);
    gd3->SetLineStyle(6);
    gd3->SetLineWidth(2);
    gd3->SetTitle("(d,d) 79DCV");

  TGraph* gd4 = new TGraph(vx.size(), &vx[0], &vdPer[0]);
    gd4->GetXaxis()->SetLimits(10.,80.);
    gd4->SetLineColor(kBlue);
    gd4->SetLineStyle(Perey);
    gd4->SetLineWidth(2);
    gd4->SetTitle("(d,d) Perey");

  TGraph* gd5 = new TGraph(vx.size(), &vx[0], &vdLH[0]);
    gd5->GetXaxis()->SetLimits(10.,80.);
    gd5->SetLineColor(kBlue);
    gd5->SetLineStyle(2);
    gd5->SetLineWidth(2);
    gd5->SetTitle("(d,d) Lohr-Haeberli");

  TGraph* gd6 = new TGraph(vx.size(), &vx[0], &vdBel[0]);
    gd6->GetXaxis()->SetLimits(10.,80.);
    gd6->SetLineColor(kBlue);
    gd6->SetLineStyle(9);
    gd6->SetLineWidth(2);
    gd6->SetTitle("(d,d) Belote 48Ca");

  TGraphErrors* expP = new TGraphErrors( 15, expPx, expPy, expPyErr2, expPyErr);
  expP->SetTitle("(p,p) Experiment");
    //expP->SetMarkerStyle(22);
  TGraphErrors* expD = new TGraphErrors(20, expDx, expDy, expDyErr2, expDyErr);
  expD->SetTitle("(d,d) Experiment");
    //expD->SetMarkerStyle(20);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(gp1);
  mg->Add(gp2);
  mg->Add(gp3);
  mg->Add(gd1);
  mg->Add(gd2);
  mg->Add(gd3);
  mg->Add(gd4);
  mg->Add(gd5);
  mg->Add(gd6);
  //mg->GetXaxis()->SetTitle("{#theta}_{CM}");
  //mg->GetYaxis()->SetTitle("Elastic counts [mb/sr]");
  mg->Draw("AC");
  mg->GetXaxis()->SetTitle("#theta_{CM} [deg]");
  mg->GetYaxis()->SetTitle("Ratio #sigma/#sigma_{Rutherford}");
  expP->Draw("same*");
  expD->Draw("same*");

  canvThick->BuildLegend();

}

void temptemp(){

chain->Draw("T_MUGAST_VAMOS>>hist(500,0,10000)","","");

}

void tempRunMe(double gamma, double width){

  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-" 
      + to_string(gamma)
      + ")<"
      + to_string(width);

  TFile *File22 = new TFile("../../../Outputs/Analysis/47K_Full_05Aug_MinWithNoMG3.root","READ");
  TTree* Tree22 = (TTree*) File22->Get("PhysicsTree");
  
  Tree22->Draw("Ex>>Ex22(120,-1,5)",gating.c_str(),"");
  TH1F* Ex22 = (TH1F*) gDirectory->Get("Ex22");
  
  //TFile *File30 = new TFile("../../../Outputs/Analysis/47K_Full_22July.root","READ");
  TFile *File30 = new TFile("../../../Outputs/Analysis/47K_Full_09Aug_retest05NoMG3.root","READ");
  TTree* Tree30 = (TTree*) File30->Get("PhysicsTree");

  Tree30->Draw("Ex>>Ex30(120,-1,5)",gating.c_str(),"");
  TH1F* Ex30 = (TH1F*) gDirectory->Get("Ex30");

  Ex22->SetLineColor(kRed);
  Ex22->Draw();
  Ex30->SetLineColor(kBlue);
  Ex30->Draw("same");

}

void ExThetaAnalysis(double gamma, double width, int version){
//  int version;
//  bool running = 1;

//  cout << "Constructing Ex:ThetaLab..." << endl;


  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && abs(AddBack_EDC-"
	        + to_string(gamma) + ") < " + to_string(width); 

  TCanvas *diagnoseTheta2 = new TCanvas("diagnoseTheta2","diagnoseTheta2",1000,1000);
  chain->Draw(
    "Ex:ThetaLab>>thetaHist(60,100,160,120,-1,5)", 
    //"abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8",
    gating.c_str(),
    "colz");
  TH2F* thetaHist = (TH2F*) gDirectory->Get("thetaHist");  


//  while(running){
//    cout << "Overlay projections (1) or candle plots (2)?" << endl;
//    cin >> version;
    cout << "Processing..." << endl;
    if(version==1){
      thetaHist->ProjectionY("tpy1",06.,15.); 
      TH1F* tpy1 = (TH1F*) gDirectory->Get("tpy1");  
      thetaHist->ProjectionY("tpy2",16.,25.); 
      TH1F* tpy2 = (TH1F*) gDirectory->Get("tpy2");  
      thetaHist->ProjectionY("tpy3",26.,35.); 
      TH1F* tpy3 = (TH1F*) gDirectory->Get("tpy3");  
      thetaHist->ProjectionY("tpy4",36.,45.); 
      TH1F* tpy4 = (TH1F*) gDirectory->Get("tpy4");  
      thetaHist->ProjectionY("tpy5",46.,55.); 
      TH1F* tpy5 = (TH1F*) gDirectory->Get("tpy5");  

      tpy1->SetLineColor(kRed);
      tpy2->SetLineColor(kOrange);
      tpy3->SetLineColor(kGreen);
      tpy4->SetLineColor(kBlue);
      tpy5->SetLineColor(kViolet);

      tpy1->Rebin(2);
      tpy2->Rebin(2);
      tpy3->Rebin(2);
      tpy4->Rebin(2);
      tpy5->Rebin(2);

      tpy1->Draw();
      tpy2->Draw("same");
      tpy3->Draw("same");
      tpy4->Draw("same");
      tpy5->Draw("same");
    }else if (version==2){
      thetaHist->GetXaxis()->SetRangeUser(105.,155.);
      thetaHist->RebinX(5);
      thetaHist->GetXaxis()->SetRangeUser(105.,155.);
      thetaHist->Draw("candlex6");

    }//else{running=0;}
//  }

}

void ForPoster_DiffCrossSec(){
  ifstream infile("DiffCrossSecInputfile.txt");
  vector<double> theta, p32, p32Pos, p32Neg, f72, f72Pos, f72Neg, zero;

  double a, b, c, d, e, f, g;
  while(infile){
    infile >> a >> b >> c >> d >> e >> f >> g;

    theta.push_back(a);
    zero.push_back(0.0);
    
    p32.push_back(b);
    p32Pos.push_back(c);
    p32Neg.push_back(d);
    
    f72.push_back(e);
    f72Pos.push_back(f);
    f72Neg.push_back(g);
  }

  TGraph *graphp32 = new TGraph(theta.size(), &theta[0], &p32[0]);
  TGraph *graphf72 = new TGraph(theta.size(), &theta[0], &f72[0]);

  int num = theta.size();

  TGraphAsymmErrors *graphp32error 
  = new TGraphAsymmErrors(num, &theta[0],  &p32[0], &zero[0], &zero[0], &p32Neg[0], &p32Pos[0]);

  TGraphAsymmErrors *graphf72error 
  = new TGraphAsymmErrors(num, &theta[0],  &f72[0], &zero[0], &zero[0], &f72Neg[0], &f72Pos[0]);


  graphp32->SetFillColor(kRed);
  graphp32->SetFillStyle(3005);
  graphf72->SetFillColor(kBlue);
  graphf72->SetFillStyle(3005);
  
  graphp32->Draw();
  graphf72->Draw("same");
  graphp32error->Draw("same");
  graphf72error->Draw("same");

/*
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(graphp32);
  mg->Add(graphf72);
  mg->Add(graphp32error);
  mg->Add(graphf72error);
  mg->Draw();
*/


}

/* MAIN FUNCTION */

void DrawPlots(){
  gStyle->SetOptStat("nemMrRi");
  LoadChainNP();
  
  cout << "==========================================" << endl;
  cout << "========== AVAILABLE FUNCTIONS ===========" << endl;
  cout << " 2D Matrices " << endl;
  cout << "\t- Draw_2DParticleGamma() "<< endl;
  cout << "\t- Draw_2DGammaGamma() "<< endl;
  cout << ""<< endl;
  cout << " Ungated histograms " << endl;
  cout << "\t- Draw_1DParticle() "<< endl;
  cout << "\t- Draw_1DParticle_MUST2() "<< endl;
  cout << "\t- Draw_1DGamma() "<< endl;
  cout << "\t- Draw_1DGamma_MG() "<< endl;
  cout << "\t- Draw_1DGamma_MM() "<< endl;
  cout << ""<< endl;
  cout << " Gated histograms " << endl;
  cout << "\t- GateParticle_SeeGamma(particle, width) "<< endl;
  cout << "\t- GateGamma_SeeParticle(gamma, width) "<< endl;
  cout << "\t- GateGamma_SeeGamma(gamma, width) "<< endl;
  cout << ""<< endl;
  cout << " Gated histograms w/ background selection" << endl;
  cout << "\t- GateParticle_SeeGamma_WithBG(particle, width, bgrnd) "<< endl;
  cout << "\t- GateGamma_SeeParticle_WithBG(gamma, width, bgrnd) "<< endl;
  cout << "\t- GateGamma_SeeGamma_WithBG(gamma, width, bgrnd) "<< endl;
  cout << ""<< endl;
  cout << " Specific functions" << endl;
  cout << "\t- CompareExsAt4MeV() "<< endl;
  cout << "\t- CompareSimExp() "<< endl;
  cout << "\t- MugastMisalignment() "<< endl;
  cout << "\t- ExPhiLab() "<< endl;
  cout << "\t- ExThetaLab() "<< endl;
  cout << "\t- ELabThetaLab() "<< endl;
  cout << "\t- MM5_ELabThetaLab() "<< endl;
  cout << "\t- MM5_ExThetaLab() "<< endl;
  cout << ""<< endl;
  cout << " Analysis functions" << endl;
  cout << "\t- FitKnownPeaks(histogram) "<< endl;
  cout << "\t- AGATA_efficiency(double Energy_kev) "<< endl;
  cout << "\t- CorrectForAGATAEffic(TH1F* hist) "<< endl;
  cout << "\t- CS(stateEnergy, stateSpin, orbital_l, orbital_j, nodes) "<< endl;
  cout << "\t---- 0.143, p3/2 = CS(0.143, 2, 1, 1.5, ?) "<< endl;
  cout << "\t---- 0.968, p1/2 = CS(0.968, 0, 1, 0.5, ?) "<< endl;
  cout << ""<< endl;
  cout << "==========================================" << endl;


/* ORIGINAL - DO NOT EDIT
  TCanvas *cG = new TCanvas("cG","cG",1000,1000);
  chain->Draw("Ex>>hcG(200,-3,7)","");
  chain->Draw("Ex>>hcG2(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600","same");
  TH1F* hcG = (TH1F*) gDirectory->Get("hcG");
  hcG->GetXaxis()->SetTitle("E_{x} [MeV]");
  hcG->GetYaxis()->SetTitle("Counts / (50 keV)");
  hcG->SetLineColor(kBlack);
  TH1F* hcG2 = (TH1F*) gDirectory->Get("hcG2");
  hcG2->SetLineColor(kRed);
*/

  /* ORIGINAL 4-PLOT KINEMATIC SCREEN */
/*
  TCanvas *c0 = new TCanvas("c0", "Kinematics", 1000, 1000);
  c0->Divide(2,2);
  c0->cd(1);
  gPad->SetLogz();
  chain->Draw("ELab:ThetaLab>>hKine(200,0,180,600,0,12)","abs(T_MUGAST_VAMOS-2777)<600","col");
  TH2F* hKine = (TH2F*) gDirectory->Get("hKine");
  hKine->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  hKine->GetYaxis()->SetTitle("E_{lab} [MeV]");
  plot_kine(Kdp, 0, kBlack, 2, 9);
  plot_kine(Kdp, 2, kRed, 2, 9);
  plot_kine(Kdp, 4, kBlue, 2, 9);
  plot_kine(Kdd, 0, kGreen+2, 2, 9);
  plot_kine(Kpp, 0, kYellow, 2, 9);
  plot_kine(Kdt, 0, kMagenta, 2, 9);
  plot_kine(K12C12C, 4, kCyan, 2, 9);
  plot_kine(Tidp, 4, 42, 2, 9);
  plot_kine(Tidt, 4, 42, 2, 9);
  plot_kine(Tidd, 4, 42, 2, 9);
  plot_kine(Ti12C12C, 4, 42, 2, 9);
  
  c0->cd(2);
  chain->Draw("Ex>>hEx(300,-3,7)","abs(T_MUGAST_VAMOS-2777)<600");
  TH1F* hEx = (TH1F*) gDirectory->Get("hEx");
  hEx->GetXaxis()->SetTitle("E_{x} [MeV]");
  hEx->GetYaxis()->SetTitle("Counts / (100 keV)");
  hEx->SetLineColor(kBlack);
  chain->Draw("Ex>>hEx_gateNoG(300,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && AddBack_EDC@.size()==0","same");
  chain->Draw("Ex>>hEx_gate(300,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && abs(AddBack_EDC-0.143)< 0.01","same");
  TH1F* hEx_gate = (TH1F*) gDirectory->Get("hEx_gate");
  hEx_gate->SetLineColor(kRed);
  auto AGATA_eff = new TF1("AGATA_eff","TMath::Exp([0]+[1]*TMath::Log(x)+[2]*pow(TMath::Log(x),2.0)+[3]*pow(TMath::Log(x),3.0)+[4]*pow(TMath::Log(x),4.0))",100,5000);
  AGATA_eff->SetParameter(0,-7.84071e+00);
  AGATA_eff->SetParameter(1, 6.44921e+00);
  AGATA_eff->SetParameter(2, -1.42899e+00);
  AGATA_eff->SetParameter(3, 1.37921e-01);
  AGATA_eff->SetParameter(4, -5.23947e-03);
  
  //TH1F* hEx_gate_scaled = (TH1F*) hEx_gate->Clone("hEx_gate_scaled");
  hEx_gate->Scale(1./(AGATA_eff->Eval(143)*0.01));
  hEx_gate->Draw("histsame");
  c0->Update();
  double ymax = gPad->GetUymax();
  plot_state(0, ymax, kBlack, 2, 9);
  plot_state(2, ymax, kRed, 2, 9);
  plot_state(3.8, ymax, kBlue, 2, 9);
  plot_state(4.644, ymax, kGreen+2, 2, 1);
  TLatex latex;
  latex.SetTextAlign(13);
  latex.SetTextSize(0.035);
  latex.SetTextAngle(90);
  latex.SetTextColor(kGreen+2);
  latex.DrawLatex(4.8,0.6*ymax,"S_{n} = 4.64 MeV");

  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);

  c0->cd(3);
  chain->Draw("AddBack_EDC:Ex>>hEgEx(300,-1,7,5000,0,5)","abs(T_MUGAST_VAMOS-2777)<600","col");
*/

/*
  chain->Draw("AddBack_EDC:Ex>>hEgEx(150,0,6,1000,0,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","col");
  TH2F* hEgEx = (TH2F*) gDirectory->Get("hEgEx");
  hEgEx->SetTitle("MUGAST#5 only");
  hEgEx->GetXaxis()->SetTitle("E_{x} [MeV]");
  hEgEx->GetYaxis()->SetTitle("E_{#gamma} [MeV]");
  TLine *line = new TLine(0,0,6,6);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->Draw();
*/
  
/*
  c0->cd(4);
  chain->Draw("AddBack_EDC>>hEg(4000,0,4)","abs(T_MUGAST_VAMOS-2777)<600");
  TH1F* hEg = (TH1F*) gDirectory->Get("hEg");
  hEg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg->GetYaxis()->SetTitle("Counts / (1 keV)");
  hEg->SetLineColor(kBlack);
  c0->Update();
  ymax = gPad->GetUymax();
  //plot_state(0.143, ymax, kMagenta, 2, 9);
  //plot_state(0.279, ymax, kCyan, 2, 9);
  //plot_state(1.863, ymax, kOrange, 2, 9);
*/

/*
  TCanvas *gammaspec = new TCanvas("gammaSpec", "gammaSpec", 1000, 1000);
  chain->Draw("AddBack_EDC>>hEg(4000,0,4)","abs(T_MUGAST_VAMOS-2777)<600");
  TH1F* hEg = (TH1F*) gDirectory->Get("hEg");
  hEg->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg->GetYaxis()->SetTitle("Counts / (1 keV)");
  hEg->SetLineColor(kBlack);
  //c0->Update();
  //ymax = gPad->GetUymax();
  //plot_state(0.143, ymax, kMagenta, 2, 9);
  //plot_state(0.279, ymax, kCyan, 2, 9);
  //plot_state(1.863, ymax, kOrange, 2, 9);
*/

/*  
  TCanvas *c1 = new TCanvas("c1", "Egamma, gated on Ex", 1000, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  chain->Draw("AddBack_EDC>>hEg_0p143(800,0,4)","abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-0.143)<0.1");
  TH1F* hEg_0p143 = (TH1F*) gDirectory->Get("hEg_0p143");
  hEg_0p143->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg_0p143->GetYaxis()->SetTitle("Counts / (5 keV)");
  hEg_0p143->SetLineColor(kBlack);

  c1->cd(2);
  chain->Draw("AddBack_EDC>>hEg_1p275(800,0,4)","abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-1.275)<0.1");
  TH1F* hEg_1p275 = (TH1F*) gDirectory->Get("hEg_1p275");
  hEg_1p275->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg_1p275->GetYaxis()->SetTitle("Counts / (5 keV)");
  hEg_1p275->SetLineColor(kBlack);

  c1->cd(3);
  chain->Draw("AddBack_EDC>>hEg_1p822(800,0,4)","abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-1.822)<0.1");
  TH1F* hEg_1p822 = (TH1F*) gDirectory->Get("hEg_1p822");
  hEg_1p822->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg_1p822->GetYaxis()->SetTitle("Counts / (5 keV)");
  hEg_1p822->SetLineColor(kBlack);

  c1->cd(4);
  chain->Draw("AddBack_EDC>>hEg_3p609(800,0,4)","abs(T_MUGAST_VAMOS-2777)<600 && abs(Ex-3.609)<0.1");
  TH1F* hEg_3p609 = (TH1F*) gDirectory->Get("hEg_3p609");
  hEg_3p609->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  hEg_3p609->GetYaxis()->SetTitle("Counts / (5 keV)");
  hEg_3p609->SetLineColor(kBlack);
  
*/

/*
  auto gr = K.GetKinematicLine3();
  gr->SetLineColor(kAzure+7);
  gr->SetLineWidth(3);
  gr->Draw("ac");
  K.SetExcitationHeavy(4);
  gr = K.GetKinematicLine3();
  gr->SetLineColor(kAzure+7);
  gr->SetLineWidth(2);
  gr->SetLineStyle(1);
  gr->Draw("c");

  AddTiStates(0); 
  AddTiStates(0.969); 
  AddTiStates(2.2292); 
  AddTiStates(2.419); 
  AddTiStates(3.223); 
  AddTiStates(3.332); 
  AddTiStates(3.622); 
  AddTiStates(4.388); 
  AddTiStates(4.458); 
  AddTiStates(4.719); 
  AddTiStates(4.852); 
  AddTiStates(5.151); 
*/   

}
