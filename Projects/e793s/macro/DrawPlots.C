#include "DoubleGausHeader.h"
#include "DrawPlots.h"


/* USE THIS SPACE TO TEST NEW FEATURES */
/* ONCE HAPPY, MOVE TO HEADER          */

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

//void ExThetaAnalysis(int version){
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












































///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void DrawPlots(){
  gStyle->SetOptStat("nei");
  LoadChainNP();
  
  cout << "==========================================" << endl;
  cout << "========== AVAILABLE FUNCTIONS ===========" << endl;
cout << " 2D Matrices " << endl;
  cout << "\t- Draw_2DParticleGamma() "<< endl;
  cout << "\t- Draw_2DGammaGamma() "<< endl;
  cout << ""<< endl;
  cout << " Ungated histograms " << endl;
  cout << "\t- Draw_1DParticle() "<< endl;
  cout << "\t- Draw_1DGamma() "<< endl;
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
  cout << " ExPhiLab_ForPoster()!!!!!!!!!!!!!!!!!"<< endl;
  cout << "==========================================" << endl;




  /************** OLD STUFF! Make these into functions, eventually ***************/

  /*
  TCanvas *cF = new TCanvas("cF","cF",1000,1000);
  gPad->SetLogz();
  chain->Draw("ELab:ThetaLab>>hcF(130,50,180,600,0,12)","abs(T_MUGAST_VAMOS-2777)<600","col");
  TH2F* hcF = (TH2F*) gDirectory->Get("hcF");
  hcF->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  hcF->GetYaxis()->SetTitle("E_{lab} [MeV]");
  plot_kine(Kdp, 0, kBlack, 2, 9);
  plot_kine(Kdp, 2, kRed, 2, 9);
  plot_kine(Kdp, 4, kBlue, 2, 9);
  plot_kine(Kdd, 0, kGreen+2, 2, 9);
  plot_kine(Kpp, 0, kYellow, 2, 9);
  */


  /* MUGAST testing phi dependance */
  //TCanvas *diagnose = new TCanvas("diagnose","diagnose",1000,1000);
//  chain->Draw("PhiLab","abs(T_MUGAST_VAMOS-2777)<600 && (Mugast.TelescopeNumber==1 ||  Mugast.TelescopeNumber==2 || Mugast.TelescopeNumber==3 || Mugast.TelescopeNumber==4 || Mugast.TelescopeNumber==5 || Mugast.TelescopeNumber==7)","");
//  chain->Draw("PhiLab>>h1","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==1","");
//  chain->Draw("PhiLab>>h2","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==2","same");
//  chain->Draw("PhiLab>>h3","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==3","same");
//  chain->Draw("PhiLab>>h4","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==4","same");
//  chain->Draw("PhiLab>>h5","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","same");
//  chain->Draw("PhiLab>>h7","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==7","same");
//  TH2F* histDiag = (TH2F*) gDirectory->Get("histDiag");
  

//  TH1F* h1 = (TH1F*) gDirectory->Get("h1");
//  h1->SetTitle("Phi of MUGAST telescopes");
//  h1->GetXaxis()->SetTitle("PhiLab [deg]");
//  h1->GetYaxis()->SetTitle("Counts");
//  h1->SetLineColor(kRed);
  //h1->SetFillStyle(3001);
  //h1->SetFillColor(kRed);
//  TH1F* h2 = (TH1F*) gDirectory->Get("h2");
//  h2->SetLineColor(kOrange);
  //h2->SetFillStyle(3001);
  //h2->SetFillColor(kOrange);
//  TH1F* h3 = (TH1F*) gDirectory->Get("h3");
//  h3->SetLineColor(kGreen);
  //h3->SetFillStyle(3001);
  //h3->SetFillColor(kGreen);
//  TH1F* h4 = (TH1F*) gDirectory->Get("h4");
//  h4->SetLineColor(kTeal);
  //h4->SetFillStyle(3001);
  //h4->SetFillColor(kTeal);
//  TH1F* h5 = (TH1F*) gDirectory->Get("h5");
//  h5->SetLineColor(kBlue);
  //h5->SetFillStyle(3001);
  //h5->SetFillColor(kBlue);
//  TH1F* h7 = (TH1F*) gDirectory->Get("h7");
//  h7->SetLineColor(kViolet);
  //h7->SetFillStyle(3001);
  //h7->SetFillColor(kViolet);
//  h7->GetYaxis()->SetRangeUser(0.,4000.);
//  h7->GetXaxis()->SetRangeUser(-180.,180.);

  //auto legend = new TLegend(0.1,0.7,0.48,0.9);
  //legend->AddEntry(h1,"MUGAST 1","f");
  //legend->AddEntry(h2,"MUGAST 2","f");
  //legend->AddEntry(h3,"MUGAST 3","f");
  //legend->AddEntry(h4,"MUGAST 4","f");
  //legend->AddEntry(h5,"MUGAST 5","f");
  //legend->AddEntry(h7,"MUGAST 7","f");
  //legend->Draw();

  /* MUGAST good Phi-Ex histogram*/
  /* MUGAST testing theta dependance by detector */
//  TCanvas *diagnose = new TCanvas("diagnose","diagnose",1000,1000);
//  chain->Draw("Ex:ThetaLab>>histDiag(36,0,180,200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600","colz");
//  chain->Draw("Ex:ThetaLab>>histDiag(36,0,180,200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==1","colz");
//  chain->Draw("Ex:ThetaLab>>histDiag(36,0,180,200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==2","colz");
//  chain->Draw("Ex:ThetaLab>>histDiag(36,0,180,200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==3","colz");
//  chain->Draw("Ex:ThetaLab>>histDiag(36,0,180,200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==4","colz");
//  chain->Draw("Ex:ThetaLab>>histDiag(36,0,180,200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","colz");
//  chain->Draw("Ex:ThetaLab>>histDiag(36,0,180,200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==7","colz");
//  TH2F* histDiag = (TH2F*) gDirectory->Get("histDiag");

  /*
  TCanvas *cG = new TCanvas("cG","cG",1000,1000);
  chain->Draw("Ex>>hcG(200,-3,7)","");
  chain->Draw("Ex>>hcG_T1(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==1","");
  chain->Draw("Ex>>hcG_T2(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==2","same");
  chain->Draw("Ex>>hcG_T3(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==3","same");
  chain->Draw("Ex>>hcG_T4(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==4","same");
  chain->Draw("Ex>>hcG_T5(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","same");
  chain->Draw("Ex>>hcG_T7(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==7","same");
  TH1F* hcG = (TH1F*) gDirectory->Get("hcG");
  hcG->GetXaxis()->SetTitle("E_{x} [MeV]");
  hcG->GetYaxis()->SetTitle("Counts / (50 keV)");
  hcG->SetLineColor(kBlack);
  TH1F* hcG_T1 = (TH1F*) gDirectory->Get("hcG_T1");
  hcG_T1->SetLineColor(kRed);
  hcG_T1->SetFillStyle(3001);
  hcG_T1->SetFillColor(kRed);
  TH1F* hcG_T2 = (TH1F*) gDirectory->Get("hcG_T2");
  hcG_T2->SetLineColor(kOrange);
  hcG_T2->SetFillStyle(3001);
  hcG_T2->SetFillColor(kOrange);
  TH1F* hcG_T3 = (TH1F*) gDirectory->Get("hcG_T3");
  hcG_T3->SetLineColor(kGreen);
  hcG_T3->SetFillStyle(3001);
  hcG_T3->SetFillColor(kGreen);
  TH1F* hcG_T4 = (TH1F*) gDirectory->Get("hcG_T4");
  hcG_T4->SetLineColor(kTeal);
  hcG_T4->SetFillStyle(3001);
  hcG_T4->SetFillColor(kTeal);
  TH1F* hcG_T5 = (TH1F*) gDirectory->Get("hcG_T5");
  hcG_T5->SetTitle("Misalignment of MUGAST telescopes");
  hcG_T5->GetXaxis()->SetTitle("E_{x} [MeV]");
  hcG_T5->GetYaxis()->SetTitle("Counts / (50 keV)");
  hcG_T5->SetLineColor(kBlue);
  hcG_T5->SetFillStyle(3001);
  hcG_T5->SetFillColor(kBlue);
  TH1F* hcG_T7 = (TH1F*) gDirectory->Get("hcG_T7");
  hcG_T7->SetLineColor(kViolet);
  hcG_T7->SetFillStyle(3001);
  hcG_T7->SetFillColor(kViolet);
  hcG_T5->GetYaxis()->SetRangeUser(0.,500.);
  hcG_T5->GetXaxis()->SetRangeUser(-3.,7.);
 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
  // MUGAST Displacement Ex histogram //
  TCanvas *cMisaligned = new TCanvas("cMisaligned","cMisaligned",1000,1000);
//  chain->Draw("Ex>>hcG(200,-3,7)","");
  chain->Draw("Ex>>hcG_T1(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==1","");
  chain->Draw("Ex>>hcG_T2(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==2","same");
  chain->Draw("Ex>>hcG_T3(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==3","same");
  chain->Draw("Ex>>hcG_T4(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==4","same");
  chain->Draw("Ex>>hcG_T5(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","same");
  chain->Draw("Ex>>hcG_T7(200,-3,7)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==7","same");
//  TH1F* hcG = (TH1F*) gDirectory->Get("hcG");
//  hcG->GetXaxis()->SetTitle("E_{x} [MeV]");
//  hcG->GetYaxis()->SetTitle("Counts / (50 keV)");
//  hcG->SetLineColor(kBlack);
  TH1F* hcG_T1 = (TH1F*) gDirectory->Get("hcG_T1");
  hcG_T1->SetTitle("Misalignment of MUGAST telescopes");
  hcG_T1->GetXaxis()->SetTitle("E_{x} [MeV]");
  hcG_T1->GetYaxis()->SetTitle("Counts / (50 keV)");
  hcG_T1->SetLineColor(kRed);
  hcG_T1->SetFillStyle(3244);
  hcG_T1->SetFillColor(kRed);
  hcG_T1->GetYaxis()->SetRangeUser(0.,200.);
  TH1F* hcG_T2 = (TH1F*) gDirectory->Get("hcG_T2");
  hcG_T2->SetLineColor(kOrange);
  hcG_T2->SetFillStyle(3344);
  hcG_T2->SetFillColor(kOrange);
  TH1F* hcG_T3 = (TH1F*) gDirectory->Get("hcG_T3");
  hcG_T3->SetLineColor(kGreen);
  hcG_T3->SetFillStyle(3444);
  hcG_T3->SetFillColor(kGreen);
  TH1F* hcG_T4 = (TH1F*) gDirectory->Get("hcG_T4");
  hcG_T4->SetLineColor(kTeal);
  hcG_T4->SetFillStyle(3544);
  hcG_T4->SetFillColor(kTeal);
  TH1F* hcG_T5 = (TH1F*) gDirectory->Get("hcG_T5");
  hcG_T5->SetLineColor(kBlue);
  hcG_T5->SetFillStyle(3644);
  hcG_T5->SetFillColor(kBlue);
  TH1F* hcG_T7 = (TH1F*) gDirectory->Get("hcG_T7");
  hcG_T7->SetLineColor(kViolet);
  hcG_T7->SetFillStyle(3644);
  hcG_T7->SetFillColor(kViolet);
//  hcG_T7->GetYaxis()->SetRangeUser(0.,750.);

  TLine *line = new TLine(0,0,0,200);
  line->SetLineColor(kBlack);
  line->Draw();

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(hcG_T1,"MUGAST 1","f");
  legend->AddEntry(hcG_T2,"MUGAST 2","f");
  legend->AddEntry(hcG_T3,"MUGAST 3","f");
  legend->AddEntry(hcG_T4,"MUGAST 4","f");
  legend->AddEntry(hcG_T5,"MUGAST 5","f");
  legend->AddEntry(hcG_T7,"MUGAST 7","f");
  legend->Draw();
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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
/*  chain->Draw("AddBack_EDC:Ex>>hEgEx(150,0,6,1000,0,5)","abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber==5","col");
  TH2F* hEgEx = (TH2F*) gDirectory->Get("hEgEx");
  hEgEx->SetTitle("MUGAST#5 only");
  hEgEx->GetXaxis()->SetTitle("E_{x} [MeV]");
  hEgEx->GetYaxis()->SetTitle("E_{#gamma} [MeV]");
  TLine *line = new TLine(0,0,6,6);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->Draw();
*/

  //  plot_state(0.143, ymax, kYellow, 2, 9);
  
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


  // K states
  /*auto gr = K.GetKinematicLine3();
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
