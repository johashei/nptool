/* Predefine functions */
vector<vector<double>> GetExpDiffCross(double Energy);
TH1F* PullThetaLabHist(int i, double minTheta, double gatesize);
TH1F* PullPhaseSpaceHist(int i, double minTheta, double gatesize);
void Scale(TGraph* g , TGraphErrors* ex);
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j);
double ToMininize(const double* parameter);
TGraph* FindNormalisation(TGraph* theory, TGraphErrors* experiment);
TList* peakFitList = new TList();

/* Global variables */
vector<Double_t> anglecentres, anglewidth;
TGraph* currentThry;
TGraphErrors* staticExp;
int indexE;
double globalS, globalSerr;

/* Output volume toggle */
bool loud = 1;

/* Scale method toggle */
bool scaleTogether = 1;

/* String for image */
string orbitalname;
string orbital;


////////////////////////////////////////////////////////////////////////////////
void canclone(TCanvas* major, int padNum, string name){
  string minorName = "./CS2_Figures/" + name + ".root";
  cout << minorName << endl;
  TFile* minorFile = new TFile(minorName.c_str(),"READ");
  TCanvas* minor = (TCanvas*) minorFile->FindObjectAny("c_peakFits");
 
  major->cd(padNum);
  TPad *pad = (TPad*)gPad;
  minor->cd();
  TObject *obj, *clone;
  pad->Range(minor->GetX1(),minor->GetY1(),minor->GetX2(),minor->GetY2());
  pad->SetTickx(minor->GetTickx());
  pad->SetTicky(minor->GetTicky());
  pad->SetGridx(minor->GetGridx());
  pad->SetGridy(minor->GetGridy());
  pad->SetLogx(minor->GetLogx());
  pad->SetLogy(minor->GetLogy());
  pad->SetLogz(minor->GetLogz());
  pad->SetBorderSize(minor->GetBorderSize());
  pad->SetBorderMode(minor->GetBorderMode());
  minor->TAttLine::Copy((TAttLine&)*pad);
  minor->TAttFill::Copy((TAttFill&)*pad);
  minor->TAttPad::Copy((TAttPad&)*pad);

  TIter next(minor->GetListOfPrimitives());
  gROOT->SetSelectedPad(pad);
  while ((obj=next())) {
     clone = obj->Clone();
     obj->Draw("SAME");
     pad->GetListOfPrimitives()->Add(clone,obj->GetDrawOption());
  }
  pad->Modified();
  pad->Update();
  major->cd(padNum);
  pad->Draw();
}

////////////////////////////////////////////////////////////////////////////////
void CS_Diagnosis(){

  auto majorCanv = new TCanvas("CompareCanv","CompareCanv",1500,1500);
  majorCanv->Divide(3,3);
  canclone(majorCanv, 1, "c_peakFits_110_112"); 
  canclone(majorCanv, 2, "c_peakFits_115_118"); 
  canclone(majorCanv, 3, "c_peakFits_120_122"); 
  canclone(majorCanv, 4, "c_peakFits_125_128"); 
  canclone(majorCanv, 5, "c_peakFits_130_132"); 
  canclone(majorCanv, 6, "c_peakFits_135_138"); 
  canclone(majorCanv, 7, "c_peakFits_140_142"); 
  canclone(majorCanv, 8, "c_peakFits_145_148"); 
  canclone(majorCanv, 9, "c_peakFits_150_152"); 

}

////////////////////////////////////////////////////////////////////////////////
void CS(){
/* Overload function */
/*
  cout << " Inputs:\n Experimental...\n\t - Energy of state\n\t - Spin of state" << endl;
  cout << " Theory...\n\t - Orbital l\n\t - Orbital j\n\t - Orbital n\n\n" << endl;

  cout << "  0f5/2 | -----  |  --?-- |   l=3, j=2.5, n=0" << endl;
  cout << "        |        |        |" << endl;
  cout << "  1p1/2 | -----  |  --?-- |   l=1, j=0.5, n=1" << endl;
  cout << "  1p3/2 | -----  |  --?-- |   l=1, j=1.5, n=1" << endl;
  cout << "  0f7/2 | -----  |  ===== |" << endl;
  cout << "        |        |        |" << endl;
  cout << "  0d3/2 | -----  |  ===== |" << endl;
  cout << "  1s1/2 | --x--  |  ===== |" << endl;
  cout << "        |        |        |" << endl;
  cout << "        |   p    |    n   |" << endl;
*/

  cout << "- CS(stateE, stateSp, orb_l, orb_j, nodes) "<< endl;
  cout << "---- 0.143, p3/2 = CS(0.143, 2, 1, 1.5) "<< endl;
  cout << "---- 0.279, p3/2 = CS(0.279, 2, 1, 1.5) "<< endl;
  cout << "---- 0.728, f7/2 = CS(0.728, 3, 3, 3.5) "<< endl;
  cout << "---- 0.968, p1/2 = CS(0.968, 0, 1, 0.5) "<< endl;
  cout << "---- 1.410, p3/2 = CS(1.410, 1, 1, 1.5) "<< endl;
  cout << "---- 1.981, p3/2 = CS(1.981, 1, 1, 0.5) "<< endl;
  cout << "---- 2.410, p3/2 = CS(2.410, 0, 1, 0.5) "<< endl;
  cout << "---- 3.2  , f7/2 = CS(3.2  , 3, 3, 3.5) "<< endl;
  cout << "---- 3.6  , f5/2 = CS(3.6  , 3, 3, 2.5) "<< endl;
  cout << "---- 3.8  , f5/2 = CS(3.8  , 3, 3, 2.5) "<< endl;
  cout << "---- 4.1  , f5/2 = CS(4.1  , 3, 3, 2.5) "<< endl;
  cout << "---- 4.4  , f5/2 = CS(4.4  , 3, 3, 2.5) "<< endl;


}

////////////////////////////////////////////////////////////////////////////////
void CS(double Energy, double Spin, double spdf, double angmom){
  // p3/2 -> spdf = 1, angmom = 1.5
  // J0 is incident spin, which is 47K g.s. therefore J0 = 1/2
  double J0 = 0.5;
  //double ElasticNorm = 5.8, ElasticNormErr = 1.3; // DeuteronNorm in elastics, 5.8 +- 1.3
  //double ElasticNorm = 2.70, ElasticNormErr = 0.46; // For thicker 
  //double ElasticNorm = 2.8, ElasticNormErr = 0.7; // DeuteronNorm in elastics, 2.8 +- 0.7
  double ElasticNorm = 2.363, ElasticNormErr = 0.000; // DeuteronNorm in elastics, reconstructed from determined thickness

  orbitalname.clear();
  orbital.clear();
  if(spdf==1){
    if(angmom==0.5){
      orbitalname="p_{1/2}";
      orbital="p12";
    } else if (angmom==1.5){
      orbitalname="p_{3/2}";
      orbital="p32";
    } else { 
      orbitalname;
      orbital="???";
    }
  } else if (spdf==3){
    if(angmom==2.5){
      orbitalname="f_{5/2}";
      orbital="f52";
    } else if (angmom==3.5) {
      orbitalname="f_{7/2}";
      orbital="f72";
    } else { 
      orbitalname="?????";
      orbital="???";
    }
  } else {
    orbitalname="?????";
    orbital="???";
  }

  /* Reduce by factor of 10,000 */
  ElasticNorm /= 10000.;
  //ElasticNorm /= 1000.;
  ElasticNormErr /= 10000.;
  //ElasticNormErr /= 1000.;
  double nodes;

  if(spdf==1){
    nodes=1;
  }
  else if(spdf==3){
    nodes=0;
  }
  else{
    cout << " INPUT NODES::" << endl;
    cin >> nodes;
  }

  /* Clean global variables, in case of multiple runs */
  indexE = 100;
  anglecentres.clear();
  anglewidth.clear();
  globalS=0.;
  globalSerr=0.;
  peakFitList->Clear();

  /* Retrieve array index of the entered peak energy */
  /* numpeaks and Energy[] defined globally in KnownPeakFitter.h */
  bool found = 0;
  for(int i=0;i<numPeaks;i++){
    if(abs(Energy-means[i])<0.01){
      cout << "========================================================" << endl;
      cout << "Identified as state #" << i << ", E = " << means[i] << endl;
      indexE = i;
      found = 1;
    }
  }
  if(!found){
    cout << "========================================================" << endl;
    cout << "NO STATE AT THAT ENERGY INDENTIFIED!! CHECK KNOWN PEAKS!!" << endl;
    return;
  }

  /* Solid Angle (from simulation) */
  /* ADD OPTION TO CHANGE SOLID ANGLE FILE DEPENDING ON PEAK!!!!*/
  //auto file = new TFile("../SolidAngle_HistFile_47Kdp_26May22_v2_0143.root");
  //auto file = new TFile("../SolidAngle_HistFile_19Jul22_47Kdp_0p000.root");
  //auto file = new TFile("../SolidAngle_HistFile_19Jul22_47Kdp_1p981.root");
  //auto file = new TFile("../SolidAngle_HistFile_19Jul22_47Kdp_4p393.root");
  auto file = new TFile("../SolidAngle_HistFile_30Jul22_47Kdp_0p000_ThetaBin0p5.root");
  //auto file = new TFile("../SolidAngle_HistFile_New.root");
  /* ADD OPTION TO CHANGE SOLID ANGLE FILE DEPENDING ON PEAK!!!!*/
  TH1F* SolidAngle = (TH1F*) file->FindObjectAny("SolidAngle_Lab_MG");
  TCanvas* c_SolidAngle = new TCanvas("c_SolidAngle","c_SolidAngle",1000,1000);
  SolidAngle->Draw("HIST");
  SolidAngle->GetXaxis()->SetRangeUser(100.,160.);
  /* (canvas deleted after Area/SA calculation) */
 
  /* Area of experimental peaks */
  TCanvas* c_PeakArea = new TCanvas("c_PeakArea","c_PeakArea",1000,1000);
  vector<vector<double>> areaArray = GetExpDiffCross(means[indexE]);
  delete c_PeakArea;

  // Array: peakenergy, peakarea, areaerror, anglemin, anglemax
  if(loud){
    for(int i=0; i<areaArray.size();i++){
      cout << i << " " 
    	   << areaArray[i][0] << " " 
	   << areaArray[i][1] << " "
	   << areaArray[i][2] << " "
	   << areaArray[i][3] << " "
	   << areaArray[i][4] << endl;
    }
  }

  /* AoSA = Area/Solid Angle [#/msr] */
  /* dSdO = Experimental Diff. Cross Sect. (Area/Solid Angle *Norm) [mb/msr] */
  vector<Double_t> AoSA, AoSAerr;
  vector<Double_t> expDCS, expDCSerr;
  for(int i=0; i<areaArray.size();i++){
    int binmin = SolidAngle->FindBin(areaArray[i][3]+0.0001);
    int binmax = SolidAngle->FindBin(areaArray[i][4]-0.0001);

//cout << "binmin " << binmin << " to binmax " << binmax << endl;

    anglecentres.push_back(((areaArray[i][4]-areaArray[i][3])/2.)+areaArray[i][3]);
    anglewidth.push_back(areaArray[i][4]-areaArray[i][3]);

    double tempsum=0, tempsumerr=0;
    for(int x=binmin;x<binmax+1;x++){
      tempsum += SolidAngle->GetBinContent(x);
      tempsumerr += SolidAngle->GetBinError(x);
      if(loud){cout << x << endl;}
    }
    if(loud){cout << "TEST CHECK " << tempsum << " +- " << tempsumerr << endl;}

    double SAerr;
    double SA = SolidAngle->IntegralAndError(binmin,binmax,SAerr);
    SA = SA*1000.;       //sr->msr
    SAerr = SAerr*1000.; //sr->msr

    /* Area over Solid angle ONLY */
    AoSA.push_back(areaArray[i][1]/SA);
    AoSAerr.push_back((areaArray[i][1]/SA) 
		    * sqrt( pow(areaArray[i][2]/areaArray[i][1],2) 
			    + pow(SAerr/SA,2)));

    /* Differential Cross Section */
    /* NOTE: DON'T INCLUDE NORM ERROR IN ERR PROPAGATION AS IT IS SYSTEMATIC! */
    double tempvalue = (areaArray[i][1]/SA)*ElasticNorm; 
    expDCS.push_back(tempvalue);
    double temperror = tempvalue
		     * sqrt( pow(areaArray[i][2]/areaArray[i][1],2)
		           + pow(SAerr/SA,2)
			   //+ pow(ElasticNormErr/ElasticNorm,2)
			   );
    expDCSerr.push_back(temperror);

    if(loud){
      cout << "Angle " << areaArray[i][3] << " to " << areaArray[i][4] 
	   << ", centre " << anglecentres[i]
	   << ": Area = " << areaArray[i][1] << " +- " << areaArray[i][2] << " cnts" 
	   << "  SA = " << SA << " +- " << SAerr << " msr" 
           << "  Area/SA = " << AoSA[i] << " +- " << AoSAerr[i] << " counts/msr"
	   << setprecision(5)
           << "  Norm = " << ElasticNorm << " +- " << ElasticNormErr
	   << " (don't include norm err in propagation)"
           << "  dSdO = " << tempvalue << " +- " << temperror  
	   << setprecision(3)
	   << endl;
    }
  }
  //delete c_SolidAngle;
  
  /* Graph of Area/Solid Angle*/
  TCanvas* c_AoSA = new TCanvas("c_AoSA","c_AoSA",1000,1000);
  c_AoSA->SetLogy();
  TGraphErrors* gAoSA = new TGraphErrors(
		  anglecentres.size(), //n
		  &(anglecentres[0]), &(AoSA[0]), //x, y
		  //&(anglewidth[0]), &(AoSAerr[0]) );
		  0, &(AoSAerr[0]) );  //errX, errY 
  gAoSA->SetTitle("Area/SolidAngle");
  gAoSA->GetXaxis()->SetTitle("ThetaLab [deg]");
  gAoSA->GetYaxis()->SetTitle("Counts/#Omega [counts/msr]");
  gAoSA->Draw();

  /* Graph of experimental diff. cross sect (dSigma/dOmega) */
  TCanvas* c_dSdO = new TCanvas("c_dSdO","c_dSdO",1000,1000);
  c_dSdO->SetLogy();
  TGraphErrors* gdSdO = new TGraphErrors(
		  anglecentres.size(),
		  &(anglecentres[0]), &(expDCS[0]),
		  //&(anglewidth[0]), &(AoSAerr[0]) );
		  0, &(expDCSerr[0]) );
  gdSdO->SetTitle("Differential Cross Section");
  gdSdO->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  gdSdO->GetXaxis()->SetTitleOffset(1.2);
  gdSdO->GetXaxis()->SetTitleSize(0.04);
  gdSdO->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/msr]");
  gdSdO->GetYaxis()->SetTitleOffset(1.2);
  gdSdO->GetYaxis()->SetTitleSize(0.04);
  gdSdO->Draw();
  c_dSdO->Update();

  /* TWOFNR diff. cross section, in mb/msr */ 
  TCanvas* c_TWOFNR = new TCanvas("c_TWOFNR","c_TWOFNR",1000,1000);
  c_TWOFNR->SetLogy();
  TGraph* TheoryDiffCross = TWOFNR(means[indexE], J0, Spin, nodes, spdf, angmom); 
  TheoryDiffCross->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/msr]"); //msr set in func above
  TheoryDiffCross->GetXaxis()->SetTitle("ThetaLab [deg]");
  TheoryDiffCross->Draw();

  /* Convert AoSA into Differential Cross Section */
  //cout << " SCALING BY NORMALIZATION = " << ElasticNorm << endl;
  //gAoSA->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/msr]");


  /* Scaled and compared on same plot */ 
/* 
  cout << "USING BASIC SCALING METHOD..." << endl;
  TGraph* ScaledTWOFNR = TheoryDiffCross;
  TCanvas* c_Compare = new TCanvas("c_Compare","c_Compare",1000,1000);
  Scale(ScaledTWOFNR,gAoSA);
  c_Compare->SetLogy();
  gAoSA->SetLineColor(kRed);
  gAoSA->SetMarkerColor(kRed);
  gAoSA->SetMarkerStyle(21);
  gAoSA->Draw("AP");
  ScaledTWOFNR->Draw("SAME");
*/

  /* Using Chi2 minimizaiton */
  cout << "USING CHI2 MINIMIZAITON..." << endl;
  TCanvas* c_Chi2Min = new TCanvas("c_Chi2Min","c_Chi2Min",1000,1000);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.03);
  //c_Chi2Min->SetLogy();

  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
  pad1->SetTopMargin(0.1);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad1->SetLogy();
  pad1->SetTickx();
  pad1->SetTicky();
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.3);
  pad2->SetBorderMode(0);
  pad2->SetTickx();
  pad2->SetTicky();
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  TGraph* Final = FindNormalisation(TheoryDiffCross,gdSdO);
  gdSdO->SetLineColor(kRed);
  gdSdO->SetMarkerColor(kRed);
  gdSdO->SetMarkerStyle(21);
  /* Construct file name string */
  /**/  ostringstream tempstream;
  /**/  if(means[indexE]<1.0){tempstream << 0;}
  /**/  tempstream << (int) (means[indexE]*1000);
  /**/  tempstream << "_" << orbital; 
  /**/  tempstream << "_spin" << Spin;
  /**/  string tempstr = tempstream.str();
  /* Construct hist title string */
  /**/  ostringstream textstream;
  /**/  textstream << std::fixed << setprecision(3);
  /**/  textstream << "   " << means[indexE];
  /**/  textstream << " MeV, ";
  /**/  textstream <<  orbitalname;
  /**/  textstream << ", spin " << (int)Spin;
  /**/	textstream << " --> S = " << globalS 
  /**/	           << " +- " << globalSerr;
  /**/  string textstring = textstream.str(); 

  gdSdO->SetTitle(textstring.c_str());
  gdSdO->GetYaxis()->SetTitleOffset(1.3);
  gdSdO->GetYaxis()->SetTitleSize(0.042);
  gdSdO->GetXaxis()->SetRangeUser(103.,157.);
  gdSdO->Draw("AP");
  Final->Draw("SAME");


  pad2->cd();
  TGraphErrors* gResid = new TGraphErrors(*gdSdO);
  for(int n=0; n < gResid->GetN(); n++){
    double x = gdSdO->GetPointX(n);
    double residual = gdSdO->GetPointY(n) - Final->Eval(x);
    gResid->SetPoint(n,x,residual);
    gResid->SetPointError(n,0,gdSdO->GetErrorY(n));
  }
  TLine* markzero = new TLine(103.,0.,157.,0.);
  gResid->SetTitle("");
  gResid->GetXaxis()->SetRangeUser(103.,157.);
  gResid->GetYaxis()->SetTitle("Residuals");
  gResid->GetYaxis()->SetTitleSize(0.15);
  gResid->GetYaxis()->SetTitleOffset(0.36);
  gResid->GetYaxis()->SetLabelSize(0.08);
  gResid->GetYaxis()->SetNdivisions(305);
  gResid->GetXaxis()->SetTitleSize(0.15);
  gResid->GetXaxis()->SetTitleOffset(0.8);
  gResid->GetXaxis()->SetLabelSize(0.1);
  gResid->GetXaxis()->SetTickLength(0.1);
  gResid->Draw();
  markzero->SetLineStyle(2);
  markzero->Draw("same");

  //pad1->cd();
  //TText* orb = new TText(0.5,0.2,"TESTING!!!!!!!!!!!");//orbitalname.c_str());
  //orb->Draw("SAME");

  string savestring1 = "./CS2_Figures/"+tempstr+".root";
  string savestring2 = "./CS2_Figures/"+tempstr+".pdf";
  c_Chi2Min->SaveAs(savestring1.c_str());
  c_Chi2Min->SaveAs(savestring2.c_str());

  cout << YELLOW
       << " -----  C2S = (2J+1)S = (2*" << (int)Spin << "+1)S = (" 
       << (int)((2.*Spin)+1.) << ")S = " 
       << ((2.*Spin)+1.) * globalS << "  ----- " 
       << RESET << endl;


  //delete c_AoSA;
  //delete c_dSdO;
}

////////////////////////////////////////////////////////////////////////////////
vector<vector<double>> GetExpDiffCross(double Energy){
  cout << "========================================================" << endl;
  vector<vector<double>> AllPeaks_OneGate;
  vector<vector<double>> OnePeak_AllGates;
  /****CHANGE ANGLE GATING****/
  int numAngleBins = 20;
  double widthAngleBins = 2.5;
  double firstAngle = 105.;
  /***************************/
  double x[numAngleBins], y[numAngleBins];
  //TList* list = new TList();

  /* Determine scaling factor for PhaseSpace */
  TCanvas* c_ExSubPSpace = new TCanvas("c_ExSubPSpace","c_ExSubPSpace",1000,1000);
  double trackScale = 0.0;
  if(scaleTogether){
    TH1F* baseEx = PullThetaLabHist(0,firstAngle,widthAngleBins);
    TH1F* basePS = PullPhaseSpaceHist(0,firstAngle,widthAngleBins);
    for(int i=1; i<numAngleBins;i++){
      TH1F* addEx = PullThetaLabHist(i,firstAngle,widthAngleBins); baseEx->Add(addEx,1.);
      TH1F* addPS = PullPhaseSpaceHist(i,firstAngle,widthAngleBins); basePS->Add(addPS,1.);
    }

    /* Subtract flat background equal to smallest bin in range */
    baseEx->GetXaxis()->SetRange(baseEx->FindBin(-1.),baseEx->FindBin(1.));
    double minValueInRange = baseEx->GetBinContent(baseEx->GetMinimumBin());
    baseEx->GetXaxis()->UnZoom();
    cout << "Subtracting background of " << minValueInRange << endl;
    for(int b=1; b<baseEx->GetNbinsX() ; b++){
      baseEx->SetBinContent(b,baseEx->GetBinContent(b)-minValueInRange);
    }

    /* Begin scaling within range, track changes */
    basePS->Scale(0.1);
    trackScale = 0.1;
    int numAngleBinsScale = baseEx->GetNbinsX();
    int nbinlow = basePS->FindBin(4.); int nbinhigh = basePS->FindBin(8.0);
    for(int b=nbinlow; b<nbinhigh; b++){
      if(baseEx->GetBinContent(b) > 0.0 && basePS->GetBinContent(b) > baseEx->GetBinContent(b)){
	while(basePS->GetBinContent(b) > baseEx->GetBinContent(b)){
          basePS->Scale(0.99999);
          trackScale *= 0.99999;
        }
      }
    }
    baseEx->Add(basePS,-1.);
    baseEx->SetName("ExSubPSpace");
    baseEx->SetTitle("ExSubPSpace");
    baseEx->Draw();
    cout << "PhaseSpace -> ExpData scaling = " << trackScale << endl;
  }

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  // TEMPORARY!!! REMOVE LAST THREE BINS ON HIGH ENERGY STATES!!!
  if(means[indexE] > 3.0){numAngleBins-=3;}
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */


  for(int i=0; i<numAngleBins;i++){
    double min = firstAngle + (i*widthAngleBins);
    double max = min + widthAngleBins;
    cout << "===================================" << endl;
    cout << "min: " << min << " max: " << max << endl;
  
    stringstream tmp; tmp << fixed << setprecision(0); 
    tmp << "c_peakFits_" << min << "_" << max; 
    string tmp2 = tmp.str();
    TCanvas* c_peakFits = new TCanvas("c_peakFits",tmp2.c_str(),1000,1000);

    /* Retrieve theta-gated Ex TH1F from file GateThetaLabHistograms.root */
    /* To change angle gates, run GateThetaLab_MultiWrite() */
    TH1F* gate = PullThetaLabHist(i,firstAngle,widthAngleBins);
    TH1F* pspace = PullPhaseSpaceHist(i,firstAngle,widthAngleBins);

    /* Scale the Phase Space at this angle... */
    /* ... for all angles together */
    if(scaleTogether){
      gate->Add(pspace,-trackScale);
    } 
    /* ... or seperately for each angular bin */
    /* NOTE THAT THIS DOES NOT ACCOUNT FOR FLAT BACKGROUND */
    else {
      if(pspace->Integral() > 50.){ // Non-garbage histogram
        pspace->Scale(0.01);
	trackScale=0.01;
        int numAngleBins = gate->GetNbinsX();
        for(int b=0; b<numAngleBins; b++){
	  if(loud){cout << " FROM " << pspace->GetBinContent(b) << 
		         " > " << gate->GetBinContent(b); 
	  }
          while(pspace->GetBinContent(b) > gate->GetBinContent(b)){
            pspace->Scale(0.9999);
	    trackScale*=0.9999;
	  }
	  if(loud){cout << " TO " << pspace->GetBinContent(b) << 
	  	      " > " << gate->GetBinContent(b) << endl;
	  }
        }
        cout << " !!! SCALE FOR THIS ANGLE = " << trackScale << endl;
        gate->Add(pspace,-1);
      }
    }

    /* Subtract flat background equal to smallest bin in range */
    /* ????? */
    /*
    gate->GetXaxis()->SetRange(gate->FindBin(-1.),gate->FindBin(1.));
    double minValueInRange = gate->GetBinContent(gate->GetMinimumBin());
    gate->GetXaxis()->UnZoom();
    cout << "Subtracting background of " << minValueInRange << endl;
    for(int b=1; b<gate->GetNbinsX() ; b++){
      gate->SetBinContent(b,gate->GetBinContent(b)-minValueInRange);
    }
    */

    /* Retrieve array containing all fits, for one angle gate. *
     * Specific peak of interest selected from the vector by   *
     * global variable indexE                                  */

    //AllPeaks_OneGate = FitKnownPeaks_RtrnArry(gate, 0.0); cout << "!!!!!!! NO SLIDING SHIFT!!!!!" << endl;
    AllPeaks_OneGate = FitKnownPeaks_RtrnArry(gate, 0.0);  cout << "!!!!!!! WITH A VARIABLE SLIDING SHIFT!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    //double slideshift = 0.00218*(((max-min)/2.)+min) - 0.29645; AllPeaks_OneGate = FitKnownPeaks_RtrnArry(gate, slideshift);  cout << "!!!!!!! WITH A FIXED SLIDING SHIFT!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    
    /* Write PS-subtracted spectrum to list */
    //list->Add(gate);
    //list->Add(c_peakFits);
    gate->GetXaxis()->SetRangeUser(-1.0,+8.0);
    string filename = "./CS2_Figures/";
    filename += tmp2 + ".root";
    c_peakFits->SaveAs(filename.c_str());
    //auto tempfile = new TFile(filename.c_str(),"UPDATE"); //Reopen newly made file
    //auto templist = new TList();
    //templist->Add();


    /* Check correct OneGate vector is selected */
    cout << "area of " << means[indexE] << " = "
	 << AllPeaks_OneGate[indexE][1] 
	 << " +- " << AllPeaks_OneGate[indexE][2] 
	 << endl;

    /* Add min and max angle to end of relevant OneGate vector */
    AllPeaks_OneGate[indexE].push_back(min);
    AllPeaks_OneGate[indexE].push_back(max);

    /* Push relevant OneGate vector to end of AllGates vector */
    OnePeak_AllGates.push_back(AllPeaks_OneGate[indexE]);
    delete c_peakFits; 
  }

  /* Write PS-subtracted spectrum to file */
  //TFile* outfile = new TFile("GateThetaLab_ExMinusGatePhaseSpace.root","RECREATE");
  //list->Write("GateThetaLab_ExMinusPhaseSpace",TObject::kSingleKey);

  return OnePeak_AllGates;
}

////////////////////////////////////////////////////////////////////////////////
TH1F* PullThetaLabHist(int i, double minTheta, double gatesize){
  //TFile* file = new TFile("GateThetaLabHistograms.root","READ");
  //TFile* file = new TFile("GateThetaLabHistograms_20May22.root","READ");
  //TFile* file = new TFile("GateThetaLabHistograms_0p05BinWidth.root","READ");
  //TFile* file = new TFile("GateThetaLabHistograms_ReadMe.root","READ");
  //TFile* file = new TFile("GateThetaLabHistograms_26May22_v2.root","READ");
  //TFile* file = new TFile("GateThetaLabHistograms_26May22_v2_bin0p02.root","READ");
  //TFile* file = new TFile("GateThetaLabHistograms_2p5degAngles_26May22v2.root","READ");
  //TFile* file = new TFile("GateThetaLabHistograms_11Apr22_20angles.root","READ");
  TFile* file = new TFile("GateThetaLabHistograms_11Jul22.root","READ");
  string histname = "cThetaLabGate_" 
	          + to_string((int) (minTheta+(i*gatesize))) + "-" 
		  + to_string((int) (minTheta+((i+1)*gatesize)));
  cout << "Loading " << histname << endl;
  TList *list = (TList*)file->Get("GateThetaLabHistograms");
  TH1F* hist = (TH1F*)list->FindObject(histname.c_str());
//  file->Close();
  return hist;
}

////////////////////////////////////////////////////////////////////////////////
TH1F* PullPhaseSpaceHist(int i, double minTheta, double gatesize){
  //TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms_ReadMe.root","READ");
  //TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms_2p5degAngles_26May22v2.root","READ");
  TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms_11Jul22.root","READ");
  string histname = "cPSpaceThetaLabGate_" 
	          + to_string((int) (minTheta+(i*gatesize))) + "-" 
		  + to_string((int) (minTheta+((i+1)*gatesize)));
  cout << "Loading " << histname << endl;
  TList *list = (TList*)file->Get("GatePhaseSpaceThetaLabHistograms");
  TH1F* hist = (TH1F*)list->FindObject(histname.c_str());
  file->Close();
  return hist;
}

////////////////////////////////////////////////////////////////////////////////
void Scale(TGraph* g , TGraphErrors* ex){
  double scale;
  double mean = 0 ;
  double* eX = ex->GetX();
  double* eY = ex->GetY();
  double totalW = 0;
  double W = 0;
  for(int i = 0 ; i < ex->GetN() ; i++){
    if(eY[i]>1 && eY[i] <1e4){
      // Incremental Error weighted average
      W = 1./ex->GetErrorY(i);
      scale = eY[i]/g->Eval(eX[i]);
      totalW +=W;
      mean = mean + W*(scale - mean)/(totalW);
    }
  }

  //scaleTWOFNR = mean;
  cout << "SCALED THEORY BY " << mean << endl;
  cout << " therefore S = " << 1/mean << " ??" << endl;  
  
  double* x = g->GetX();
  double* y = g->GetY();
  for(unsigned int i = 0 ; i < g->GetN() ; i++)
    g->SetPoint(i,x[i],y[i]*mean);
}

////////////////////////////////////////////////////////////////////////////////
//TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j, const char* model){
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j){
  /* This function mved between directories in order to run TWOFNR in proper *
   * location. This is, weirdly, the least tempremental way of doing this.   */

  cout << "========================================================" << endl;
  int johnson, tandyval;
  cout << "Using Johnson-Soper ..."; johnson=5; tandyval=0;
  //cout << "Using Johnson-Tandy 1 ..."; johnson=6; tandyval=1;
  //cout << "Using Johnson-Tandy 2 ..."; johnson=6; tandyval=2;
  //cout << "Using Johnson-Tandy 3 ..."; johnson=6; tandyval=3;
  //cout << "Using Johnson-Tandy 4 ..."; johnson=6; tandyval=4;

  int modelA,modelB;
//  switch (model):{
//    case 'K': case 'k':{
//      cout << " ... Koning-Delaroche." << endl; modelA=6; modelB=4;
//    }
//    case 'C': case 'c':{
      cout << " ... and Chapel-Hill." << endl; modelA=2; modelB=2;
//    }      
//    case 'B': case 'b':{
//      cout << " ... Bechetti-Greenlees." << endl; modelA=1; modelB=1;
//    }
//  }


  char origDirchar[200];
  getcwd(origDirchar,200);
  string origDir{origDirchar};
  string twofnrDir = "/home/charlottepaxman/Programs/TWOFNR";
  cout << "Current directory    " << origDir << endl;
  cout << "Moving to directory  " << twofnrDir << endl;
  chdir(twofnrDir.c_str());
  //Check
  system("pwd"); 
  cout << "===================================" << endl;

  /* Delete existing tran.jjj & 24.jjj files */
  remove("tran.jjj");
  remove("24.jjj");

  double BeamEnergy =  7.7;
  double QValue = 2.274 - E;

  std::ofstream Front_Input("in.front");
  Front_Input << "jjj" << std::endl;
  Front_Input << "pipo" << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << BeamEnergy << std::endl;
  Front_Input << 47 << " " << 19 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << "0 0 0" << std::endl;
  Front_Input << l << " " << j << std::endl;
  Front_Input << n << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << QValue << std::endl; 
  Front_Input << 1 << std::endl;
  Front_Input << J0 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << johnson << std::endl;
  if(johnson==6){//JTandy selected, give version
    Front_Input << tandyval << std::endl;
  }
  Front_Input << 1 << std::endl;
  Front_Input << J << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << modelA << std::endl;
  Front_Input << modelB << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1.25 << " " << 0.65 << std::endl;
  Front_Input << 6.0 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << 0 << std::endl;

  Front_Input.close();

  cout << "Filled front20 input file." << endl;
  cout << "Executing front20..." << endl;
  system("(exec ~/Programs/TWOFNR/front20 < in.front > /dev/null)"); 
    ifstream checkfront("tran.jjj");
    if(!checkfront){
      //front20 fail!
      cout << " !! ERROR !! \n front20 failed to complete" << endl;
      return 0;
    } else {
      cout << "-> tran.jjj generated (success not guaranteed)" << endl;
      checkfront.close();
    }

  /* in.twofnr instructs twofnr20 to evaluate tran.jjj */
  cout << "Executing twofnr20..." << endl;
  system("(exec ~/Programs/TWOFNR/twofnr20 < in.twofnr > /dev/null)");
    ifstream checktwofnr("24.jjj");
    if(!checktwofnr){
      //twofnr20 fail!
      cout << " !! ERROR !! \n twofnr20 failed to complete" << endl;
      terminate();
//      return 0;
    } else {
      cout << "-> twofnr20 complete!" << endl;
      checktwofnr.close();
    }

  TGraph* CS = new TGraph("24.jjj");

  //mb/sr->mb/msr is x1/1000
  for(int i=0; i<CS->GetN(); i++){
    double x, newy;
    CS->GetPoint(i,x,newy);
    newy=newy/1000;
    CS->SetPoint(i,x,newy);
  }


  cout << "===================================" << endl;
  cout << "Current directory    " << twofnrDir << endl;
  cout << "Moving to directory  " << origDir << endl;
  chdir(origDir.c_str());
  system("pwd"); 
  cout << "========================================================" << endl;

  return CS;
}

////////////////////////////////////////////////////////////////////////////////
double Chi2(TGraph* theory, TGraphErrors* exper){
  double Chi2 = 0;
  double chi = 0;

  //cout << setprecision(8);
  //for(int i = 1 ; i < exper->GetN() ; i++){
  for(int i = 0 ; i < exper->GetN() ; i++){
    if(exper->GetPointY(i)>1.0e-8){ //0){
      //chi=(exper->Eval(anglecentres[i])-theory->Eval(anglecentres[i]) ) / (exper->GetErrorY(i));
      chi=(exper->GetPointY(i) - theory->Eval(anglecentres[i]) ) / (exper->GetErrorY(i));
      //cout << "COMPARE::::: " << exper->Eval(anglecentres[i]) << " TO " << exper->GetPointY(i) << endl;
      Chi2 +=chi*chi;
    }
  }
  if(loud){cout << "Chi2 = " << Chi2 << endl;}
  return Chi2;
  //cout << setprecision(3);
}

////////////////////////////////////////////////////////////////////////////////
double ToMininize(const double* parameter){
  static int f = 0 ;
  TGraph* g = new TGraph();
  double* X = currentThry->GetX(); // gets valies from global g1 = tgraph passed to find norm
  double* Y = currentThry->GetY();
  for(int i = 0 ; i < currentThry->GetN() ; i++){
    g->SetPoint(g->GetN(),X[i],parameter[0]*Y[i]); // scales this tgraph by parameter
  }
  double chi2  = Chi2(g,staticExp);  //compares theory tgraph to experimental tgrapherrors by chi2
  return chi2;
}

////////////////////////////////////////////////////////////////////////////////
TGraph* FindNormalisation(TGraph* theory, TGraphErrors* experiment){
  /* (dSdO)meas = S * (dSdO)calc */

  // Set global variable
  currentThry = theory;
  staticExp = experiment;

  // Construct minimiser
  const char* minName ="Minuit";const char* algoName="Migrad";
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);

  // Number of parameters (should only be 1 for me)
  int mysize = 1;

  // Create funciton wrapper for minimizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(&ToMininize,mysize);
  min->SetFunction(f);
  
  // Set range of parameter(??)
  double* parameter = new double[mysize];
  for(unsigned int i = 0 ; i < mysize ; i++){
    parameter[i] = 0.8;
    char name[4];
    sprintf(name,"S%d",i+1);
    min->SetLimitedVariable(i,name,parameter[i],0.01,0,10);
  }
 

  ///// TO IMPROVE: FIND WAY OF OBTAINING NDF AND PRINT CHI2/NDF /////

  // Minimise
  min->Minimize();
  const double *xs = min->X();
  const double *err = min->Errors(); 

  // Write out
  for(int i = 0  ; i < mysize ; i++){
    cout << Form("S%d=",i+1) << xs[i] << "(" << err[i] << ")" << endl;
  }

  /* Store S value in global variable, to access for drawing on plots */
  globalS = xs[0];
  globalSerr = err[0];

  // Return the Fitted CS
  TGraph* g = new TGraph(); 
  double* X = theory->GetX();
  double* Y = theory->GetY();
  if(loud){
    cout << setprecision(8);
    cout << "Start: X[0] = " << theory->GetPointX(4) << " Y[0] = " << theory->GetPointY(4) << endl;
    cout << "multip by " << xs[0] << endl;
  }
  
  for(int i=0; i<theory->GetN(); i++){ g->SetPoint(g->GetN(),X[i],xs[0]*Y[i]); }

  //for(int i=0; i<theory->GetN(); i++){ g->SetPoint(i,X[i],xs[0]*Y[i]); }

  if(loud){
    cout << "End:   X[0] = " << g->GetPointX(4) << " Y[0] = " << g->GetPointY(4) << endl;
    cout << setprecision(3);
  }

  return g;
}
