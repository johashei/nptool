vector<vector<double>> ExpDiffCross(double Energy);
TH1F* PullThetaLabHist(int i, double minTheta, double gatesize);
void Scale(TGraph* g , TGraphErrors* ex);
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j);
double ToMininize(const double* parameter);
TGraph* FindNormalisation(TGraph* theory, TGraphErrors* experiment);

vector<Double_t> anglecentres, anglewidth;
TGraph* currentThry;
TGraphErrors* staticExp;

//BASIC RUN: CS(0.143,2,1,1.5,1)

////////////////////////////////////////////////////////////////////////////////
void CS(double Energy, double Spin, double spdf, double angmom, double nodes){
  // p3/5 -> spdf = 1, angmom = 1.5
  // J0 is incident spin, which is 47K g.s. therefore J0 = 1/2
  double J0 = 0.5;

  /* Solid Angle (from simulation) */
  auto file = new TFile("../SolidAngle_HistFile_06Dec_47Kdp.root");
  TH1F* SolidAngle = (TH1F*) file->FindObjectAny("SolidAngle_Lab_MG");
  TCanvas* c_SolidAngle = new TCanvas("c_SolidAngle","c_SolidAngle",1000,1000);
  SolidAngle->Draw();
 
  /* Area of experimental peaks */
  TCanvas* c_PeakArea = new TCanvas("c_PeakArea","c_PeakArea",1000,1000);
  vector<vector<double>> areaArray = ExpDiffCross(0.143);
  //cout << std::setprecision(3) << std::fixed;
  delete c_PeakArea;

  // Array: i, peakenergy, peakarea, areaerror, anglemin, anglemax
//  for(int i=0; i<areaArray.size();i++){
//    cout << i << " " 
//	 << areaArray[i][0] << " " 
//	 << areaArray[i][1] << " "
//	 << areaArray[i][2] << " "
//	 << areaArray[i][3] << " "
//	 << areaArray[i][4] << endl;
//  }

  /* Area over Solid Angle */
  vector<Double_t> AoSA, AoSAerr;
  for(int i=0; i<areaArray.size();i++){
    int binmin = SolidAngle->FindBin(areaArray[i][3]+0.0001);
    int binmax = SolidAngle->FindBin(areaArray[i][4]-0.0001);

    anglecentres.push_back(((areaArray[i][4]-areaArray[i][3])/2.)+areaArray[i][3]);
    anglewidth.push_back(areaArray[i][4]-areaArray[i][3]);

    double SAerr; // IS THIS IN MSR OR SR????
    double SA = 1000. * SolidAngle->IntegralAndError(binmin,binmax,SAerr);
    //SAerr = SAerr*1000.; //????

    AoSA.push_back(areaArray[i][1] / SA);
    AoSAerr.push_back((areaArray[i][1]/SA) 
		    * sqrt( pow(areaArray[i][2]/areaArray[i][1],2) + pow(SAerr/SA,2) ) );

    cout << "Angle " << areaArray[i][3] << " to " << areaArray[i][4] 
	 << ", centre " << anglecentres[i]
	 << ": SA = " << SA << " +- " << SAerr 
         << "  Area/SA = " << AoSA[i] << " +- " << AoSAerr[i] << " counts/msr"<< endl;
  }
  
  TCanvas* c_AoSA = new TCanvas("c_AoSA","c_AoSA",1000,1000);
  TGraphErrors* gAoSA = new TGraphErrors(
		  anglecentres.size(),
		  &(anglecentres[0]), &(AoSA[0]),
		  //&(anglewidth[0]), &(AoSAerr[0]) );
		  0, &(AoSAerr[0]) );
  gAoSA->SetTitle("Area/SolidAngle");
  gAoSA->GetXaxis()->SetTitle("ThetaLab [deg]");
  gAoSA->GetYaxis()->SetTitle("___ [counts/msr]");
  gAoSA->Draw();

  /* TWOFNR diff. cross section */ 
  TCanvas* c_TWOFNR = new TCanvas("c_TWOFNR","c_TWOFNR",1000,1000);
  TGraph* DiffCross = TWOFNR(Energy, J0, Spin, nodes, spdf, angmom); 
  DiffCross->Draw();

  /* Scaled and compared on same plot */ 
/* 
  cout << "USING BASIC SCALING METHOD..." << endl;
  TGraph* ScaledTWOFNR = DiffCross;
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
  TGraph* Final = FindNormalisation(DiffCross,gAoSA);
  gAoSA->SetLineColor(kRed);
  gAoSA->SetMarkerColor(kRed);
  gAoSA->SetMarkerStyle(21);
  gAoSA->Draw("AP");
  Final->Draw("SAME");
}



////////////////////////////////////////////////////////////////////////////////
vector<vector<double>> ExpDiffCross(double Energy){
  cout << "========================================================" << endl;
  // Implement some kind of energy selection later, for now, hard code that we are selecting 0.143, array index 1
  vector<vector<double>> AllPeaks_OneGate;
  vector<vector<double>> OnePeak_AllGates;
  int numbins = 10;
  double x[numbins], y[numbins];
  for(int i=0; i<numbins;i++){
    double bin = 5.;
    double min = 105. + (i*bin);
    double max = min + bin;
    cout << "===================================" << endl;
    cout << "min: " << min << " max: " << max << endl;

    TH1F* gate = PullThetaLabHist(i,105.,5.);
    AllPeaks_OneGate = FitKnownPeaks_RtrnArry(gate);
    //cout << "area of " << peakAreas[1][0] << " = " << peakAreas[1][1] 
    cout << "area of 0.143 = " << AllPeaks_OneGate[1][1] 
	 << " +- " << AllPeaks_OneGate[1][2] 
	 << endl;
    AllPeaks_OneGate[1].push_back(min);
    AllPeaks_OneGate[1].push_back(max);
    OnePeak_AllGates.push_back(AllPeaks_OneGate[1]);
  }
  return OnePeak_AllGates;
}

////////////////////////////////////////////////////////////////////////////////
TH1F* PullThetaLabHist(int i, double minTheta, double gatesize){
  TFile* file = new TFile("GateThetaLabHistograms.root","READ");
  string histname = "cThetaLabGate_" 
	          + to_string(minTheta+(i*gatesize)) + "-" 
		  + to_string(minTheta+((i+1)*gatesize));
  TList *list = (TList*)file->Get("GateThetaLabHistograms");
  TH1F* hist = (TH1F*)list->FindObject(histname.c_str());
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
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j){

  cout << "========================================================" << endl;
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
  Front_Input << 5 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << J << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 6 << std::endl;
  Front_Input << 4 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1.25 << " " << 0.65 << std::endl;
  Front_Input << 6.0 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << 0 << std::endl;

  Front_Input.close();

  cout << "Filled front input file." << endl;
  cout << "Executing front20..." << endl;
  system("(exec ~/Programs/TWOFNR/front20 < in.front > /dev/null)"); 
  cout << "Executing twofnr20..." << endl;
  system("(exec ~/Programs/TWOFNR/twofnr20 < in.twofnr > /dev/null)");

  cout << "twofnr20 complete." << endl;


  TGraph* CS = new TGraph("24.jjj");

  cout << "===================================" << endl;
  cout << "Current directory    " << twofnrDir << endl;

  cout << "Moving to directory  " << origDir << endl;
  chdir(origDir.c_str());
  system("pwd"); 
  cout << "========================================================" << endl;


  return CS;
}

////////////////////////////////////////////////////////////////////////////////
double Chi2(TGraph* theory , TGraphErrors* exper){
  double Chi2 = 0 ;
  double chi;
  //for(int i = 1 ; i < exper->GetN() ; i++){
  for(int i = 0 ; i < exper->GetN() ; i++){
    if(exper->Eval(anglecentres[i])>0)  {
      chi=(exper->Eval(anglecentres[i])-theory->Eval(anglecentres[i]) ) / (exper->GetErrorY(i));
      Chi2 +=chi*chi;
    }
  }
  cout << "Chi2 = " << Chi2 << endl;
  return Chi2;
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
  // Set global variable
  currentThry = theory;
  staticExp = experiment;

  // Construct minimiser
  const char* minName ="Minuit";const char* algoName="Migrad";
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);

  // Number of parameter
  int mysize = 1; // Originally 2 - why??

  // Create funciton wrapper for minimizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(&ToMininize,mysize);
  min->SetFunction(f);
  
  // Set range of paramenter(??)
  double* parameter = new double[mysize];
  for(unsigned int i = 0 ; i < mysize ; i++){
    parameter[i] = 0.5;
    char name[4];
    sprintf(name,"S%d",i+1);
    min->SetLimitedVariable(i,name,parameter[i],0.1,0,1000);
  }
  
  // Minimise
  min->Minimize();
  const double *xs = min->X();
  const double *err = min->Errors(); 

  // Write out
  for(int i = 0  ; i < mysize ; i++){
    cout << Form("S%d=",i+1) << xs[i] <<"("<<err[i] << ") "; 
  }
  cout << endl;

  // Return the Fitted CS
  TGraph* g = new TGraph(); 
  double* X = theory->GetX();
  double* Y = theory->GetY();
  for(int i=0; i<theory->GetN(); i++){ g->SetPoint(g->GetN(),X[i],xs[0]*Y[i]); }

  return g;
}
