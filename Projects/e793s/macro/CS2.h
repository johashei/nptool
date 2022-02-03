//TGraph* ExpDiffCross(double Energy);
double* ExpDiffCross(double Energy);
TH1F* GateThetaLab_RtrnHist(double minTheta, double maxTheta, double binsize);
TH1F* PullThetaLabHist(int i, double minTheta, double gatesize);
void Scale(TGraph* g , TGraphErrors* ex);
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j);
double ToMininize(const double* parameter);
//TGraph* FindNormalisation(TGraph* cs1, TH1* hexp);
TGraph* FindNormalisation(TGraph* cs1, double* hexp);

TGraph* g1;
//TH1* current;
double* current;
TH1F* SolidAngle;
//auto diffcrossexp = new TGraph("diffcrossexp","diffcrossexp",180,0,180);

////////////////////////////////////////////////////////////////////////////////
void CS(double Energy, double Spin, double spdf, double angmom, double nodes){
  auto file = new TFile("../SolidAngle_HistFile_06Dec_47Kdp.root");
  SolidAngle = (TH1F*) file->FindObjectAny("SolidAngle_Lab_MG");

  // p3/5 -> spdf = 1, angmom = 1.5

  // J0 is incident spin, which is 47K g.s. therefore J0 = 1/2
  double J0 = 0.5;

  // Solid Angle
  TCanvas* c_SolidAngle = new TCanvas("c_SolidAngle","c_SolidAngle",1000,1000);
  SolidAngle->Draw();
 
  // 
  TCanvas* c_TWOFNR = new TCanvas("c_TWOFNR","c_TWOFNR",1000,1000);

  // some work around here.
  // TWOFNR is calling the ADWA calculation
  auto diffcross = TWOFNR(Energy, J0, Spin, nodes, spdf, angmom); 
  diffcross->Draw();
  // need to define hexp, best to have to have some function that take Ex and Edc 
  // min and max, make the angular distrib, divide by solid angle and return the 
  // normalised CS.
  //
  //

  //TCanvas* c_ExpDiffCross = new TCanvas("c_ExpDiffCross","c_ExpDiffCross",1000,1000);
  auto temp = ExpDiffCross(0.143);
  //TCanvas* c_ExpDiffCross = new TCanvas("c_ExpDiffCross","c_ExpDiffCross",1000,1000);
  //temp->Draw();  

//  TCanvas* c_Final = new TCanvas("c_Final","c_Final",1000,1000);
//  temp->Divide(SolidAngle);
//  temp->Draw();

  //GateThetaLab(105,110,0.1);FitKnownPeaks(Ex_ThetaLabGate)
  //
  FindNormalisation(diffcross,temp);
}

//TGraph* ExpDiffCross(double Energy){
double* ExpDiffCross(double Energy){
  // Implement some kind of energy selection later, for now, hard code that we are selecting 0.143, array index 1
  int numbins = 10;
  double x[numbins], y[numbins];
  //for(int i=0; i<10;i++){
  for(int i=0; i<numbins;i++){
    double* peakAreas;
    double bin = 5.;
    double min = 105. + (i*bin);
    double max = min + bin;
    cout << "min: " << min << " max: " << max << endl;

    TH1F* gate = PullThetaLabHist(i,105.,5.);
    peakAreas = FitKnownPeaks_RtrnArry(gate);
    cout << "area of 0.143 = " << peakAreas[1] << endl;
//    delete gate;
    //diffcrossexp->Fill(min+(bin/2.),peakAreas[1]);    

//    x[i] = min+(bin/2.);
    y[i] = peakAreas[1];
  }

//  TGraph* diffcrossexp = new TGraph(numbins,x,y);

//  return diffcrossexp;
  return y;
}

/*
TH1F* GateThetaLab_RtrnHist(double minTheta, double maxTheta, double binsize){
  string gating = "abs(T_MUGAST_VAMOS-2777)<600 && Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8 && ThetaLab > " 
      + to_string(minTheta)
      + " && ThetaLab < "
      + to_string(maxTheta);

  string title = to_string(minTheta)+" < ThetaLab < "+to_string(maxTheta);
  string ytitle = "Counts / " + to_string(binsize) + " MeV";
  string draw = "Ex>>Ex_ThetaLabGate(" + to_string(8.0/binsize) + ",-1,7)";

  TCanvas *cEx_ThetaLabGate = new TCanvas("cEx_ThetaLabGate","cEx_ThetaLabGate",1000,1000);
  chain->Draw(draw.c_str(),gating.c_str(),"colz");
  TH1F* Ex_ThetaLabGate = (TH1F*) gDirectory->Get("Ex_ThetaLabGate");
  Ex_ThetaLabGate->GetXaxis()->SetTitle("Ex [MeV]");
  Ex_ThetaLabGate->GetYaxis()->SetTitle(ytitle.c_str());
  Ex_ThetaLabGate->SetTitle(title.c_str());

  return Ex_ThetaLabGate;
}
*/

TH1F* PullThetaLabHist(int i, double minTheta, double gatesize){
  TFile* file = new TFile("GateThetaLabHistograms.root","READ");

  string histname = "cThetaLabGate_" + to_string(minTheta+(i*gatesize)) + "-" + to_string(minTheta+((i+1)*gatesize));
  TList *list = (TList*)file->Get("GateThetaLabHistograms");
  TH1F* hist = (TH1F*)list->FindObject(histname.c_str());
  return hist;
}


////////////////////////////////////////////////////////////////////////////////
/*
 double Chi2(TGraph* g , TH1* h){
 double Chi2 = 0 ;
 double chi;

  for(int i = 1 ; i < h->GetNbinsX() ; i++){
    if(h->GetBinContent(i)>0)  {
      chi = (h->GetBinContent(i) - g->Eval(h->GetBinCenter(i)) ) / (h->GetBinError(i));
      Chi2 +=chi*chi;
    }
  }

 return Chi2;
}
*/

double Chi2(TGraph* g , double* h){
 double Chi2 = 0 ;
 double chi;

  for(int i = 1 ; i < 10 ; i++){
    if(h[i]>0)  {
      chi = (h[i] - g->Eval(h[i]) ); // DIVIDE BY ERROR!!!   /(h->GetBinError(i));
      Chi2 +=chi*chi;
    }
  }

 return Chi2;
}

////////////////////////////////////////////////////////////////////////////////
double ToMininize(const double* parameter){
static int f = 0 ;
  TGraph* g = new TGraph();
  double* X = g1->GetX();
  double* Y = g1->GetY();
  for(int i = 0 ; i < g1->GetN() ; i++){
    g->SetPoint(g->GetN(),X[i],parameter[0]*Y[i]);
  }
  double chi2  = Chi2(g,current);  
  return chi2;
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

  double* x = g->GetX();
  double* y = g->GetY();
  for(unsigned int i = 0 ; i < g->GetN() ; i++)
    g->SetPoint(i,x[i],y[i]*mean);
}
////////////////////////////////////////////////////////////////////////////////
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j){
  cout << " ==================================== " << endl;
  char origDirchar[200];
  getcwd(origDirchar,200);
  string origDir{origDirchar};
  string twofnrDir = "/home/charlottepaxman/Programs/TWOFNR";
  cout << "Current directory    " << origDir << endl;

  cout << "Moving to directory  " << twofnrDir << endl;
  chdir(twofnrDir.c_str());
  //Check
  system("pwd"); 
  cout << " ==================================== " << endl;

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

  cout << " ==================================== " << endl;
  cout << "Current directory    " << twofnrDir << endl;

  cout << "Moving to directory  " << origDir << endl;
//  string line = "cd ";
//  line += twofnrDir;
//  system(line.c_str()); 
  chdir(origDir.c_str());
  //Check
  system("pwd"); 
  cout << " ==================================== " << endl;


  return CS;
  //return 0;
}

////////////////////////////////////////////////////////////////////////////////
//TGraph* FindNormalisation(TGraph* cs1, TH1* hexp){
TGraph* FindNormalisation(TGraph* cs1, double* hexp){
  // Set global variable
  g1 = cs1;
  current = hexp;

  const char* minName ="Minuit";
  const char* algoName="Migrad";
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);

  // Number of parameter
  int mysize = 2;

  // create funciton wrapper for minimizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(&ToMininize,mysize);
  min->SetFunction(f);
  
  double* parameter = new double[mysize];
  for(unsigned int i = 0 ; i < mysize ; i++){
    parameter[i] = 0.5;
    char name[4];
    sprintf(name,"S%d",i+1);
    min->SetLimitedVariable(i,name,parameter[i],0.1,0,1000);
  }
  
  // do the minimization
  min->Minimize();
  const double *xs = min->X();
  const double *err = min->Errors(); 

  for(int i = 0  ; i < mysize ; i++){
    cout << Form("S%d=",i+1) << xs[i] <<"("<<err[i] << ") "; 
  }
  cout << endl;
  // Return the Fitted CS
  TGraph* g = new TGraph(); 
  double* X = cs1->GetX();
  double* Y = cs1->GetY();
  for(int i = 0 ; i < cs1->GetN() ; i++){
    g->SetPoint(g->GetN(),X[i],xs[0]*Y[i]);
  }
  return g;
}

