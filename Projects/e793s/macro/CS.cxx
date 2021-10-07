void Scale(TGraph* g , TGraphErrors* ex);
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j);
double ToMininize(const double* parameter);
TGraph* FindNormalisation(TGraph* cs1, TH1* hexp);

TGraph* g1;
TH1* current;

////////////////////////////////////////////////////////////////////////////////
void CS(){
  file = new TFile("Efficiency.root");
  TH1* Eff = (TH1*) file->FindObjectAny("SolidAngle");

  new TCanvas();
  h->Draw("colz");
  // CS
  TCanvas* c = new TCanvas("CS","CS",1000,1000);
  // some work around here.
  // TWOFNR is calling the ADWA calculation
  auto g = TWOFNR(double E, double J0, double J, double n, double l, double j); 
  // need to define hexp, best to have to have some function that take Ex and Edc 
  // min and max, make the angular distrib, divide by solid angle and return the 
  // normalised CS.
  FindNormation(g,hexp);
}

////////////////////////////////////////////////////////////////////////////////
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
  double BeamEnergy =  8;

  NPL::Reaction r("28Mg(d,p)29Mg@224");
  r.SetExcitationHeavy(E);
  double QValue = r.GetQValue();

  std::ofstream Front_Input("in.front");
  Front_Input << "jjj" << std::endl;
  Front_Input << "pipo" << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << BeamEnergy << std::endl;
  Front_Input << 28 << " " << 12 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << "0 0 0" << std::endl;
  Front_Input << l << " " << j << std::endl;
  Front_Input << n << std::endl;
  Front_Input << 2 << std::endl;

  // unbound case:
//  if( QValue < 0 )
//    Front_Input << 0.01 << std::endl;
//  else
    Front_Input << QValue << std::endl; 

  Front_Input << 1 << std::endl;
  Front_Input << J0 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 5 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << J << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1.25 << " " << 0.65 << std::endl;
  Front_Input << 6.0 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << 0 << std::endl;

  Front_Input.close() ;

  system("(exec FRONT < in.front &> /dev/null)"); 
  system("(exec echo tran.jjj | TWOFNR &> /dev/null)");
 // system("exec FRONT < in.front"); 
 // system("(exec echo tran.jjj | TWOFNR)");


  TGraph* CS = new TGraph("24.jjj");
  return CS;
}

////////////////////////////////////////////////////////////////////////////////
TGraph* FindNormalisation(TGraph* cs1, TH1* hexp){
  // Set global variable
  g1 = cs1;
  current = hexp;

  const char* minName ="Minuit";const char* algoName="Migrad";
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);

  // Number of parameter
  mysize = 2;

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

