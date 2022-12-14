/* #include "TGraphErrors.h" */

TGraphErrors* GetT(int ring);
TGraphErrors* GetZ(int ring);
TGraphErrors* Calibrate(const TGraphErrors* g, double offset, double vdrift, double base_line);
void Scale(const TGraphErrors* ref,TGraphErrors* mod);
double chi2(TGraphErrors* g1,TGraphErrors* g2);
double Chi2(const double* parameter);
void NumericalMinimization(const char * minName ,const char *algoName);
TGraphErrors* z;// simulated Z distribution (ref)
TGraphErrors* t;// measured T distribution
TGraphErrors* c;// calibrated Z distribution from t
vector<double> vdrift;
vector<double> evdrift;
vector<double> offset;
vector<double> eoffset;
unsigned int current_ring;
double percent = 0.02;
////////////////////////////////////////////////////////////////////////////////

void cal(){
  TCanvas* cv = new TCanvas();
  cv->Divide(4,5);
  
  ofstream outfile; 
  outfile.open("CalibVDrift.txt"); 
  for(unsigned int i = 0 ; i < 18; i++){
    cout << "Processing Ring " << i << endl;
    cv->cd(i+1);
    current_ring=i+1;
    t = GetT(i+1);
    z = GetZ(i+1);
    NumericalMinimization("Minuit2","Fumili2"); //SCAN, MIGRAD, COMBINED
    z->Draw("ap");

    z->SetLineColor(kRed);
    c->SetLineColor(kBlack);
    c->Draw("p");
    outfile << "Minos_R" << i+1 << "_VDRIFT " << vdrift[i] << endl;
    outfile << "Minos_R" << i+1 << "_OFFSET " << offset[i] << endl;
    }
   outfile.close(); 

  TCanvas* c2 = new TCanvas();
  c2->Divide(2,1);
  c2->cd(1);
  
  vector<double> x = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
  vector<double> ex;// = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  ex.resize(18,0);
  TGraphErrors* VDrift = new TGraphErrors(vdrift.size(),&x[0],&vdrift[0],&ex[0],&evdrift[0]);

  TGraphErrors* Offset= new TGraphErrors(offset.size(),&x[0],&offset[0],&ex[0],&eoffset[0]);
  VDrift->Draw("ap");
  c2->cd(2);
  Offset->Draw("ap");
}


////////////////////////////////////////////////////////////////////////////////
TGraphErrors* GetT(int ring){
  static TFile* Tfile = new TFile("TPad_distrib_new.root"); 
  TString name = Form("TPad_ring%d",ring);
  TH1F* h = (TH1F*) Tfile->FindObjectAny(name);
  h->Rebin(8);
  int maxBin=h->GetMaximumBin();
  double scale = h->Integral(maxBin-10,maxBin+10)/20.;
  h->Scale(1./scale);
  cout << h->GetBinCenter(maxBin) << endl;
  unsigned int size = h->GetNbinsX();
  vector<double> x ;
  vector<double> y ;
  vector<double> ex;
  vector<double> ey;
  for(unsigned int i = 0 ; i < size ; i++){
    x.push_back(h->GetBinCenter(i)); 
    y.push_back(h->GetBinContent(i));
    ex.push_back(0);
    ey.push_back(h->GetBinError(i));
    }

  TGraphErrors* g = new TGraphErrors(size,&x[0],&y[0],&ex[0],&ey[0]);
  return g;
}  

////////////////////////////////////////////////////////////////////////////////
TGraphErrors* GetZ(int ring){
  static TFile* Tfile = new TFile("ZPad_distrib.root"); 
  TString name = Form("ZPad_Ring%d",ring);
  TH1F* h = (TH1F*) Tfile->FindObjectAny(name);
  h->Rebin(8);
  int maxBin=h->GetMaximumBin();
  double scale = h->Integral(maxBin-10,maxBin+10)/20.;
  h->Scale(1./scale);
  unsigned int size = h->GetNbinsX();
  vector<double> x ;
  vector<double> y ;
  vector<double> ex;
  vector<double> ey;
  for(unsigned int i = 0 ; i < size ; i++){
    x.push_back(h->GetBinCenter(i)); 
    y.push_back(h->GetBinContent(i));
    ex.push_back(0);
    ey.push_back(h->GetBinError(i));
  }

  TGraphErrors* g = new TGraphErrors(size,&x[0],&y[0],&ex[0],&ey[0]);
  g->SetMarkerColor(kAzure+7);
  g->SetLineColor(kAzure+7);
  return g;
} 

////////////////////////////////////////////////////////////////////////////////
TGraphErrors* Calibrate(const TGraphErrors* g, double offset, double vdrift,double base_line){
  TGraphErrors* res = new TGraphErrors(*g); 
  unsigned int size = g->GetN();
  double x, y; 
  for(unsigned int i = 0 ; i < size ; i++){
    g->GetPoint(i, x, y);
    res->SetPoint(i,(x*30+offset)*vdrift*(0.98-percent*current_ring/18.),y+base_line);
  }
  res->SetLineColor(kGreen-3);
  res->SetMarkerColor(kGreen-3);
  return res;
}

////////////////////////////////////////////////////////////////////////////////
void Scale(const TGraphErrors* ref,TGraphErrors* mod){

  unsigned int size = ref->GetN();
  double scale=0;
  int point=0;
  double x, y; 
  for(unsigned int i = 0 ; i < size ; i++){
    ref->GetPoint(i,x,y);
    double xref = x;
    double yref = y;
    double ymod = mod->Eval(xref);
    if(yref>0.9){
      scale+= ymod/yref; 
      point++;
    }
  } 

  double val = scale/point;
  size = mod->GetN(); 
  for(unsigned int i = 0 ; i < size ; i++){
    mod->GetPoint(i,x,y);
    mod->SetPoint(i,x,y/val);  
  }

}
////////////////////////////////////////////////////////////////////////////////
double Chi2(const double* parameter){
  c = Calibrate(t,parameter[0],parameter[1],parameter[2]);
  Scale(z,c);
  return chi2(z,c);
} 


////////////////////////////////////////////////////////////////////////////////
double chi2(TGraphErrors* g1,TGraphErrors* g2){

  double chi2 = 0; 
  unsigned int size = g1->GetN();
  double x, y;
  for(unsigned int i = 0 ; i < size ; i++){
    g1->GetPoint(i, x, y);
    double xg1 = x;
    double yg1 = y;
    double yg2 = g2->Eval(xg1);
    if(yg2 && yg1)
      chi2+=(yg1-yg2)*(yg1-yg2)/(g1->GetErrorY(i));
  } 
  return chi2; 

}
////////////////////////////////////////////////////////////////////////////////
void NumericalMinimization(const char * minName ,const char *algoName){
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);

  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  /* min->SetPrecision(1e-4); */
  //min->SetTolerance(0.001);
  min->SetValidError(true);
  //  min->SetPrintLevel(1);

  // create funciton wrapper for minimizer
  // a IMultiGenFunction type 
  ROOT::Math::Functor f(&Chi2,3); 

  min->SetFunction(f);
  min->Clear();
  double* parameter = new double[3];
  parameter[0] = -1000;
  //parameter[1] =0.0291536; //p5mm
 // parameter[1] =0.030712; //m5mm
//  parameter[1] =0.0304835; // d5mm

//  parameter[1] =0.0306181; // d5mm

//  parameter[1] =0.0314412; // speed of ring 18

  parameter[1] = 0.0345; // speed for right target length

//  parameter[1] =0.0291703; // speed of ring 1

  parameter[2] =0;// base line
//  parameter[1] =0.0300424 ;
  //parameter[1] =0.029 ;

  
  // Set the free variables to be minimized!
  min->SetLimitedVariable(0,"Offset",parameter[0],10,-2000,2000);
  //min->SetLimitedVariable(1,"VDrift",parameter[1],1e-2,0.02,0.05); 
  min->SetFixedVariable(1,"VDrift",parameter[1]);
  //min->SetFixedVariable(2,"BaseLine",parameter[2]); 
  min->SetLimitedVariable(2,"BaseLine",parameter[2],1,0,20); 
  // Set a fixed VDrift

  // do the minimization
  min->Minimize(); 

  const double *xs = min->X();
  std::cout << "Offset = " << xs[0] <<  " +/- " << min->Errors()[0] << std::endl;
  std::cout << "VDrift = " << xs[1] <<  " +/- " << min->Errors()[1] << std::endl;
  std::cout << "Chi2= " << Chi2(xs) << std::endl;

  vdrift.push_back(xs[1]*(0.98+percent*current_ring/18.));
  evdrift.push_back(min->Errors()[1]);
  offset.push_back(xs[0]);
  eoffset.push_back(min->Errors()[0]);
  return 0;
}
