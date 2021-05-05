void CalibFDC(){
  TChain* c = new TChain("RawTree");
  c->Add("~/mrdc/FDC02.root");
  c->SetBranchStatus("*",false);

  // FDC
  int FDC = 0;
  int NLayer = 8;
  c->SetBranchStatus("SamuraiFDC0",true);
  ofstream fout("FDC0_Calib.txt");  
  double inter = 2.5;// inter wire distance
  double low=1720;
  double high=1820;


/*  // FDC2
  int FDC=2;
  int NLayer = 14;
  c->SetBranchStatus("SamuraiFDC2",true);
  ofstream fout("FDC2_Calib.txt");  
  double inter = 10;// inter wire distance
  double low=1200;
  double high=1700;
*/

  string cond;

  TF1* f = new TF1("sigmoid","[0]/(exp([1]*([2]-x))+1)");

 //TF1* f = new TF1("sigmoid","0*(x<[0])+(x-[0])*[1]*(x>[0]&&x<[2])+[3]*(x>[2])");

 TH1D* h = new TH1D("h","h",1500,0,3000);
 for(unsigned int i = 0 ; i < NLayer; i++){
  cond = Form("fFDC%d_LayerNbr==%d&&fFDC%d_Edge==0&&fFDC%d_Time>%f&&fFDC%d_Time<%f&&fFDC%d_Time@.size()>10&&fFDC%d_Time@.size()<100",FDC,i,FDC,FDC,low,FDC,high,FDC,FDC);
  c->Draw(Form("fFDC%d_Time>>h",FDC),cond.c_str());
  TH1D* g = new TH1D(*h);
  unsigned int size = h->GetNbinsX();
  for(unsigned int i = 0 ; i < size ; i++){
    g->SetBinContent(i,h->Integral(0,i)); 
  }
  //new TCanvas();
  //g->Draw();
  cout << g->GetMaximum() << endl;
  f->SetParameter(0,g->GetMaximum());
  f->SetParameter(1,8e-3);
  f->SetParameter(2,1500);
  
  g->Fit(f);
  // Renormalize the distribution to 10 mm
  f->SetParameter(0,inter);
  double p0 = f->GetParameter(0);
  double p1 = f->GetParameter(1);
  double p2 = f->GetParameter(2);
  //new TCanvas();
  //c->Draw(Form("%f/(exp(%f*(%f-(Time+ToT)))+1)>>hh(200,0,2.5)",p0,p1,p2),cond.c_str()); 
  fout << Form("FDC%d_L%d ",FDC,i) << p0 << " " << p1 << " " << p2 << endl;
   } 
   fout.close();
  }
