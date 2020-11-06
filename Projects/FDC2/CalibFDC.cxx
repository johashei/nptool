void CalibFDC(){
  TChain* c = new TChain("RawTree");
  c->Add("FDC2.root");

 string cond;
 ofstream fout("FDC2_Calib.txt");  
 TF1* f = new TF1("sigmoid","[0]/(exp([1]*([2]-x))+1)");

 //TF1* f = new TF1("sigmoid","0*(x<[0])+(x-[0])*[1]*(x>[0]&&x<[2])+[3]*(x>[2])");

 TH1D* h = new TH1D("h","h",1500,0,3000);
 for(unsigned int i = 0 ; i < 14 ; i++){
  cond = Form("fDC_LayerNbr==%d&&fDC_Edge==0&&fDC_Time>1200&&fDC_Time<1700&&fDC_Time@.size()>10&&fDC_Time@.size()<100",i);
  c->Draw("fDC_Time>>h",cond.c_str());
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
/*
  f->SetParameter(0,1200);
  f->SetParameter(1,1);
  f->SetParameter(2,1700);  
  f->SetParameter(3,g->GetMaximum());
*/
  g->Fit(f);
  // Renormalize the distribution to 10 mm
  f->SetParameter(0,10);
  double p0 = f->GetParameter(0);
  double p1 = f->GetParameter(1);
  double p2 = f->GetParameter(2);
  //new TCanvas();
  //c->Draw(Form("%f/(exp(%f*(%f-(Time+ToT)))+1)>>hh(200,0,2.5)",p0,p1,p2),cond.c_str()); 
  fout << Form("FDC2_L%d ",i) << p0 << " " << p1 << " " << p2 << endl;
   } 
   fout.close();
  }
