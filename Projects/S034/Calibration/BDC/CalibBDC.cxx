void CalibBDC(){
  
  for(int run12 =824; run12 <= 824 ; run12++){
  
  TChain* c = new TChain("RawTree");
  c->Add(Form("../../root/mrdc/run%d/run%d_SAMURAIBDC.root",run12,run12));
  c->SetBranchStatus("*",false);

  // BDC2
  int BDC=1;
  int NLayer = 8;
  c->SetBranchStatus("SamuraiBDC",true);
  ofstream fout(Form("BDC%d_Calib_%d.txt",BDC, run12));  
  double inter = 2.5;// inter wire distance / 2
  double low=1600;
  double high=1900;

  cout <<"Calibrating BDC" << BDC <<  " run number : " << run12 << endl;
  
  string cond;

  TF1* f = new TF1("sigmoid","[0]/(exp([1]*([2]-x))+1)");
 //TF1* f = new TF1("sigmoid","0*(x<[0])+(x-[0])*[1]*(x>[0]&&x<[2])+[3]*(x>[2])");

 TH1D* h = new TH1D("h","h",1500,0,3000);
 for(unsigned int i = 0 ; i < NLayer; i++){
  cond = Form("fBDC_LayerNbr==%d&&fBDC_Edge==0&&fBDC_Time>%f&&fBDC_Time<%f&&fBDC_Time@.size()>10&&fBDC_Time@.size()<100 &&fBDC_DetectorNbr==%d ",i,low,high,BDC);
  c->Draw("fBDC_Time>>h",cond.c_str());
  TH1D* g = new TH1D(*h);
  unsigned int size = h->GetNbinsX();
  for(unsigned int i = 0 ; i < size ; i++){
    g->SetBinContent(i,h->Integral(0,i)); 
  }

  //new TCanvas();
  //g->Draw();
  cout << g->GetMaximum() << endl;
  f->SetParameter(0,g->GetMaximum());
  f->SetParameter(1,0.1);
  f->SetParameter(2,1700);
  f->SetParLimits(1,0,1);
  f->SetParLimits(2,1600,1900);

  g->Fit(f);
  // Renormalize the distribution to the maximum driftlength 
  f->SetParameter(0,inter);
  double p0 = f->GetParameter(0);
  double p1 = f->GetParameter(1);
  double p2 = f->GetParameter(2);
  //new TCanvas();
  //c->Draw(Form("%f/(exp(%f*(%f-(Time+ToT)))+1)>>hh(200,0,2.5)",p0,p1,p2),cond.c_str()); 
  fout << Form("BDC%d_L%d ",BDC,i) << p0 << " " << p1 << " " << p2 << endl;
//  gPad->Update();
//  gPad->WaitPrimitive();
   } 
   fout.close();

/*   ////////////////////////////////////////// */ 
/*   // BDC2 */
/*   ////////////////////////////////////////// */
  
/*   BDC=2; */
/*   NLayer = 14; */
/*   c->SetBranchStatus("SamuraiBDC2",true); */
/*   ofstream fout2(Form("CalibFiles/BDC2_Calib_%d.txt", run12)); */  
/*   inter = 10;// inter wire distance */
/*   low=1200; */
/*   high=1700; */
/*   h->Clear(); */
  
/*   for(unsigned int i = 0 ; i < NLayer; i++){ */
/*     cond = Form("fBDC%d_LayerNbr==%d&&fBDC%d_Edge==0&&fBDC%d_Time>%f&&fBDC%d_Time<%f&&fBDC%d_Time@.size()>10&&fBDC%d_Time@.size()<100",BDC,i,BDC,BDC,low,BDC,high,BDC,BDC); */
/*     c->Draw(Form("fBDC%d_Time>>h",BDC),cond.c_str()); */
/*     TH1D* g = new TH1D(*h); */
/*     unsigned int size = h->GetNbinsX(); */
/*     for(unsigned int i = 0 ; i < size ; i++){ */
/*       g->SetBinContent(i,h->Integral(0,i)); */ 
/*     } */
/*     //new TCanvas(); */
/*     //g->Draw(); */
/*     cout << g->GetMaximum() << endl; */
/*     f->SetParameter(0,g->GetMaximum()); */
/*     f->SetParameter(1,8e-3); */
/*     f->SetParameter(2,1500); */

/*     g->Fit(f); */
/*     // Renormalize the distribution to 10 mm */
/*     f->SetParameter(0,inter); */
/*     double p0 = f->GetParameter(0); */
/*     double p1 = f->GetParameter(1); */
/*     double p2 = f->GetParameter(2); */
/*     //new TCanvas(); */
/*     //c->Draw(Form("%f/(exp(%f*(%f-(Time+ToT)))+1)>>hh(200,0,2.5)",p0,p1,p2),cond.c_str()); */ 
/*     fout2 << Form("BDC%d_L%d ",BDC,i) << p0 << " " << p1 << " " << p2 << endl; */
/*   } */ 
/*   fout2.close(); */
  
  }
}
