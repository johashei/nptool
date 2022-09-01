TFile* ifile;

int NumberOfDetectors=72;
int NumberOfAnodes=1;
int m_NumberOfGammaPeak=1;

double PosGammaPeak = 3.2;

bool Finder(TH1F* h, Double_t *mean, Double_t *sigma);

/////////////////////////////////////////////////////
void OpenRootFile(){
  ifile = new TFile("histo_tof_file_q1_80ns.root");
}

/////////////////////////////////////////////////////
void FitTofGammaPeak(){
  OpenRootFile();

  TCanvas *cLG[11];
  TCanvas *cHG[11];

  ofstream ofile;
  ofile.open("Vendeta_Time.cal");

  for(int i=0; i<11; i++){
    TString canvas_name = Form("LG_Anode%i",i+1);
    cLG[i] = new TCanvas(canvas_name,canvas_name,1800,1800);
    cLG[i]->Divide(10,8);

    canvas_name = Form("HG_Anode%i",i+1);
    cHG[i] = new TCanvas(canvas_name,canvas_name,1800,1800);
    cHG[i]->Divide(10,8);

  }

  Double_t* mean = new Double_t[1];
  Double_t* sigma = new Double_t[1];
		TGraph* gSigma_LG = new TGraph();
		TGraph* gSigma_HG = new TGraph();

  for(int i=0; i<NumberOfDetectors; i++){
    for(int j=0; j<NumberOfAnodes; j++){
						j=5;
      // LG //
      TString histo_name = Form("hLG_Det%i_Anode%i",i+1,j+1);
      TH1F* h = (TH1F*) ifile->FindObjectAny(histo_name);

      cLG[j]->cd(i+1);
      h->Draw();

      mean[0] = 0;
      sigma[0] = 0;
      TString LG_token = Form("Vendeta_DET%i_LG_ANODE%i_TIMEOFFSET",i+1,j+1);
      if(Finder(h, mean, sigma)){
        ofile << LG_token << " " << -mean[0]+PosGammaPeak << endl;
								gSigma_LG->SetPoint(i,i+1,sigma[0]);
      }
      else{
        ofile << LG_token << " 0" << endl;
      }

      // HG //
      histo_name = Form("hHG_Det%i_Anode%i",i+1,j+1);
      h = (TH1F*) ifile->FindObjectAny(histo_name);

      cHG[j]->cd(i+1);
      h->Draw();

      mean[0] = 0;
      sigma[0] = 0; 
      TString HG_token = Form("Vendeta_DET%i_HG_ANODE%i_TIMEOFFSET",i+1,j+1);
      if(Finder(h, mean, sigma)){
        ofile << HG_token << " " << -mean[0]+PosGammaPeak << endl;
								gSigma_HG->SetPoint(i,i+1,sigma[0]);
      }
      else{
        ofile << HG_token << " 0" << endl;
      }
    }
  }

		TCanvas* c1 = new TCanvas("Sigma","Sigma",1200,600);
		c1->Divide(2,1);

		gSigma_LG->SetMarkerStyle(8);
		gSigma_HG->SetMarkerStyle(8);

		c1->cd(1);
		gSigma_LG->Draw();
		c1->cd(2);
		gSigma_HG->Draw();



}


/////////////////////////////////////////////////////
bool Finder(TH1F* h, Double_t *mean, Double_t *sigma){
  Double_t resolsig = 2;
  Double_t resolsigTSpec = 1;
  Double_t seuil = 0.2;

  TSpectrum* s = new TSpectrum(m_NumberOfGammaPeak,resolsigTSpec);

  Int_t nfound = s->Search(h,resolsig,"new",seuil);
  Double_t *xpeaks = s->GetPositionX();

  Double_t linf=0;
  Double_t lsup=0;
  if(nfound == m_NumberOfGammaPeak){
    cout << "Gamma Peak Found" << endl;
    for(int i=0; i<nfound; i++){
      linf = xpeaks[i]-2;
      lsup = xpeaks[i]+0.4;

      TF1* gauss = new TF1("gaus","gaus",linf,lsup);
      h->Fit(gauss,"RQ");
      mean[i] = gauss->GetParameter(1);
      sigma[i] = gauss->GetParameter(2);
    }
    return true;
  }

  if(nfound != m_NumberOfGammaPeak){
    cout << "Warning. Number of peak different of " << m_NumberOfGammaPeak << " !! / nfound = " << nfound << endl; 
    return false;
  }


}
