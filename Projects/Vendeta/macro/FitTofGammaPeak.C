TFile* ifile;

int NumberOfDetectors=72;
int NumberOfAnodes=11;
int m_NumberOfGammaPeak=1;
int run=86;
double PosGammaPeak = 3.2;

double sigma_anode[11];

bool Finder(TH1F* h, Double_t *mean, Double_t *sigma);

/////////////////////////////////////////////////////
void OpenRootFile(){
  ifile = new TFile(Form("histo_tof_file_run%i.root",run));
}

/////////////////////////////////////////////////////
void FitTofGammaPeak(){
  OpenRootFile();

  TCanvas *cLG[11];
  TCanvas *cHG[11];

  ofstream ofile;
  ofile.open(Form("Vendeta_Time_run%i.cal",run));

		TFile* orootfile = new TFile(Form("histo_tof_fitted_run%i.root",run),"recreate");

  for(int i=0; i<NumberOfAnodes; i++){
			sigma_anode[i]=0;
  }

  Double_t* mean = new Double_t[1];
  Double_t* sigma = new Double_t[1];
		TGraph* gSigma_LG = new TGraph();
		TGraph* gSigma_HG = new TGraph();

		int countLG=0;
		int countHG=0;
  for(int i=0; i<NumberOfDetectors; i++){
    for(int j=0; j<NumberOfAnodes; j++){
      // LG //
      TString histo_name = Form("hLG_Det%i_Anode%i",i+1,j+1);
      TH1F* h = (TH1F*) ifile->FindObjectAny(histo_name);

      //cLG[j]->cd(i+1);
      h->Draw();

      mean[0] = 0;
      sigma[0] = 0;
      TString LG_token = Form("Vendeta_DET%i_LG_ANODE%i_TIMEOFFSET",i+1,j+1);
      if(Finder(h, mean, sigma)){
        ofile << LG_token << " " << -mean[0]+PosGammaPeak << endl;
								gSigma_LG->SetPoint(countLG,countLG+1,sigma[0]);

								sigma_anode[j] += sigma[0];

								countLG++;
								h->Write();
      }
      else{
        ofile << LG_token << " 0" << endl;
      }

      // HG //
      histo_name = Form("hHG_Det%i_Anode%i",i+1,j+1);
      h = (TH1F*) ifile->FindObjectAny(histo_name);

      //cHG[j]->cd(i+1);
      h->Draw();

      mean[0] = 0;
      sigma[0] = 0; 
      TString HG_token = Form("Vendeta_DET%i_HG_ANODE%i_TIMEOFFSET",i+1,j+1);
      if(Finder(h, mean, sigma)){
        ofile << HG_token << " " << -mean[0]+PosGammaPeak << endl;
								gSigma_HG->SetPoint(countHG,countHG+1,sigma[0]);
								countHG++;
								h->Write();
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

		gSigma_LG->SetName("sigma_LG");
		gSigma_HG->SetName("sigma_HG");
		
		gSigma_LG->Write();
		gSigma_HG->Write();

		orootfile->Write();
		orootfile->Close();

		for(int i=0; i<NumberOfAnodes; i++){
			cout << "Anode= " << i+1 << " / " << sigma_anode[i]/NumberOfDetectors << endl; 
		}

}


/////////////////////////////////////////////////////
bool Finder(TH1F* h, Double_t *mean, Double_t *sigma){
  Double_t resolsig = 2;
  Double_t resolsigTSpec = 1;
  Double_t seuil = 0.3;

  TSpectrum* s = new TSpectrum(m_NumberOfGammaPeak,resolsigTSpec);

  Int_t nfound = s->Search(h,resolsig,"new",seuil);
  Double_t *xpeaks = s->GetPositionX();

  Double_t linf=0;
  Double_t lsup=0;
  if(nfound == m_NumberOfGammaPeak){
    cout << "Gamma Peak Found" << endl;
    for(int i=0; i<nfound; i++){
      linf = xpeaks[i]-2;
      lsup = xpeaks[i]+1;

      TF1* gauss = new TF1("gaus","gaus",linf,lsup);
      h->Fit(gauss,"RQ");
      mean[i] = gauss->GetParameter(1);
      sigma[i] = gauss->GetParameter(2);
    }
    return true;
  }

		else if(nfound != m_NumberOfGammaPeak){
    cout << "Warning. Number of peak different of " << m_NumberOfGammaPeak << " !! / nfound = " << nfound << endl; 
    return false;
  }
}
