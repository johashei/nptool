const int MUGAST_TELESCOPES = 6; /* Mugast: 1, 2, 3, 4, 5, 7 */
const int MUST2_TELESCOPES = 5; 

void DrawTimeRawSpectra(const char* input_file){

  TFile f(input_file);
  TDirectory *dir = (TDirectory*)(f.Get("ControlSpectra"));
  dir->cd();
  TDirectoryFile *dir_mugast = (TDirectoryFile*)(dir->Get("Mugast"));
  dir_mugast->cd();
  TH2F *h_mugast[MUGAST_TELESCOPES*2];
  h_mugast[0] = (TH2F*) dir_mugast->Get("MG1_STRX_T_RAW");
  h_mugast[1] = (TH2F*) dir_mugast->Get("MG2_STRX_T_RAW");
  h_mugast[2] = (TH2F*) dir_mugast->Get("MG3_STRX_T_RAW");
  h_mugast[3] = (TH2F*) dir_mugast->Get("MG4_STRX_T_RAW");
  h_mugast[4] = (TH2F*) dir_mugast->Get("MG5_STRX_T_RAW");
//  h_mugast[5] = (TH2F*) dir_mugast->Get("MG6_STRX_T_RAW");
  h_mugast[5] = (TH2F*) dir_mugast->Get("MG7_STRX_T_RAW");
//  h_mugast[7] = (TH2F*) dir_mugast->Get("MG9_STRX_T_RAW");
//  h_mugast[8] = (TH2F*) dir_mugast->Get("MG11_STRX_T_RAW");
  h_mugast[6] = (TH2F*) dir_mugast->Get("MG1_STRY_T_RAW");
  h_mugast[7] = (TH2F*) dir_mugast->Get("MG2_STRY_T_RAW");
  h_mugast[8] = (TH2F*) dir_mugast->Get("MG3_STRY_T_RAW");
  h_mugast[9] = (TH2F*) dir_mugast->Get("MG4_STRY_T_RAW");
  h_mugast[10] = (TH2F*) dir_mugast->Get("MG5_STRY_T_RAW");
//  h_mugast[14] = (TH2F*) dir_mugast->Get("MG6_STRY_T_RAW");
  h_mugast[11] = (TH2F*) dir_mugast->Get("MG7_STRY_T_RAW");
//  h_mugast[16] = (TH2F*) dir_mugast->Get("MG9_STRY_T_RAW");
//  h_mugast[17] = (TH2F*) dir_mugast->Get("MG11_STRY_T_RAW");

  TCanvas *c0 = new TCanvas("c0", "Raw time spectra of Mugast (X e Y)", 1200, 900);
  c0->Divide(4,3);
  for(int i=0; i<MUGAST_TELESCOPES*2; i++){
    h_mugast[i]->SetDirectory(0);
    c0->cd(i+1); 
    h_mugast[i]->Draw("colz");
  }
  dir_mugast->Close();
  dir->cd();

  TDirectoryFile *dir_must2 = (TDirectoryFile*)(dir->Get("Must2"));
  dir_must2->cd();
  TH2F *h_must[MUST2_TELESCOPES*2];
  h_must[0] = (TH2F*) dir_must2->Get("MM1_STRX_T_RAW");
  h_must[1] = (TH2F*) dir_must2->Get("MM2_STRX_T_RAW");
  h_must[2] = (TH2F*) dir_must2->Get("MM3_STRX_T_RAW");
  h_must[3] = (TH2F*) dir_must2->Get("MM4_STRX_T_RAW");
  h_must[4] = (TH2F*) dir_must2->Get("MM5_STRX_T_RAW");
  h_must[5] = (TH2F*) dir_must2->Get("MM1_STRY_T_RAW");
  h_must[6] = (TH2F*) dir_must2->Get("MM2_STRY_T_RAW");
  h_must[7] = (TH2F*) dir_must2->Get("MM3_STRY_T_RAW");
  h_must[8] = (TH2F*) dir_must2->Get("MM4_STRY_T_RAW");
  h_must[9] = (TH2F*) dir_must2->Get("MM5_STRY_T_RAW");

  TCanvas *c1 = new TCanvas("c1", "Raw time spectra of MUST2 (X e Y)", 1200, 600);
  c1->Divide(5,2);
  for(int i=0; i<MUST2_TELESCOPES*2; i++){
    h_must[i]->SetDirectory(0);
    c1->cd(i+1); 
    h_must[i]->Draw("colz");
  }
  dir_must2->Close();
}

void DrawTimeCalSpectra(const char* input_file){

  TFile f(input_file);
  TDirectory *dir = (TDirectory*)(f.Get("ControlSpectra"));
  dir->cd();
  TDirectoryFile *dir_mugast = (TDirectoryFile*)(dir->Get("Mugast"));
  dir_mugast->cd();
  TH2F *h_mugast[MUGAST_TELESCOPES*2];
  h_mugast[0] = (TH2F*) dir_mugast->Get("MG1_STRX_T_CAL");
  h_mugast[1] = (TH2F*) dir_mugast->Get("MG2_STRX_T_CAL");
  h_mugast[2] = (TH2F*) dir_mugast->Get("MG3_STRX_T_CAL");
  h_mugast[3] = (TH2F*) dir_mugast->Get("MG4_STRX_T_CAL");
  h_mugast[4] = (TH2F*) dir_mugast->Get("MG5_STRX_T_CAL");
//  h_mugast[5] = (TH2F*) dir_mugast->Get("MG6_STRX_T_CAL");
  h_mugast[5] = (TH2F*) dir_mugast->Get("MG7_STRX_T_CAL");
//  h_mugast[7] = (TH2F*) dir_mugast->Get("MG9_STRX_T_CAL");
//  h_mugast[8] = (TH2F*) dir_mugast->Get("MG11_STRX_T_CAL");
  h_mugast[6] = (TH2F*) dir_mugast->Get("MG1_STRY_T_CAL");
  h_mugast[7] = (TH2F*) dir_mugast->Get("MG2_STRY_T_CAL");
  h_mugast[8] = (TH2F*) dir_mugast->Get("MG3_STRY_T_CAL");
  h_mugast[9] = (TH2F*) dir_mugast->Get("MG4_STRY_T_CAL");
  h_mugast[10] = (TH2F*) dir_mugast->Get("MG5_STRY_T_CAL");
//  h_mugast[14] = (TH2F*) dir_mugast->Get("MG6_STRY_T_CAL");
  h_mugast[11] = (TH2F*) dir_mugast->Get("MG7_STRY_T_CAL");
//  h_mugast[16] = (TH2F*) dir_mugast->Get("MG9_STRY_T_CAL");
//  h_mugast[17] = (TH2F*) dir_mugast->Get("MG11_STRY_T_CAL");

  TCanvas *c0 = new TCanvas("c0", "Calibrated time spectra of Mugast (X e Y)", 1200, 900);
  c0->Divide(4,3);
  for(int i=0; i<MUGAST_TELESCOPES*2; i++){
    h_mugast[i]->SetDirectory(0);
    c0->cd(i+1); 
    h_mugast[i]->Draw("colz");
  }
  dir_mugast->Close();
  dir->cd();

  TDirectoryFile *dir_must2 = (TDirectoryFile*)(dir->Get("Must2"));
  dir_must2->cd();
  TH2F *h_must[MUST2_TELESCOPES*2];
  h_must[0] = (TH2F*) dir_must2->Get("MM1_STRX_T_CAL");
  h_must[1] = (TH2F*) dir_must2->Get("MM2_STRX_T_CAL");
  h_must[2] = (TH2F*) dir_must2->Get("MM3_STRX_T_CAL");
  h_must[3] = (TH2F*) dir_must2->Get("MM4_STRX_T_CAL");
  h_must[4] = (TH2F*) dir_must2->Get("MM5_STRX_T_CAL");
  h_must[5] = (TH2F*) dir_must2->Get("MM1_STRY_T_CAL");
  h_must[6] = (TH2F*) dir_must2->Get("MM2_STRY_T_CAL");
  h_must[7] = (TH2F*) dir_must2->Get("MM3_STRY_T_CAL");
  h_must[8] = (TH2F*) dir_must2->Get("MM4_STRY_T_CAL");
  h_must[9] = (TH2F*) dir_must2->Get("MM5_STRY_T_CAL");

  TCanvas *c1 = new TCanvas("c1", "Calibrated time spectra of MUST2 (X e Y)", 1200, 600);
  c1->Divide(5,2);
  for(int i=0; i<MUST2_TELESCOPES*2; i++){
    h_must[i]->SetDirectory(0);
    c1->cd(i+1); 
    h_must[i]->Draw("colz");
  }
  dir_must2->Close();
}

void DrawEnergyRawSpectra(const char* input_file){

  TFile f(input_file);
  TDirectory *dir = (TDirectory*)(f.Get("ControlSpectra"));
  dir->cd();
  TDirectoryFile *dir_mugast = (TDirectoryFile*)(dir->Get("Mugast"));
  dir_mugast->cd();
  int binning = 200;
  TH2F *h_mugast[MUGAST_TELESCOPES*2];
  h_mugast[0] = (TH2F*) dir_mugast->Get("MG1_STRX_E_RAW");
  h_mugast[1] = (TH2F*) dir_mugast->Get("MG2_STRX_E_RAW");
  h_mugast[2] = (TH2F*) dir_mugast->Get("MG3_STRX_E_RAW");
  h_mugast[3] = (TH2F*) dir_mugast->Get("MG4_STRX_E_RAW");
  h_mugast[4] = (TH2F*) dir_mugast->Get("MG5_STRX_E_RAW");
//  h_mugast[5] = (TH2F*) dir_mugast->Get("MG6_STRX_E_RAW");
  h_mugast[5] = (TH2F*) dir_mugast->Get("MG7_STRX_E_RAW");
//  h_mugast[7] = (TH2F*) dir_mugast->Get("MG9_STRX_E_RAW");
//  h_mugast[8] = (TH2F*) dir_mugast->Get("MG11_STRX_E_RAW");
  h_mugast[6] = (TH2F*) dir_mugast->Get("MG1_STRY_E_RAW");
  h_mugast[7] = (TH2F*) dir_mugast->Get("MG2_STRY_E_RAW");
  h_mugast[8] = (TH2F*) dir_mugast->Get("MG3_STRY_E_RAW");
  h_mugast[9] = (TH2F*) dir_mugast->Get("MG4_STRY_E_RAW");
  h_mugast[10] = (TH2F*) dir_mugast->Get("MG5_STRY_E_RAW");
//  h_mugast[14] = (TH2F*) dir_mugast->Get("MG6_STRY_E_RAW");
  h_mugast[11] = (TH2F*) dir_mugast->Get("MG7_STRY_E_RAW");
//  h_mugast[16] = (TH2F*) dir_mugast->Get("MG9_STRY_E_RAW");
//  h_mugast[17] = (TH2F*) dir_mugast->Get("MG11_STRY_E_RAW");

  TCanvas *c0 = new TCanvas("c0", "Raw energy spectra of Mugast (X e Y)", 1200, 600);
  c0->Divide(4,3);
  for(int i=0; i<MUGAST_TELESCOPES*2; i++){
    h_mugast[i]->SetDirectory(0);
    c0->cd(i+1);
    if(i<MUGAST_TELESCOPES)
      h_mugast[i]->GetYaxis()->SetRangeUser(8200,9000);
    else
      h_mugast[i]->GetYaxis()->SetRangeUser(7200,8200);
    
    h_mugast[i]->Draw("colz");
    gPad->SetLogz();
  }
  dir_mugast->Close();
  dir->cd();

  TDirectoryFile *dir_must2 = (TDirectoryFile*)(dir->Get("Must2"));
  dir_must2->cd();
  binning = 10;
  TH2F *h_must[MUST2_TELESCOPES*2];
  h_must[0] = (TH2F*) dir_must2->Get("MM1_STRX_E_RAW");
  h_must[1] = (TH2F*) dir_must2->Get("MM2_STRX_E_RAW");
  h_must[2] = (TH2F*) dir_must2->Get("MM3_STRX_E_RAW");
  h_must[3] = (TH2F*) dir_must2->Get("MM4_STRX_E_RAW");
  h_must[4] = (TH2F*) dir_must2->Get("MM5_STRX_E_RAW");
  h_must[5] = (TH2F*) dir_must2->Get("MM1_STRY_E_RAW");
  h_must[6] = (TH2F*) dir_must2->Get("MM2_STRY_E_RAW");
  h_must[7] = (TH2F*) dir_must2->Get("MM3_STRY_E_RAW");
  h_must[8] = (TH2F*) dir_must2->Get("MM4_STRY_E_RAW");
  h_must[9] = (TH2F*) dir_must2->Get("MM5_STRY_E_RAW");

  TCanvas *c1 = new TCanvas("c1", "Raw energy spectra of MUST2 (X e Y)", 1200, 600);
  c1->Divide(5,2);
  for(int i=0; i<MUST2_TELESCOPES*2; i++){
    h_must[i]->SetDirectory(0);
    c1->cd(i+1); 
    if(i<MUST2_TELESCOPES)
      h_must[i]->GetYaxis()->SetRangeUser(8200,9000);
    else
      h_must[i]->GetYaxis()->SetRangeUser(7200,8200);
 
    
    h_must[i]->Draw("colz");
    gPad->SetLogz();
  }
  dir_must2->Close();
}

void DrawEnergyCalSpectra(const char* input_file){

  TFile f(input_file);
  TDirectory *dir = (TDirectory*)(f.Get("ControlSpectra"));
  dir->cd();
  TDirectoryFile *dir_mugast = (TDirectoryFile*)(dir->Get("Mugast"));
  dir_mugast->cd();
  int binning = 200;
  TH2F *h_mugast[MUGAST_TELESCOPES*2];

  h_mugast[0] = (TH2F*) dir_mugast->Get("MG1_STRX_E_CAL");
  h_mugast[1] = (TH2F*) dir_mugast->Get("MG2_STRX_E_CAL");
  h_mugast[2] = (TH2F*) dir_mugast->Get("MG3_STRX_E_CAL");
  h_mugast[3] = (TH2F*) dir_mugast->Get("MG4_STRX_E_CAL");
  h_mugast[4] = (TH2F*) dir_mugast->Get("MG5_STRX_E_CAL");
//  h_mugast[5] = (TH2F*) dir_mugast->Get("MG6_STRX_E_CAL");
  h_mugast[5] = (TH2F*) dir_mugast->Get("MG7_STRX_E_CAL");
//  h_mugast[7] = (TH2F*) dir_mugast->Get("MG9_STRX_E_CAL");
//  h_mugast[8] = (TH2F*) dir_mugast->Get("MG11_STRX_E_CAL");
  h_mugast[6] = (TH2F*) dir_mugast->Get("MG1_STRY_E_CAL");
  h_mugast[7] = (TH2F*) dir_mugast->Get("MG2_STRY_E_CAL");
  h_mugast[8] = (TH2F*) dir_mugast->Get("MG3_STRY_E_CAL");
  h_mugast[9] = (TH2F*) dir_mugast->Get("MG4_STRY_E_CAL");
  h_mugast[10] = (TH2F*) dir_mugast->Get("MG5_STRY_E_CAL");
//  h_mugast[14] = (TH2F*) dir_mugast->Get("MG6_STRY_E_CAL");
  h_mugast[11] = (TH2F*) dir_mugast->Get("MG7_STRY_E_CAL");
//  h_mugast[16] = (TH2F*) dir_mugast->Get("MG9_STRY_E_CAL");
//  h_mugast[17] = (TH2F*) dir_mugast->Get("MG11_STRY_E_CAL");

  TCanvas *c0 = new TCanvas("c0", "Calibrated energy spectra of Mugast (X e Y)", 1200, 900);
  c0->Divide(4,3);
  for(int i=0; i<MUGAST_TELESCOPES*2; i++){
    h_mugast[i]->SetDirectory(0);
    h_mugast[i]->GetYaxis()->SetRange(4*binning, 7*binning);
    c0->cd(i+1); 
    h_mugast[i]->Draw("colz");
    gPad->SetLogz();
  }
  dir_mugast->Close();
  dir->cd();

  TDirectoryFile *dir_must2 = (TDirectoryFile*)(dir->Get("Must2"));
  dir_must2->cd();
  binning = 10;
  TH2F *h_must[MUST2_TELESCOPES*2];
  h_must[0] = (TH2F*) dir_must2->Get("MM1_STRX_E_CAL");
  h_must[1] = (TH2F*) dir_must2->Get("MM2_STRX_E_CAL");
  h_must[2] = (TH2F*) dir_must2->Get("MM3_STRX_E_CAL");
  h_must[3] = (TH2F*) dir_must2->Get("MM4_STRX_E_CAL");
  h_must[4] = (TH2F*) dir_must2->Get("MM5_STRX_E_CAL");
  h_must[5] = (TH2F*) dir_must2->Get("MM1_STRY_E_CAL");
  h_must[6] = (TH2F*) dir_must2->Get("MM2_STRY_E_CAL");
  h_must[7] = (TH2F*) dir_must2->Get("MM3_STRY_E_CAL");
  h_must[8] = (TH2F*) dir_must2->Get("MM4_STRY_E_CAL");
  h_must[9] = (TH2F*) dir_must2->Get("MM5_STRY_E_CAL");

  TCanvas *c1 = new TCanvas("c1", "Calibrated energy spectra of MUST2 (X e Y)", 1200, 600);
  c1->Divide(5,2);
  for(int i=0; i<MUST2_TELESCOPES*2; i++){
    h_must[i]->SetDirectory(0);
    h_must[i]->GetYaxis()->SetRange(4*binning, 7*binning);
    c1->cd(i+1); 
    h_must[i]->Draw("colz");
    gPad->SetLogz();
  }
  dir_must2->Close();
}
