/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace contact address: pierre.morfouace2@cea.fr   *
 *                                                                           *
 * Creation Date  : June 2021                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Sofia analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
#include<algorithm>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPPhysicalConstants.h"
#include"NPGlobalSystemOfUnits.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  SofBeamID = new TSofBeamID();
  SofFF = new TSofFissionFragment();
  SofSci= (TSofSciPhysics*) m_DetectorManager->GetDetector("SofSci");
  SofTrim= (TSofTrimPhysics*) m_DetectorManager->GetDetector("SofTrim");
  SofTwim= (TSofTwimPhysics*) m_DetectorManager->GetDetector("SofTwim");
  SofTofW= (TSofTofWPhysics*) m_DetectorManager->GetDetector("SofTofW");
  SofAt= (TSofAtPhysics*) m_DetectorManager->GetDetector("SofAt");
  SofMwpc= (TSofMwpcPhysics*) m_DetectorManager->GetDetector("SofMwpc");

  InitParameter();
  InitOutputBranch();
  LoadCut();
  LoadSpline();

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  //cout << "************" << endl;
  RunID = fRunID;
  BeamAnalysis();

  unsigned int sofsci_size = SofSci->DetectorNbr.size();
  if(sofsci_size==2){
    double start_time = SofSci->TimeNs[1];
    SofTofW->SetTofAlignedValue(36);
    SofTofW->SetStartTime(start_time);
    SofTofW->BuildPhysicalEvent();

    FissionFragmentAnalysis();
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::FissionFragmentAnalysis(){
  unsigned int softofw_size = SofTofW->PlasticNbr.size();
  unsigned int softwim_size = SofTwim->SectionNbr.size();

  double TOF_CC[2];
  double Plastic[2];
  double Plastic_left = -1;
  double Plastic_right = -1;
  double TOF_left = -1;
  double TOF_right = -1;
  double TOF_up = -1;
  double TOF_down = -1;
  double E1 = -1;
  double E2 = -1;
  double E3 = -1;
  double E4 = -1;
  double DT1 = -1000;
  double DT2 = -1000;
  double DT3 = -1000;
  double DT4 = -1000;
  double L_CC = 8.45;
  double Beta_left = -1;
  double Beta_right = -1;
  double Beta_up = -1;
  double Beta_down = -1;
  double Beta1 = -1;
  double Beta2 = -1;
  double Beta3 = -1;
  double Beta4 = -1;
  double Beta_norm = 0.745;
  double Gamma1 = -1;
  double Gamma2 = -1;
  double Beta_Z1 = -1;
  double Beta_Z2 = -1;
  double Brho1 = -1;
  double Brho2 = -1;
  double AoQ1 = -1;
  double AoQ2 = -1;
  double A1 = -1;
  double A2 = -1;
  double Z1 = -1;
  double Z2 = -1;
  double Zsum = -1;
  int iZ1 = -1;
  int iZ2 = -1;
  int iZsum = -1;
  double ThetaIn1 = -1000;
  double ThetaIn2 = -1000;

  for(int i = 0; i<2; i++){
    TOF_CC[i] = -1;
    Plastic[i] = -1;
  }

  vector<double> PosY;
  vector<double> PosX;
  if(softofw_size==2){
    for(unsigned int i=0; i<softofw_size; i++){
      TOF_CC[i] = SofTofW->CalTof[i];
      Plastic[i] = SofTofW->PlasticNbr[i];
      PosX.push_back(SofTofW->CalPosX[i]);
      PosY.push_back(SofTofW->CalPosY[i]);

      //SofFF->SetTofPosX(SofTofW->CalPosX[i]);
      //SofFF->SetTofPosY(SofTofW->CalPosY[i]);
    }
  }

  vector<double> Xmwpc4;
  vector<double> Ymwpc4;
  for(unsigned int i=0; i<SofMwpc->DetectorNbr.size(); i++){
    if(SofMwpc->DetectorNbr[i]==4){
      Xmwpc4.push_back(SofMwpc->PositionX1[i]);
      Ymwpc4.push_back(SofMwpc->PositionY[i]);

      //SofFF->SetMwpcPosX(SofMwpc->PositionX1[i]);
      //SofFF->SetMwpcPosY(SofMwpc->PositionY[i]);
    }
  }

  vector<double> good_posx;
  vector<double> good_posy;
  for(unsigned int i=0; i<PosX.size(); i++){
    double tofx = PosX[i];
    double tofy = PosY[i];
    for(unsigned int k=0; k<Xmwpc4.size(); k++){
      double posx = Xmwpc4[k];
      if(abs(posx-tofx) < 100){
        good_posx.push_back(posx);
        SofFF->SetMwpcPosX(posx);
        SofFF->SetTofPosX(tofx);
      }
    }
    for(unsigned int p=0; p<Ymwpc4.size(); p++){
      double posy = Ymwpc4[p];
      if(abs(posy-tofy) < 20){
        good_posy.push_back(posy);
        SofFF->SetMwpcPosY(posy);
        SofFF->SetTofPosY(tofy);
      }
    }
  }

  int mult1 = SofTwim->mult1;
  int mult2 = SofTwim->mult2;
  int mult3 = SofTwim->mult3;
  int mult4 = SofTwim->mult4;

  int multL = mult1 + mult2;
  int multR = mult3 + mult4;

  if(softwim_size>1){
    if( (mult1>1 && mult1<17) || (mult2>1 && mult2<17) || (mult3>1 && mult3<17) || (mult4>1 && mult4<17)){
      for(unsigned int i=0; i< softwim_size; i++){
        int sec = SofTwim->SectionNbr[i];

        if(sec==1){
          E1 = SofTwim->EnergySection[i];
          DT1 = SofTwim->DriftTime[i];
        }
        if(sec==2){
          E2 = SofTwim->EnergySection[i];
          DT2 = SofTwim->DriftTime[i];
        }
        if(sec==3){
          E3 = SofTwim->EnergySection[i]; 
          DT3 = SofTwim->DriftTime[i];
        }
        if(sec==4){
          E4 = SofTwim->EnergySection[i];     
          DT4 = SofTwim->DriftTime[i];
        }
      }

      if(softwim_size>2){
        if(E1>0 && E2>0){
          E1 = E1+E2;
          E2 = -1;
        }
        if(E3>0 && E4>0){
          E3 = E3+E4;
          E4 = -1;
        }
      }

      if(E1>0)
        E1 = E1/16;
      if(E2>0)
        E2 = E2/16;
      if(E3>0)
        E3 = E3/16;
      if(E4>0)
        E4 = E4/16;

      if(Plastic[0]<Plastic[1]){
        Plastic_left = Plastic[0];
        Plastic_right = Plastic[1];
        TOF_left = TOF_CC[0];
        TOF_right = TOF_CC[1];
      }
      if(Plastic[0]>Plastic[1]){
        Plastic_left = Plastic[1];
        Plastic_right = Plastic[0];
        TOF_left = TOF_CC[1];
        TOF_right = TOF_CC[0];
      }

      if(PosY.size()==2){
        if(PosY[0]>PosY[1]){
          TOF_up = TOF_CC[0];
          TOF_down = TOF_CC[1];
        }
        if(PosY[0]<PosY[1]){
          TOF_up = TOF_CC[1];
          TOF_down = TOF_CC[0];
        }
      }

      double velocity_left = L_CC/TOF_left;
      double velocity_right = L_CC/TOF_right;
      Beta_left = velocity_left * m/ns / NPUNITS::c_light;
      Beta_right = velocity_right * m/ns / NPUNITS::c_light;

      double velocity_down = L_CC/TOF_down;
      double velocity_up = L_CC/TOF_up;
      Beta_down = velocity_down * m/ns / NPUNITS::c_light;
      Beta_up = velocity_up * m/ns / NPUNITS::c_light;

      if(E1!=-1 && E2!=-1){
        Beta1 = Beta_down;
        Beta2 = Beta_up;
      }
      if(E1!=-1 && E3!=-1){
        Beta1 = Beta_down;
        Beta3 = Beta_up;
      }
      if(E1!=-1 && E4!=-1){
        Beta1 = Beta_left;
        Beta4 = Beta_right;
      }
      if(E2!=-1 && E3!=-1){
        Beta2 = Beta_left;
        Beta3 = Beta_right;
      }
      if(E2!=-1 && E4!=-1){
        Beta2 = Beta_up;
        Beta4 = Beta_down;
      }
      if(E3!=-1 && E4!=-1){
        Beta3 = Beta_up;
        Beta4 = Beta_down;
      }


      if(Beta1!=-1 && E1>0){ 
        E1 = E1 / fcorr_z_beta[0]->Eval(Beta1) * fcorr_z_beta[0]->Eval(Beta_norm);
        E1 = E1 / fcorr_z_dt[0]->Eval(DT1) * fcorr_z_dt[0]->Eval(55);
      }
      if(Beta2!=-1 && E2>0){ 
        E2 = E2 / fcorr_z_beta[1]->Eval(Beta2) * fcorr_z_beta[1]->Eval(Beta_norm);
        E2 = E2 / fcorr_z_dt[1]->Eval(DT2) * fcorr_z_dt[1]->Eval(55);
      }
      if(Beta3!=-1 && E3>0){ 
        E3 = E3 / fcorr_z_beta[2]->Eval(Beta3) * fcorr_z_beta[2]->Eval(Beta_norm);
        E3 = E3 / fcorr_z_dt[2]->Eval(DT3) * fcorr_z_dt[2]->Eval(-55);
      }
      if(Beta4!=-1 && E4>0){ 
        E4 = E4 / fcorr_z_beta[3]->Eval(Beta4) * fcorr_z_beta[3]->Eval(Beta_norm);
        E4 = E4 / fcorr_z_dt[3]->Eval(DT4) * fcorr_z_dt[3]->Eval(-55);
      }


      // Z calibration //
      if(E1>0 && E2>0 && E3==-1 && E4==-1){
        Z1 = E1;
        Z2 = E2;
        Beta_Z1 = Beta1;
        Beta_Z2 = Beta2;
      }
      if(E1>0 && E2==-1 && E3>0 && E4==-1){
        Z1 = E1;
        Z2 = E3;
        Beta_Z1 = Beta1;
        Beta_Z2 = Beta3;
      }
      if(E1>0 && E2==-1 && E3==-1 && E4>0){
        Z1 = E1;
        Z2 = E4;
        Beta_Z1 = Beta1;
        Beta_Z2 = Beta4;
      }
      if(E1==-1 && E2>0 && E3>0 && E4==-1){
        Z1 = E2;
        Z2 = E3;
        Beta_Z1 = Beta2;
        Beta_Z2 = Beta3;
      }
      if(E1==-1 && E2>0 && E3==-1 && E4>0){
        Z1 = E2;
        Z2 = E4;
        Beta_Z1 = Beta2;
        Beta_Z2 = Beta4;
      }
      if(E1==-1 && E2==-1 && E3>0 && E4>0){
        Z1 = E3;
        Z2 = E4;
        Beta_Z1 = Beta3;
        Beta_Z2 = Beta4;
      }

      if(Z1>0 && Z2>0){
        Z1 = fZff_p0 + fZff_p1*Z1 + fZff_p2*Z1*Z1;
        Z2 = fZff_p0 + fZff_p1*Z2 + fZff_p2*Z2*Z2;

        Z1 = sqrt(Z1);
        Z2 = sqrt(Z2);

        Zsum = Z1+Z2;

        iZ1 = (int) round(Z1);
        iZ2 = (int) round(Z2);
        iZsum = iZ1 + iZ2;
      }

      Gamma1 = 1. / sqrt(1 - Beta_Z1 * Beta_Z1);
      Gamma2 = 1. / sqrt(1 - Beta_Z2 * Beta_Z2);

      AoQ1 = Brho1 / (3.10761 * Beta_Z1 * Gamma1);
      AoQ2 = Brho2 / (3.10761 * Beta_Z2 * Gamma2);

      A1 = AoQ1 * iZ1;
      A2 = AoQ2 * iZ2;

      //*** Filling the Fission Fragment Tree ***//
      SofFF->SetTOF(TOF_left);
      SofFF->SetTOF(TOF_right);
      SofFF->SetBeta(Beta_Z1);
      SofFF->SetBeta(Beta_Z2);
      SofFF->SetGamma(Gamma1);
      SofFF->SetGamma(Gamma2);
      SofFF->SetiZ(iZ1);
      SofFF->SetiZ(iZ2);
      SofFF->SetZ(Z1);
      SofFF->SetZ(Z2);
      SofFF->SetAoQ(AoQ1);
      SofFF->SetAoQ(AoQ1);
      SofFF->SetA(A1);
      SofFF->SetA(A2);
      SofFF->SetBrho(Brho1);
      SofFF->SetBrho(Brho2);

      SofFF->SetDT(DT1);
      SofFF->SetDT(DT2);
      SofFF->SetDT(DT3);
      SofFF->SetDT(DT4);

      SofFF->SetZsum(Zsum);
      SofFF->SetIntZsum(iZsum);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::BeamAnalysis(){
  unsigned int sofsci_size = SofSci->DetectorNbr.size();
  if(sofsci_size==2){
    double beta = SofSci->Beta[0];
    //cout << "Set beta to " << beta << endl;
    SofTrim->SetBeta(beta);
    SofTrim->BuildSimplePhysicalEvent();
    double Zbeam,Qmax,Theta;
    if(SofTrim->EnergySection.size()>0){
      double Anode1 = SofTrim->EnergySection[0];
      double Anode2 = SofTrim->EnergySection[1];
      double Anode3 = SofTrim->EnergySection[2];
      if(fRunID==13)
        Zbeam = max(Anode1, Anode2);
      else 
        Zbeam = SofTrim->GetMaxEnergySection();

      //Qmax = DetermineQmax();
      Theta = SofTrim->Theta[0];

      double TofFromS2    = SofSci->CalTof[0];
      double velocity_mns = SofSci->VelocityMNs[0];
      double Beta         = SofSci->Beta[0];
      double XS2          = SofSci->PosMm[0];
      //double XCC          = SofSci->PosMm[1];
      double XCC=0;
      double YCC=0;
      for(unsigned int i=0; i<SofMwpc->DetectorNbr.size(); i++){
        if(SofMwpc->DetectorNbr[i]==1){
          XCC = SofMwpc->PositionX1[i];
          YCC = SofMwpc->PositionY[i];
        }
      }

      double LS2;
      LS2 = fLS2_0;//*(1 + fK_LS2*Theta);
      velocity_mns = LS2/TofFromS2;
      Beta = velocity_mns * m/ns / NPUNITS::c_light;
      double Gamma        = 1./(TMath::Sqrt(1 - TMath::Power(Beta,2)));
      double Brho = fBrho0 * (1 - XS2/fDS2 - XCC/fDCC);
      double AoQ  = Brho / (3.10716*Gamma*Beta);
      double A    = AoQ * Qmax;

      // Y dependence correction //
      double Y_p0 = 23943.8;
      double Y_p1 = 12.362;
      Zbeam = Zbeam/(Y_p0 + Y_p1*YCC)*Y_p0;

      // Z calibration //
      Zbeam = fZbeam_p0 + fZbeam_p1*Zbeam + fZbeam_p2*Zbeam*Zbeam;
      Zbeam = sqrt(Zbeam);

      // Last beta correction //
      double Beta_norm = 0.8355;
      Zbeam = Zbeam/(fZBeta_p0 + fZBeta_p1*Beta)*(fZBeta_p0 + fZBeta_p1*Beta_norm);

      // Filling Beam tree
      SofBeamID->SetZbeam(Zbeam);
      //SofBeamID->SetQmax(rand.Gaus(Qmax,0.15));
      SofBeamID->SetAoQ(AoQ);
      SofBeamID->SetAbeam(A);
      SofBeamID->SetBeta(Beta);
      SofBeamID->SetGamma(Gamma);
      SofBeamID->SetBrho(Brho);
      SofBeamID->SetXS2(XS2);
      SofBeamID->SetXCC(XCC);
      SofBeamID->SetYCC(YCC);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadCut(){
  TString input_path = "./calibration/SofTrim/cut/";

  TString rootfile;
  TString cutfile;
  TFile* file;
  for(int i=0; i<3; i++){
    // Q=77
    rootfile = Form("cutsec%iQ77.root",i+1);
    cutfile = input_path + rootfile;
    file = new TFile(cutfile);
    cutQ77[i] = (TCutG*) file->Get(Form("cutsec%iQ77",i+1));
    // Q=78
    rootfile = Form("cutsec%iQ78.root",i+1);
    cutfile = input_path + rootfile;
    file = new TFile(cutfile);
    cutQ78[i] = (TCutG*) file->Get(Form("cutsec%iQ78",i+1));
    // Q=79
    rootfile = Form("cutsec%iQ79.root",i+1);
    cutfile = input_path + rootfile;
    file = new TFile(cutfile);
    cutQ79[i] = (TCutG*) file->Get(Form("cutsec%iQ79",i+1));
    // Q=80
    rootfile = Form("cutsec%iQ80.root",i+1);
    cutfile = input_path + rootfile;
    file = new TFile(cutfile);
    cutQ80[i] = (TCutG*) file->Get(Form("cutsec%iQ80",i+1));
  }
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadSpline(){
  TString input_path = "./calibration/SofTwim/spline/";

  TString rootfile = input_path + "spline_beta.root";
  TFile* ifile = new TFile(rootfile,"read");

  TString splinename;
  if(ifile->IsOpen()){
    cout << "Loading Beta spline for fission fragment analysis..." << endl;
    for(int i=0; i<4; i++){
      splinename = Form("spline_beta_sec%i",i+1);
      fcorr_z_beta[i] = (TSpline3*) ifile->FindObjectAny(splinename);
    }
    ifile->Close();
  }
  else
    cout << "File " << rootfile << " not found!" << endl;

  //*** ***//
  rootfile = input_path + "spline_dt.root";
  ifile = new TFile(rootfile,"read");

  if(ifile->IsOpen()){
    cout << "Loading DT spline for fission fragment analysis..." << endl;
    for(int i=0; i<4; i++){
      splinename = Form("spline_dt_sec%i",i+1);
      fcorr_z_dt[i] = (TSpline3*) ifile->FindObjectAny(splinename);
    }
    ifile->Close();
  }
  else
    cout << "File " << rootfile << " not found!" << endl;

}

////////////////////////////////////////////////////////////////////////////////
int Analysis::DetermineQmax(){
  int Qmax;
  int Qsec[3];

  unsigned int size = SofTrim->EnergySection.size();
  for(unsigned int i=0; i<size; i++){
    int SecNbr   = SofTrim->SectionNbr[i];
    double Theta = SofTrim->Theta[i];
    double Esec  = SofTrim->EnergySection[i];


    if(cutQ77[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 77;
    else if(cutQ78[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 78;
    else if(cutQ79[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 79;
    else if(cutQ80[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 80;
    //else if(cutQ81[SecNbr-1]->IsInside(Theta,Esec))
    //Qsec[SecNbr-1] = 81;
  }

  Qmax = max(Qsec[0],Qsec[1]);
  Qmax = max(Qsec[2],Qmax);

  return Qmax;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitParameter(){
  fLS2_0 = 135.614;
  fDS2   = 8000;
  fDCC   = -10000;
  fK_LS2 = -30e-8;

  fRunID = 6;

  // Beam parameter //
  fZBeta_p0 = 1;
  fZBeta_p1 = 0;

  // FF parameter //
  fZff_p0 = 2.80063;
  fZff_p1 = 6.91985e-2;
  fZff_p2 = 1.01598e-7;

  if(fRunID==1 || fRunID==2){
    //fBrho0 = 10.6813; // 180Hg
    fBrho0 = 10.6955; // 180Hg
    fZbeam_p0 = -5303.06;
    fZbeam_p1 = 0.674945;
    fZbeam_p2 = -8.32085e-6;

    fZBeta_p0 = 72.946;
    fZBeta_p1 = 6.0644;
  }
  if(fRunID==3){
    fBrho0 = 10.8183; // 182Hg
    fZbeam_p0 = -2737.25;
    fZbeam_p1 = 0.452017;
    fZbeam_p2 = -3.48831e-6;

    fZBeta_p0 = 76.6738;
    fZBeta_p1 = 1.60128;
  }
  if(fRunID==4){
    fBrho0 = 10.9558; // 184Hg
    fZbeam_p0 = -5044.61;
    fZbeam_p1 = 0.639986;
    fZbeam_p2 = -7.3077e-6;
  }
  if(fRunID==5){
    fBrho0 = 10.8138; // 187Pb
    fZbeam_p0 = -2858.72;
    fZbeam_p1 = 0.454064;
    fZbeam_p2 = -3.36443e-6;

    fZBeta_p0 = 71.0975;
    fZBeta_p1 = 10.7007;
  }
  if(fRunID==6){
    fBrho0 = 10.9476; // 189Pb
    fZbeam_p0 = 1590.66;
    fZbeam_p1 = 0.0956455;
    fZbeam_p2 = 3.84585e-6;

    fZBeta_p0 = 74.6063;
    fZBeta_p1 = 6.4635;
  }
  if(fRunID==7){
    fBrho0 = 10.6814; // 175Pt
    fZbeam_p0 = 459.68;
    fZbeam_p1 = 0.162277;
    fZbeam_p2 = 3.10164e-6;

    fZBeta_p0 = 66.9433;
    fZBeta_p1 = 10.8664;
  }
  if(fRunID==8){
    fBrho0 = 11.0864; // 204Fr
    fZbeam_p0 = 4122.94;
    fZbeam_p1 = -0.119867;
    fZbeam_p2 = 8.29115e-6;

    fZBeta_p0 = 63.9575;
    fZBeta_p1 = 25.1988;
  }
  if(fRunID==9){
    fBrho0 = 11.2712; // 207Fr
    fZbeam_p0 = -1752.27;
    fZbeam_p1 = 0.346018;
    fZbeam_p2 = -8.64673e-7;

    fZBeta_p0 = 63.9575;
    fZBeta_p1 = 25.1988;
  }
  if(fRunID==10){
    fBrho0 = 11.0955; // 199At run 423 & 424
    fZbeam_p0 = -116.425;
    fZbeam_p1 = 0.218256;
    fZbeam_p2 = 1.62399e-6;

    fZBeta_p0 = 61.3889;
    fZBeta_p1 = 25.8908;
  }
  if(fRunID==11){
    fBrho0 = 10.9970; // 199At run 425 & 426
    fZbeam_p0 = -116.425;
    fZbeam_p1 = 0.218256;
    fZbeam_p2 = 1.62399e-6;

    fZBeta_p0 = 61.3889;
    fZBeta_p1 = 25.8908;
  }
  if(fRunID==12){
    fBrho0 = 10.8697; //197At
    fZbeam_p0 = -2683.52;
    fZbeam_p1 = 0.422551;
    fZbeam_p2 = -2.44471e-6;

    fZBeta_p0 = 62.9188;
    fZBeta_p1 = 22.8827;
  }
  if(fRunID==13){
    fBrho0 = 11.3418; // 216Th
    fZbeam_p0 = 1651.57;
    fZbeam_p1 = 0.0876127;
    fZbeam_p2 = 4.02563e-6;

    fZBeta_p0 = 38.879;
    fZBeta_p1 = 58.7667;
  }
  if(fRunID==14){
    fBrho0 = 11.5067; // 221Pa
    fZbeam_p0 = 186.892;
    fZbeam_p1 = 0.20739;
    fZbeam_p2 = 1.61797e-6;

  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  SofBeamID->Clear();
  SofFF->Clear();
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("RunID",&RunID,"RunID/I");
  RootOutput::getInstance()->GetTree()->Branch("SofBeamID","TSofBeamID",&SofBeamID);
  RootOutput::getInstance()->GetTree()->Branch("SofFissionFragment","TSofFissionFragment",&SofFF);

}
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy{
    public:
      proxy(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
      }
  };

  proxy p;
}


