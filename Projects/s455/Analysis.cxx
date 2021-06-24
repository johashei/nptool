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

  InitParameter();
  InitOutputBranch();
  LoadCut();
  LoadSpline();

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  //cout << "************" << endl;
  BeamAnalysis();

  unsigned int sofsci_size = SofSci->DetectorNbr.size();
  if(sofsci_size==2){
    double start_time = SofSci->TimeNs[1];
    SofTofW->SetTofAlignedValue(36);
    SofTofW->SetStartTime(start_time);
    SofTofW->BuildPhysicalEvent();
  }

  FissionFragmentAnalysis();

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::FissionFragmentAnalysis(){
  unsigned int softofw_size = SofTofW->PlasticNbr.size();
  unsigned int softwim_size = SofTwim->SectionNbr.size();

  double TOF_CC[2];
  double Plastic[2];
  double TOF_left = -1;
  double TOF_right = -1;
  double Esec[2];
  double Section[2];
  double E_left = -1;
  double E_right = -1;
  double E1 = -1;
  double E2 = -1;
  double E3 = -1;
  double E4 = -1;
  double L_CC = 8.;
  double Beta_left;
  double Beta_right;
  double Beta_norm = 0.7;

  for(int i = 0; i<2; i++){
    TOF_CC[i] = -1;
    Plastic[i] = -1;
    Esec[i] = -1;
    Section[i] = -1;
  }


  if(softofw_size==2 && softwim_size==2){ 
    for(unsigned int i=0; i< softofw_size; i++){
      TOF_CC[i] = SofTofW->CalTof[i];
      Plastic[i] = SofTofW->PlasticNbr[i];

      Esec[i] = SofTwim->EnergySection[i];
      int sec = SofTwim->SectionNbr[i];
      Section[i] = sec;

      if(sec==1)
        E1 = SofTwim->EnergySection[i];
      else if(sec==2)
        E2 = SofTwim->EnergySection[i];
      else if(sec==3)
        E3 = SofTwim->EnergySection[i]; 
      else if(sec==4)
        E4 = SofTwim->EnergySection[i];     
    }
  }

  if(Plastic[0]<Plastic[1]){
    TOF_left = TOF_CC[0];
    TOF_right = TOF_CC[1];
  }
  else{
    TOF_left = TOF_CC[1];
    TOF_right = TOF_CC[0];
  }

  if(TOF_left != -1 && TOF_right != -1 && abs(Plastic[0]-Plastic[1]) != 1){
    double velocity_left = L_CC/TOF_left;
    double velocity_right = L_CC/TOF_right;

    Beta_left = velocity_left * m/ns / NPUNITS::c_light;
    Beta_right = velocity_right * m/ns / NPUNITS::c_light;

    /*E1 = E1 / fcorr_z_beta[0]->Eval(Beta_left) * fcorr_z_beta[0]->Eval(Beta_norm);
    E2 = E2 / fcorr_z_beta[1]->Eval(Beta_left) * fcorr_z_beta[1]->Eval(Beta_norm);
    E3 = E3 / fcorr_z_beta[2]->Eval(Beta_right) * fcorr_z_beta[2]->Eval(Beta_norm);
    E4 = E4 / fcorr_z_beta[3]->Eval(Beta_right) * fcorr_z_beta[3]->Eval(Beta_norm);
*/
    double Zsum = E_left + E_right;

    SofFF->SetBeta(Beta_left);
    SofFF->SetBeta(Beta_right);
    SofFF->SetZ(E1);
    SofFF->SetZ(E2);
    SofFF->SetZ(E3);
    SofFF->SetZ(E4);
    SofFF->SetZsum(Zsum);
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
      Zbeam = SofTrim->GetMaxEnergySection();
      Qmax = DetermineQmax();
      Theta = SofTrim->Theta[0];

      double TofFromS2    = SofSci->CalTof[0];
      double velocity_mns = SofSci->VelocityMNs[0];
      double Beta         = SofSci->Beta[0];
      double XS2          = SofSci->PosMm[0];
      double XCC          = SofSci->PosMm[1];
      double LS2;

      LS2 = fLS2_0*(1 + fK_LS2*Theta);
      velocity_mns = LS2/TofFromS2;
      Beta = velocity_mns * m/ns / NPUNITS::c_light;
      double Gamma        = 1./(TMath::Sqrt(1 - TMath::Power(Beta,2)));
      double Brho = fBrho0 * (1 - XS2/fDS2 - XCC/fDCC);
      double AoQ  = Brho / (3.10716*Gamma*Beta);
      double A    = AoQ * Qmax;

      // Filling Beam tree
      SofBeamID->SetZbeam(Zbeam);
      SofBeamID->SetQmax(rand.Gaus(Qmax,0.15));
      SofBeamID->SetAoQ(AoQ);
      SofBeamID->SetAbeam(A);
      SofBeamID->SetBeta(Beta);
      SofBeamID->SetGamma(Gamma);
      SofBeamID->SetBrho(Brho);
      SofBeamID->SetXS2(XS2);
      SofBeamID->SetXCC(XCC);
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
    // Q=81
    rootfile = Form("cutsec%iQ81.root",i+1);
    cutfile = input_path + rootfile;
    file = new TFile(cutfile);
    cutQ81[i] = (TCutG*) file->Get(Form("cutsec%iQ81",i+1));
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

    if(cutQ78[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 78;
    else if(cutQ79[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 79;
    else if(cutQ80[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 80;
    else if(cutQ81[SecNbr-1]->IsInside(Theta,Esec))
      Qsec[SecNbr-1] = 81;
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
  fLS2_0 = 136.3706933;
  fDS2   = 9500;
  //fDCC   = -30000;
  fDCC   = -10000;
  fK_LS2 = -2.5e-8;
  fBrho0 = 10.8183; // run401 -> 182Hg
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  SofBeamID->Clear();
  SofFF->Clear();
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  //RootOutput::getInstance()->GetTree()->Branch("Zbeam",&Zbeam,"Zbeam/D");
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

