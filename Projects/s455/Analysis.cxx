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
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  SofBeamID = new TSofBeamID();
  SofSci= (TSofSciPhysics*) m_DetectorManager->GetDetector("SofSci");
  SofTrim= (TSofTrimPhysics*) m_DetectorManager->GetDetector("SofTrim");
  //SofTofW= (TSofTofWPhysics*) m_DetectorManager->GetDetector("SofTofW");

  InitParameter();
  InitOutputBranch();
  LoadCut();


}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  //cout << "************" << endl;
  BeamAnalysis();

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
    }

    double TofFromS2    = SofSci->CalTof[0];
    double velocity_mns = SofSci->VelocityMNs[0];
    double Beta         = SofSci->Beta[0];
    double Gamma        = 1./(TMath::Sqrt(1 - TMath::Power(Beta,2)));
    double XS2          = SofSci->PosMm[0];
    double XCC          = SofSci->PosMm[1];
    double LS2;
    LS2 = fLS2_0*(1 + fK_LS2*Theta);
    velocity_mns = LS2/TofFromS2;
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
  fDCC   = -30000;
  fK_LS2 = -2.5e-8;
  fBrho0 = 10.8183; // run401 -> 182Hg

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  SofBeamID->Clear();
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  //RootOutput::getInstance()->GetTree()->Branch("Zbeam",&Zbeam,"Zbeam/D");
  RootOutput::getInstance()->GetTree()->Branch("SofBeamID","TSofBeamID",&SofBeamID);

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

