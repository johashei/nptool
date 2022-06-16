/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Vendeta analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPOptionManager.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  Vendeta= (TVendetaPhysics*) m_DetectorManager->GetDetector("Vendeta");
  FC= (TFissionChamberPhysics*) m_DetectorManager->GetDetector("FissionChamber");

  InitOutputBranch();

  neutron = new NPL::Particle("1n");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();

  unsigned int FC_mult = FC->AnodeNumber.size();
  if(FC_mult==1){
    int anode = FC->AnodeNumber[0];
    double Time_FC = FC->Time[0];

    Vendeta->SetAnodeNumber(anode);
    Vendeta->BuildPhysicalEvent();
    unsigned int Vendeta_mult = Vendeta->DetectorNumber.size();
    for(unsigned int i=0; i<Vendeta_mult; i++){
      int DetNbr          = Vendeta->DetectorNumber[i];
      double Time_Vendeta = Vendeta->Time[i];
      double Rdet         = Vendeta->GetDistanceFromTarget(DetNbr);
      TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);

      double DT = Time_Vendeta - Time_FC;

      double DeltaTheta = atan(63.5/Rdet);
      double Theta_Vendeta = DetPos.Theta();
      double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);

      neutron->SetTimeOfFlight(DT/(Rdet));
      double En = neutron->GetEnergy();

      // Filling output tree
      Tof.push_back(DT);
      ELab.push_back(En);
      ThetaLab.push_back(Theta_random);
      Q1.push_back(Vendeta->Q1[i]);
      Q2.push_back(Vendeta->Q2[i]);
      HG_status.push_back(Vendeta->isHG[i]);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab);
  RootOutput::getInstance()->GetTree()->Branch("Tof",&Tof);
  RootOutput::getInstance()->GetTree()->Branch("Q1",&Q1);
  RootOutput::getInstance()->GetTree()->Branch("Q2",&Q2);
  RootOutput::getInstance()->GetTree()->Branch("HG_status",&HG_status);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  ThetaLab.clear();
  ELab.clear();
  Tof.clear();
  Q1.clear();
  Q2.clear();
  HG_status.clear();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
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

