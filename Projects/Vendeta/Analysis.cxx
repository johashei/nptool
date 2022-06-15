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

  unsigned int FC_mult = FC->Energy.size();
  if(FC_mult==1){
    int anode = FC->AnodeNumber[0];

    Vendeta->SetAnodeNumber(anode);
    Vendeta->BuildPhysicalEvent();
  }

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab);
  RootOutput::getInstance()->GetTree()->Branch("Tof",&Tof);
  RootOutput::getInstance()->GetTree()->Branch("Charge",&Charge);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  ThetaLab.clear();
  ELab.clear();
  Tof.clear();
  Charge.clear();
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

