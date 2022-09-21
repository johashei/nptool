/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2025                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPOptionManager.h"
#include"RootOutput.h"
#include"RootInput.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  m_ChiNu= (TChiNuPhysics*) m_DetectorManager->GetDetector("ChiNu");
  InitialConditions = new TInitialConditions();
  InteractionCoordinates = new TInteractionCoordinates();
  ReactionConditions = new TReactionConditions();

  InitInputBranch();
  InitOutputBranch();

  my_Reaction = new NPL::Reaction("1n(238U,1n)238U@0.75");

  neutron = new NPL::Nucleus("1n");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  Einit = InitialConditions->GetKineticEnergy(0);
  double init_ThetaLab = ReactionConditions->GetTheta(0)*deg;
  double init_BeamEnergy = ReactionConditions->GetBeamEnergy();
  //neutron->SetKineticEnergy(init_BeamEnergy);
  neutron->SetKineticEnergy(Einit);
  double beam_TOF = neutron->GetTimeOfFlight();

  double Xtarget = InitialConditions->GetIncidentPositionX();
  double Ytarget = InitialConditions->GetIncidentPositionY();
  double Ztarget = 0;//InitialConditions->GetIncidentPositionZ();
  TVector3 TargetPos = TVector3(Xtarget,Ytarget,Ztarget);

  for(int i=0; i<m_ChiNu->Energy.size(); i++){
    if(m_ChiNu->Energy.size()>0){
      double Rdet, R;
      Rdet = m_ChiNu->GetDetectorPosition(m_ChiNu->DetectorNumber[i]);
      TVector3 DetPos = m_ChiNu->GetVectorDetectorPosition(m_ChiNu->DetectorNumber[i]);
      TVector3 HitPos = DetPos-TargetPos;
      //R= HitPos.Mag()*1e-3;
      R= Rdet*mm;
      Distance.push_back(R);	
      Det.push_back(m_ChiNu->DetectorNumber[i]); 
      T.push_back(m_ChiNu->Time[i]);
      double T_stop = (m_ChiNu->Time[i])*1e-9;
      //neutron->SetTimeOfFlight((T_stop-beam_TOF)/(Rdet*1e-3));
      neutron->SetTimeOfFlight((T_stop)/(Rdet*1e-3-8e-3));
      E.push_back(m_ChiNu->Energy[i]);
      Elab.push_back(neutron->GetEnergy());


      double DeltaTheta = atan(89.0/Rdet);
      double exp_ThetaLab = m_ChiNu->GetVectorDetectorPosition(m_ChiNu->DetectorNumber[i]).Theta();
      double random_ThetaLab = ra.Uniform(exp_ThetaLab-DeltaTheta, exp_ThetaLab+DeltaTheta);
      double dEx = my_Reaction->ReconstructRelativistic(Elab[i], random_ThetaLab);
      
      ThetaLab.push_back(random_ThetaLab/deg);
      //ThetaLab.push_back(exp_ThetaLab/deg);
      Ex.push_back(dEx);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Einit",&Einit,"Einit/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("Elab",&Elab);   
  RootOutput::getInstance()->GetTree()->Branch("E",&E);   
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex);   
  RootOutput::getInstance()->GetTree()->Branch("T",&T);   
  RootOutput::getInstance()->GetTree()->Branch("Distance",&Distance);   
  RootOutput::getInstance()->GetTree()->Branch("Det",&Det);   
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
  RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
  RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);

  RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
  RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
  RootInput:: getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&ReactionConditions);
}

////////////////////////////////////////////////////////////////////////////////     
void Analysis::ReInitValue(){
  Einit      = -100;
  Ex.clear();
  ThetaLab.clear();
  Elab.clear();
  E.clear();
  T.clear();
  Distance.clear();
  Det.clear();
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

