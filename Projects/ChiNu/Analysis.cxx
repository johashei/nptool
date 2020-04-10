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
	
	InitInputBranch();
	InitOutputBranch();
    
	neutron = new NPL::Nucleus("1n");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
	ReInitValue();
	Einit = InitialConditions->GetKineticEnergy(0);
	
	double Xtarget = InitialConditions->GetIncidentPositionX();
	double Ytarget = InitialConditions->GetIncidentPositionY();
	double Ztarget = InitialConditions->GetIncidentPositionZ();
        TVector3 TargetPos = TVector3(Xtarget,Ytarget,Ztarget);
	    
   	for(int i=0; i<m_ChiNu->Energy.size(); i++){
		if(m_ChiNu->Energy.size()>0){
			double Rdet, R;
			Rdet = m_ChiNu->GetDetectorPosition(m_ChiNu->DetectorNumber[i]);
			TVector3 DetPos = m_ChiNu->GetVectorDetectorPosition(m_ChiNu->DetectorNumber[i]);
			TVector3 HitPos = DetPos-TargetPos;
			R= HitPos.Mag()*1e-3;
        		Distance.push_back(R);	
			Det.push_back(m_ChiNu->DetectorNumber[i]); 
			T.push_back(m_ChiNu->Time[i]);
        		neutron->SetTimeOfFlight(m_ChiNu->Time[i]*1e-9/R);
			E.push_back(m_ChiNu->Energy[i]);
        		Elab.push_back(neutron->GetEnergy());
		}
    	}



}

///////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
	RootOutput::getInstance()->GetTree()->Branch("Einit",&Einit,"Einit/D");
	RootOutput::getInstance()->GetTree()->Branch("Elab",&Elab);   
	RootOutput::getInstance()->GetTree()->Branch("E",&E);   
	RootOutput::getInstance()->GetTree()->Branch("T",&T);   
	RootOutput::getInstance()->GetTree()->Branch("Distance",&Distance);   
	RootOutput::getInstance()->GetTree()->Branch("Det",&Det);   
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  	RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
  	RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
  	RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);
}

////////////////////////////////////////////////////////////////////////////////     
void Analysis::ReInitValue(){
	Einit      = -100;
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

