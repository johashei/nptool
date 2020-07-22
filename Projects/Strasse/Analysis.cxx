/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2012                                               *
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
#include"NPFunction.h"
#include"NPTrackingUtility.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  IC= new TInitialConditions;
  DC= new TInteractionCoordinates;
  RC= new TReactionConditions;
  

  InitOutputBranch();
  InitInputBranch();
  
  Strasse = (TStrassePhysics*)  m_DetectorManager -> GetDetector("Strasse");
  myReaction = new NPL::Reaction();
  myReaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
} 

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  // Reinitiate calculated variable
  ReInitValue();
  unsigned int size = Strasse->GetEventMultiplicity();
  if(size==2){ // 2 proton detected
    // Proton 1
    TVector3 InnerPos1 = Strasse->GetInnerPositionOfInteraction(0);
    TVector3 OuterPos1 = Strasse->GetOuterPositionOfInteraction(0);
    TVector3 Proton1 = OuterPos1-InnerPos1;
    // Proton 2
    TVector3 InnerPos2 = Strasse->GetInnerPositionOfInteraction(1);
    TVector3 OuterPos2 = Strasse->GetOuterPositionOfInteraction(1);
    TVector3 Proton2 = OuterPos2-InnerPos2;

    double deltaPhi = abs(Proton1.Phi()/deg-Proton2.Phi()/deg);
    double sumTheta = Proton1.Theta()/deg+Proton2.Theta()/deg;
    double OpeningAngle = Proton1.Angle(Proton2)/deg;
   cout << OpeningAngle << endl;  
    // reject event that make no physical sense
    if(deltaPhi<170 && sumTheta<80){
      return;
      }
    
    // computing minimum distance of the two lines
    TVector3 Vertex;
    Distance = NPL::MinimumDistance(InnerPos1,OuterPos1,InnerPos2,OuterPos2,Vertex);
    VertexX=Vertex.X();
    VertexY=Vertex.Y();
    VertexZ=Vertex.Z();
    }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab,"ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("VerteX",&VertexX,"VertexX/D");
  RootOutput::getInstance()->GetTree()->Branch("VerteY",&VertexY,"VertexY/D");
  RootOutput::getInstance()->GetTree()->Branch("VerteZ",&VertexZ,"VertexZ/D");
  RootOutput::getInstance()->GetTree()->Branch("Distance",&Distance,"Distance/D");
  RootOutput::getInstance()->GetTree()->Branch("InteractionCoordinates","TInteractionCoordinates",&DC);
  RootOutput::getInstance()->GetTree()->Branch("ReactionConditions","TReactionConditions",&RC);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
    RootInput:: getInstance()->GetChain()->SetBranchAddress("InteractionCoordinates",&DC);
    RootInput:: getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&RC);
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  Ex = -1000 ;
  ELab = -1000;
  ThetaLab = -1000;
  ThetaCM = -1000;
  VertexX=-1000;
  VertexY=-1000;
  VertexZ=-1000;
  Distance=-1000;
}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the AnalysisFactory              //
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

