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
  // reaction properties
  myQFS = new NPL::QFS();
  myQFS->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  // reaction properties
  myBeam = new NPL::Beam();
  myBeam->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  InitialBeamEnergy = myBeam->GetEnergy() * myBeam->GetA();
  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  // EnergyLoss Tables
  string BeamName = NPL::ChangeNameToG4Standard(myBeam->GetName());
  BeamTarget = NPL::EnergyLoss(BeamName+"_"+TargetMaterial+".G4table","G4Table",10000);
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
    Theta12  = Proton1.Angle(Proton2)/deg;

    // reject event that make no physical sense
    /*if(deltaPhi<170 && sumTheta<80){
      return;
      }
    */
    // computing minimum distance of the two lines
    TVector3 Vertex;
    TVector3 delta;
    Distance = NPL::MinimumDistance(InnerPos1,OuterPos1,InnerPos2,OuterPos2,Vertex,delta);
    VertexX=Vertex.X();
    VertexY=Vertex.Y();
    VertexZ=Vertex.Z();
    deltaX=delta.X();
    deltaY=delta.Y();
    deltaZ=delta.Z();
  }

    //double thickness_before = 0;
    //double EA_vertex = BeamTarget.Slow(InitialBeamEnergy,thickness_before,0);

    // setting up Lorentz Vector from measured trajectories and energies
    //LV_A.SetVect(PA); LV_p1.SetE(EA_vertex); 
    //LV_p1.SetVect(P1); LV_p1.SetE(E1); 
    //LV_p2.SetVect(P2); LV_p1.SetE(E2); 

    // computing Ex from Missing Mass
    //double EB = LV_A.E() + LV_T.E() - LV_p1.E() - LV_p2.E();   
    //TVector3 PB = LV_A.Vect() + LV_p1.Vect() - LV_p2.Vect();   
    //Ex = TMath::Sqrt( EB*EB - PB.Mag2() ) - myQFS->GetNucleusB()->Mass();
}


////////////////////////////////////////////////////////////////////////////////
TVector3 Analysis::InterpolateInPlaneZ(TVector3 V0, TVector3 V1, double Zproj){
    TVector3 Vproj(-999,-999,-999);
    double t = (Zproj - V1.Z()) / (V1.Z()-V0.Z());
    double Xproj= V1.X() + (V1.X()-V0.X()) * t;
    double Yproj= V1.Y() + (V1.Y()-V0.Y()) * t; 
    Vproj.SetXYZ(Xproj,Yproj,Zproj);
    return Vproj;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab,"ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("Theta12",&Theta12,"Theta12/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("VertexX",&VertexX,"VertexX/D");
  RootOutput::getInstance()->GetTree()->Branch("VertexY",&VertexY,"VertexY/D");
  RootOutput::getInstance()->GetTree()->Branch("VertexZ",&VertexZ,"VertexZ/D");
  RootOutput::getInstance()->GetTree()->Branch("deltaX",&deltaX,"deltaX/D");
  RootOutput::getInstance()->GetTree()->Branch("deltaY",&deltaY,"deltaY/D");
  RootOutput::getInstance()->GetTree()->Branch("deltaZ",&deltaZ,"deltaZ/D");

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
  Theta12 = -1000;
  ThetaCM = -1000;
  VertexX=-1000;
  VertexY=-1000;
  VertexZ=-1000;
  deltaX=-1000;
  deltaY=-1000;
  deltaZ=-1000;
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

