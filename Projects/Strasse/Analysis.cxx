/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : july  2020                                               *
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
#include"NPPhysicalConstants.h"
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
  Catana = (TCatanaPhysics*)  m_DetectorManager -> GetDetector("Catana");
  // reaction properties
  m_QFS = new NPL::QFS();
  m_QFS->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  // reaction properties
  myBeam = new NPL::Beam();
  myBeam->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  InitialBeamEnergy = myBeam->GetEnergy()/* * myBeam->GetA()*/;
  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  // EnergyLoss Tables
  string BeamName = NPL::ChangeNameToG4Standard(myBeam->GetName());
  BeamTarget = NPL::EnergyLoss(BeamName+"_"+TargetMaterial+".G4table","G4Table",100);
  protonTarget = NPL::EnergyLoss("proton_"+TargetMaterial+".G4table","G4Table",100);
  protonAl = NPL::EnergyLoss("proton_Al.G4table","G4Table",100);
  protonSi = NPL::EnergyLoss("proton_Si.G4table","G4Table",100);
  LV_T.SetVectM(TVector3(0,0,0),NPUNITS::proton_mass_c2);
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

    deltaPhi = abs(Proton1.Phi()/deg-Proton2.Phi()/deg);
    sumTheta = Proton1.Theta()/deg+Proton2.Theta()/deg;
    Theta12  = Proton1.Angle(Proton2)/deg;

    // reject event that make no physical sense
    /*if(deltaPhi<170 && sumTheta<80){
      return;
      }
    */
    // computing minimum distance of the two lines
    TVector3 Vertex;
    TVector3 delta;
    Distance = NPL::MinimumDistanceTwoLines(InnerPos1,OuterPos1,InnerPos2,OuterPos2,Vertex,delta);
    VertexX=Vertex.X();
    VertexY=Vertex.Y();
    VertexZ=Vertex.Z();
    deltaX=delta.X();
    deltaY=delta.Y();
    deltaZ=delta.Z();
    
    // Look for associated Catana event
    double d1,d2;
    unsigned int i1,i2;
    i1 = Catana->FindClosestHitToLine(InnerPos1,OuterPos1,d1);
    i2 = Catana->FindClosestHitToLine(InnerPos2,OuterPos2,d2);
    if(i1!=i2){
      E1 = ReconstructProtonEnergy(Vertex,Proton1,Catana->Energy[i1]); 
      E2 = ReconstructProtonEnergy(Vertex,Proton2,Catana->Energy[i2]);
      double TA = BeamTarget.Slow(InitialBeamEnergy,abs(VertexZ-75),0);
      // setting up Lorentz Vector from measured trajectories and energies
      TVector3 PA(0,0,sqrt(TA*(TA+2*m_QFS->GetNucleusA()->Mass()))); // for like there is no BDC
      Proton1=E1*Proton1.Unit();
      Proton2=E2*Proton2.Unit();
      
      LV_A.SetVectM(PA,m_QFS->GetNucleusA()->Mass());
      double P1= sqrt(E1*(E1+2*NPUNITS::proton_mass_c2));
      double P2= sqrt(E2*(E2+2*NPUNITS::proton_mass_c2));

      LV_p1.SetVectM(Proton1.Unit()*P1,NPUNITS::proton_mass_c2); 
      LV_p2.SetVectM(Proton2.Unit()*P2,NPUNITS::proton_mass_c2); 

      // computing Ex from Missing Mass
      LV_B = LV_A + LV_T - LV_p1 - LV_p2;
      //LV_B = RC->GetParticleMomentum(2);
      Ex = LV_B.M() - m_QFS->GetNucleusB()->Mass();
    }
    
  }
}

////////////////////////////////////////////////////////////////////////////////
double Analysis::ReconstructProtonEnergy(const TVector3& x0, const TVector3& dir,const double& Ecatana){
    TVector3 Normal = TVector3(0,0,1);
    Normal.SetPhi(dir.Phi());
    double Theta = dir.Angle(Normal);  
    // Catana Al housing 
    double E = protonAl.EvaluateInitialEnergy(Ecatana,0.5*mm,Theta);
    // Strasse Chamber
    E = protonAl.EvaluateInitialEnergy(E,3*mm,Theta);
    // Outer Barrel
    E = protonSi.EvaluateInitialEnergy(E,300*micrometer,Theta);
    // Inner Barrel
    E = protonSi.EvaluateInitialEnergy(E,200*micrometer,Theta);
    // LH2 target
    static TVector3 x1;
    x1= x0+dir;
    TVector3 T1(0,30,0);
    TVector3 T2(0,30,1);
    T1.SetPhi(dir.Phi());
    T2.SetPhi(dir.Phi());
    TVector3 Vertex,delta;
    double d = NPL::MinimumDistanceTwoLines(x0,x1,T1,T2,Vertex,delta);
    E = protonTarget.EvaluateInitialEnergy(E,d,Theta);
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
  RootOutput::getInstance()->GetTree()->Branch("E1",&E1,"E1/D");
  RootOutput::getInstance()->GetTree()->Branch("E2",&E2,"E2/D");
  RootOutput::getInstance()->GetTree()->Branch("Theta12",&Theta12,"Theta12/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("VertexX",&VertexX,"VertexX/D");
  RootOutput::getInstance()->GetTree()->Branch("VertexY",&VertexY,"VertexY/D");
  RootOutput::getInstance()->GetTree()->Branch("VertexZ",&VertexZ,"VertexZ/D");
  RootOutput::getInstance()->GetTree()->Branch("deltaX",&deltaX,"deltaX/D");
  RootOutput::getInstance()->GetTree()->Branch("deltaY",&deltaY,"deltaY/D");
  RootOutput::getInstance()->GetTree()->Branch("deltaZ",&deltaZ,"deltaZ/D");
  RootOutput::getInstance()->GetTree()->Branch("deltaPhi",&deltaPhi,"deltaPhi/D");
  RootOutput::getInstance()->GetTree()->Branch("sumTheta",&sumTheta,"sumTheta/D");


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
  E1= -1000;
  E2 = -1000;
  Theta12 = -1000;
  ThetaCM = -1000;
  VertexX=-1000;
  VertexY=-1000;
  VertexZ=-1000;
  deltaX=-1000;
  deltaY=-1000;
  deltaZ=-1000;
  Distance=-1000;
  sumTheta=-1000;
  deltaPhi=-1000;
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

