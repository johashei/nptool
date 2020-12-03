/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adriment Matta contact address: matta@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : Octover 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe eAGanil analysis project                             *
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
#include"NPPhysicalConstants.h"
#include"NPSystemOfUnits.h"
#include"NPPhysicalConstants.h"
using namespace NPUNITS;
#include"RootInput.h"
#include"RootOutput.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
   eAGanil= (TeAGanilPhysics*) m_DetectorManager->GetDetector("eAGanil");
   Inter = new TInteractionCoordinates();
   RootInput:: getInstance()->GetChain()->SetBranchAddress("InteractionCoordinates",&Inter);
   RootInput:: getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&Initial ); 
   RootInput:: getInstance()->GetChain()->SetBranchAddress("Run",&Run ); 
   RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex);
   RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab);
   RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab);
   RootOutput::getInstance()->GetTree()->Branch("Detected",&Detected);
   RootOutput::getInstance()->GetTree()->Branch("Resolution",&Resolution);
   RootOutput::getInstance()->GetTree()->Branch("Run",&Run);
   m_reaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
   Resolution={5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6};
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  Ex.clear(); ELab.clear(); ThetaLab.clear(); Detected.clear();
  // Get the generated values
  double Energy = Initial->GetKineticEnergy(0); 
  double Theta = Initial->GetParticleDirection(0).Angle(TVector3(0,0,1));
  unsigned int sizeR = Resolution.size();
 for(unsigned int i = 0 ; i < sizeR ; i++){
  
  //double E = Rand.Gaus(Energy,Energy*Resolution[i]);
  //double E = Rand.Gaus(Energy,Energy*1e-4);
//  double T = Rand.Gaus(Theta,Theta*Resolution[i]);
  //double ExO = m_reaction.ReconstructRelativistic(E,Theta);
//  double ExO = m_reaction.ReconstructRelativistic(Energy,T);

  double ExO = m_reaction.ReconstructRelativistic(Energy,Theta);
  Ex.push_back(ExO);
  ThetaLab.push_back(Theta/deg);
  //ELab.push_back(E);
  if(Inter->GetDetectedMultiplicity())
    Detected.push_back(1);
  else
    Detected.push_back(0);
  } 
 /* static double m2 = electron_mass_c2*electron_mass_c2;
  static double m = electron_mass_c2;
    // Momentum in MeV.c
    double q = sqrt((Energy+m)*(Energy+m)-m2);
    q = Rand.Gaus(q,q*1e-4);
    ELab = sqrt(q*q+m2)-m;
    */
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

