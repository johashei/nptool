/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto, Omar Nasr                                  *
 * contact address: pilottoelia@gmail.com, omar.nasr@etu.unicaen.fr          *
 *                                                                           *
 * Creation Date  : September 2021                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Use to kill the beam track and replace it with the reaction product       *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "SamuraiFieldPropagation.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4IonTable.hh"
#include "G4UserLimits.hh"
#include "NPFunction.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "TLorentzVector.h"
#include "NPSFunction.hh"
#include <Randomize.hh>
#include <iostream>
#include <string>
#include "G4UserLimits.hh"

////////////////////////////////////////////////////////////////////////////////
NPS::SamuraiFieldPropagation::SamuraiFieldPropagation(G4String modelName, G4Region* envelope)
  : G4VFastSimulationModel(modelName, envelope) {
    ReadConfiguration();
    // m_shoot=false;
    //m_rand=0;
    //m_Z=0;

    ABLA = new G4AblaInterface();
  }

////////////////////////////////////////////////////////////////////////////////
NPS::SamuraiFieldPropagation::SamuraiFieldPropagation(G4String modelName)
  : G4VFastSimulationModel(modelName) {}

////////////////////////////////////////////////////////////////////////////////
NPS::SamuraiFieldPropagation::~SamuraiFieldPropagation() {}

////////////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::AttachReactionConditions() {
  // Reasssigned the branch address
  /*
  if (RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
    RootOutput::getInstance()->GetTree()->SetBranchAddress(
		  "ReactionConditions", &m_ReactionConditions);*/
}

////////////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::ReadConfiguration() {
  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());
  /*
  if(input.GetAllBlocksWithToken("TwoBodyReaction").size()>0) m_ReactionType =TwoBody;
  else if(input.GetAllBlocksWithToken("QFSReaction").size()>0) m_ReactionType =QFS;
  else if(input.GetAllBlocksWithToken("FusionReaction").size()>0) m_ReactionType =Fusion;

  // Two body
  if (m_ReactionType==TwoBody ) {
    m_Reaction.ReadConfigurationFile(input);
    m_BeamName = NPL::ChangeNameToG4Standard(m_Reaction.GetParticle1()->GetName());
    if(m_Reaction.GetParticle3()->GetName() != ""){
      m_active = true;
      m_ReactionConditions = new TReactionConditions();
      AttachReactionConditions();
      if (!RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
        RootOutput::getInstance()->GetTree()->Branch(
            "ReactionConditions", "TReactionConditions", &m_ReactionConditions);
    }
  } 

  // QFS
  else  if (m_ReactionType==QFS) {
    m_QFS.ReadConfigurationFile(input);
    m_BeamName = NPL::ChangeNameToG4Standard(m_QFS.GetParticleA()->GetName());
    m_active = true;
    m_ReactionConditions = new TReactionConditions();
    AttachReactionConditions();
    if (!RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
      RootOutput::getInstance()->GetTree()->Branch(
          "ReactionConditions", "TReactionConditions", &m_ReactionConditions);
  }

  // Fusion
  else  if (m_ReactionType==Fusion) {
    vector<InputBlock*> blocks=input.GetAllBlocksWithToken("FusionReaction");
    m_BeamName = NPL::ChangeNameToG4Standard(blocks[0]->GetString("Beam"));
    m_TargetNuclei = blocks[0]->GetString("Target");
    m_FusionProduct = blocks[0]->GetString("Product");
    m_FusionExcitation = blocks[0]->GetDouble("ExcitationEnergy","MeV");
    m_active = true;
    // not used
    m_ReactionConditions = new TReactionConditions();
  }
  else {
    m_active = false;
  }
  */
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::SamuraiFieldPropagation::IsApplicable(const G4ParticleDefinition& particleType) {
  if (particleType.GetPDGCharge() == 0) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::SamuraiFieldPropagation::ModelTrigger(const G4FastTrack& fastTrack) {
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::DoIt(const G4FastTrack& fastTrack,
    G4FastStep&        fastStep) {
  //std::cout << "DOIT" << std::endl;

  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();

  m_length = 100*cm;
  m_B = G4ThreeVector( 0. , 1. , 0. ) * tesla;
  
  //G4double energy = PrimaryTrack->GetKineticEnergy();
  G4double speed = PrimaryTrack->GetVelocity();
  G4double time  = PrimaryTrack->GetGlobalTime()+m_length/speed;

  
  //Get Track info
  G4ThreeVector localDir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();
  G4ThreeVector newDir;
  G4ThreeVector newPosition;
  G4ThreeVector newMomentum;



  
  double B = m_B.mag() / tesla;
  double q = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;
  double ConF = (1.e6 * (1.602176634e-19) ) / 2.99792458e8 ;
  G4ThreeVector rho = (localMomentum*ConF).cross(m_B/tesla);
 
  //cout << "Rho : " << Cart(rho) << "\t" << rho.mag()<< endl;
  rho = rho / (q*B*B) * meter;
  /* cout << "Rho : " << Cart(rho) << "\t" << (rho/m).mag()<< endl;

  cout << "q : " << q << endl;
  cout << "B mag : " << B << endl;
  cout << "B vec : " << Cart(m_B/tesla) << endl;
  cout << "Mom : " << Cart(localMomentum * ConF*1e20) << "\t" << localMomentum.mag() << endl;
  cout << "Pos : " << Cart(localPosition/cm) << endl;
  cout << "conversion: " << ConF << endl;
  cout << "tesla " << tesla << endl;
  cout << "meter " << meter << endl;*/


  double angle = m_length / rho.mag() * rad;
  cout << "Rho: " << Prt(-rho) << endl;
  G4ThreeVector rho2 = (-rho).rotateY(-angle);
  newPosition = rho + rho2 + localPosition;
  cout << "newPos: " << Cart(newPosition) << endl;
  cout << "locPos: " << Cart(localPosition) << endl;
  cout << "Rho: " << Prt(-rho) << endl;
  cout << "Rho2: " << Prt(rho2) << endl;
  


  cout << "localDir: " << Prt(localDir) << endl;

  newDir = localDir.rotateY(-angle);
  newMomentum = localMomentum.getR() * newDir;
  
  cout << "localDir: " << Prt(localDir) << endl;
  cout << "newDir: " << Prt(newDir) << endl;
  cout << "localMom: " << Prt(localMomentum) << endl;
  cout << "newMom: " << Prt(newMomentum) << endl;



  cout << "angle " << angle << endl;



  
  /*
  newDir = localDir.rotateY(-1*deg);
  newPosition = localPosition + m_length * newDir;
  newMomentum = localMomentum.getR() * newDir;
  */
  //cout << " New Theta: " << newTheta << endl;
  //cout << " Local Pos: " << Prt(localPosition) << endl;
  //cout << " Local Mom: " << Prt(localMomentum) << endl;
  //cout << " Local Dir: " << Prt(localDir) << endl;
  //cout << " New Dir: "      << Prt(newDir) << endl;
  //cout << " New Position: " << Prt(newPosition) << endl;
  //cout << " New Momentum: " << Prt(newMomentum) << endl;
  //cout << endl;

  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );//FIXME
  fastStep.ProposePrimaryTrackFinalTime( time );

  return;
}
  

  
  
  
    /*; 
  m_shoot=false;

  // Set the end of the step conditions
  fastStep.SetPrimaryTrackFinalKineticEnergyAndDirection(0, pdirection);
  fastStep.SetPrimaryTrackFinalPosition(worldPosition);
  fastStep.SetTotalEnergyDeposited(0);
  fastStep.SetPrimaryTrackFinalTime(time);// FIXME
  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);

}

    */
