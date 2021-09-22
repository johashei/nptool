/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : Octobre 2017                                             *
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
  //G4ThreeVector  pdirection   = PrimaryTrack->GetMomentum().unit();

  m_length = 300*cm;
  
  //G4double energy = PrimaryTrack->GetKineticEnergy();
  G4double speed = PrimaryTrack->GetVelocity();
  G4double time  = PrimaryTrack->GetGlobalTime()+m_length/speed;

  
  //Get Track info
  G4ThreeVector localDir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();
  /*
  G4ThreeVector globalDir = PrimaryTrack->GetMomentumDirection();
  G4ThreeVector globalPosition = PrimaryTrack->GetPosition();
  G4ThreeVector globalMomentum = PrimaryTrack->GetMomentum();
  */
  G4ThreeVector newDir;
  G4ThreeVector newPosition;
  G4ThreeVector newMomentum;
  
  newDir = localDir.rotateY(-40*deg);/*
  double newTheta = localDir.getTheta() - 2*deg;
  if (newTheta < 0) newTheta += 2*M_PI;
  newDir.setTheta( newTheta );*/
  newPosition = localPosition + m_length * newDir;
  newMomentum = localMomentum.getR() * newDir;

  //cout << " New Theta: " << newTheta << endl;
  //cout << " Local Pos: " << Prt(localPosition) << endl;
  //cout << " Local Mom: " << Prt(localMomentum) << endl;
  //cout << " Local Dir: " << Prt(localDir) << endl;
  //cout << " New Dir: "      << Prt(newDir) << endl;
  //cout << " New Position: " << Prt(newPosition) << endl;
  //cout << " New Momentum: " << Prt(newMomentum) << endl;
  //cout << endl;
  
  /*
  newDir = globalDir;
  double newTheta = globalDir.getTheta() - 3*deg;
  if (newTheta < 0) newTheta += 2*M_PI;
  newDir.setTheta( newTheta );
  newPosition = globalPosition + m_length * newDir;
  newMomentum = globalMomentum.mag() * newDir;*/
  /*
  cout << " Global Dir: " << Prt(globalDir) << endl;
  cout << " New Theta: " << newTheta << endl;
  cout << " New Dir: "      << Prt(newDir) << endl;
  cout << endl;
  */
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
