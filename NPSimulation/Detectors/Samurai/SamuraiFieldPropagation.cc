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

//C++ libraries
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>//FIXME
#include <Randomize.hh>//FIXME
//G4 libraries
#include "G4Electron.hh"//FIXME
#include "G4Gamma.hh"//FIXME
#include "G4EmCalculator.hh"//FIXME
#include "G4VPhysicalVolume.hh"//FIXME
#include "G4IonTable.hh"//FIXME
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
//nptool libraries
#include "NPFunction.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSFunction.hh"
//other
#include "SamuraiFieldPropagation.hh"
#include "SamuraiFieldMap.h"
#include "RootOutput.h"
#include "TLorentzVector.h"

////////////////////////////////////////////////////////////////////////////////
NPS::SamuraiFieldPropagation::SamuraiFieldPropagation(G4String modelName, G4Region* envelope)
  : G4VFastSimulationModel(modelName, envelope) {
  //ReadConfiguration();
    //m_B = G4ThreeVector( 0. , 1., 0. ) * tesla;

    m_Map = NULL;
    m_Initialized = false;


    //ABLA = new G4AblaInterface();
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
/*
void NPS::SamuraiFieldPropagation::ReadConfiguration() {
  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());
  
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
  
}
*/
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
  //std::cout << "FIELDMAP " << m_FieldMap << std::endl;
  
  if(!m_Initialized){
    m_Map = new SamuraiFieldMap();
    //m_Map->LoadMap(m_Angle, "/local/pilotto/nptool/Projects/Samurai/" + m_FieldMap, 1);
    m_Map->LoadMap(0, m_FieldMap, 1);
    m_Initialized = true;
  }

  switch (m_Meth){
    case RungeKutta:
      RungeKuttaPropagation(fastTrack, fastStep);
      break;
    case EliaOmar:
      EliaOmarPropagation(fastTrack, fastStep);
      break;
    default:
      cout << endl;
      cout << "//////WARNING///////" << endl;
      cout << "In SamuraiFieldPropagation.cc:\n";
      cout << "Propagation method not defined - THIS MESSAGE SHOULD NEVER APPEAR\n";
      cout << endl;

  }
  
  return;
}

/////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::EliaOmarPropagation (const G4FastTrack& fastTrack,
							G4FastStep& fastStep){
  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  
  
  //G4double energy = PrimaryTrack->GetKineticEnergy();
  G4double speed = PrimaryTrack->GetVelocity();
  G4double time  = PrimaryTrack->GetGlobalTime()+m_StepSize/speed;

  
  //Get Track info
  G4ThreeVector localDir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();


  //Calculate the curvature radius
  
  //G4ThreeVector m_B = MagField(localPosition);
  /*
  G4ThreeVector m_B = VtoG4( m_Map->InterpolateB(G4toV(localPosition)) );
  double B = m_B.mag() / tesla;
  double q = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;

  //units conversion factor MeV/c to kg*m/s (SI units)
  double ConF = (1.e6 * e_SI ) / (c_light / (m/s)) ;
  
  G4ThreeVector rho = (localMomentum*ConF).cross(m_B/tesla);
  rho = rho / (q*B*B) * meter;

  
  //Calculate new position
  double angle = m_StepSize / rho.mag() * rad;
  if( -q*m_B.y() <= 0) angle *= -1; 
  G4ThreeVector rho2 = (-rho).rotateY(angle);
  newPosition = rho + rho2 + localPosition;

  
  //Calculate new momentum
  newDir = localDir.rotateY(angle);
  newMomentum = localMomentum.getR() * newDir;

  
  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );//FIXME
  fastStep.ProposePrimaryTrackFinalTime( time );


  static int count = 0;
  count++;
  
  if(newPosition.getZ() > 1750*mm){
    cout << setprecision(15) << "\nFINAL STEP" << endl;
    cout << setprecision(15) << "S (um) " << m_StepSize/um << endl;
    cout << setprecision(15) << "N " << count << endl;
    cout << setprecision(15) << "X " << newPosition.getX() << endl;
    cout << setprecision(15) << "Y " << newPosition.getY() << endl;
    cout << setprecision(15) << "Theta " << newDir.getTheta() << endl;
    }*/

  G4ThreeVector B = VtoG4( m_Map->InterpolateB(G4toV(localPosition)) );
  double magB = B.mag() / tesla;
  double charge = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;
  double mass = PrimaryTrack->GetParticleDefinition()->GetPDGMass() / kg;

  //units conversion factor MeV/c to kg*m/s (SI units)
  double ConF = (1.e6 * e_SI ) / (c_light / (m/s)) ;
  
  G4ThreeVector rho = -(localMomentum*ConF).cross(B/tesla);
  rho = rho / (charge*magB*magB) * meter;
  G4ThreeVector omega = ( -(charge/mass) * B / tesla ) / s;

  
  //Calculate new position
  double angle = m_StepSize / rho.mag() * rad;
  G4ThreeVector rho2 = rho;
  rho2.rotate(angle, omega);

  //Setting new kinematic properties
  G4ThreeVector newPosition = localPosition - rho + rho2;
  G4ThreeVector newDir = localDir;
  newDir.rotate(angle, omega);
  G4ThreeVector newMomentum = localMomentum.mag() * newDir;

  cout<< "Newpos = "<< Cart(newPosition) << endl;
  cout<< "NewDir = "<< Cart(newDir) << endl;
  cout<< "NewMom = "<< Cart(newMomentum) << endl;
  cout<< "B = "<< Cart(B) << endl;
  cout<< "omega = "<< Cart(omega) << endl;
  cout<< "rho = "<< Cart(rho) << endl;
  cout<< "ConF = "<< ConF << endl;
  
  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );//FIXME
  fastStep.ProposePrimaryTrackFinalTime( time );

  return;
  
}

/////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::RungeKuttaPropagation (const G4FastTrack& fastTrack,
							  G4FastStep& fastStep){
  cout << "Method not defined yet\n";
  return;
}

/////////////////////////////////////////////////////////////////////////
G4ThreeVector NPS::SamuraiFieldPropagation::MagField (G4ThreeVector pos){
  double x = pos.x()/meter;
  double z = pos.z()/meter;
  double a = 5;
  double c = 2;
  double By = c * exp( -x*x / a ) * tesla;
  //double By = - c * z * tesla;
  //double By;
  //if (z <= 0) By = -c * tesla;
  //else By = c * tesla;
  G4ThreeVector B (0,By,0);
  return B;
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
