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
  //std::cout << "\nDOIT" << std::endl;
  //std::cout << "FIELDMAP " << m_FieldMap << std::endl;
  
  if(!m_Initialized){
    m_Map = new SamuraiFieldMap();
    //cout << "before load map\n";
    m_Map->LoadMap(0, m_FieldMap, 10);//FIXME

    //Needed for the RungeKutta method
    m_Map->SetStepLimit( 10.*m / m_StepSize );
    m_Map->SetRmax( m_Rmax );

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
  
  G4ThreeVector localVel = localDir * speed;



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

  static int count = 0;
  count++;

  G4ThreeVector B = VtoG4( m_Map->InterpolateB(G4toV(localPosition/mm)) );
  //G4ThreeVector B (100.*tesla, 200.*tesla, 300.*tesla);
  //G4ThreeVector B = MagField(localPosition);


  bool debugging = false;
  if (debugging && localPosition.y() <= 400*mm && localPosition.y() >= -400*mm){
    vector <double> dummy = m_Map->InterpolateB(G4toV(localPosition/mm));
    for (auto x:dummy) cout << x << endl;
    cout << "loc Pos " << Cart(localPosition/millimeter) << endl;
    cout << count << " B = " << Cart(B/tesla) << "\t" << Prt(B/tesla) << endl;
  }
  
  double magB = B.mag() / tesla;

  
  G4ThreeVector newPosition;
  G4ThreeVector newDir;
  G4ThreeVector newMomentum;  

  int N_print = -1;
  
  if(count == 1 && count < N_print){
      cout << "\nIteration " << count << endl;
    
      if (magB == 0) cout << "//////////////// ZERO ///////////////////// "  << endl;
      cout<< "localPos = "<< Cart(localPosition/meter) << "\t" << Prt(localPosition/meter) << endl;
      cout<< "localDir = "<< Cart(localDir) << "\t" << Prt(localDir) << endl;
      cout<< "localMom = "<< Cart(localMomentum) << "\t" << Prt(localMomentum) << endl;
  }
    
  if(magB != 0){

    //unit conversion factors 
    double ConF_p = 1 / (joule * c_light / (m/s)) ; // MeV/c to kg*m/s (SI units)
    double ConF_m = 1 / ( joule * (c_light/(m/s)) * (c_light/(m/s)) ) ; // MeV/c^2 to kg (SI units)
  
    double charge = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;
    double mass = PrimaryTrack->GetParticleDefinition()->GetPDGMass() * ConF_m;

    //Calculate curvature radius and pulsation
    G4ThreeVector rho = -(localMomentum * ConF_p).cross(B/tesla);
    rho = rho / (charge * magB * magB) * m;//CLEANME
    G4ThreeVector omega = ( -(charge/mass) * B / tesla ) * hertz;


    double StepTime = m_StepSize / speed;
    
    double L_B = B.unit().dot(localDir) * m_StepSize;
    double L_perp = sqrt(m_StepSize*m_StepSize - L_B*L_B);
    
    
    //Calculate new position
    double angle = L_perp / rho.mag() * rad;
    // double angle = omega.mag() * StepTime;
    G4ThreeVector rho2 = rho;
    rho2.rotate(angle, omega);

    //motion along the B direction
    G4ThreeVector B_motion = L_B * B.unit();
    //G4ThreeVector B_motion = B.unit().dot(localVel) * StepTime * B.unit();

    //Setting new kinematic properties
    newPosition = localPosition - rho + rho2 + B_motion;
    newDir = localDir;
    newDir.rotate(angle, omega);
    newMomentum = localMomentum.mag() * newDir;

    
    if (count < N_print){
      cout << "\nIteration " << count << endl;
    
      if (magB == 0) cout << "//////////////// ZERO ///////////////////// "  << endl;
      cout << "mass(kg) = " << mass << " " << endl;
      cout << "charge (C) = " << charge << endl;
      cout << "speed (m/s) = " << speed / (m/s)  << endl;
      cout << "StepTime (s) = " << StepTime / s << endl;
      cout << "speed_B = " << B.unit().dot(localVel) / (m/s)  << endl;
      cout<< "B_motion = "<< Cart(B_motion/meter) << "\t" << Prt(B_motion/meter) << endl;
      cout<< "Newpos = "<< Cart(newPosition/meter) << "\t" << Prt(newPosition/meter) << endl;
      cout<< "NewDir = "<< Cart(newDir) << "\t" << Prt(newDir) << endl;
      cout<< "NewMom = "<< Cart(newMomentum) << "\t" << Prt(newMomentum) << endl;
      cout<< "B = "<< Cart(B/tesla) << "\t" << Prt(B/tesla) << endl;
      cout<< "omega = "<< Cart(omega/hertz) << "\t" << Prt(omega/hertz) << endl;
      cout<< "rho = "<< Cart(rho/meter) << "\t" << Prt(rho/meter) << endl;
      cout<< "rho2 = "<< Cart(rho2/meter) << "\t" << Prt(rho2/meter) << endl;
      cout << "angle " << angle/rad << endl;
      cout<< "ConF_p = "<< ConF_p << endl;
      cout<< "ConF_m = "<< ConF_m << endl;
    }
  }
  else{
    //Setting new kinematic properties with no magnetic field
    newPosition = localPosition + localDir * m_StepSize;
    newDir = localDir;
    newMomentum = localMomentum;
    
    if (count < N_print){
      cout << "\nIteration " << count << endl;
    
      if (magB == 0) cout << "//////////////// ZERO ///////////////////// " << endl;
      cout<< "Newpos = "<< Cart(newPosition/meter) << "\t" << Prt(newPosition/meter) << endl;
      cout<< "NewDir = "<< Cart(newDir) << "\t" << Prt(newDir) << endl;
      cout<< "NewMom = "<< Cart(newMomentum) << "\t" << Prt(newMomentum) << endl;
      cout<< "B = "<< Cart(B/tesla) << "\t" << Prt(B/tesla) << endl;
    }
  }


  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );//FIXME
  fastStep.ProposePrimaryTrackFinalTime( time );

  return;
  
}

/////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::RungeKuttaPropagation (const G4FastTrack& fastTrack,
							  G4FastStep& fastStep){
  //cout << "RKO " << fastTrack.GetPrimaryTrack()->GetKineticEnergy()/GeV << " \n";

  static int counter = 0;//debugging purposes
  static bool inside = false;//true if previous step is inside the volume
  static vector<TVector3> trajectory;
  static int count;//keeps track of the step reached inside trajectory

  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  
  //G4double energy = PrimaryTrack->GetKineticEnergy();
  G4double speed = PrimaryTrack->GetVelocity();
  //G4double time  = PrimaryTrack->GetGlobalTime()+m_StepSize/speed;
  
  //Get Track info
  G4ThreeVector localDir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();
  
  //G4ThreeVector localVel = localDir * speed;

  double charge = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;
  double ConF_p = 1 / (joule * c_light / (m/s)) ; // MeV/c to kg*m/s (SI units)

  if(PrimaryTrack->GetTouchable()->GetSolid()->Inside(localPosition) == kInside) cout << "INSIDE " << PrimaryTrack->GetTouchable()->GetSolid()->GetName() << " \n";
  else if(PrimaryTrack->GetTouchable()->GetSolid()->Inside(localPosition) == kSurface) cout << "SURFACE\n";
  else cout << "OUTSIDE\n";
  cout << "locPos " << Cart(localPosition) << endl;
  if (inside) cout << "inside is true\n";
  else cout << "inside is false" << endl;

  if (!inside){
    count = 2;//skip first two positions as they are the same as the current position
    double Brho = localMomentum.mag() * ConF_p / charge;
    TVector3 pos (localPosition.x(), localPosition.y(), localPosition.z());
    TVector3 dir (localDir.x(), localDir.y(), localDir.z());

    m_Map->SetTimeIntervalSize( m_StepSize / speed );

    trajectory.clear();
    //cout << "before propagate\n";
    trajectory = m_Map->Propagate(Brho, pos, dir);
    //cout << "\n Propagate \n";

    inside = true;
  }

  G4ThreeVector newPosition (trajectory[count].x(), trajectory[count].y(), trajectory[count].z());


  if (inside){
    G4VSolid* solid = PrimaryTrack->GetTouchable()->GetSolid();

    if (solid->Inside(newPosition) != kInside || count >= trajectory.size()-2){
      inside = false;
      //cout << "///////// NOT INSIDE //////////" << endl;
      G4ThreeVector toOut = solid->DistanceToOut(localPosition, localDir) * localDir;
      newPosition = localPosition + toOut;
      //if (solid->Inside(newPosition) == kSurface) cout << "I'm surfing baby\n";
      //else if (solid->Inside(newPosition) == kOutside) cout << "I'm outside baby\n";
      //else cout << "I'm inside baby\n";
      //cout << "newPos " << Cart(newPosition) << endl;
      counter++;
      cout << "\n" << counter << " inside = false\n";
    }
  }

/*
  if(counter > 0){
    //cout << "newP " << Cart(newPosition) << endl;
    //cout << "speed " << speed/(m/s) << endl;
    //G4ThreeVector B = VtoG4( m_Map->InterpolateB(G4toV(localPosition/mm)) );
    //cout << "B " << Cart(B) << endl;
    if(!inside) cout << "Size " << trajectory.size() << endl;
    if(!inside) for(unsigned int i=0; i<trajectory.size(); i++) cout << i << " " << trajectory[i].x() << " " << trajectory[i].y() << " " << trajectory[i].z() << endl;
  }
*/
  //cout << Cart(newPosition) << endl;

  G4ThreeVector newDir = (newPosition - localPosition).unit();
  G4ThreeVector newMomentum = newDir * localMomentum.mag(); 
  G4double time  = PrimaryTrack->GetGlobalTime()+(newPosition - localPosition).mag()/speed;

  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );//FIXME
  fastStep.ProposePrimaryTrackFinalTime( time );

  count++;
  
  return;
}

/////////////////////////////////////////////////////////////////////////
G4ThreeVector NPS::SamuraiFieldPropagation::MagField (G4ThreeVector pos){
  double x = pos.x()/meter;
  double y = pos.y()/meter;
  double z = pos.z()/meter;
  double a = 5;
  double c = 1;
  double Bx = 0.1 * c * exp( -x*x / a ) * tesla;
  double By = c * exp( -y*y / a ) * tesla;
  double Bz = 0.3 * c * exp( -z*z / a ) * tesla;

  //double By = - c * z * tesla;
  //double By;
  //if (z <= 0) By = -c * tesla;
  //else By = c * tesla;
  G4ThreeVector B (Bx,By,Bz);
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
