/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2021                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  A quite Standard Geant4 EventAction class.                               *
 *  Call the Fill method of the output tree.                                 *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "SteppingAction.hh"
#include "NPOptionManager.h"
#include "RootOutput.h"

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(){
   m_record_track=NPOptionManager::getInstance()->GetRecordTrack();
   m_tree =  RootOutput::getInstance()->GetTree();
   // Attach track info class to m_tree
  }

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void 	SteppingAction::UserSteppingAction (const G4Step* /*step*/){
 if(m_record_track){
    //TrackRecording(step);
   }
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::TrackRecording(const G4Step* /*step*/){
// Clear the tracks class
// m_Tracks->Clear();
//  unsigned int size = traj->size();
  //    for(unsigned int i = 0 ; i < size ; i++){
       // FILL 
       // Particle name
       // Interaction points
       // Momentum
       // Kinetic energy
       // Volume Name
      //  }

  }
