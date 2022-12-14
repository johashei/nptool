#ifndef Sweeper_h
#define Sweeper_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: B. Monteagudo  contact address: monteagu@frib.msu.edu                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Sweeper simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TSweeperData.h"
#include "NPInputParser.h"
//#include "MagField.hh"

class Sweeper : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Sweeper() ;
    virtual ~Sweeper() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:

    // Spherical
    void AddDetector(G4ThreeVector POS,double Theta,double Brho, double* Dist);  
  
    //Logical Volumes
    G4LogicalVolume* BuildMotherVolume();
    G4LogicalVolume* BuildSweeper(double theta);
    G4LogicalVolume* BuildSweeperMagField(double theta);
    G4LogicalVolume* BuildCRDC();
    G4LogicalVolume* BuildIonChamber();
    G4LogicalVolume* BuildThinScint();
    G4LogicalVolume* BuildOldHodo();
    G4LogicalVolume* BuildNewHodo();

    //Magnetic Field
  void SetSweeperField(bool kMap, double bfield);
  
  private:
    //Logical volumes
    G4LogicalVolume* m_MotherLog;
    G4LogicalVolume* m_SweeperLog;
    G4LogicalVolume* m_SweeperMagFieldLog;
    G4LogicalVolume* m_CRDCLog;
    G4LogicalVolume* m_IonChamberLog;
    G4LogicalVolume* m_ThinScintLog;
    G4LogicalVolume* m_OldHodoLog;
    G4LogicalVolume* m_NewHodoLog;
  
  
    //Physical Volumes
    G4VPhysicalVolume *m_SweeperPhys;
    G4VPhysicalVolume *m_CRDCPhys;

    //Magnetic Field
    G4MagneticField *fSweeperMagField;
  
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetextorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;

    // Construct detector and inialise sensitive part.
    // Called After DetecorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Detector branch to the EventTree.
    // Called After DetecorConstruction::AddDetector Method
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event) ;

  public:   // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_SweeperScorer;
  
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TSweeperData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_R; 
    vector<double>  m_Theta;
    vector<double>  m_Phi;
    vector<double>  m_Brho; 
    vector<G4ThreeVector> m_Pos;
    //Detectors distances
    vector<double> m_DistToExit;
    vector<double> m_DistToDC1;
    vector<double> m_DistToDC2;
    vector<double> m_DistToIC;
    vector<double> m_DistToHodo;
  
   
    // Visualisation Attributes
    G4VisAttributes* m_VisCRDC;
    G4VisAttributes* m_VisSweeper;
    G4VisAttributes* m_VisIonChamber;
    G4VisAttributes* m_VisThinScint;
    G4VisAttributes* m_VisHodo;

  
  
  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
