#ifndef SAMURAI_h
#define SAMURAI_h 1
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto  contact address: pilottoelia@gmail.com     *
 *                                                                           *
 * Creation Date  : septembre 2021                                           *
 * Last update    : septembre 2021                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  SAMURAI simulation                                  *
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
#include "G4VFastSimulationModel.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TSAMURAIData.h"
#include "NPInputParser.h"

class SAMURAI : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    SAMURAI() ;
    virtual ~SAMURAI() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddMagnet(G4ThreeVector POS, double Angle);
    // Spherical
    void AddMagnet(double R, double Theta, double Phi, double Angle);

    G4LogicalVolume* BuildMagnet();//FIXME
    G4LogicalVolume* BuildYoke();//FIXME
  
  private:
  
    G4LogicalVolume* m_Magnet;
    G4LogicalVolume* m_Yoke;
    
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetectorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;

    // Construct detector and inialise sensitive part.
    // Called After DetecorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Detector branch to the EventTree.
    // Called After DetecorConstruction::AddDetector Method
    //void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event) ;

  public:
    // Scorer
    //   Initialize all Scorer used by the MUST2Array
    //void InitializeScorers() ;

    // Set region were magnetic field is active:
    void SetPropagationRegion();

    //   Associated Scorer
    //G4MultiFunctionalDetector* m_SAMURAIScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TSAMURAIData* m_Event ; //FIXME

    // Region were magnetic field is active:
    G4Region* m_PropagationRegion;
    vector<G4VFastSimulationModel*> m_PropagationModel;
  

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private:
    // Geometry
    // Detector Coordinate
    double m_R;
    double m_Theta;
    double m_Phi;

    // Angle of Rotation
    double m_Angle;
    
    // Visualisation Attributes
    G4VisAttributes* m_VisMagnet;
    G4VisAttributes* m_VisYokes;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
