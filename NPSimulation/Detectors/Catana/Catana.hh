#ifndef Catana_h
#define Catana_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Catana simulation                                   *
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
#include "TCatanaData.h"
#include "NPInputParser.h"

class Catana : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Catana() ;
    virtual ~Catana() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDummyDetector(double Z);
    void AddDetectorType1(double R, double Theta, double Phi);
    void AddDetectorType2(double R, double Theta, double Phi);
    void AddDetectorType3(double R, double Theta, double Phi);

    G4LogicalVolume* BuildDummyDetector();
    G4LogicalVolume* BuildDetectorType1();
    G4LogicalVolume* BuildDetectorType2();
    G4LogicalVolume* BuildDetectorType3();
  private:
    G4LogicalVolume* m_DummyDetector;
    G4LogicalVolume* m_DetectorType1;
    G4LogicalVolume* m_DetectorType2;
    G4LogicalVolume* m_DetectorType3;
    
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
    G4MultiFunctionalDetector* m_CatanaScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TCatanaData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_Z; 
    vector<double>  m_R1; 
    vector<double>  m_R2; 
    vector<double>  m_R3; 
    vector<double>  m_Theta1; 
    vector<double>  m_Theta2; 
    vector<double>  m_Theta3; 
    vector<double>  m_Phi1; 
    vector<double>  m_Phi2; 
    vector<double>  m_Phi3; 

    // Visualisation Attribute
    G4VisAttributes* m_VisCrystal;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
