#ifndef Samurai_h
#define Samurai_h 1
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Omar Nasr  contact address: omar.h.nasr@outlook.com      *
 *                                                                           *
 * Creation Date  : septembre 2021                                           *
 * Last update    : septembre 2021                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Samurai simulation                                  *
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
#include "G4VSolid.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "NPInputParser.h"

#include "SamuraiFieldPropagation.hh"
#include "TSamuraiIdealData.h"


class SamuraiFDC2 : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    SamuraiFDC2() ;
    virtual ~SamuraiFDC2() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
  
    // Cartesian fdc0
    void AddDetector(G4ThreeVector Mag_Pos, double Mag_Angle, G4ThreeVector Offset, double Off_Angle);
    // Spherical fdc0
    /*void AddFDC2(double R, double Theta, double Phi, double Angle, 
                  double r, double theta, double phi);*/
    //G4VSolid* BuildRYokeSolid();

    G4LogicalVolume* BuildFDC2();
  private:
  
    //Solids used for the magnet construction
    //G4VSolid* m_RYokeSolid;
    G4LogicalVolume* m_FDC2;

    // Visualisation Attributes
    G4VisAttributes* m_VisFDC2;

    
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
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAction
    void ReadSensitive(const G4Event* event) ;

  public:
    // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;


    //   Associated Scorer
    G4MultiFunctionalDetector* m_FDC2Scorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    //TSamuraiFDC2Data* m_Event;
    //////////////////////////////////////////////////////////////////

    TSamuraiIdealData* m_Event;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private:
    //Detector coordinates
    G4ThreeVector m_Pos;
    // Angle of Rotation
    double m_Angle;
  

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif







