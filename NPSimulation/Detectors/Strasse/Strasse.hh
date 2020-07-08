#ifndef Strasse_h
#define Strasse_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny  contact address: flavigny@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Strasse simulation                                  *
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
#include "TStrasseData.h"
#include "NPInputParser.h"

class Strasse : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Strasse() ;
    virtual ~Strasse() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cylindrical coordinate
    void AddInnerDetector(double R,double Z,double Phi);  
    void AddOuterDetector(double R,double Z,double Phi);  

    G4LogicalVolume* BuildInnerDetector();
    G4LogicalVolume* BuildOuterDetector();
    G4LogicalVolume* BuildElectronic();
    G4LogicalVolume* BuildFrame();
    G4LogicalVolume* BuildChamber();

  private:
    G4LogicalVolume* m_InnerDetector;
    G4LogicalVolume* m_OuterDetector;
    G4LogicalVolume* m_Electronic;
    G4LogicalVolume* m_Frame;
    G4LogicalVolume* m_Chamber;


  private:
    //    Initialize material used in detector definition
    void InitializeMaterial();


    //   List of material
    G4Material* m_MaterialSilicon ;
    G4Material* m_MaterialAl      ;
    G4Material* m_MaterialVacuum  ;
    G4Material* m_MaterialPCB     ;



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
    G4MultiFunctionalDetector* m_InnerScorer ;
    G4MultiFunctionalDetector* m_OuterScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TStrasseData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_Inner_R; 
    vector<double>  m_Inner_Z;
    vector<double>  m_Inner_Phi; 
    vector<double>  m_Outer_R; 
    vector<double>  m_Outer_Z;
    vector<double>  m_Outer_Phi; 


    // Visualisation Attribute
    //G4VisAttributes* m_VisTrap;

    // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();


  private: // Visualisation
    // Dark Grey
    G4VisAttributes* SiliconVisAtt  ;
    // Green
    G4VisAttributes* PCBVisAtt;
    // Gold Yellow
    G4VisAttributes* PADVisAtt  ;
    // Light Grey
    G4VisAttributes* FrameVisAtt ;
    // Light Blue
    G4VisAttributes* GuardRingVisAtt ;


};
#endif
