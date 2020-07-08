/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F.Flavigny  contact address: flavigny@lpccaen.in2p3.fr   *
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

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4TwoVector.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// NPTool header
#include "Strasse.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Strasse_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoEnergy = 0.015*MeV ;

  ////////////////////
  // Inner Detector //
  ////////////////////
  // Wafer parameter
  double Inner_Wafer_Length=100*mm;
  double Inner_Wafer_Width=50*mm;
  double Inner_Wafer_Thickness=300*micrometer;
  double Inner_Wafer_AlThickness=0.4*micrometer;
  double Inner_Wafer_PADExternal=1*cm;
  double Inner_Wafer_PADInternal=1*mm;
  double Inner_Wafer_GuardRing=0.5*mm;

  // PCB parameter
  double Inner_PCB_PortWidth=1*cm;
  double Inner_PCB_StarboardWidth=2*mm;
  double Inner_PCB_BevelAngle= 60*deg;
  double Inner_PCB_UpstreamWidth=1*cm;
  double Inner_PCB_DownstreamWidth=2*mm;
  double Inner_PCB_MidWidth=2*mm;
  double Inner_PCB_Thickness=3*mm;

  ////////////////////
  // Outer Detector //
  ////////////////////
  // Wafer parameter
  double Outer_Wafer_Length=150*mm;
  double Outer_Wafer_Width=75*mm;
  double Outer_Wafer_Thickness=300*micrometer;
  double Outer_Wafer_AlThickness=0.4*micrometer;
  double Outer_Wafer_PADExternal=1*cm;
  double Outer_Wafer_PADInternal=1*mm;
  double Outer_Wafer_GuardRing=0.5*mm;

  // PCB parameter
  double Outer_PCB_PortWidth=1*cm;
  double Outer_PCB_StarboardWidth=2*mm;
  double Outer_PCB_BevelAngle= 60*deg;
  double Outer_PCB_UpstreamWidth=1*cm;
  double Outer_PCB_DownstreamWidth=2*mm;
  double Outer_PCB_MidWidth=2*mm;
  double Outer_PCB_Thickness=3*mm;

}

using namespace Strasse_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Strasse Specific Method
Strasse::Strasse(){
  InitializeMaterial();
  m_Event = new TStrasseData() ;
  m_InnerScorer = 0;
  m_OuterScorer = 0;
  m_InnerDetector=0;
  m_OuterDetector=0;
  m_Chamber=0;
  m_Frame=0;
  m_Electronic=0;
  // Dark Grey
  SiliconVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
  // Green
  PCBVisAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2)) ;
  // Gold Yellow
  PADVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.2)) ;
  // Light Grey
  FrameVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
  // Light Blue
  GuardRingVisAtt = new G4VisAttributes(G4Colour(0, 0, 0,0.5)) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Strasse::~Strasse(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::AddInnerDetector(double  R, double  Z, double  Phi){
  m_Inner_R.push_back(R);
  m_Inner_Z.push_back(Z);
  m_Inner_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::AddOuterDetector(double  R, double  Z, double  Phi){
  m_Outer_R.push_back(R);
  m_Outer_Z.push_back(Z);
  m_Outer_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildInnerDetector(){
  if(!m_InnerDetector){
    // Compute the needed full length of the PCB
    // along beam axis
    double Inner_PCB_Length= 2*Inner_Wafer_Length
      +Inner_PCB_UpstreamWidth
      +Inner_PCB_MidWidth
      +Inner_PCB_DownstreamWidth;

    // perpendicular to beam axis
    double Inner_PCB_Width= Inner_Wafer_Width
      +Inner_PCB_StarboardWidth
      +Inner_PCB_PortWidth;


    vector<G4TwoVector> PCBCrossSection;
    double l1 = Inner_PCB_Thickness*0.5/tan(Inner_PCB_BevelAngle);

    PCBCrossSection.push_back(G4TwoVector(Inner_PCB_Width*0.5-l1,-Inner_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(Inner_PCB_Width*0.5,Inner_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Inner_PCB_Width*0.5-l1,Inner_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Inner_PCB_Width*0.5,-Inner_PCB_Thickness*0.5));

    G4ExtrudedSolid* PCBFull =
      new G4ExtrudedSolid("PCBFull",
          PCBCrossSection,
          Inner_PCB_Length*0.5,// half length
          G4TwoVector(0,0),1,// offset, scale
          G4TwoVector(0,0),1);// offset, scale

    // Master Volume that encompass everything else
    m_InnerDetector =
      new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector", 0, 0, 0);
    m_InnerDetector->SetVisAttributes(G4VisAttributes::Invisible);

    // Build the PCB
    // Calculate the hole shift within the PCB
    double Width_Shift= -0.5*Inner_PCB_Width + 0.5*Inner_Wafer_Width // Flush to border
      +Inner_PCB_PortWidth; // add the port side shift

    double Length_Shift1 = -0.5*Inner_PCB_Length + 0.5*Inner_Wafer_Length // Flush to border
      + Inner_PCB_UpstreamWidth;// add Upstream side shift

    double Length_Shift2 = Length_Shift1 // overlap detector 1
      + Inner_Wafer_Length // at opposing edge
      + Inner_PCB_MidWidth; // after mid width

    G4ThreeVector HoleShift1 = G4ThreeVector(Width_Shift, 0, Length_Shift1);
    G4ThreeVector HoleShift2 = G4ThreeVector(Width_Shift, 0, Length_Shift2);

    G4Box*  HoleShape = new G4Box("HoleShape",
        Inner_Wafer_Width*0.5,
        Inner_PCB_Thickness*0.5+0.1*mm,
        Inner_Wafer_Length*0.5);

    // Substracting the hole Shape from the Stock PCB
    G4SubtractionSolid* PCB_1 = new G4SubtractionSolid("PCB_1", PCBFull, HoleShape,
        new G4RotationMatrix,HoleShift1);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_1, HoleShape,
        new G4RotationMatrix,HoleShift2);

    // Sub Volume PCB
    G4LogicalVolume* logicPCB =
      new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(PCBVisAtt);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0),
        logicPCB,"Strasse_Inner_PCB",m_InnerDetector,
        false,0);

    // Sub volume Wafer
    G4Box*  WaferShape = new G4Box("WaferShape",
        Inner_Wafer_Width*0.5,
        Inner_Wafer_Thickness*0.5+Inner_Wafer_AlThickness,
        Inner_Wafer_Length*0.5);

    G4LogicalVolume* logicWafer1 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer1", 0, 0, 0);
    logicWafer1->SetVisAttributes(GuardRingVisAtt);

    G4LogicalVolume* logicWafer2 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer2", 0, 0, 0);
    logicWafer2->SetVisAttributes(GuardRingVisAtt);


    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0.5*Inner_Wafer_Thickness
          +Inner_Wafer_AlThickness
          -0.5*Inner_PCB_Thickness,0)// flush the wafer to the pcb on one side
        +HoleShift1, // Shift wafer in the hole 
        logicWafer1,"Strasse_Inner_Wafer1",m_InnerDetector,
        false,0);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0.5*Inner_Wafer_Thickness
          +Inner_Wafer_AlThickness
          -0.5*Inner_PCB_Thickness,0)// flush the wafer to the pcb on one side
        +HoleShift2, // Shift wafer in the hole 
        logicWafer2,"Strasse_Inner_Wafer2",m_InnerDetector,
        false,0);

    // Sub volume Active Wafer
    G4Box*  ActiveWaferShape = new G4Box("ActiveWaferShape",
        Inner_Wafer_Width*0.5-Inner_Wafer_GuardRing,
        Inner_Wafer_Thickness*0.5,
        0.5*(Inner_Wafer_Length-Inner_Wafer_PADExternal-Inner_Wafer_PADInternal-Inner_Wafer_GuardRing));

    G4LogicalVolume* logicActiveWafer =
      new G4LogicalVolume(ActiveWaferShape,m_MaterialSilicon,"logicActiveWafer", 0, 0, 0);
    logicActiveWafer->SetVisAttributes(SiliconVisAtt);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0.5*(Inner_Wafer_PADExternal-Inner_Wafer_PADInternal)), // assymetric pading for bounding
        logicActiveWafer,"Strasse_Inner_ActiveWafer1",logicWafer1,
        false,1);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,-0.5*(Inner_Wafer_PADExternal-Inner_Wafer_PADInternal)), // assymetric pading for bounding
        logicActiveWafer,"Strasse_Inner_ActiveWafer2",logicWafer2,
        false,1);
  }
  return m_InnerDetector;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildOuterDetector(){
  if(!m_OuterDetector){
    // Compute the needed full length of the PCB
    // along beam axis
    double Outer_PCB_Length= 2*Outer_Wafer_Length
      +Outer_PCB_UpstreamWidth
      +Outer_PCB_MidWidth
      +Outer_PCB_DownstreamWidth;

    // perpendicular to beam axis
    double Outer_PCB_Width= Outer_Wafer_Width
      +Outer_PCB_StarboardWidth
      +Outer_PCB_PortWidth;


    vector<G4TwoVector> PCBCrossSection;
    double l1 = Outer_PCB_Thickness*0.5/tan(Outer_PCB_BevelAngle);

    PCBCrossSection.push_back(G4TwoVector(Outer_PCB_Width*0.5-l1,-Outer_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(Outer_PCB_Width*0.5,Outer_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Outer_PCB_Width*0.5-l1,Outer_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Outer_PCB_Width*0.5,-Outer_PCB_Thickness*0.5));

    G4ExtrudedSolid* PCBFull =
      new G4ExtrudedSolid("PCBFull",
          PCBCrossSection,
          Outer_PCB_Length*0.5,// half length
          G4TwoVector(0,0),1,// offset, scale
          G4TwoVector(0,0),1);// offset, scale

    // Master Volume that encompass everything else
    m_OuterDetector =
      new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector", 0, 0, 0);
    m_OuterDetector->SetVisAttributes(G4VisAttributes::Invisible);

    // Build the PCB
    // Calculate the hole shift within the PCB
    double Width_Shift= -0.5*Outer_PCB_Width + 0.5*Outer_Wafer_Width // Flush to border
      +Outer_PCB_PortWidth; // add the port side shift

    double Length_Shift1 = -0.5*Outer_PCB_Length + 0.5*Outer_Wafer_Length // Flush to border
      + Outer_PCB_UpstreamWidth;// add Upstream side shift

    double Length_Shift2 = Length_Shift1 // overlap detector 1
      + Outer_Wafer_Length // at opposing edge
      + Outer_PCB_MidWidth; // after mid width

    G4ThreeVector HoleShift1 = G4ThreeVector(Width_Shift, 0, Length_Shift1);
    G4ThreeVector HoleShift2 = G4ThreeVector(Width_Shift, 0, Length_Shift2);

    G4Box*  HoleShape = new G4Box("HoleShape",
        Outer_Wafer_Width*0.5,
        Outer_PCB_Thickness*0.5+0.1*mm,
        Outer_Wafer_Length*0.5);

    // Substracting the hole Shape from the Stock PCB
    G4SubtractionSolid* PCB_1 = new G4SubtractionSolid("PCB_1", PCBFull, HoleShape,
        new G4RotationMatrix,HoleShift1);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_1, HoleShape,
        new G4RotationMatrix,HoleShift2);

    // Sub Volume PCB
    G4LogicalVolume* logicPCB =
      new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(PCBVisAtt);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0),
        logicPCB,"Strasse_Outer_PCB",m_OuterDetector,
        false,0);

    // Sub volume Wafer
    G4Box*  WaferShape = new G4Box("WaferShape",
        Outer_Wafer_Width*0.5,
        Outer_Wafer_Thickness*0.5+Outer_Wafer_AlThickness,
        Outer_Wafer_Length*0.5);

    G4LogicalVolume* logicWafer1 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer1", 0, 0, 0);
    logicWafer1->SetVisAttributes(GuardRingVisAtt);

    G4LogicalVolume* logicWafer2 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer2", 0, 0, 0);
    logicWafer2->SetVisAttributes(GuardRingVisAtt);


    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0.5*Outer_Wafer_Thickness
          +Outer_Wafer_AlThickness
          -0.5*Outer_PCB_Thickness,0)// flush the wafer to the pcb on one side
        +HoleShift1, // Shift wafer in the hole 
        logicWafer1,"Strasse_Outer_Wafer1",m_OuterDetector,
        false,0);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0.5*Outer_Wafer_Thickness
          +Outer_Wafer_AlThickness
          -0.5*Outer_PCB_Thickness,0)// flush the wafer to the pcb on one side
        +HoleShift2, // Shift wafer in the hole 
        logicWafer2,"Strasse_Outer_Wafer2",m_OuterDetector,
        false,0);

    // Sub volume Active Wafer
    G4Box*  ActiveWaferShape = new G4Box("ActiveWaferShape",
        Outer_Wafer_Width*0.5-Outer_Wafer_GuardRing,
        Outer_Wafer_Thickness*0.5,
        0.5*(Outer_Wafer_Length-Outer_Wafer_PADExternal-Outer_Wafer_PADInternal-Outer_Wafer_GuardRing));

    G4LogicalVolume* logicActiveWafer =
      new G4LogicalVolume(ActiveWaferShape,m_MaterialSilicon,"logicActiveWafer", 0, 0, 0);
    logicActiveWafer->SetVisAttributes(SiliconVisAtt);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0.5*(Outer_Wafer_PADExternal-Outer_Wafer_PADInternal)), // assymetric pading for bounding
        logicActiveWafer,"Strasse_Outer_ActiveWafer1",logicWafer1,
        false,1);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,-0.5*(Outer_Wafer_PADExternal-Outer_Wafer_PADInternal)), // assymetric pading for bounding
        logicActiveWafer,"Strasse_Outer_ActiveWafer2",logicWafer2,
        false,1);
  }
  return m_OuterDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Strasse::ReadConfiguration(NPL::InputParser parser){
  // Info block
  vector<NPL::InputBlock*> blocks_info = parser.GetAllBlocksWithTokenAndValue("Strasse","Info");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_info.size() << " info block founds " << endl; 

  if(blocks_info.size()>1){
    cout << "ERROR: can only accepte one info block, " << blocks_info.size() << " info block founds." << endl; 
    exit(1); 
  }

  vector<string> info = {
    "Inner_Wafer_Length",         
    "Inner_Wafer_Width",          
    "Inner_Wafer_Thickness",     
    "Inner_Wafer_AlThickness",    
    "Inner_Wafer_PADExternal",    
    "Inner_Wafer_PADInternal",  
    "Inner_Wafer_GuardRing",    
    "Inner_PCB_PortWidth",      
    "Inner_PCB_StarboardWidth", 
    "Inner_PCB_BevelAngle",     
    "Inner_PCB_UpstreamWidth",  
    "Inner_PCB_DownstreamWidth",
    "Inner_PCB_MidWidth",       
    "Inner_PCB_Thickness",      
    "Outer_Wafer_Length",       
    "Outer_Wafer_Width",        
    "Outer_Wafer_Thickness",    
    "Outer_Wafer_AlThickness",  
    "Outer_Wafer_PADExternal",  
    "Outer_Wafer_PADInternal",  
    "Outer_Wafer_GuardRing",    
    "Outer_PCB_PortWidth",      
    "Outer_PCB_StarboardWidth", 
    "Outer_PCB_BevelAngle",     
    "Outer_PCB_UpstreamWidth",  
    "Outer_PCB_DownstreamWidth",
    "Outer_PCB_MidWidth",       
    "Outer_PCB_Thickness",      
  };

  if(blocks_info[0]->HasTokenList(info)){
    cout << endl << "////  Strasse info block" <<  endl;
    Inner_Wafer_Length = blocks_info[0]->GetDouble("Inner_Wafer_Length","mm");
    Inner_Wafer_Width = blocks_info[0]->GetDouble("Inner_Wafer_Width","mm");          
    Inner_Wafer_Thickness = blocks_info[0]->GetDouble("Inner_Wafer_Thickness","micrometer");      
    Inner_Wafer_AlThickness = blocks_info[0]->GetDouble("Inner_Wafer_AlThickness","micrometer");     
    Inner_Wafer_PADExternal = blocks_info[0]->GetDouble("Inner_Wafer_PADExternal","mm");     
    Inner_Wafer_PADInternal = blocks_info[0]->GetDouble("Inner_Wafer_PADInternal","mm");   
    Inner_Wafer_GuardRing = blocks_info[0]->GetDouble("Inner_Wafer_GuardRing","mm");     
    Inner_PCB_PortWidth = blocks_info[0]->GetDouble("Inner_PCB_PortWidth","mm");       
    Inner_PCB_StarboardWidth = blocks_info[0]->GetDouble("Inner_PCB_StarboardWidth","mm");  
    Inner_PCB_BevelAngle = blocks_info[0]->GetDouble("Inner_PCB_BevelAngle","mm");      
    Inner_PCB_UpstreamWidth = blocks_info[0]->GetDouble("Inner_PCB_UpstreamWidth","mm");   
    Inner_PCB_DownstreamWidth = blocks_info[0]->GetDouble("Inner_PCB_DownstreamWidth","mm"); 
    Inner_PCB_MidWidth = blocks_info[0]->GetDouble("Inner_PCB_MidWidth","mm");        
    Inner_PCB_Thickness = blocks_info[0]->GetDouble("Inner_PCB_Thickness","mm");       
    Outer_Wafer_Length = blocks_info[0]->GetDouble("Outer_Wafer_Length","mm");        
    Outer_Wafer_Width = blocks_info[0]->GetDouble("Outer_Wafer_Width","mm");         
    Outer_Wafer_Thickness = blocks_info[0]->GetDouble("Outer_Wafer_Thickness","mm");     
    Outer_Wafer_AlThickness = blocks_info[0]->GetDouble("Outer_Wafer_AlThickness","micrometer");   
    Outer_Wafer_PADExternal = blocks_info[0]->GetDouble("Outer_Wafer_PADExternal","mm");   
    Outer_Wafer_PADInternal = blocks_info[0]->GetDouble("Outer_Wafer_PADInternal","mm");   
    Outer_Wafer_GuardRing = blocks_info[0]->GetDouble("Outer_Wafer_GuardRing","mm");     
    Outer_PCB_PortWidth = blocks_info[0]->GetDouble("Outer_PCB_PortWidth","mm");       
    Outer_PCB_StarboardWidth = blocks_info[0]->GetDouble("Outer_PCB_StarboardWidth","mm");  
    Outer_PCB_BevelAngle = blocks_info[0]->GetDouble("Outer_PCB_BevelAngle","deg");      
    Outer_PCB_UpstreamWidth = blocks_info[0]->GetDouble("Outer_PCB_UpstreamWidth","mm");   
    Outer_PCB_DownstreamWidth = blocks_info[0]->GetDouble("Outer_PCB_DownstreamWidth","mm"); 
    Outer_PCB_MidWidth = blocks_info[0]->GetDouble("Outer_PCB_MidWidth","mm");        
    Outer_PCB_Thickness = blocks_info[0]->GetDouble("Outer_PCB_Thickness","mm");       

  }

  else{
    cout << "ERROR: check your input file formatting " << endl;
    exit(1);
  }


  // Inner Barrel
  vector<NPL::InputBlock*> blocks_inner = parser.GetAllBlocksWithTokenAndValue("Strasse","Inner");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_inner.size() << " inner detectors found " << endl; 

  vector<string> coord = {"Radius","Z","Phi"};

  for(unsigned int i = 0 ; i < blocks_inner.size() ; i++){
    if(blocks_inner[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse inner detector" << i+1 <<  endl;

      double R = blocks_inner[i]->GetDouble("Radius","mm");
      double Z= blocks_inner[i]->GetDouble("Z","mm");
      double Phi = blocks_inner[i]->GetDouble("Phi","deg");
      AddInnerDetector(R,Z,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  // Outer barrel
  vector<NPL::InputBlock*> blocks_outer = parser.GetAllBlocksWithTokenAndValue("Strasse","Outer");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_outer.size() << " outer detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_outer.size() ; i++){
    if(blocks_outer[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse outer detector" << i+1 <<  endl;

      double R = blocks_outer[i]->GetDouble("Radius","mm");
      double Z= blocks_outer[i]->GetDouble("Z","mm");
      double Phi = blocks_outer[i]->GetDouble("Phi","deg");
      AddOuterDetector(R,Z,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Strasse::ConstructDetector(G4LogicalVolume* world){
  // Inner Barrel
  for (unsigned short i = 0 ; i < m_Inner_R.size() ; i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,m_Inner_R[i],m_Inner_Z[i]) ;
    Det_pos.rotate(-m_Inner_Phi[i],G4ThreeVector(0,0,1));
    G4RotationMatrix* Rot =  new G4RotationMatrix(0*deg,0*deg,m_Inner_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
        BuildInnerDetector(),
        "Strasse",world,false,i+1);

  }

  // Outer Barrel 
  for (unsigned short i = 0 ; i < m_Outer_R.size() ; i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,m_Outer_R[i],m_Outer_Z[i]) ;
    Det_pos.rotate(-m_Outer_Phi[i],G4ThreeVector(0,0,1));
    G4RotationMatrix* Rot =  new G4RotationMatrix(0*deg,0*deg,m_Outer_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
        BuildOuterDetector(),
        "Strasse",world,false,i+1);
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Strasse::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Strasse")){
    pTree->Branch("Strasse", "TStrasseData", &m_Event) ;
  }
  pTree->SetBranchAddress("Strasse", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Strasse::ReadSensitive(const G4Event* ){
  m_Event->Clear();
  /*
  ///////////
  // First Stage scorer
  DSSDScorers::PS_Rectangle* InnerScorer= (DSSDScorers::PS_Rectangle*) m_InnerScorer->GetPrimitive(0);

  unsigned int sizeFront = InnerScorer->GetLengthMult(); 
  for(unsigned int i = 0 ; i < sizeFront ; i++){
  double Energy = RandGauss::shoot(InnerScorer->GetEnergyLength(i), ResoEnergy);   
  if(Energy>EnergyThreshold){
  int DetNbr  = InnerScorer->GetDetectorLength(i);
  int StripFront = InnerScorer->GetStripLength(i);
  m_Event->SetInnerXE(DetNbr, StripFront, Energy);
  }
  }
  unsigned int sizeBack = InnerScorer->GetWidthMult(); 
  for(unsigned int i = 0 ; i < sizeBack ; i++){
  double Energy = RandGauss::shoot(InnerScorer->GetEnergyWidth(i), ResoEnergy);   
  if(Energy>EnergyThreshold){
  int DetNbr  = InnerScorer->GetDetectorWidth(i);
  int StripFront = InnerScorer->GetStripWidth(i);
  m_Event->SetInnerYE(DetNbr, StripFront, Energy);
  }
  }
  InnerScorer->clear();

  ///////////
  // Second Stage scorer
  DSSDScorers::PS_Rectangle* OuterScorer= (DSSDScorers::PS_Rectangle*) m_OuterScorer->GetPrimitive(0);

  unsigned int sizeFrontOuter = OuterScorer->GetLengthMult(); 
  for(unsigned int i = 0 ; i < sizeFrontOuter ; i++){
  double Energy = RandGauss::shoot(OuterScorer->GetEnergyLength(i), ResoEnergy);   
  if(Energy>EnergyThreshold){
  int DetNbr  = OuterScorer->GetDetectorLength(i);
  int StripFront = OuterScorer->GetStripLength(i);
  m_Event->SetOuterXE(DetNbr, StripFront, Energy);
  }
  }
  unsigned int sizeBackOuter = OuterScorer->GetWidthMult(); 
  for(unsigned int i = 0 ; i < sizeBackOuter ; i++){
  double Energy = RandGauss::shoot(OuterScorer->GetEnergyWidth(i), ResoEnergy);   
  if(Energy>EnergyThreshold){
  int DetNbr  = OuterScorer->GetDetectorWidth(i);
  int StripFront = OuterScorer->GetStripWidth(i);
  m_Event->SetOuterYE(DetNbr, StripFront, Energy);
  }
  }
  OuterScorer->clear();
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Strasse::InitializeScorers() { 
  /*  // This check is necessary in case the geometry is reloaded
      bool already_exist = false; 
      m_InnerScorer = CheckScorer("InnerScorer",already_exist) ;
      m_OuterScorer = CheckScorer("OuterScorer",already_exist) ;

      if(already_exist) 
      return ;

  // Otherwise the scorer is initialised
  G4VPrimitiveScorer* InnerScorer = new DSSDScorers::PS_Rectangle("InnerScorer",1,
  TrapezoidBaseLarge,
  TrapezoidHeight,
  128,128);
  G4VPrimitiveScorer* OuterScorer = new DSSDScorers::PS_Rectangle("OuterScorer",1,
  TrapezoidBaseLarge,
  TrapezoidHeight,
  16,16);

  G4VPrimitiveScorer* InteractionInner = new InteractionScorers::PS_Interactions("InteractionInner",ms_InterCoord,0);
  G4VPrimitiveScorer* InteractionOuter = new InteractionScorers::PS_Interactions("InteractionOuter",ms_InterCoord,0);

  // Register it to the multifunctionnal detector
  m_InnerScorer->RegisterPrimitive(InnerScorer);
  m_InnerScorer->RegisterPrimitive(InteractionInner);
  m_OuterScorer->RegisterPrimitive(OuterScorer);
  m_OuterScorer->RegisterPrimitive(InteractionOuter);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_InnerScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterScorer);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Strasse::Construct(){
  return  (NPS::VDetector*) new Strasse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::InitializeMaterial(){
  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialAl = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Strasse{
    public:
      proxy_nps_Strasse(){
        NPS::DetectorFactory::getInstance()->AddToken("Strasse","Strasse");
        NPS::DetectorFactory::getInstance()->AddDetector("Strasse",Strasse::Construct);
      }
  };

  proxy_nps_Strasse p_nps_Strasse;
}


