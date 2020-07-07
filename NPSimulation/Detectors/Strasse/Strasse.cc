/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Strasse simulation                             *
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

  // Trapezoid dimension
  //const double TrapezoidBaseLarge = 95*mm;
  const double TrapezoidBaseLarge = 78.1*mm;
  //const double TrapezoidBaseSmall = 45*mm;
  const double TrapezoidBaseSmall = 43.3*mm;
  //const double TrapezoidHeight = 118*mm;
  const double TrapezoidHeight = 61.8*mm;
  const double TrapezoidLength = 1*cm;
  const double InnerThickness = 100*um;
  const double OutterThickness = 1*mm;
  const double DistanceBetweenSi = 7*mm;
  //const double InnerNbrOfStrips = 128;
  //const double OutterNbrOfStrips = 16;
}
using namespace Strasse_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Strasse Specific Method
Strasse::Strasse(){
  m_Event = new TStrasseData() ;
  m_InnerScorer = 0;
  m_OutterScorer = 0;
  m_TrapezoidDetector = 0;
}

Strasse::~Strasse(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildTrapezoidDetector(){
  if(!m_TrapezoidDetector){
    // Definittion of the volume containing the sensitive detectors
    G4Trap* solidTrapezoid = new G4Trap("Strasse",
        TrapezoidLength*0.5, 0*deg, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg, 
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg);
    G4LogicalVolume* logicTrapezoid = new G4LogicalVolume(solidTrapezoid,
        MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"),
        "Strasse",
        0,0,0);

    G4VisAttributes* TrapezoidVisAtt = new G4VisAttributes(G4Colour(0.90, 0.90, 0.90));
    TrapezoidVisAtt->SetForceWireframe(true);
    logicTrapezoid->SetVisAttributes(TrapezoidVisAtt);

    // First stage silicon detector
    G4ThreeVector positionInner = G4ThreeVector(0,0,-4*mm);

    G4Trap* solidInner = new G4Trap("solidFirstSatge",
        InnerThickness*0.5, 0*deg, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg);
    G4LogicalVolume* logicInner = new G4LogicalVolume(solidInner,
        MaterialManager::getInstance()->GetMaterialFromLibrary("Si"),
        "logicInner",
        0,0,0);
    new G4PVPlacement(0,
        positionInner,
        logicInner,
        "Strasse_Inner",
        logicTrapezoid,
        false,
        0);
    // Set First Stage sensitive
    logicInner->SetSensitiveDetector(m_InnerScorer);

    // Visualisation of First Stage strips
    G4VisAttributes* InnerVisAtt = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
    logicInner->SetVisAttributes(InnerVisAtt);

    //////
    // Second stage silicon detector
    G4ThreeVector positionOutter = G4ThreeVector(0,0,-0.5*TrapezoidLength+DistanceBetweenSi);

    G4Trap* solidOutter = new G4Trap("solidSecondSatge",
        OutterThickness*0.5, 0*deg, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg);
    G4LogicalVolume* logicOutter = new G4LogicalVolume(solidOutter,
        MaterialManager::getInstance()->GetMaterialFromLibrary("Si"),
        "logicOutter",
        0,0,0);
    new G4PVPlacement(0,
        positionOutter,
        logicOutter,
        "Strasse_Outter",
        logicTrapezoid,
        false,
        0);
    // Set Second Stage sensitive
    logicOutter->SetSensitiveDetector(m_OutterScorer);

    // Visualisation of Second Stage strips
    G4VisAttributes* OutterVisAtt = new G4VisAttributes(G4Colour(0.4,0.5,0.5));
    logicOutter->SetVisAttributes(OutterVisAtt);


    m_TrapezoidDetector = logicTrapezoid;
  }
  return m_TrapezoidDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Strasse::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Strasse");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetector(R,Theta,Phi);
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
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = TrapezoidHeight*0.5 + m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*Strasse_NS::TrapezoidLength*0.5;
    // Building Detector reference frame
    G4double ii = cos(m_Theta[i]) * cos(m_Phi[i]);
    G4double jj = cos(m_Theta[i]) * sin(m_Phi[i]);
    G4double kk = -sin(m_Theta[i]);
    G4ThreeVector Y(ii,jj,kk);
    G4ThreeVector w = Det_pos.unit();
    G4ThreeVector u = w.cross(Y);
    G4ThreeVector v = w.cross(u);
    v = v.unit();
    u = u.unit();

    G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
        BuildTrapezoidDetector(),
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
  DSSDScorers::PS_Rectangle* OutterScorer= (DSSDScorers::PS_Rectangle*) m_OutterScorer->GetPrimitive(0);

  unsigned int sizeFrontOutter = OutterScorer->GetLengthMult(); 
  for(unsigned int i = 0 ; i < sizeFrontOutter ; i++){
    double Energy = RandGauss::shoot(OutterScorer->GetEnergyLength(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = OutterScorer->GetDetectorLength(i);
      int StripFront = OutterScorer->GetStripLength(i);
      m_Event->SetOutterXE(DetNbr, StripFront, Energy);
    }
  }
  unsigned int sizeBackOutter = OutterScorer->GetWidthMult(); 
  for(unsigned int i = 0 ; i < sizeBackOutter ; i++){
    double Energy = RandGauss::shoot(OutterScorer->GetEnergyWidth(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = OutterScorer->GetDetectorWidth(i);
      int StripFront = OutterScorer->GetStripWidth(i);
      m_Event->SetOutterYE(DetNbr, StripFront, Energy);
    }
  }
  OutterScorer->clear();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Strasse::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_InnerScorer = CheckScorer("InnerScorer",already_exist) ;
  m_OutterScorer = CheckScorer("OutterScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  G4VPrimitiveScorer* InnerScorer = new DSSDScorers::PS_Rectangle("InnerScorer",1,
      TrapezoidBaseLarge,
      TrapezoidHeight,
      128,128);
  G4VPrimitiveScorer* OutterScorer = new DSSDScorers::PS_Rectangle("OutterScorer",1,
      TrapezoidBaseLarge,
      TrapezoidHeight,
      16,16);

  G4VPrimitiveScorer* InteractionInner = new InteractionScorers::PS_Interactions("InteractionInner",ms_InterCoord,0);
  G4VPrimitiveScorer* InteractionOutter = new InteractionScorers::PS_Interactions("InteractionOutter",ms_InterCoord,0);

  // Register it to the multifunctionnal detector
  m_InnerScorer->RegisterPrimitive(InnerScorer);
  m_InnerScorer->RegisterPrimitive(InteractionInner);
  m_OutterScorer->RegisterPrimitive(OutterScorer);
  m_OutterScorer->RegisterPrimitive(InteractionOutter);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_InnerScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_OutterScorer);
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
