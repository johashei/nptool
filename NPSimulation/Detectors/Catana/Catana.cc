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

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
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
#include "Catana.hh"
#include "CalorimeterScorers.hh"
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
namespace Catana_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 4.5*ns ;
  const double ResoEnergy = 1.0*MeV ;
  double DummyInnerRadius = 80*mm ; 
  double DummyOuterRadius = 160*mm ; 
  double Length = 200*mm ;
  string Material = "CsI";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Catana Specific Method
Catana::Catana(){
  m_Event = new TCatanaData() ;
  m_CatanaScorer = 0;
  m_DummyDetector = 0;

  // RGB Color + Transparency
  m_VisCrystal = new G4VisAttributes(G4Colour(1, 0.8, 0, 0.5));   

}

Catana::~Catana(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Catana::AddDummyDetector(double Z){
  m_Z.push_back(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Catana::BuildDummyDetector(){
  if(!m_DummyDetector){
    G4Tubs* tub = new G4Tubs("Catana_Dummy",Catana_NS::DummyInnerRadius,Catana_NS::DummyOuterRadius,Catana_NS::Length*0.5,0,360*deg);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Catana_NS::Material);
    m_DummyDetector = new G4LogicalVolume(tub,DetectorMaterial,"logic_Catana_tub",0,0,0);
    m_DummyDetector->SetVisAttributes(m_VisCrystal);
    m_DummyDetector->SetSensitiveDetector(m_CatanaScorer);

  }
  return m_DummyDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Catana::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("Catana","Dummy");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> token = {"Z"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      double Z = blocks[i]->GetDouble("Z","mm");
      AddDummyDetector(Z);
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
void Catana::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_Z.size() ; i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,0,m_Z[i]) ;
    G4RotationMatrix* Rot = new G4RotationMatrix();
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildDummyDetector(),
          "CatanaDummy",world,false,i+1);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Catana::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Catana")){
    pTree->Branch("Catana", "TCatanaData", &m_Event) ;
  }
  pTree->SetBranchAddress("Catana", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Catana::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_CatanaScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Catana_NS::ResoEnergy);
    if(Energy>Catana_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Catana_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Catana::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_CatanaScorer = CheckScorer("CatanaScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_CatanaScorer->RegisterPrimitive(Calorimeter);
  m_CatanaScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_CatanaScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Catana::Construct(){
  return  (NPS::VDetector*) new Catana();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Catana{
    public:
      proxy_nps_Catana(){
        NPS::DetectorFactory::getInstance()->AddToken("Catana","Catana");
        NPS::DetectorFactory::getInstance()->AddDetector("Catana",Catana::Construct);
      }
  };

  proxy_nps_Catana p_nps_Catana;
}
