/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : March 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Scone simulation                             *
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
#include "G4SubtractionSolid.hh"

// NPTool header
#include "Scone.hh"
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
namespace Scone_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 4.5*ns ;
  const double ResoEnergy = 1.0*MeV ;
  const double XSection = 25.1*mm ; 
  const double YSection = 25.6*mm ;
  const double LengthR1 = 1000*mm ;
  const double LengthR2 = 500*mm ;
  const double Length2x2 = 400*mm ;
  const string Material = "CH2";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Scone Specific Method
Scone::Scone(){
  m_Event = new TSconeData() ;
  m_SconeScorer = 0;
  
  m_SquareDetector = 0;

  m_BuildRing1 = 1;
  m_BuildRing2 = 1;

  m_Assembly = 0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));   
  m_Vis2x2 = new G4VisAttributes(G4Colour(0, 0.5, 0.6, 1.0));   
  m_Vis6x6 = new G4VisAttributes(G4Colour(0.2, 1, 1, 1));   
  m_Vis6x6R2 = new G4VisAttributes(G4Colour(0.2, 0.8, 0.6, 0.8));   

}

Scone::~Scone(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Scone::Build2x2Assembly(){
  if(!m_SquareDetector){
    m_2x2Assembly = new G4AssemblyVolume();
    G4RotationMatrix *Rv = new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0); Tv.setY(0); Tv.setZ(0);

    // One bar definitation
    G4Box* box = new G4Box("Scone_Box",Scone_NS::XSection*0.5,
        Scone_NS::YSection*0.5,Scone_NS::Length2x2*0.5);
    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Scone_NS::Material);
    m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_Scone_Box",0,0,0);
    m_SquareDetector->SetVisAttributes(m_Vis2x2);
    m_SquareDetector->SetSensitiveDetector(m_SconeScorer);
    
    double posX = 0;
    double posY = 0;
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        posX = i*Scone_NS::XSection;
        posY = j*Scone_NS::YSection;
        Tv.setX(posX);
        Tv.setY(posY);
        m_2x2Assembly->AddPlacedVolume(m_SquareDetector, Tv, Rv);
      }
    }
  }

  return m_2x2Assembly;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Scone::Build6x6Assembly(double plastic_length){
  m_6x6Assembly = new G4AssemblyVolume();
  G4RotationMatrix *Rv = new G4RotationMatrix(0,0,0);
  G4ThreeVector Tv;
  Tv.setX(0); Tv.setY(0); Tv.setZ(0);

  // One bar definitation
  G4Box* box = new G4Box("Scone_Box",Scone_NS::XSection*0.5,
       Scone_NS::YSection*0.5,plastic_length*0.5);
  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Scone_NS::Material);
  m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_Scone_Box",0,0,0);
  if(plastic_length == Scone_NS::LengthR2)
    m_SquareDetector->SetVisAttributes(m_Vis6x6R2);
  else 
    m_SquareDetector->SetVisAttributes(m_Vis6x6);
  m_SquareDetector->SetSensitiveDetector(m_SconeScorer);
    
  double posX = 0;
  double posY = 0;
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      posX = i*Scone_NS::XSection;
      posY = j*Scone_NS::YSection;
      Tv.setX(posX);
      Tv.setY(posY);
      m_6x6Assembly->AddPlacedVolume(m_SquareDetector, Tv, Rv);
    }
  }

  return m_6x6Assembly;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Scone::BuildSquareDetector(){
  if(!m_SquareDetector){
    G4Box* box = new G4Box("Scone_Box",Scone_NS::XSection*0.5,
        Scone_NS::YSection*0.5,Scone_NS::LengthR1*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Scone_NS::Material);
    m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_Scone_Box",0,0,0);
    m_SquareDetector->SetVisAttributes(m_VisSquare);
    m_SquareDetector->SetSensitiveDetector(m_SconeScorer);
  }
  return m_SquareDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in Detecor
void Scone::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Scone");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Ring1","Ring2"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Scone " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      m_BuildRing1 = blocks[i]->GetInt("Ring1");
      m_BuildRing2 = blocks[i]->GetInt("Ring2");
      AddDetector(Pos);
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
void Scone::ConstructDetector(G4LogicalVolume* world){
  Build2x2Block(world);

  if(m_BuildRing1==1)
    BuildRing1(world);

  if(m_BuildRing2==1)
    BuildRing2(world);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::Build2x2Block(G4LogicalVolume* world){

    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    double posX_orig = -0.5*Scone_NS::XSection;
    double posX_neg = posX_orig - 2*Scone_NS::XSection;
    double posX_pos = posX_orig + 2*Scone_NS::XSection;
    double posY_orig = -0.5*Scone_NS::YSection;
    double posY_neg = posY_orig - 2*Scone_NS::YSection;
    double posY_pos = posY_orig + 2*Scone_NS::YSection;
    
    double posX[16] = {posX_orig, posX_orig, posX_neg, posX_neg, posX_neg, posX_pos, posX_pos, posX_pos, posX_orig, posX_orig, posX_neg, posX_neg, posX_neg, posX_pos, posX_pos, posX_pos};
    double posY[16] = {posY_neg, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_orig, posY_pos};
    double posZ[16] = {-300,-300,-300,-300,-300,-300,-300,-300,300,300,300,300,300,300,300,300};

    G4ThreeVector Det_pos;

    for(int i=0; i<16; i++)
    {
      Det_pos = G4ThreeVector(posX[i], posY[i], posZ[i]);
      Build2x2Assembly()->MakeImprint(world, Det_pos, Rot, m_Assembly);
      m_Assembly++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::BuildRing1(G4LogicalVolume* world){
    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    double posX_orig = -2.5*Scone_NS::XSection;
    double posX_neg = posX_orig - 6*Scone_NS::XSection;
    double posX_pos = posX_orig + 6*Scone_NS::XSection;
    double posY_orig = -2.5*Scone_NS::YSection;
    double posY_neg = posY_orig - 6*Scone_NS::YSection;
    double posY_pos = posY_orig + 6*Scone_NS::YSection;
    
    double posX[8] = {posX_orig, posX_orig, posX_neg, posX_neg, posX_neg, posX_pos, posX_pos, posX_pos};
    double posY[8] = {posY_neg, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_orig, posY_pos};
    double posZ[8] = {0,0,0,0,0,0,0,0};

    G4ThreeVector Det_pos;

    for(int i=0; i<8; i++)
    {
      Det_pos = G4ThreeVector(posX[i], posY[i], posZ[i]);
      Build6x6Assembly(Scone_NS::LengthR1)->MakeImprint(world, Det_pos, Rot, m_Assembly);
      m_Assembly++;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::BuildRing2(G4LogicalVolume* world){
    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    double posX_orig = -2.5*Scone_NS::XSection;
    double posX_left = posX_orig - 12*Scone_NS::XSection;
    double posX_right = posX_orig + 12*Scone_NS::XSection;
    double posY_orig = -2.5*Scone_NS::YSection;
    double posY_down = posY_orig - 12*Scone_NS::YSection;
    double posY_up = posY_orig + 12*Scone_NS::YSection;
  
    //double posX[16] = {posX_left, posX_left, posX_left, posX_left, posX_left};
    //double posY[16] = {posY_down};
    //double posZ[16] = {0};

    G4ThreeVector Det_pos;

    for(int i=0; i<5; i++)
    {
      double posX = posX_left;
      double posY = posY_down + 6*i*Scone_NS::YSection;
      Det_pos = G4ThreeVector(posX, posY, 0);
      Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
      m_Assembly++;
    }
    for(int i=0; i<5; i++)
    {
      double posX = posX_right;
      double posY = posY_down + 6*i*Scone_NS::YSection;
      Det_pos = G4ThreeVector(posX, posY, 0);
      Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
      m_Assembly++;
    }

    Det_pos = G4ThreeVector(posX_left + 6*Scone_NS::XSection, posY_down, 0);
    Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
    m_Assembly++;
    Det_pos = G4ThreeVector(posX_left + 12*Scone_NS::XSection, posY_down, 0);
    Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
    m_Assembly++;
    Det_pos = G4ThreeVector(posX_left + 18*Scone_NS::XSection, posY_down, 0);
    Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
    m_Assembly++;

    Det_pos = G4ThreeVector(posX_left + 6*Scone_NS::XSection, posY_up, 0);
    Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
    m_Assembly++;
    Det_pos = G4ThreeVector(posX_left + 12*Scone_NS::XSection, posY_up, 0);
    Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
    m_Assembly++;
    Det_pos = G4ThreeVector(posX_left + 18*Scone_NS::XSection, posY_up, 0);
    Build6x6Assembly(Scone_NS::LengthR2)->MakeImprint(world, Det_pos, Rot, m_Assembly);
    m_Assembly++;




}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Scone::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Scone")){
    pTree->Branch("Scone", "TSconeData", &m_Event) ;
  }
  pTree->SetBranchAddress("Scone", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Scone::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_SconeScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Scone_NS::ResoEnergy);
    if(Energy>Scone_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Scone_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Scone::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_SconeScorer = CheckScorer("SconeScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_SconeScorer->RegisterPrimitive(Calorimeter);
  m_SconeScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SconeScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Scone::Construct(){
  return  (NPS::VDetector*) new Scone();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Scone{
    public:
      proxy_nps_Scone(){
        NPS::DetectorFactory::getInstance()->AddToken("Scone","Scone");
        NPS::DetectorFactory::getInstance()->AddDetector("Scone",Scone::Construct);
      }
  };

  proxy_nps_Scone p_nps_Scone;
}
