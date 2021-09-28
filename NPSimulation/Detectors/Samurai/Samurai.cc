/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
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
 *  This class describe  Samurai simulation                                 *
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
#include "G4RegionStore.hh"
#include "G4FastSimulationManager.hh"
#include "G4UserLimits.hh"

// NPTool header
#include "Samurai.hh"
#include "SamuraiFieldPropagation.hh"
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
namespace Samurai_NS{
  // Samurai magnet construction paramethers
  const double Magnet_Width = 6700*mm;
  const double Magnet_Height = 4640*mm;
  const double Magnet_Depth = 3500*mm;
  const string Magnet_Material = "G4_Galactic"; //where the main beam will travel
  const double Yoke_Height = 1880*mm; //(4640-880)/2
  const string Yoke_Material = "G4_Fe";

  const double StepSize = 1*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Samurai Specific Method
Samurai::Samurai(){

  //Visualization attributes
  m_VisMagnet = new G4VisAttributes(G4Colour(0,0,1,0.5));
  m_VisYokes = new G4VisAttributes(G4Colour(0,1,0,0.5));

  //Logical volumes
  m_Magnet = NULL;
  m_Yoke = NULL;

  //Propagation region
  m_PropagationRegion = NULL;
}

Samurai::~Samurai(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Samurai::AddMagnet(G4ThreeVector POS, double Angle){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R = POS.mag();
  m_Theta = POS.theta();
  m_Phi = POS.phi();
  m_Angle = Angle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Samurai::AddMagnet(double R, double Theta, double Phi, double Angle){
  m_R = R;
  m_Theta = Theta;
  m_Phi = Phi;
  m_Angle = Angle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Samurai::BuildMagnet(){
  if(!m_Magnet){
    //Shape - G4Box
    G4Box* box = new G4Box("Samurai_Box",Samurai_NS::Magnet_Width*0.5,
			   Samurai_NS::Magnet_Height*0.5,Samurai_NS::Magnet_Depth*0.5);
  
    //Material - vacuum
    G4Material* VacuumMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(Samurai_NS::Magnet_Material);

    //Logical Volume
    m_Magnet = new G4LogicalVolume(box, VacuumMaterial, "logic_Samurai_box",0,0,0);
    m_Magnet->SetVisAttributes(m_VisMagnet);
  }
  return m_Magnet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Samurai::BuildYoke(){
  if(!m_Yoke){
    //Shape - G4Box
    G4Box* yoke = new G4Box("Samurai_yoke_1",Samurai_NS::Magnet_Width*0.5,
			    Samurai_NS::Yoke_Height*0.5,Samurai_NS::Magnet_Depth*0.5);

    //Material
    G4Material* YokeMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(Samurai_NS::Yoke_Material);

    //Logical Volume
    m_Yoke = new G4LogicalVolume(yoke, YokeMaterial, "logic_Samurai_yoke",0,0,0);
    m_Yoke->SetVisAttributes(m_VisYokes);
  }
  return m_Yoke;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetectorConfiguration Method
void Samurai::ReadConfiguration(NPL::InputParser parser){
  
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Samurai");

  if(blocks.size()==1){
    if(NPOptionManager::getInstance()->GetVerboseLevel()) {
      cout << "//// Samurai magnet found " << endl;
    }
    vector<string> cart = {"POS","ANGLE"};
    vector<string> sphe = {"R","Theta","Phi","ANGLE"};

    if(blocks[0]->HasTokenList(cart)){
      G4ThreeVector Pos = NPS::ConvertVector(blocks[0]->GetTVector3("POS", "cm"));
      double Angle = blocks[0]->GetDouble("ANGLE","deg");
      AddMagnet(Pos,Angle);
    }
    else if(blocks[0]->HasTokenList(sphe)){
      double R = blocks[0]->GetDouble("R","mm");
      double Theta = blocks[0]->GetDouble("Theta","deg");
      double Phi = blocks[0]->GetDouble("Phi","deg");
      double Angle = blocks[0]->GetDouble("ANGLE","deg");
      AddMagnet(R,Theta,Phi,Angle);
    }
  }
  else{
    cout << "ERROR: there should be only one Samurai magnet, check your input file" << endl;
    exit(1);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetectorConstruction::AddDetector Method
void Samurai::ConstructDetector(G4LogicalVolume* world){

  //Yokes placement
  double distance = ( Samurai_NS::Magnet_Height - Samurai_NS::Yoke_Height ) * 0.5;
  G4ThreeVector upper (0, distance,0);
  G4ThreeVector lower (0,-distance,0);
  
  new G4PVPlacement(0, upper, BuildYoke(), "Upper_Yoke", BuildMagnet(), false, 0);
  new G4PVPlacement(0, lower, BuildYoke(), "Lower_Yoke", BuildMagnet(), false, 0);
  

  //Magnet placement
  G4double wX = m_R * sin(m_Theta) * cos(m_Phi);
  G4double wY = m_R * sin(m_Theta) * sin(m_Phi);
  G4double wZ = m_R * cos(m_Theta);
  G4ThreeVector Mag_pos = G4ThreeVector(wX, wY, wZ);
  G4ThreeVector u ( cos(m_Angle) , 0, -sin(m_Angle) );
  G4ThreeVector v (0,1,0);
  G4ThreeVector w ( sin(m_Angle) , 0,  cos(m_Angle) );
  
  G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
  
  //new G4PVPlacement(G4Transform3D(*Rot,Mag_pos),
  //        m_Magnet, "Samurai",world, false, 0);
  new G4PVPlacement(Rot, Mag_pos,
          BuildMagnet(), "Samurai", world, false, 0);

  // Set region were magnetic field is active
  SetPropagationRegion();

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Set region were magnetic field is active
// Called in Samurai::ConstructDetector
void Samurai::SetPropagationRegion(){
  
  if(!m_PropagationRegion){
    m_PropagationRegion= new G4Region("NPSamuraiFieldPropagation");
    m_PropagationRegion -> AddRootLogicalVolume(m_Magnet);
    m_PropagationRegion->SetUserLimits(new G4UserLimits(Samurai_NS::StepSize));
  }
  
  G4FastSimulationManager* mng = m_PropagationRegion->GetFastSimulationManager();
  //To make sure no other models are present in the region
  unsigned int size = m_PropagationModel.size();
  for(unsigned int i = 0 ; i < size ; i++)
      mng->RemoveFastSimulationModel(m_PropagationModel[i]);
  m_PropagationModel.clear();
 
  G4VFastSimulationModel* fsm;
  fsm = new NPS::SamuraiFieldPropagation("SamuraiFieldPropagation", m_PropagationRegion);
  ((NPS::SamuraiFieldPropagation*) fsm)->SetStepSize(Samurai_NS::StepSize);
  ((NPS::SamuraiFieldPropagation*) fsm)->SetAngle(m_Angle);
  m_PropagationModel.push_back(fsm);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
//FIXME
/*
void Samurai::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Samurai")){
    pTree->Branch("Samurai", "TSamuraiData", &m_Event) ;
  }
  pTree->SetBranchAddress("Samurai", &m_Event) ;
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
//FIXME
void Samurai::ReadSensitive(const G4Event* ){
}
  
/*
//FIXME
void Samurai::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_SamuraiScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Samurai_NS::ResoEnergy);
    if(Energy>Samurai_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Samurai_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
/*
//FIXME
void Samurai::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_SamuraiScorer = CheckScorer("SamuraiScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_SamuraiScorer->RegisterPrimitive(Calorimeter);
  m_SamuraiScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SamuraiScorer) ;
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Samurai::Construct(){
  return  (NPS::VDetector*) new Samurai();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Samurai{
    public:
      proxy_nps_Samurai(){
        NPS::DetectorFactory::getInstance()->AddToken("Samurai","Samurai");
        NPS::DetectorFactory::getInstance()->AddDetector("Samurai",Samurai::Construct);
      }
  };
  
  proxy_nps_Samurai p_nps_Samurai;
}


