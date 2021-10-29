/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
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

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
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
#include "SamuraiFDC2.hh"
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
namespace SamuraiFDC2_NS{
  // Samurai magnet construction paramethers

  //Main outer box
  const double FDC2_Width = 6700*mm; //(x)
  const double FDC2_Height = 4640*mm;//(y)
  const double FDC2_Depth = 500*mm;//(z)
  const string FDC2_Material = "G4_Fe"; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Samurai Specific Method
SamuraiFDC2::SamuraiFDC2(){

  //Visualization attributes
  m_VisFDC2 = new G4VisAttributes(G4Colour(1,0,1,0.5));
  //Logical volumes
  m_FDC2 = NULL;
  //Scorer
  m_FDC2Scorer = NULL;

}

SamuraiFDC2::~SamuraiFDC2(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiFDC2::AddDetector(G4ThreeVector Mag_Pos, double Angle, G4ThreeVector Offset){

  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_Angle = Angle;
  Offset.rotateY(-m_Angle);
  m_Pos = Mag_Pos + Offset;



  return;
}
/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiFDC2::AddFDC2(double r, double theta, double phi, 
                          double R, double Theta, double Phi, double Angle){

  m_R_d = r;
  m_Theta_d = theta;
  m_Phi_d = phi;

  m_R = R;
  m_Theta = Theta;
  m_Phi = Phi;
  m_Angle = Angle;

  return;
}
//Spherical Magnet
void SamuraiFDC2::AddMagnet(double R, double Theta, double Phi, double Angle){

  
  return;
}
//Cartezian Magnet
void SamuraiFDC2::AddMagnet(G4ThreeVector POS, double Angle){

  m_R = POS.mag();
  m_Theta = POS.theta();
  m_Phi = POS.phi();
  m_Angle = Angle;

  return;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SamuraiFDC2::BuildFDC2(){
  if(!m_FDC2){
    //Shape - G4Box
    G4Box* box = new G4Box("FDC2_Box",SamuraiFDC2_NS::FDC2_Width*0.5,
			   SamuraiFDC2_NS::FDC2_Height*0.5,SamuraiFDC2_NS::FDC2_Depth*0.5);
  
    //Material - vacuum
    G4Material* VacuumMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(SamuraiFDC2_NS::FDC2_Material);

    //Logical Volume
    m_FDC2 = new G4LogicalVolume(box, VacuumMaterial, "logic_SamuraiFDC2_box",0,0,0);
    m_FDC2->SetVisAttributes(m_VisFDC2);
    m_FDC2->SetSensitiveDetector(m_FDC2Scorer);
  }
  return m_FDC2;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetectorConfiguration Method
void SamuraiFDC2::ReadConfiguration(NPL::InputParser parser){

  cout << "I am inside SAMURAIFDC2 read" << endl;
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Samurai");
  vector<NPL::InputBlock*> blocks2 = parser.GetAllBlocksWithToken("SAMURAIFDC2");

  if(blocks.size()==1 && blocks2.size()==1){
    if(NPOptionManager::getInstance()->GetVerboseLevel()) {
      cout << "/////// Samurai FDC2 found ///////" << endl;
    }
    vector<string> cart = {"POS","ANGLE"};
    vector<string> sphe = {"R","Theta","Phi","ANGLE"};

    G4ThreeVector Mag_Pos;
    double Angle;

    if(blocks[0]->HasTokenList(cart)){
      Mag_Pos = NPS::ConvertVector(blocks[0]->GetTVector3("POS", "cm"));
      Angle = blocks[0]->GetDouble("ANGLE","deg");
    }
    else if(blocks[0]->HasTokenList(sphe)){
      double R = blocks[0]->GetDouble("R","mm");
      double Theta = blocks[0]->GetDouble("Theta","deg");
      double Phi = blocks[0]->GetDouble("Phi","deg");
      Mag_Pos.setMag(R);
      Mag_Pos.setTheta(Theta);
      Mag_Pos.setPhi(Phi);
      Angle = blocks[0]->GetDouble("ANGLE","deg");
    }
    
    G4ThreeVector Offset = NPS::ConvertVector(blocks2[0]->GetTVector3("Offset", "mm"));
    //string xml = blocks2[0]->GetString("XML");
    //bool invert_x = blocks2[0]->GetBool("InvertX");
    //bool invert_y = blocks2[0]->GetBool("InvertY");
    //bool invert_z = blocks2[0]->GetBool("InvertD");
  
    AddDetector(Mag_Pos, Angle, Offset);

  }
  else{
    cout << "ERROR: there should be only one Samurai magnet, check your input file" << endl;
    exit(1);
  }
    
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetectorConstruction::AddDetector Method
void SamuraiFDC2::ConstructDetector(G4LogicalVolume* world){
  //magnet placement
  /*G4double wX = m_R * sin(m_Theta) * cos(m_Phi);
  G4double wY = m_R * sin(m_Theta) * sin(m_Phi);
  G4double wZ = m_R * cos(m_Theta);

  G4double mX = m_R * sin(m_Theta) * cos(m_Phi);
  G4double mY = m_R * sin(m_Theta) * sin(m_Phi);
  G4double mZ = m_R * cos(m_Theta);

  G4ThreeVector DFC0_pos = G4ThreeVector(wX+mX, wY+mY, wZ+mZ);*/

  

  G4ThreeVector u ( cos(m_Angle) , 0, -sin(m_Angle) );
  G4ThreeVector v (0,1,0);
  G4ThreeVector w ( sin(m_Angle) , 0,  cos(m_Angle) );
  
  G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
  
  new G4PVPlacement(Rot, m_Pos,
          BuildFDC2(), "SamuraiFDC2", world, false, 0);


  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
//FIXME

void SamuraiFDC2::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("SamuraiFDC2")){
    pTree->Branch("SamuraiFDC2", "TSamuraiData", &m_Event) ; //WATCH OUT !!!!!!
  }
  pTree->SetBranchAddress("SamuraiFDC2", &m_Event) ; //WATCH OUT !!!!!!
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion

//FIXME
void SamuraiFDC2::ReadSensitive(const G4Event* event){
  
  m_Event->Clear();
  //Interaction Scorer
  InteractionScorers::PS_Interactions* Scorer= (InteractionScorers::PS_Interactions*) m_FDC2Scorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    //vector<unsigned int> level = Scorer->GetLevel(i); 
    //double Energy = RandGauss::shoot(Scorer->GetEnergy(i),SamuraiFDC2_NS::ResoEnergy);
    //double Energy = Scorer->GetEnergy(i);
    m_Event->SetDetectorNumber( 2 );//2 is for FDC2
    m_Event->SetEnergy(Scorer->GetEnergy(i));
  }
  /*
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  //CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_FDC2Scorer->GetPrimitive(0);
  
  //Interaction Scorer
  InteractionScorers::PS_Interactions* Scorer= (InteractionScorers::PS_Interactions*) m_FDC2Scorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    //double Energy = RandGauss::shoot(Scorer->GetEnergy(i),SamuraiFDC2_NS::ResoEnergy);
    double Energy = Scorer->GetEnergy();

    if(Energy>SamuraiFDC2_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),SamuraiFDC2_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }*/
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   

//FIXME
void SamuraiFDC2::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_FDC2Scorer = CheckScorer("FDC2Scorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  //G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ; //WATCH OUT 
  //and register it to the multifunctionnal detector
  //m_FDC2Scorer->RegisterPrimitive(Calorimeter);
  m_FDC2Scorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_FDC2Scorer) ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* SamuraiFDC2::Construct(){
  return  (NPS::VDetector*) new SamuraiFDC2();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {////// WHAT IS THIS??
  class proxy_nps_samuraiFDC2{
    public:
      proxy_nps_samuraiFDC2(){
        NPS::DetectorFactory::getInstance()->AddToken("SAMURAIFDC2","SAMURAIFDC2");
        NPS::DetectorFactory::getInstance()->AddDetector("SAMURAIFDC2",SamuraiFDC2::Construct);
      }
  };
  
  proxy_nps_samuraiFDC2 p_nps_samuraiFDC2;
}


