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

// NPTool header
#include "SamuraiFDC2.hh"
//#include "CalorimeterScorers.hh"
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
  const double FDC2_Width = 2616*mm; //(x)
  const double FDC2_Height = 1156*mm;//(y)
  const double FDC2_Depth = 876*mm;//(z)
  const string FDC2_Material = "G4_Fe"; 

  //Detector Number
  const short int FDC2_DetectorNumber = 2;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Samurai Specific Method
SamuraiFDC2::SamuraiFDC2(){

  //Visualization attributes
  m_VisFDC2 = new G4VisAttributes(G4Colour(1,0,1,0.5));
  //Logical volumes
  m_FDC2 = NULL;
  //Scorer
  m_FDC2Scorer = NULL;
  //Ideal Data event
  m_Event = new TSamuraiIdealData;

}

SamuraiFDC2::~SamuraiFDC2(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiFDC2::AddDetector(G4ThreeVector Mag_Pos, double Mag_Angle, G4ThreeVector Offset, double Off_Angle){

  m_Angle = Mag_Angle + (90.*deg - Off_Angle);

  Offset.rotateY(-(m_Angle));
  m_Pos = Mag_Pos + Offset;

  return;
}

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

// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// (Called in DetecorConstruction::ReadDetectorConfiguration Method)
void SamuraiFDC2::ReadConfiguration(NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Samurai");
  vector<NPL::InputBlock*> blocks2 = parser.GetAllBlocksWithToken("SAMURAIFDC2");

  if(blocks.size()==1 && blocks2.size()==1){
    if(NPOptionManager::getInstance()->GetVerboseLevel()) {
      cout << "/////// Samurai FDC2 found ///////" << endl;
    }
    vector<string> cart = {"POS","ANGLE"};
    vector<string> sphe = {"R","Theta","Phi","ANGLE"};

    G4ThreeVector Mag_Pos;
    double Mag_Angle;

    if(blocks[0]->HasTokenList(cart)){
      Mag_Pos = NPS::ConvertVector(blocks[0]->GetTVector3("POS", "cm"));
      Mag_Angle = blocks[0]->GetDouble("ANGLE","deg");
    }
    else if(blocks[0]->HasTokenList(sphe)){
      double R = blocks[0]->GetDouble("R","mm");
      double Theta = blocks[0]->GetDouble("Theta","deg");
      double Phi = blocks[0]->GetDouble("Phi","deg");
      Mag_Pos.setMag(R);
      Mag_Pos.setTheta(Theta);
      Mag_Pos.setPhi(Phi);
      Mag_Angle = blocks[0]->GetDouble("ANGLE","deg");
    }
    
    G4ThreeVector Offset = NPS::ConvertVector(blocks2[0]->GetTVector3("Offset", "mm"));
    double Off_Angle = blocks2[0]->GetDouble("OffAngle","deg");

    //string xml = blocks2[0]->GetString("XML");
    //bool invert_x = blocks2[0]->GetBool("InvertX");
    //bool invert_y = blocks2[0]->GetBool("InvertY");
    //bool invert_z = blocks2[0]->GetBool("InvertD");
  
    AddDetector(Mag_Pos, Mag_Angle, Offset, Off_Angle);

  }
  else{
    cout << "ERROR: there should be only one Samurai magnet, check your input file" << endl;
    exit(1);
  }
    
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// (Called After DetectorConstruction::AddDetector Method)
void SamuraiFDC2::ConstructDetector(G4LogicalVolume* world){

  G4RotationMatrix* Rot = new G4RotationMatrix();
  Rot->rotateY(m_Angle);
  new G4PVPlacement(Rot, m_Pos,
          BuildFDC2(), "SamuraiFDC2", world, false, 0);

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SamuraiFDC2::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("IdealDataFDC2")){
    pTree->Branch("IdealDataFDC2", "TSamuraiIdealData", &m_Event) ;
  }
  pTree->SetBranchAddress("IdealDataFDC2", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Read sensitive part and fill the Root tree.
// (Called at in the EventAction::EndOfEventAvtion)

void SamuraiFDC2::ReadSensitive(const G4Event* event){
  
  m_Event->Clear();
  //Interaction Scorer
  InteractionScorers::PS_Interactions* Scorer= (InteractionScorers::PS_Interactions*) m_FDC2Scorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    //vector<unsigned int> level = Scorer->GetLevel(i); 
    //double Energy = RandGauss::shoot(Scorer->GetEnergy(i),SamuraiFDC2_NS::ResoEnergy);
    //double Energy = Scorer->GetEnergy(i);

    double energy = Scorer->GetEnergy(i);
    double brho = Scorer->GetBrho(i);
    double posx = Scorer->GetPositionX(i);
    double posy = Scorer->GetPositionY(i);
    double posz = Scorer->GetPositionZ(i);
    double mom_mag = brho*Scorer->GetCharge(i);
    double theta = Scorer->GetTheta(i);
    double phi = Scorer->GetPhi(i);
    m_Event->SetData(SamuraiFDC2_NS::FDC2_DetectorNumber, energy, posx, posy, posz, mom_mag, theta, phi, brho);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  

void SamuraiFDC2::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_FDC2Scorer = CheckScorer("FDC2Scorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  //Interaction 
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ; 
  m_FDC2Scorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_FDC2Scorer) ;
}

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
extern"C" {
  class proxy_nps_samuraiFDC2{
    public:
      proxy_nps_samuraiFDC2(){
        NPS::DetectorFactory::getInstance()->AddToken("SAMURAIFDC2","SAMURAIFDC2");
        NPS::DetectorFactory::getInstance()->AddDetector("SAMURAIFDC2",SamuraiFDC2::Construct);
      }
  };
  
  proxy_nps_samuraiFDC2 p_nps_samuraiFDC2;
}


