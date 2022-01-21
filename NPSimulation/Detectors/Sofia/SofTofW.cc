/******************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project          *
 *                                                                            *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                 *
 * For the list of contributors see $NPTOOL/Licence/Contributors              *
 ******************************************************************************/

/******************************************************************************
 * Original Author: Pierre Morfouace contact address: pierre.morfouace2@cea.fr*
 *                                                                            *
 * Creation Date  : November 2020                                             *
 * Last update    :                                                           *
 *----------------------------------------------------------------------------*
 * Decription:                                                                *
 *  This class describe a simple SofTofW setup for simulation                   *
 *                                                                            *
 *----------------------------------------------------------------------------*
 * Comment:                                                                   *
 *                                                                            *
 ******************************************************************************/

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
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

// NPTool header
#include "SofTofW.hh"
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
namespace SofTofW_NS{
  // Energy and time Resolution
  const double EnergyThreshold       = 0.1*MeV;
  const double ResoTime              = 0.007*ns;
  const double ResoEnergy            = 1.0*MeV;
  //const double TwinResoEnergy        = 0*keV;
  const double tof_plastic_height    = 660*mm;
  const double tof_plastic_width     = 32*mm;
  const double tof_plastic_thickness = 0.5*mm;
  const string Material              = "BC400";

  const double GLAD_height           = 1.5*m;
  const double GLAD_width            = 5*m;
  const double GLAD_Leff             = 2*m;

  //const double twin_anode_width = 10.*cm;
  //const double twin_anode_height = 10.*cm;
  //const double twin_anode_thickness = 3.1*cm;

  //const double twin_cathode_width = 30*um;
  //const double twin_cathode_height = 20*cm;
  //const double twin_cathode_thickness = 50*cm;


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// SofTofW Specific Method
SofTofW::SofTofW(){
  m_Event = new TSofTofWData() ;
  m_TofScorer = 0;
  m_PlasticTof = 0;
  m_GLAD= 0;
  m_VacuumPipe= 0;
  m_TofWall = 0;

  m_Build_GLAD= 0;
  m_Build_VacuumPipe= 0;
  m_VacuumPipeX= 0;
  m_VacuumPipeY= 0;
  m_VacuumPipeZ= 0;
  m_GLAD_MagField = 0;
  m_GLAD_DistanceFromTarget = 0;

  //m_Build_Twin_Music= 0;
  //m_Twin_Music_DistanceFromTarget= 0;
  //m_Twin_Music_Gas= "P10_1atm";

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0.53, 0.81, 0.98, 0.5));   
  m_VisGLAD = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5));   

}

SofTofW::~SofTofW(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SofTofW::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SofTofW::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* SofTofW::BuildTOFDetector(){
  m_TofWall = new G4AssemblyVolume();

  if(!m_PlasticTof){
    G4Box* box = new G4Box("SofTofW_Box",SofTofW_NS::tof_plastic_height*0.5,
        SofTofW_NS::tof_plastic_width*0.5,SofTofW_NS::tof_plastic_thickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(SofTofW_NS::Material);
    m_PlasticTof = new G4LogicalVolume(box,DetectorMaterial,"logic_SofTofW_Box",0,0,0);
    m_PlasticTof->SetVisAttributes(m_VisSquare);
    m_PlasticTof->SetSensitiveDetector(m_TofScorer);

    G4RotationMatrix* Rv = new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0);
    Tv.setY(0);
    Tv.setZ(0);
    for(unsigned int i=0; i<28; i++){
      int k = -14+i;
      Tv.setY(k*(SofTofW_NS::tof_plastic_width+0.5));
      m_TofWall->AddPlacedVolume(m_PlasticTof, Tv, Rv);
    }

  }
  return m_TofWall;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*G4LogicalVolume* SofTofW::BuildTwinMusic(){

  double twin_width = 21*cm;
  double twin_height = 21*cm;
  double twin_thickness = 51*cm;

  double sector_width = 10*cm;
  double sector_height = 10*cm;
  double sector_thickness = 50*cm;

  if(m_Build_Twin_Music==1){
  if(!m_TwinMusic){
// Full Twin volume
G4Box* Twinbox = new G4Box("Twin_Box", twin_width*0.5, twin_height*0.5, twin_thickness*0.5);
G4Material* TwinMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
m_TwinMusic = new G4LogicalVolume(Twinbox, TwinMaterial, "logic_twin", 0,0,0);

G4VisAttributes* m_VisDet = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.3));
m_VisDet->SetForceWireframe(1);
m_TwinMusic->SetVisAttributes(m_VisDet);

// Sector Twin volume
G4Box* Sectorbox = new G4Box("Sector_Box", sector_width*0.5, sector_height*0.5, sector_thickness*0.5);
G4LogicalVolume* LogicalSector = new G4LogicalVolume(Sectorbox, TwinMaterial, "logic_twin", 0,0,0);
LogicalSector->SetVisAttributes(G4VisAttributes::GetInvisible());
// Drift Anode Area //
G4Box* Anodebox = new G4Box("Anode_Box", SofTofW_NS::twin_anode_width*0.5, SofTofW_NS::twin_anode_height*0.5, SofTofW_NS::twin_anode_thickness*0.5);

G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_Twin_Music_Gas);
m_AnodeDriftArea = new G4LogicalVolume(Anodebox, DetectorMaterial, "logic_twin_anode", 0, 0, 0);
m_AnodeDriftArea->SetVisAttributes(m_VisTwin);
m_AnodeDriftArea->SetSensitiveDetector(m_TwinScorer);

// Cathode plane in the middle //
G4Box* cathode_box = new G4Box("Cathode_box", SofTofW_NS::twin_cathode_width*0.5, SofTofW_NS::twin_cathode_height*0.5, SofTofW_NS::twin_cathode_thickness*0.5);
G4Material* CathodeMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
G4VisAttributes* m_VisCathode = new G4VisAttributes(G4Colour(0.7,0.4,0,1));
G4LogicalVolume* LogicalCathode = new G4LogicalVolume(cathode_box, CathodeMaterial, "logic_cathode", 0,0,0);
LogicalCathode->SetVisAttributes(m_VisCathode);

//G4RotationMatrix* Rv = new G4RotationMatrix(0,0,0);
G4ThreeVector Tv;
Tv.setX(0);
Tv.setY(0);
Tv.setZ(0);


new G4PVPlacement(0, Tv,
LogicalCathode,
"Cathode",
m_TwinMusic, false, 0);

int anode_nbr= 0;
for(unsigned int j=0; j<16; j++){
anode_nbr++;
Tv.setZ(j*SofTofW_NS::twin_anode_thickness -7.5*SofTofW_NS::twin_anode_thickness);
new G4PVPlacement(0,Tv,
m_AnodeDriftArea,
"Anode",
LogicalSector, false, anode_nbr);
}

for(unsigned int i=0; i<4; i++){
if(i==0){
Tv.setX(0.5*SofTofW_NS::twin_anode_width+SofTofW_NS::twin_cathode_width);
Tv.setY(0.5*SofTofW_NS::twin_anode_height);
}
if(i==1){
Tv.setX(-0.5*SofTofW_NS::twin_anode_width-SofTofW_NS::twin_cathode_width);
Tv.setY(0.5*SofTofW_NS::twin_anode_height);
}
if(i==2){
  Tv.setX(-0.5*SofTofW_NS::twin_anode_width-SofTofW_NS::twin_cathode_width);
  Tv.setY(-0.5*SofTofW_NS::twin_anode_height);
}
if(i==3){
  Tv.setX(0.5*SofTofW_NS::twin_anode_width+SofTofW_NS::twin_cathode_width);
  Tv.setY(-0.5*SofTofW_NS::twin_anode_height);
}

Tv.setZ(0);
new G4PVPlacement(0,Tv,
    LogicalSector,
    "Sector",
    m_TwinMusic,false,i+1);
}
}
}

return m_TwinMusic;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* SofTofW::BuildVacuumPipe(){
  if(!m_VacuumPipe){
    m_VacuumPipe = new G4AssemblyVolume;

    G4Tubs* tube = new G4Tubs("tube",8.*cm,15*cm,155./2*cm,0,360*deg);

    G4Material* tube_mat = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");

    G4LogicalVolume* tube_vol = new G4LogicalVolume(tube,tube_mat,"logic_tube",0,0,0);
    G4ThreeVector Pos = G4ThreeVector(0,0,0);
    G4RotationMatrix* Rot = new G4RotationMatrix();
    m_VacuumPipe->AddPlacedVolume(tube_vol,Pos,Rot);
  }
  
  return m_VacuumPipe;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* SofTofW::BuildGLAD()
{
  if(!m_GLAD){
    m_GLAD = new G4AssemblyVolume;
    string basepath = getenv("NPTOOL");
    //string path = basepath + "/NPSimulation/Detectors/Sofia/gdml/glad.gdml";
    //m_gdmlparser.Read(path);

    //G4LogicalVolume* vol1 = m_gdmlparser.GetVolume("GEcrans");
    //G4LogicalVolume* vol2 = m_gdmlparser.GetVolume("G2202001_Demi_Ecran_thermique_interne");
    
    //G4LogicalVolume* vol3 = m_gdmlparser.GetVolume("G2402001_Enceinte_interne");
    
    //G4LogicalVolume* vol4 = m_gdmlparser.GetVolume("GEnceinte_externe");
    //G4LogicalVolume* vol5 = m_gdmlparser.GetVolume("G2403002_Fonf_cote_sortie");
    //G4LogicalVolume* vol6 = m_gdmlparser.GetVolume("G2403001_Fond_cote_entree");
    //G4LogicalVolume* vol7 = m_gdmlparser.GetVolume("GToles");

    // *** GLAD field *** //    
    G4Box* box = new G4Box("glad_Box",SofTofW_NS::GLAD_width*0.5,SofTofW_NS::GLAD_height*0.5,SofTofW_NS::GLAD_Leff*0.5);
    G4Material* vac = MaterialManager::getInstance()->GetMaterialFromLibrary("Vaccuum");
    G4LogicalVolume* vol_field = new G4LogicalVolume(box,vac,"logic_GLAD_field",0,0,0);
    vol_field->SetVisAttributes(m_VisGLAD);

    G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0,m_GLAD_MagField,0));
    //G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    G4FieldManager* fieldMgr = new G4FieldManager(magField);

    //fieldMgr->SetDetectorField(magField);

    fieldMgr->CreateChordFinder(magField);

    vol_field->SetFieldManager(fieldMgr,true);

    // *** vol1 *** //
    //G4ThreeVector Pos1 = G4ThreeVector(0,0,0);
    //G4RotationMatrix* Rot1 = new G4RotationMatrix();
    //m_GLAD->AddPlacedVolume(vol1,Pos1,Rot1);

    // *** vol2 *** //
    //G4ThreeVector Pos2 = G4ThreeVector(0,0,0);
    //G4RotationMatrix* Rot2 = new G4RotationMatrix();
    //m_GLAD->AddPlacedVolume(vol2,Pos2,Rot2);
   
    // *** vol3 *** //
    /*G4ThreeVector Pos3 = G4ThreeVector(0,0,0);
    G4RotationMatrix* Rot3 = new G4RotationMatrix();
    Rot3->rotateX(90*deg);
    Rot3->rotateY(90*deg);
    m_GLAD->AddPlacedVolume(vol3,Pos3,Rot3);*/
   
    // *** vol4 *** //
    //G4ThreeVector Pos4 = G4ThreeVector(0*cm,0,0*cm);
    //G4RotationMatrix* Rot4 = new G4RotationMatrix();
    //Rot4->rotateY(180*deg);
    //Rot4->rotateZ(90*deg); 
    //m_GLAD->AddPlacedVolume(vol4,Pos4,Rot4);
    
    // *** vol5 *** //
    //G4VisAttributes* Vis_vol5 = new G4VisAttributes(G4Colour(0,0,1,0.5));
    //vol5->SetVisAttributes(Vis_vol5);
    //G4ThreeVector Pos5 = G4ThreeVector(0,0,0);
    //G4RotationMatrix* Rot5 = new G4RotationMatrix(); 
    //Rot5->rotateY(90*deg);
    //Rot5->rotateZ(90*deg);
    //m_GLAD->AddPlacedVolume(vol5,Pos5,Rot5);

    // *** vol6 ***//
    //G4ThreeVector Pos6 = G4ThreeVector(0,0,0);
    //G4RotationMatrix* Rot6 = new G4RotationMatrix(); 
    //m_GLAD->AddPlacedVolume(vol6,Pos6,Rot6);
   
    // *** vol7 *** //
    //G4ThreeVector Pos7 = G4ThreeVector(0,0,0);
    //G4RotationMatrix* Rot7 = new G4RotationMatrix(); 
    //m_GLAD->AddPlacedVolume(vol7,Pos7,Rot7);


    G4ThreeVector Pos_field = G4ThreeVector(0,0,1.5*m);
    G4RotationMatrix* Rot_field = new G4RotationMatrix();
    m_GLAD->AddPlacedVolume(vol_field,Pos_field,Rot_field);

  }
  return m_GLAD;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*G4AssemblyVolume* SofTofW::BuildGLAD(){
  if(!m_GLAD){
  m_GLAD = new G4AssemblyVolume;
  G4Box* box = new G4Box("glad_Box",SofTofW_NS::GLAD_width*0.5,
  SofTofW_NS::GLAD_height*0.5,SofTofW_NS::GLAD_Leff*0.5);

  G4Material* GLADMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vaccuum");
  G4LogicalVolume* vol1 = new G4LogicalVolume(box,GLADMaterial,"logic_GLAD_Box",0,0,0);
  vol1->SetVisAttributes(m_VisGLAD);

  G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0,m_GLAD_MagField,0));
//G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
G4FieldManager* fieldMgr = new G4FieldManager(magField);

//fieldMgr->SetDetectorField(magField);

fieldMgr->CreateChordFinder(magField);

vol1->SetFieldManager(fieldMgr,true);

G4ThreeVector Pos1 = G4ThreeVector(0,0,0);
G4RotationMatrix* Rot1 = new G4RotationMatrix(0,0,0);
m_GLAD->AddPlacedVolume(vol1,Pos1,Rot1);
}
return m_GLAD;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void SofTofW::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofTofW");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Build_GLAD","Build_VacuumPipe"};
  vector<string> sphe = {"R","Theta","Phi","Build_GLAD","Build_VacuumPipe"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTofW " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      m_Build_GLAD = blocks[i]->GetInt("Build_GLAD");
      m_Build_VacuumPipe = blocks[i]->GetInt("Build_VacuumPipe");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTofW " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      m_Build_GLAD = blocks[i]->GetInt("Build_GLAD");
      m_Build_VacuumPipe = blocks[i]->GetInt("Build_VacuumPipe");
      m_GLAD_MagField = blocks[i]->GetDouble("GLAD_MagField","T");
      m_GLAD_DistanceFromTarget = blocks[i]->GetDouble("GLAD_DistanceFromTarget", "m");
      m_VacuumPipeX = blocks[i]->GetDouble("VacuumPipeX","m");
      m_VacuumPipeY = blocks[i]->GetDouble("VacuumPipeY","m");
      m_VacuumPipeZ = blocks[i]->GetDouble("VacuumPipeZ","m");

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
void SofTofW::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*SofTofW_NS::tof_plastic_thickness*0.5;
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

    BuildTOFDetector()->MakeImprint(world,Det_pos,Rot);
  }

  if(m_Build_GLAD==1){
    G4ThreeVector GLAD_pos = G4ThreeVector(0,0,m_GLAD_DistanceFromTarget);
    //G4ThreeVector GLAD_pos = G4ThreeVector(0,0,0);
    BuildGLAD()->MakeImprint(world,GLAD_pos,0);
    /*new G4PVPlacement(0, GLAD_pos,
      BuildGLAD(),
      "GLAD",
      world, false, 0);
     */
  }
  if(m_Build_VacuumPipe==1){
    G4ThreeVector Tube_Pos = G4ThreeVector(m_VacuumPipeX,m_VacuumPipeY,m_VacuumPipeZ);
    BuildVacuumPipe()->MakeImprint(world,Tube_Pos,0);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SofTofW::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("SofTofW")){
    pTree->Branch("SofTofW", "TSofTofWData", &m_Event) ;
  }
  pTree->SetBranchAddress("SofTofW", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void SofTofW::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // TOF scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_TofScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),SofTofW_NS::ResoEnergy);
    if(Energy>SofTofW_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),SofTofW_NS::ResoTime);
      //int DetectorNbr = level[0];
      int PlasticNbr = level[1]-1;
      //m_Event->SetDetectorNbr(DetectorNbr);
      m_Event->SetPlasticNbr(PlasticNbr);
      m_Event->SetEnergy(Energy);
      m_Event->SetCoarseTime(Time); 
    }
  }
  Scorer->clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void SofTofW::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_TofScorer = CheckScorer("TofScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; 
  level.push_back(1);
  level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  //G4VPrimitiveScorer* TwinCalorimeter= new CalorimeterScorers::PS_Calorimeter("TwinCalorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_TofScorer->RegisterPrimitive(Calorimeter);
  m_TofScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_TofScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* SofTofW::Construct(){
  return  (NPS::VDetector*) new SofTofW();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_SofTofW{
    public:
      proxy_nps_SofTofW(){
        NPS::DetectorFactory::getInstance()->AddToken("SofTofW","SofTofW");
        NPS::DetectorFactory::getInstance()->AddDetector("SofTofW",SofTofW::Construct);
      }
  };

  proxy_nps_SofTofW p_nps_SofTofW;
}
