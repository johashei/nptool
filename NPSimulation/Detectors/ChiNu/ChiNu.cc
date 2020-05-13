/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2019                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  ChiNu simulation                             *
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
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"

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
#include "ChiNu.hh"
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
namespace ChiNu_NS{
  // EJ309 Scintillator - Energy and time Resolution
  const double EnergyThreshold = 0.*MeV;
  const double ResoTime = 1*ns ;
  const double ResoEnergy = 0.1*MeV ;
  const double Radius = 8.90*cm ; 
  const double Thickness = 5.08*cm ; 
  const string Material = "EJ309";

  // PMT
  const double PMT_Height = 392*mm; 
  const string PMT_Material = "Al"; 
 
  // Light guide
  const double LG_Rmin1 = 0*mm;
  const double LG_Rmax1 = Radius;
  const double LG_Rmin2 = 0*mm; 
  const double LG_Rmax2 = 50*mm; 
  const double LG_Thickness = 30*mm; 
  const string LG_Material = "PMMA";

  // Pyrex Window
  const double Pyrex_radius = Radius;
  const double Pyrex_thickness = 6.4*mm;
  const string Pyrex_material = "Pyrex"; 

  // Lead shield
  const double Lead_Radius = 9*cm;
  const double Lead_Thickness = 2*mm;

  // Fission Chamber
  const string FCWall_Material = "CH2";
  const double Cu_Thickness = 17*micrometer;
 // const double Al_Thickness = 12*micrometer;
  const double Kapton_Thickness = 50*micrometer;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ChiNu Specific Method
ChiNu::ChiNu(){
  m_Event = new TChiNuData() ;
  m_ChiNuScorer = 0;
  m_PMT = 0;
  m_LightGuide = 0;
  m_CylindricalDetector = 0;
  m_LeadShield = 0;
  m_AssemblyVolume = 0;

  m_BuildLeadShield = 0;

  m_FissionChamberWall = 0;
  m_FissionChamberVolume = 0;

  // RGB Color + Transparency
  m_VisCylinder = new G4VisAttributes(G4Colour(0.0, 0.5, 1, 1));   
  m_VisPMT = new G4VisAttributes(G4Colour(0.3, 0.1, 0.1, 0.3));   
  m_VisLightGuide = new G4VisAttributes(G4Colour(0.1,0.5,0.7,1));
  m_VisPyrex = new G4VisAttributes(G4Colour(0.1,0.5,0.7,0.7));
  m_VisLeadShield = new G4VisAttributes(G4Colour(0.2,0.2,0.2,1));
  m_VisFCWall = new G4VisAttributes(G4Colour(0.1,0.5,0.7,1));
  m_VisAl = new G4VisAttributes(G4Colour(0.839,0.803,0.803,1));
  m_VisTi = new G4VisAttributes(G4Colour(0.776,0.662,0.662,0.5));
  m_VisCu = new G4VisAttributes(G4Colour(0.70, 0.40, 0. ,1));
  m_VisRogers4003C = new G4VisAttributes(G4Colour(0.60, 0.60, 0.2 ,1));
}

ChiNu::~ChiNu(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChiNu::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChiNu::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* ChiNu::BuildFissionChamber(){
  if(!m_FissionChamberVolume){
    m_FissionChamberVolume = new G4AssemblyVolume();
    
    G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0); Tv.setY(0); Tv.setZ(0);

    // Bottom PCB plate //
    double PCB_width = 18.*cm;
    double PCB_length = 33.*cm;
    double PCB_Rogers_height = 1.6*mm;
    double PCB_Cu_height = 6*35.*um;
    double PCB_PosY = -8.5*cm; 
    // Cu layers
    G4Box* PCB_Cu_solid = new G4Box("PCB_Cu_solid",0.5*PCB_width,0.5*PCB_Cu_height,0.5*PCB_length); 
    G4Material* Cu_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
    G4LogicalVolume* PCB_Cu_vol = new G4LogicalVolume(PCB_Cu_solid, Cu_material,"PCB_Cu_logic",0,0,0);
    PCB_Cu_vol->SetVisAttributes(m_VisCu);
    Tv.setY(PCB_PosY);
    m_FissionChamberVolume->AddPlacedVolume(PCB_Cu_vol, Tv, Rv);

    // Rogers 4003C layers
    G4Box* PCB_Rogers_solid = new G4Box("PCB_Rogers_solid",0.5*PCB_width,0.5*PCB_Rogers_height,0.5*PCB_length); 
    G4Material* Rogers_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Rogers4003C");
    G4LogicalVolume* PCB_Rogers_vol = new G4LogicalVolume(PCB_Rogers_solid, Rogers_material,"PCB_Rogers_logic",0,0,0);
    PCB_Rogers_vol->SetVisAttributes(m_VisRogers4003C);
    Tv.setY(PCB_PosY + 0.5*PCB_Cu_height + 0.5*PCB_Rogers_height);
    m_FissionChamberVolume->AddPlacedVolume(PCB_Rogers_vol, Tv, Rv);

    //Al frame //
    double frame1_width = 18.*cm;
    double frame1_height = 5.*mm;
    double frame1_length = 33.*cm;
    double frame2_width = 15.2*cm;
    double frame2_height = 5.1*mm;
    double frame2_length = 30.2*cm;

    G4Box* frame1 = new G4Box("frame1", 0.5*frame1_width, 0.5*frame1_height, 0.5*frame1_length);
    G4Box* frame2 = new G4Box("frame2", 0.5*frame2_width, 0.5*frame2_height, 0.5*frame2_length);
    G4VSolid* Al_frame = (G4VSolid*) new G4SubtractionSolid("Al_frame",frame1,frame2,0,G4ThreeVector(0,0,0));
    G4Material* Al_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4LogicalVolume* Al_frame_vol = new G4LogicalVolume(Al_frame,Al_material,"Al_frame_logic",0,0,0);
    Al_frame_vol->SetVisAttributes(m_VisAl);
    Tv.setY(PCB_PosY+ 0.5*PCB_Cu_height + PCB_Rogers_height + 0.5*frame1_height);
    m_FissionChamberVolume->AddPlacedVolume(Al_frame_vol, Tv, Rv);
    Tv.setY(PCB_PosY- 0.5*PCB_Cu_height - 0.5*frame1_height);
    m_FissionChamberVolume->AddPlacedVolume(Al_frame_vol, Tv, Rv);
    
    double box1_width = 15.*cm;
    double box1_height = 16.*cm;
    double box1_length = 30.*cm;
    double box2_width = 14.8*cm;
    double box2_height = 15.8*cm;
    double box2_length = 29.8*cm;
    double box3_width = 12.5*cm;
    double box3_height = 11.7*cm;
    double box3_length = 30.1*cm;
    double box4_width = 15.1*cm;
    double box4_height = 11.8*cm;
    double box4_length = 27.4*cm;
    double box5_width = 12.4*cm;
    double box5_height = 16.1*cm;
    double box5_length = 27.4*cm;

    G4Box* box1 = new G4Box("box1", 0.5*box1_width, 0.5*box1_height, 0.5*box1_length);
    G4Box* box2 = new G4Box("box2", 0.5*box2_width, 0.5*box2_height, 0.5*box2_length);
    G4Box* box3 = new G4Box("box3", 0.5*box3_width, 0.5*box3_height, 0.5*box3_length);
    G4Box* box4 = new G4Box("box4", 0.5*box4_width, 0.5*box4_height, 0.5*box4_length);
    G4Box* box5 = new G4Box("box5", 0.5*box5_width, 0.5*box5_height, 0.5*box5_length);

    G4VSolid* box_int1 = (G4VSolid*) new G4SubtractionSolid("box_int1",box1,box2,0,G4ThreeVector(0,0,0));
    G4VSolid* box_int2 = (G4VSolid*) new G4SubtractionSolid("box_int2",box_int1,box3,0,G4ThreeVector(0,0,0));
    G4VSolid* box_int3 = (G4VSolid*) new G4SubtractionSolid("box_int3",box_int2,box4,0,G4ThreeVector(0,0,0));
    G4VSolid* box_int4 = (G4VSolid*) new G4SubtractionSolid("box_int4",box_int3,box5,0,G4ThreeVector(0,0,0));

    G4LogicalVolume* full_box_vol = new G4LogicalVolume(box_int4, Al_material, "full_box_logic", 0,0,0);
    full_box_vol->SetVisAttributes(m_VisAl);
    Tv.setY(0);
    m_FissionChamberVolume->AddPlacedVolume(full_box_vol, Tv, Rv);

    // Ti foils //
    double foil1_width = 13*cm;
    double foil1_length = 29*cm;
    double foil1_thickness = 100*um;
    double foil2_width = 13*cm;
    double foil2_height = 14*cm;
    double foil2_thickness = 50*um;

    G4Material* Ti_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Ti");
    G4Box* foil1_solid = new G4Box("foil1", 0.5*foil1_width, 0.5*foil1_thickness, 0.5*foil1_length);
    G4LogicalVolume* foil1_vol = new G4LogicalVolume(foil1_solid, Ti_material, "foil1_logic", 0, 0, 0);
    foil1_vol->SetVisAttributes(m_VisTi);
    Tv.setY(0.5*box2_height);
    m_FissionChamberVolume->AddPlacedVolume(foil1_vol, Tv, Rv);
    Tv.setY(0);
    Tv.setX(-0.5*box2_width);
    Rv->rotateZ(90*deg);
    m_FissionChamberVolume->AddPlacedVolume(foil1_vol, Tv, Rv);
    Tv.setX(0.5*box2_width);
    m_FissionChamberVolume->AddPlacedVolume(foil1_vol, Tv, Rv);

    G4Box* foil2_solid = new G4Box("foil2", 0.5*foil2_width, 0.5*foil2_height, 0.5*foil2_thickness);
    G4LogicalVolume* foil2_vol = new G4LogicalVolume(foil2_solid, Ti_material, "foil2_logic", 0, 0, 0); 
    foil2_vol->SetVisAttributes(m_VisTi);
    Tv.setX(0);Tv.setY(0);Tv.setZ(-0.5*box2_length);
    m_FissionChamberVolume->AddPlacedVolume(foil2_vol, Tv, Rv);
    Tv.setZ(0.5*box2_length);
    m_FissionChamberVolume->AddPlacedVolume(foil2_vol, Tv, Rv);

    // Cathode and Anode //
    BuildCathode(-27.5);
    double origine_anode = -25*mm;
    double origine_cathode = -22.5*mm;
    for(int i=0; i<11; i++){
	BuildAnode(origine_anode+i*5*mm); 
    }
    for(int i=0; i<10; i++){
	BuildCathode(origine_cathode+i*5*mm);
    }
    BuildCathode(27.5);


  }
  return m_FissionChamberVolume;
} 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChiNu::BuildCathode(double Zpos){
    // Al plate: 12 um
    G4Tubs* Al_plate_solid = new G4Tubs("Al_plate",0,40*mm,12*micrometer,0,360*deg);
    G4Material* Al_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4LogicalVolume* Al_vol = new G4LogicalVolume(Al_plate_solid, Al_material,"logic_Al",0,0,0);
    Al_vol->SetVisAttributes(m_VisAl);
     
    G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0); Tv.setY(0); Tv.setZ(0);

    Tv.setZ(Zpos);
    m_FissionChamberVolume->AddPlacedVolume(Al_vol, Tv, Rv);
    
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChiNu::BuildAnode(double Zpos){
    // Cu plate: 17 um
    G4Tubs* Cu_plate_solid = new G4Tubs("Cu_plate",0,40*mm,0.5*ChiNu_NS::Cu_Thickness,0,360*deg);
    G4Material* Cu_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
    G4LogicalVolume* Cu_vol = new G4LogicalVolume(Cu_plate_solid, Cu_material,"logic_Cu",0,0,0);
    Cu_vol->SetVisAttributes(m_VisCu);
   
    // Kapton: 50 um
    G4Tubs* Kapton_solid = new G4Tubs("Kapton",0,40*mm,0.5*ChiNu_NS::Kapton_Thickness,0,360*deg);
    G4Material* Kapton_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Kapton");
    G4LogicalVolume* Kapton_vol = new G4LogicalVolume(Kapton_solid, Kapton_material,"logic_Kapton",0,0,0);
    Kapton_vol->SetVisAttributes(m_VisFCWall);
    
    G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0); Tv.setY(0); Tv.setZ(0);

    Tv.setZ(Zpos);
    m_FissionChamberVolume->AddPlacedVolume(Kapton_vol, Tv, Rv);
    Tv.setZ(Zpos-0.5*ChiNu_NS::Kapton_Thickness-0.5*ChiNu_NS::Cu_Thickness);
    m_FissionChamberVolume->AddPlacedVolume(Cu_vol, Tv, Rv);
    Tv.setZ(Zpos+0.5*ChiNu_NS::Kapton_Thickness+0.5*ChiNu_NS::Cu_Thickness);
    m_FissionChamberVolume->AddPlacedVolume(Cu_vol, Tv, Rv);


} 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* ChiNu::BuildDetector(){
  if(!m_CylindricalDetector){
    m_AssemblyVolume = new G4AssemblyVolume();
    G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0); Tv.setY(0); Tv.setZ(0);

    // Scintillator
    G4Tubs* tub = new G4Tubs("ChiNu_Cyl",0,ChiNu_NS::Radius,ChiNu_NS::Thickness*0.5,0,360*deg);
    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::Material);
    m_CylindricalDetector = new G4LogicalVolume(tub,DetectorMaterial,"logic_ChiNu_tub",0,0,0);
    m_CylindricalDetector->SetVisAttributes(m_VisCylinder);
    m_CylindricalDetector->SetSensitiveDetector(m_ChiNuScorer);
    m_AssemblyVolume->AddPlacedVolume(m_CylindricalDetector, Tv, Rv);

    // Pyrex Window
    G4Tubs* Pyrex_tub = new G4Tubs("Pyrex_tub",0,ChiNu_NS::Pyrex_radius, ChiNu_NS::Pyrex_thickness*0.5, 0 , 360*deg);
    G4Material* PyrexMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::Pyrex_material);
    G4LogicalVolume* LogPyrex = new G4LogicalVolume(Pyrex_tub, PyrexMaterial,"logic_pyrex",0,0,0);
    LogPyrex->SetVisAttributes(m_VisPyrex);
    Tv.setZ(ChiNu_NS::Thickness*0.5 + ChiNu_NS::Pyrex_thickness*0.5);
    m_AssemblyVolume->AddPlacedVolume(LogPyrex, Tv, Rv);


    // Light guide
    G4Cons* LGsolid = new G4Cons("light_guide", ChiNu_NS::LG_Rmin1, ChiNu_NS::LG_Rmax1, ChiNu_NS::LG_Rmin2, ChiNu_NS::LG_Rmax2, ChiNu_NS::LG_Thickness*0.5,0,360*deg); 
    G4Material* LGMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::LG_Material);
    m_LightGuide = new G4LogicalVolume(LGsolid,LGMaterial,"logic_light_guide",0,0,0);
    m_LightGuide->SetVisAttributes(m_VisLightGuide);
    Tv.setZ(ChiNu_NS::Thickness*0.5 + ChiNu_NS::Pyrex_thickness + ChiNu_NS::LG_Thickness*0.5);
    m_AssemblyVolume->AddPlacedVolume(m_LightGuide, Tv, Rv);

    // PMT
    //G4Tubs* pmt = new G4Tubs("ChiNu_pmt",ChiNu_NS::PMT_InnerDiameter*0.5,ChiNu_NS::PMT_OuterDiameter*0.5,ChiNu_NS::PMT_Thickness*0.5,0,360*deg);
    
    double zplane[4] ={0, 18*cm, 24*cm, ChiNu_NS::PMT_Height};
    double rin[4] = {ChiNu_NS::Radius+0.2*mm, ChiNu_NS::Radius+0.2*mm, 40*mm, 40*mm};
    double rout[4] = {ChiNu_NS::Radius+2.2*mm, ChiNu_NS::Radius+2.2*mm, 42*mm, 42*mm};
    G4Polycone* pmt = new G4Polycone("ChiNu_PMT", 0, 360*deg, 4, zplane, rin, rout);
    G4Material* pmtMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::PMT_Material);
    m_PMT = new G4LogicalVolume(pmt,pmtMaterial,"logic_pmt_tub",0,0,0);
    m_PMT->SetVisAttributes(m_VisPMT);
    //Tv.setZ(ChiNu_NS::Thickness*0.5 + ChiNu_NS::LG_Thickness + ChiNu_NS::Pyrex_thickness + ChiNu_NS::PMT_Height*0.5);
    Tv.setZ(-ChiNu_NS::Thickness*0.5);
    m_AssemblyVolume->AddPlacedVolume(m_PMT, Tv, Rv);

    // Lead shield
    if(m_BuildLeadShield){
    	G4Tubs* lead = new G4Tubs("lead_shield", 0, ChiNu_NS::Lead_Radius, ChiNu_NS::Lead_Thickness*0.5, 0, 360*deg);
    	G4Material* LeadMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Pb");
	m_LeadShield = new G4LogicalVolume(lead, LeadMaterial, "logic_lead_shield",0,0,0);
	m_LeadShield->SetVisAttributes(m_VisLeadShield);
	Tv.setZ(-ChiNu_NS::Thickness*0.5 - ChiNu_NS::Lead_Thickness*0.5-10*mm);
        m_AssemblyVolume->AddPlacedVolume(m_LeadShield, Tv, Rv);
    }

  }
  return m_AssemblyVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void ChiNu::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ChiNu");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  ChiNu " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  ChiNu " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      m_BuildLeadShield = blocks[i]->GetInt("LeadShield");
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
void ChiNu::ConstructDetector(G4LogicalVolume* world){
  
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*ChiNu_NS::Thickness*0.5;
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
    BuildDetector()->MakeImprint(world,Det_pos,Rot,i);
  }

  G4RotationMatrix* Rot_FC = new G4RotationMatrix(0,0,0);
  G4ThreeVector Pos_FC = G4ThreeVector(0,0,0) ;
  BuildFissionChamber()->MakeImprint(world,Pos_FC,Rot_FC,0);


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void ChiNu::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("ChiNu")){
    pTree->Branch("ChiNu", "TChiNuData", &m_Event) ;
  }
  pTree->SetBranchAddress("ChiNu", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void ChiNu::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_ChiNuScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),ChiNu_NS::ResoEnergy);
    if(Energy>ChiNu_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),ChiNu_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void ChiNu::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_ChiNuScorer = CheckScorer("ChiNuScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_ChiNuScorer->RegisterPrimitive(Calorimeter);
  m_ChiNuScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ChiNuScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* ChiNu::Construct(){
  return  (NPS::VDetector*) new ChiNu();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_ChiNu{
    public:
      proxy_nps_ChiNu(){
        NPS::DetectorFactory::getInstance()->AddToken("ChiNu","ChiNu");
        NPS::DetectorFactory::getInstance()->AddDetector("ChiNu",ChiNu::Construct);
      }
  };

  proxy_nps_ChiNu p_nps_ChiNu;
}
