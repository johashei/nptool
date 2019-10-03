/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2019                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Mugast simulation                                   *
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
#include "Mugast.hh"
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
namespace Mugast_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 4.5*ns ;
  const double ResoEnergy = 1.0*MeV ;
  const double Radius = 50*mm ; 
  const double Width = 100*mm ;
  const double Thickness = 300*mm ;
  const string Material = "BC400";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Mugast Specific Method
Mugast::Mugast(){
  m_Event = new TMugastData() ;
  m_MugastScorer = 0;
  m_SquareDetector = 0;
  m_TrapezoidDetector = 0;
  m_AnnularDetector = 0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));   
  m_VisCylinder = new G4VisAttributes(G4Colour(0, 0, 1, 0.5));   

}

Mugast::~Mugast(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Mugast::AddDetector(string Shape,G4ThreeVector PX1_Y1 ,G4ThreeVector PX1_Y128 ,G4ThreeVector PX128_Y1,G4ThreeVector PX128_Y128){
   m_X1_Y1.push_back(PX1_Y1); // Top Left Corner Position Vector
   m_X1_Y128.push_back(PX1_Y128); // Bottom Left Corner Position Vector
   m_X128_Y1.push_back(PX128_Y1); // Bottom Right Corner Position Vector
   m_X128_Y128.push_back(PX128_Y128); // Center Corner Position Vector
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Mugast::BuildSquareDetector(){
  if(!m_SquareDetector){
    G4Box* box = new G4Box("Mugast_Box",Mugast_NS::Width*0.5,
        Mugast_NS::Width*0.5,Mugast_NS::Thickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Mugast_NS::Material);
    m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_Mugast_Box",0,0,0);
    m_SquareDetector->SetVisAttributes(m_VisSquare);
    m_SquareDetector->SetSensitiveDetector(m_MugastScorer);
  }
  return m_SquareDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Mugast::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Mugast");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector found" << endl;
    
  // Cartesian Case
  vector<string> cart = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128"};
  
  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (NPOptionManager::getInstance()->GetVerboseLevel())
      cout << endl << "//// Mugast detector " << i + 1 << endl;
    
    string shape = blocks[i]->GetMainValue();

    if (blocks[i]->HasTokenList(cart)&& (shape=="Square"|| shape=="Trapezoid")) {
      G4ThreeVector A
        = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1", "mm"));
      G4ThreeVector B
        = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y1", "mm"));
      G4ThreeVector C
        = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y128", "mm"));
      G4ThreeVector D
        = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y128", "mm"));
      AddDetector(shape,A,B,C,D);
    }
    else if (blocks[i]->HasTokenList(cart)&& (shape=="Annular")) {
      G4ThreeVector A
        = NPS::ConvertVector(blocks[i]->GetTVector3("Center", "mm"));
      AddDetector(shape,A);
    }

    else {
      cout << "WARNING: Missing token for Mugast blocks, check your input "
        "file"
        << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Mugast::ConstructDetector(G4LogicalVolume* world){
/*  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*Mugast_NS::Thickness*0.5;
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

    if(m_Shape[i] == "Square"){
      new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildSquareDetector(),
          "Mugast",world,false,i+1);
    }
  }

  */
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Mugast::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Mugast")){
    pTree->Branch("Mugast", "TMugastData", &m_Event) ;
  }
  pTree->SetBranchAddress("Mugast", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Mugast::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_MugastScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Mugast_NS::ResoEnergy);
    if(Energy>Mugast_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Mugast_NS::ResoTime);
      int DetectorNbr = level[0];
   //   m_Event->SetEnergy(DetectorNbr,Energy);
   //   m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Mugast::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_MugastScorer = CheckScorer("MugastScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_MugastScorer->RegisterPrimitive(Calorimeter);
  m_MugastScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_MugastScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Mugast::Construct(){
  return  (NPS::VDetector*) new Mugast();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Mugast{
    public:
      proxy_nps_Mugast(){
        NPS::DetectorFactory::getInstance()->AddToken("Mugast","Mugast");
        NPS::DetectorFactory::getInstance()->AddDetector("Mugast",Mugast::Construct);
      }
  };

  proxy_nps_Mugast p_nps_Mugast;
}
