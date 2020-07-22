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
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"

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
  double DummyInnerRadius = 200*mm ; 
  double DummyOuterRadius = 600*mm ; 
  double Length = 600*mm ;
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
  m_DetectorType1 = 0;
  m_DetectorType2 = 0;
  m_DetectorType3 = 0;

  // RGB Color + Transparency
  m_VisCrystal = new G4VisAttributes(G4Colour(1, 0.8, 0, 0.2));   

}

Catana::~Catana(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Catana::AddDummyDetector(double Z){
  m_Z.push_back(Z);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Catana::AddDetectorType1(double R, double Theta, double Phi){
  m_R1.push_back(R);
  m_Theta1.push_back(Theta);
  m_Phi1.push_back(Phi);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Catana::AddDetectorType2(double R, double Theta, double Phi){
  m_R2.push_back(R);
  m_Theta2.push_back(Theta);
  m_Phi2.push_back(Phi);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Catana::AddDetectorType3(double R, double Theta, double Phi){
  m_R3.push_back(R);
  m_Theta3.push_back(Theta);
  m_Phi3.push_back(Phi);
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
G4LogicalVolume* Catana::BuildDetectorType1(){

  if(!m_DetectorType1){
   
  G4Material* CsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
  G4Material* Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4Material* Vacuum= MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4Material* Teflon= MaterialManager::getInstance()->GetMaterialFromLibrary("G4_TEFLON");
  G4RotationMatrix* Rot= new G4RotationMatrix();
  // -- Al housing outside --
  double Dz = 98/2.*mm;
  double Theta = 0.*deg;
  double Phi = 0.*deg;
  double Dy1 = 36.6/2*mm;
  double Dx1 = 62.3/2*mm;
  double Dx2 = Dx1;
  double Alp1 = 0*deg;
  double Dy2 = 54.0/2*mm;
  double Dx3 = 92.1/2*mm;
  double Dx4 = Dx3;
  double Alp2 = 0*deg;

  G4Trap *solidDetectorType1 = new G4Trap("CatanaDetectorType1",Dz, Theta, Phi,
					Dy1, Dx1,Dx2, Alp1, Dy2, Dx3,
					Dx4, Alp2);

  m_DetectorType1  = new G4LogicalVolume(solidDetectorType1,
					   Vacuum,
					   "logicDetectorType1",
					   0,0,0);
  
  G4Trap *solidHousingOUT1 = new G4Trap("solidHousingOUT1",Dz, Theta, Phi,
					Dy1, Dx1,Dx2, Alp1, Dy2, Dx3,
					Dx4, Alp2);

  // -- Al housing inside --
  Dz = 97./2.*mm;
  Theta = 0*deg;
  Phi = 0*deg;
  Dy1 = 35.6/2.*mm;
  Dx1 = 61.3/2.*mm;
  Dx2 = Dx1;
  Alp1 = 0*deg;
  Dy2 = 53.0/2.*mm;
  Dx3 = 91.1/2.*mm;
  Dx4 = Dx3;
  Alp2 = 0*deg;

  G4Trap* solidHousingIN1 = new G4Trap("solidHousingIN1",Dz, Theta, Phi,
				       Dy1, Dx1, Dx2, Alp1, Dy2, Dx3,
				       Dx4, Alp2);
  
  G4SubtractionSolid* solidHousing1 =
    new G4SubtractionSolid("solidHousing1",solidHousingOUT1,// mother
			   solidHousingIN1, Rot,
			   G4ThreeVector(0.,0.,0.));

  G4LogicalVolume* LogicType1Housing = new G4LogicalVolume(solidHousing1,
					   Al,
					   "logicHousing1",
					   0,0,0);
  
   new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          LogicType1Housing,
          "CatanaHousingType1",m_DetectorType1,false,0);

  // -- Crystal -- 
  Dz = 93./2.*mm;
  Theta = 0*deg;
  Phi = 0*deg;
  Dy1 = 31.6/2.*mm;
  Dx1 = 57.3/2.*mm;
  Dx2 = Dx1;
  Alp1 = 0*deg;
  Dy2 = 49.0/2.*mm;
  Dx3 = 87.1/2.*mm;
  Dx4 = Dx3;
  Alp2 = 0*deg;

  G4Trap* solidCrystal1 = new G4Trap("solidCrystal1",Dz, Theta, Phi,
				     Dy1, Dx1, Dx2, Alp1, Dy2, Dx3,
				     Dx4, Alp2);

  G4LogicalVolume* logicCrystal1 = new G4LogicalVolume(solidCrystal1,// solid
					   CsI, // Material
					   "logicCrystal1", // name
					   0,0,0);

    m_DetectorType1->SetVisAttributes(m_VisCrystal);
    m_DetectorType1->SetSensitiveDetector(m_CatanaScorer);

  new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          logicCrystal1,
          "CatanaCrystalType1",m_DetectorType1,false,0);


  // -- Teflon reflector: thickness is 0.25mm
  Dz += 0.25*mm;
  Dy1 += 0.25*mm;
  Dx1 += 0.25*mm;
  Dx2 = Dx1;
  Dy2 += 0.25*mm;
  Dx3 += 0.25*mm;
  Dx4 = Dx3;

  G4Trap* solidReflectorOUT1 = new G4Trap("solidReflectorOUT1",Dz, Theta,
					  Phi, Dy1, Dx1, Dx2, Alp1,
					  Dy2, Dx3, Dx4, Alp2);

  G4RotationMatrix rotmat(0,0,0);
  G4SubtractionSolid* solidReflector1 =
    new G4SubtractionSolid("solidReflector1",solidReflectorOUT1,// mother
			   solidCrystal1, &rotmat,
			   G4ThreeVector(0.,0.,0.));

  G4LogicalVolume* logicReflectorType1 = new G4LogicalVolume(solidReflector1,
					     Teflon,
					     "logicReflector1",
					     0,0,0);
  
  new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          logicReflectorType1,
          "CatanaReflectorType1",m_DetectorType1,false,0);


  }
  return m_DetectorType1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Catana::BuildDetectorType2(){

  if(!m_DetectorType2){
   
  G4Material* CsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
  G4Material* Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4Material* Vacuum= MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4Material* Teflon= MaterialManager::getInstance()->GetMaterialFromLibrary("G4_TEFLON");
  G4RotationMatrix* Rot= new G4RotationMatrix();
  // -- Al housing outside --
  //
  double Dz = 107/2.*mm;
  double Dy1 = 34.9/2.*mm;
  double Dy2 = 55.4/2.*mm;//measured by ruler
  double Dx1 = 57.1/2.*mm;
  double Dx2 = 63.6/2.*mm;
  double Dx3 = 84.5/2.*mm;//measured by ruler
  double Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition
  double Theta = 0.*deg;
  double Phi = 0.*deg;
  double Alp1 = 0*deg;
  double Alp2 = 0*deg;


  G4Trap *solidDetectorType2 = new G4Trap("CatanaDetectorType2",Dz, Theta, Phi,
					Dy1, Dx1,Dx2, Alp1, Dy2, Dx3,
					Dx4, Alp2);

  m_DetectorType2  = new G4LogicalVolume(solidDetectorType2,
					   Vacuum,
					   "logicDetectorType2",
					   0,0,0);
  
  G4Trap *solidHousingOUT2 = new G4Trap("solidHousingOUT2",Dz, Theta, Phi,
					Dy1, Dx1,Dx2, Alp1, Dy2, Dx3,
					Dx4, Alp2);

  // -- Al housing inside --
  Dz  -= 0.5*mm;
  Dy1 -= 0.5*mm;
  Dy2 -= 0.5*mm;
  Dx1 -= 0.5*mm;
  Dx2 -= 0.5*mm;
  Dx3 -= 0.5*mm;

  Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition
  G4Trap* solidHousingIN2 = new G4Trap("solidHousingIN2",Dz, Theta, Phi,
				       Dy1, Dx1, Dx2, Alp1, Dy2, Dx3,
				       Dx4, Alp2);
  
  G4SubtractionSolid* solidHousing2 =
    new G4SubtractionSolid("solidHousing2",solidHousingOUT2,// mother
			   solidHousingIN2, Rot,
			   G4ThreeVector(0.,0.,0.));

  G4LogicalVolume* LogicType2Housing = new G4LogicalVolume(solidHousing2,
					   Al,
					   "logicHousing2",
					   0,0,0);
  
   new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          LogicType2Housing,
          "CatanaHousingType2",m_DetectorType2,false,0);

  // -- Crystal -- 
  Dz = 103./2.*mm;
  Dy1 -= 2.*mm;
  Dy2 -= 2.*mm;
  Dx1 -= 2.*mm;
  Dx2 -= 2.*mm;
  Dx3 -= 2.*mm;
  Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition

  G4Trap* solidCrystal2= new G4Trap("solidCrystal2",Dz, Theta, Phi,
				     Dy1, Dx1, Dx2, Alp1, Dy2, Dx3,
				     Dx4, Alp2);

  G4LogicalVolume* logicCrystal2 = new G4LogicalVolume(solidCrystal2,// solid
					   CsI, // Material
					   "logicCrystal1", // name
					   0,0,0);

    m_DetectorType2->SetVisAttributes(m_VisCrystal);
    m_DetectorType2->SetSensitiveDetector(m_CatanaScorer);

  new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          logicCrystal2,
          "CatanaCrystalType2",m_DetectorType2,false,0);


  // -- Teflon reflector: thickness is 0.25mm
  Dz += 0.25*mm;
  Dy1 += 0.25*mm;
  Dy2 += 0.25*mm;
  Dx1 += 0.25*mm;
  Dx2 += 0.25*mm;
  Dx3 += 0.25*mm;
  Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition

  G4Trap* solidReflectorOUT2 = new G4Trap("solidReflectorOUT2",Dz, Theta,
					  Phi, Dy1, Dx1, Dx2, Alp1,
					  Dy2, Dx3, Dx4, Alp2);

  G4RotationMatrix rotmat(0,0,0);
  G4SubtractionSolid* solidReflector2 =
    new G4SubtractionSolid("solidReflector2",solidReflectorOUT2,// mother
			   solidCrystal2, &rotmat,
			   G4ThreeVector(0.,0.,0.));

  G4LogicalVolume* logicReflectorType2 = new G4LogicalVolume(solidReflector2,
					     Teflon,
					     "logicReflector2",
					     0,0,0);
  
  new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          logicReflectorType2,
          "CatanaReflectorType2",m_DetectorType2,false,0);


  }
  return m_DetectorType2;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Catana::BuildDetectorType3(){

  if(!m_DetectorType3){
   
  G4Material* CsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
  G4Material* Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4Material* Vacuum= MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4Material* Teflon= MaterialManager::getInstance()->GetMaterialFromLibrary("G4_TEFLON");
  G4RotationMatrix* Rot= new G4RotationMatrix();
  // -- Al housing outside --
  //
  double Dz = 127/2.*mm;
  double Dy1 = 38.3/2.*mm;
  double Dy2 = 64.7/2.*mm;
  double Dx1 = 49.7/2.*mm;
  double Dx2 = 58.5/2.*mm;
  double Dx3 = 74.9/2.*mm;
  double Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition
  double Theta = 0.*deg;
  double Phi = 0.*deg;
  double Alp1 = 0*deg;
  double Alp2 = 0*deg;


  G4Trap *solidDetectorType3 = new G4Trap("CatanaDetectorType3",Dz, Theta, Phi,
					Dy1, Dx1,Dx2, Alp1, Dy2, Dx3,
					Dx4, Alp2);

  m_DetectorType3  = new G4LogicalVolume(solidDetectorType3,
					   Vacuum,
					   "logicDetectorType3",
					   0,0,0);
  
  G4Trap *solidHousingOUT3 = new G4Trap("solidHousingOUT3",Dz, Theta, Phi,
					Dy1, Dx1,Dx2, Alp1, Dy2, Dx3,
					Dx4, Alp2);

  // -- Al housing inside --
  Dz = 127./2.*mm;
  Dy1 -= 0.5*mm;
  Dy2 -= 0.5*mm;
  Dx1 -= 0.5*mm;
  Dx2 -= 0.5*mm;
  Dx3 -= 0.5*mm;
  Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition

  G4Trap* solidHousingIN3 = new G4Trap("solidHousingIN3",Dz, Theta, Phi,
				       Dy1, Dx1, Dx2, Alp1, Dy2, Dx3,
				       Dx4, Alp2);
  
  G4SubtractionSolid* solidHousing3 =
    new G4SubtractionSolid("solidHousing3",solidHousingOUT3,// mother
			   solidHousingIN3, Rot,
			   G4ThreeVector(0.,0.,0.));

  G4LogicalVolume* LogicType3Housing = new G4LogicalVolume(solidHousing3,
					   Al,
					   "logicHousing3",
					   0,0,0);
  
   new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          LogicType3Housing,
          "CatanaHousingType3",m_DetectorType3,false,0);

  // -- Crystal -- 
  Dz = 123./2.*mm;
  Dy1 -= 2.*mm;
  Dy2 -= 2.*mm;
  Dx1 -= 2.*mm;
  Dx2 -= 2.*mm;
  Dx3 -= 2.*mm;
  Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition

  G4Trap* solidCrystal3= new G4Trap("solidCrystal3",Dz, Theta, Phi,
				     Dy1, Dx1, Dx2, Alp1, Dy2, Dx3,
				     Dx4, Alp2);

  G4LogicalVolume* logicCrystal3 = new G4LogicalVolume(solidCrystal3,// solid
					   CsI, // Material
					   "logicCrystal3", // name
					   0,0,0);

    m_DetectorType3->SetVisAttributes(m_VisCrystal);
    m_DetectorType3->SetSensitiveDetector(m_CatanaScorer);

  new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          logicCrystal3,
          "CatanaCrystalType3",m_DetectorType3,false,0);


  // -- Teflon reflector: thickness is 0.25mm
  Dz += 0.25*mm;
  Dy1 += 0.25*mm;
  Dy2 += 0.25*mm;
  Dx1 += 0.25*mm;
  Dx2 += 0.25*mm;
  Dx3 += 0.25*mm;
  Dx4 = Dx3 + (Dy2/Dy1)*(Dx2-Dx1);//planarity condition

  G4Trap* solidReflectorOUT3 = new G4Trap("solidReflectorOUT3",Dz, Theta,
					  Phi, Dy1, Dx1, Dx2, Alp1,
					  Dy2, Dx3, Dx4, Alp2);

  G4RotationMatrix rotmat(0,0,0);
  G4SubtractionSolid* solidReflector3 =
    new G4SubtractionSolid("solidReflector3",solidReflectorOUT3,// mother
			   solidCrystal3, &rotmat,
			   G4ThreeVector(0.,0.,0.));

  G4LogicalVolume* logicReflectorType3 = new G4LogicalVolume(solidReflector3,
					     Teflon,
					     "logicReflector3",
					     0,0,0);
  
  new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
          logicReflectorType3,
          "CatanaReflectorType3",m_DetectorType3,false,0);


  }
  return m_DetectorType3;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Catana::ReadConfiguration(NPL::InputParser parser){
  // Dummy
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

  // Type 1
  blocks = parser.GetAllBlocksWithTokenAndValue("Catana","Type1");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  token = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetectorType1(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
 // Type 2
  blocks = parser.GetAllBlocksWithTokenAndValue("Catana","Type2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetectorType2(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
   // Type 2
  blocks = parser.GetAllBlocksWithTokenAndValue("Catana","Type3");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetectorType3(R,Theta,Phi);

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
  for (unsigned short i = 0 ; i < m_R1.size() ; i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,0,m_R1[i]) ;
    Det_pos.setTheta(m_Theta1[i]);
    Det_pos.setPhi(m_Phi1[i]);
    G4RotationMatrix* Rot = new G4RotationMatrix();
    Rot->rotateY(m_Theta1[i]);
    Rot->rotateZ(m_Phi1[i]);
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildDetectorType1(),
          "CatanaType1",world,false,i+1);
  }
  
  for (unsigned short i = 0 ; i < m_R2.size() ; i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,0,m_R2[i]) ;
    Det_pos.setTheta(m_Theta2[i]);
    Det_pos.setPhi(m_Phi2[i]);
    G4RotationMatrix* Rot = new G4RotationMatrix();
    Rot->rotateY(m_Theta2[i]);
    Rot->rotateZ(m_Phi2[i]);
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildDetectorType2(),
          "CatanaType2",world,false,i+2);
  }
  
  for (unsigned short i = 0 ; i < m_R3.size() ; i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,0,m_R3[i]) ;
    Det_pos.setTheta(m_Theta3[i]);
    Det_pos.setPhi(m_Phi3[i]);
    G4RotationMatrix* Rot = new G4RotationMatrix();
    Rot->rotateY(m_Theta3[i]);
    Rot->rotateZ(m_Phi3[i]);
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildDetectorType3(),
          "CatanaType3",world,false,i+3);
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
