/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: morfouac@nscl.msu.edu *
 *                                                                           *
 * Creation Date  : June 2016                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  NeutronWall simulation                              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
#include <cstring>
#include <string>

//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
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
#include "NeutronWall.hh"
#include "CalorimeterScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace NeutronWall_NS{
    // Energy and time Resolution
    const double EnergyThreshold = 0.1*MeV;
    const double ResoTime = 4.5*ns ;
    const double ResoEnergy = 5.0*MeV ;
    //The size of NS should depend on the distance between NeutronWall and plastic Bar right now
    double NS_X = 2020.0*mm;
    double NS_Y = 2020.0*mm;
    //the front and back aluminum sheet are both 0.8 thick whereas 143.5 is user assumed heigh in z
    double NS_Z = (143.5+0.8+0.8)*mm;
    //using Alouter minus Alinner, one get an Al frame (including front and back sheets)
    const double Alinner_X = 2000.0*mm;
    const double Alinner_Y = 2000.0*mm;
    const double Alinner_Z = 143.5*mm;
    const double Alouter_X = 2020.0*mm;
    const double Alouter_Y = 2020.0*mm;
    const double Alouter_Z = (143.5+0.8+0.8)*mm;
    
    const double frame_thickness = 10*mm;
    
    const double Scintillator_X = 1994.0*mm;
    const double Scintillator_Y = 70.2*mm;
    const double Scintillator_Z = 57.5*mm;
    //do the same thing to Pyrex tube to create a volume to store Scintillator
    
    const double Py_Xinner = 1994.0*mm;
    const double Py_Yinner = 70.2*mm;
    const double Py_Zinner = 57.5*mm;
    const double Py_Xouter = 2000.0*mm;
    const double Py_Youter = 76.2*mm;
    const double Py_Zouter = 63.5*mm;
    
    const double seperation_between_pyrex = 3.0*mm;
    const double upper_gap = 10.0*mm;
    
    //Add elements about the plastic bars
    double PlasticBar_X = 0.0*mm;
    double PlasticBar_Y = 0.0*mm;
    double PlasticBar_Z = 10.0*mm;
    
    //Add total height of neutronwall and vetowall for comparision
    double TotalHeightOfNeutronWall = 0.0*mm;
    double TotalHeightOfVetoWall = 0.0*mm;
    
    //Scale down factor
    double ScaleDownFactor = 0;
    
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// NeutronWall Specific Method
NeutronWall::NeutronWall(){
    m_Event = new TNeutronWallData() ;
    m_NeutronWallScorer = 0;
    m_VetoWallScorer = 0;
    m_NeutronWall_log = 0;
    m_AlCase_log = 0;
    m_Quartz_log = 0;
    m_QuartzCap_log = 0;
    m_PMTube_log = 0;
    m_Scintillator_log = 0;
    m_ShadowBar_log = 0;
    m_PlasticBar_log = 0;
    
    
    // RGB Color + Transparency
    m_VisScintillator = new G4VisAttributes(G4Colour(1, 0.843137, 0, 0.3)); //gold
    m_VisQuartz = new G4VisAttributes(G4Colour(0,1, 0, 0.1)); //green
    m_VisAl = new G4VisAttributes(G4Colour(173.0/255.0,178.0/255.0,189.0/255.0,0.3)); //Al
    m_VisNW = new G4VisAttributes(G4Colour(0.972549,0.972549,1,0.1)); //ghostwhite
    m_VisPlasticBar = new G4VisAttributes(G4Colour(0.9,0,0.9,1.0)); //pink
}

NeutronWall::~NeutronWall(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NeutronWall::AddNeutronWall(double  R, double  Theta, double  Phi, double X, double Y, double Z, double Rotation, int Bars, string NWMaterial, double VWDistance, int VetoWall, string VWMaterial, double Overlap){
    m_R.push_back(R);
    m_Theta.push_back(Theta);
    m_Phi.push_back(Phi);
    m_X.push_back(X);
    m_Y.push_back(Y);
    m_Z.push_back(Z);
    m_Rot.push_back(Rotation);
    m_Bars.push_back(Bars);
    m_NWMaterial.push_back(NWMaterial);
    m_VWDistance.push_back(VWDistance);
    m_VetoWall.push_back(VetoWall);
    m_VWMaterial.push_back(VWMaterial);
    m_Overlap.push_back(Overlap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NeutronWall::BuildDetector(){
    
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void NeutronWall::ReadConfiguration(string Path){
    ifstream ConfigFile           ;
    ConfigFile.open(Path.c_str()) ;
    string LineBuffer          ;
    string DataBuffer          ;
    
    double Theta = 0 , Phi = 0 , R = 0 ;
    double X = 0 , Y = 0 , Z = 0 ;
    double Rot =0;
    int Bars = 0;
    string NWMaterial = "NE213";
    double VWDistance = 0.0;
    int VetoWall = 0;
    string VWMaterial = "BC400";
    double Overlap = 3;
    
    bool check_Theta = false ;
    bool check_Phi = false ;
    bool check_R = false ;
    bool check_rotation = false ;
    bool check_X = false ;
    bool check_Y = false ;
    bool check_Z = false ;
    bool ReadingStatus = false ;
    bool check_Bars = false ;
    bool check_NWMaterial = false ;
    
    while (!ConfigFile.eof()) {
        getline(ConfigFile, LineBuffer);
        
        //   If line is a Start Up NeutronWall bloc, Reading toggle to true
        string name = "NeutronWall";
        
        if (LineBuffer.compare(0, name.length(), name) == 0) {
            G4cout << "///" << G4endl           ;
            G4cout << "NeutronWall found: " << G4endl   ;
            ReadingStatus = true ;
        }
        
        //   Else don't toggle to Reading Block Status
        else ReadingStatus = false ;
        
        //   Reading Block
        while(ReadingStatus){
            // Pickup Next Word
            ConfigFile >> DataBuffer ;
            
            //   Comment Line
            if (DataBuffer.compare(0, 1, "%") == 0) {
                ConfigFile.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' );
            }
            
            //   Finding another telescope (safety), toggle out
            else if (DataBuffer.compare(0, name.length(),name) == 0) {
                G4cout << "WARNING: Another Detector is find before standard sequence of Token, Error may occured in Telecope definition" << G4endl ;
                ReadingStatus = false ;
            }
            
            //Angle method
            else if (DataBuffer.compare(0, 6, "THETA=") == 0) {
                check_Theta = true;
                ConfigFile >> DataBuffer ;
                Theta = atof(DataBuffer.c_str()) ;
                Theta = Theta * deg;
                G4cout << "Theta:  " << Theta / deg << G4endl;
            }
            
            else if (DataBuffer.compare(0, 4, "PHI=") == 0) {
                check_Phi = true;
                ConfigFile >> DataBuffer ;
                Phi = atof(DataBuffer.c_str()) ;
                Phi = Phi * deg;
                G4cout << "Phi:  " << Phi / deg << G4endl;
            }
            
            else if (DataBuffer.compare(0, 2, "R=") == 0) {
                check_R = true;
                ConfigFile >> DataBuffer ;
                R = atof(DataBuffer.c_str()) ;
                R = R * mm;
                G4cout << "R:  " << R/mm << G4endl;
            }
            
            //Position method
            else if (DataBuffer.compare(0, 2, "X=") == 0) {
                check_X = true;
                ConfigFile >> DataBuffer ;
                X = atof(DataBuffer.c_str()) ;
                X = X * cm;
                G4cout << "X:  " << X / cm << G4endl;
            }
            
            else if (DataBuffer.compare(0, 2, "Y=") == 0) {
                check_Y = true;
                ConfigFile >> DataBuffer ;
                Y = atof(DataBuffer.c_str()) ;
                Y = Y * cm;
                G4cout << "Y:  " << Y / cm << G4endl;
            }
            
            else if (DataBuffer.compare(0, 2, "Z=") == 0) {
                check_Z = true;
                ConfigFile >> DataBuffer ;
                Z = atof(DataBuffer.c_str()) ;
                Z = Z * cm;
                G4cout << "Z:  " << Z / cm << G4endl;
            }
            
            
            //General
            else if (DataBuffer.compare(0, 4, "Rot=") == 0) {
                check_rotation = true;
                ConfigFile >> DataBuffer ;
                Rot = atof(DataBuffer.c_str());
                Rot = Rot*deg ;
                G4cout << "Rotation:  " << Rot/deg << G4endl;
            }
            
            //Bar number
            else if (DataBuffer.compare(0, 5, "BARS=") == 0){
                check_Bars = true;
                ConfigFile >> DataBuffer ;
                Bars = atoi(DataBuffer.c_str()) ;
                G4cout << "Bars:  " << Bars << G4endl;
            }
            
            
            //Material type
            else if (DataBuffer.compare(0, 11, "NWMATERIAL=") == 0){
                check_NWMaterial = true;
                ConfigFile >> DataBuffer ;
                NWMaterial = DataBuffer;
                G4cout << "NWMaterials:  " << NWMaterial << G4endl;
            }
            
            //Distance
            else if (DataBuffer.compare(0, 11, "VWDISTANCE=") == 0){
                //check_VWDistance = true;
                ConfigFile >> DataBuffer ;
                VWDistance = atof(DataBuffer.c_str());
                VWDistance = VWDistance * mm;
                G4cout << "VWDistance: " << VWDistance << G4endl;
            }
            
            //Decide whether to add the vetowall or not, 1 means yes, 0 means no
            else if (DataBuffer.compare(0, 9, "VETOWALL=") == 0){
                //check_VetoWall = true;
                ConfigFile >> DataBuffer ;
                VetoWall = atoi(DataBuffer.c_str());
                G4cout << "VetoWall:  " << VetoWall << G4endl;
            }
            
            //VetoWall Material
            else if (DataBuffer.compare(0, 11, "VWMATERIAL=") == 0){
                //check_VWMaterial = true;
                ConfigFile >> DataBuffer ;
                VWMaterial = DataBuffer ;
                G4cout << "VWMaterial:  " << VWMaterial << G4endl;
            }
            
            //Overlap
            else if (DataBuffer.compare(0, 8, "OVERLAP=") == 0){
                ConfigFile >> DataBuffer ;
                Overlap = atof(DataBuffer.c_str());
                Overlap = Overlap*mm;
                G4cout << "Overlap:  " << Overlap << G4endl;
            }
            
            ///////////////////////////////////////////////////
            //   If no Detector Token and no comment, toggle out
            else{
                ReadingStatus = false;
                G4cout << "Wrong Token Sequence: Getting out " << DataBuffer << G4endl ;
            }
            
            /////////////////////////////////////////////////
            //   If All necessary information there, toggle out
            
            if (( check_Theta && check_Phi && check_R && check_Bars && check_NWMaterial)
                ||
                ( check_X && check_Y && check_Z && check_rotation && check_Bars && check_NWMaterial)){
                
                
                // Convert Cartesian to Spherical (detector always face the target)
                if (check_X){
                    R = sqrt (X*X+Y*Y+Z*Z);
                    Theta = acos(Z / (R) );
                    Phi = atan2(Y,X);
                }
                
                AddNeutronWall(R,Theta,Phi,X,Y,Z,Rot, Bars, NWMaterial, VWDistance, VetoWall, VWMaterial, Overlap);
                
                //   Reinitialisation of Check Boolean
                check_Theta = false ;
                check_Phi = false ;
                check_R = false ;
                check_rotation = false ;
                check_X = false ;
                check_Y = false ;
                check_Z = false ;
                ReadingStatus = false ;
                check_Bars = false ;
                check_NWMaterial = false ;
                G4cout << "///"<< G4endl ;
            }
        }
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void NeutronWall::ConstructDetector(G4LogicalVolume* world){
    
    for (unsigned short i = 0 ; i < m_R.size() ; i++) {
        
        G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
        G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
        G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
        G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
        // So the face of the Scintillator Bar is at R instead of the middle
        Det_pos+=Det_pos.unit()*NeutronWall_NS::Scintillator_Z*0.5;
        // Building Detector reference frame
        G4double ii = cos(m_Theta[i]) * cos(m_Phi[i]);
        G4double jj = cos(m_Theta[i]) * sin(m_Phi[i]);
        G4double kk = -sin(m_Theta[i]);
        G4ThreeVector Y(ii,jj,kk);
        G4ThreeVector w = Det_pos.unit();
        G4ThreeVector u = w.cross(Y);
        G4ThreeVector v = u.cross(w);
        v = v.unit();
        u = u.unit();
        
        
        
        //Initialize scale down factor , measured from face to face but VWDistance is from center to center
        NeutronWall_NS::ScaleDownFactor = (m_R[i]-(m_VWDistance[i]-NeutronWall_NS::Scintillator_Z*0.5+0.5*NeutronWall_NS::PlasticBar_Z))/(m_R[i]);
	
	G4cout << "//////" << "ScaleDownFactor: " << NeutronWall_NS::ScaleDownFactor << "//////" << endl;
						//one layer case 
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       
        //if (m_VetoWall[i] == 1){
            //Initialize property about the plastic bar, if exists
           // NeutronWall_NS::PlasticBar_X = NeutronWall_NS::Py_Xinner*NeutronWall_NS::ScaleDownFactor;
            //NeutronWall_NS::PlasticBar_Y = NeutronWall_NS::Py_Yinner*NeutronWall_NS::ScaleDownFactor+4*mm;
	    //m_Overlap[i] = -(NeutronWall_NS::seperation_between_pyrex+6*mm)*NeutronWall_NS::ScaleDownFactor;
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
						//Overlap case
	double ExtrudeY = 36*mm;
	double ExtrudeX = 8*mm;
	if (m_VetoWall[i] == 1){
            //Initialize property about the plastic bar, if exists
            NeutronWall_NS::PlasticBar_X = NeutronWall_NS::Py_Xinner+ExtrudeX;
            NeutronWall_NS::PlasticBar_Y = NeutronWall_NS::Py_Yinner*NeutronWall_NS::ScaleDownFactor+ExtrudeY;
	    m_Overlap[i] = -(NeutronWall_NS::seperation_between_pyrex+6*mm)*NeutronWall_NS::ScaleDownFactor;
	
					
            
            // If VetoWall exists, then extend the space of NeutronWall_NS add 3mm to Z in order to house the 1mm seperation between front and back layer of the veto wall
            NeutronWall_NS::NS_Z = 2.0*(m_VWDistance[i]+1.5*NeutronWall_NS::PlasticBar_Z+3.0*mm);
        }
        
        G4RotationMatrix* Rot = new G4RotationMatrix(v,u,w);
        
        G4Material* ScintMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_NWMaterial[i]);
        G4Material* vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
        G4Material* Aluminum = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
        G4Material* Pyrex = MaterialManager::getInstance()->GetMaterialFromLibrary("Pyrex");
        
        //Neutron Wall Box
        G4Box* NeutronWall_box = new G4Box("NeutronWall_Box",NeutronWall_NS::NS_X*0.5,
                                           NeutronWall_NS::NS_Y*0.5,NeutronWall_NS::NS_Z*0.5);
        m_NeutronWall_log = new G4LogicalVolume(NeutronWall_box,vacuum,"NeutronWall_Log",0,0,0);
        m_NeutronWall_log->SetVisAttributes(m_VisNW);
        
        //Aluminum inner box (subtractee)
        G4Box* Alinner_box = new G4Box("Alinner_box",NeutronWall_NS::Alinner_X*0.5,NeutronWall_NS::Alinner_Y*0.5,NeutronWall_NS::Alinner_Z*0.5);
        //Aluminum outer box (subtracter)
        G4Box* Alouter_box = new G4Box("Alouter_box",NeutronWall_NS::Alouter_X*0.5,NeutronWall_NS::Alouter_Y*0.5,NeutronWall_NS::Alouter_Z*0.5);
        
        G4SubtractionSolid* AlCase_frame = new G4SubtractionSolid("AlCase_frame", Alouter_box, Alinner_box);
        
        m_AlCase_log = new G4LogicalVolume(AlCase_frame,Aluminum,"AlCase_Log",0,0,0);
        
        m_AlCase_log->SetVisAttributes(m_VisAl);
        
        //Quartz tube
        G4Box* Quartz_boxinner = new G4Box("Quartz_Boxinner",NeutronWall_NS::Py_Xinner*0.5,NeutronWall_NS::Py_Yinner*0.5,NeutronWall_NS::Py_Zinner*0.5);
        G4Box* Quartz_boxouter = new G4Box("Quartz_Box",NeutronWall_NS::Py_Xouter*0.5,NeutronWall_NS::Py_Youter*0.5,NeutronWall_NS::Py_Zouter*0.5);
        G4SubtractionSolid* Quartz_box = new G4SubtractionSolid("Quartz_box", Quartz_boxouter, Quartz_boxinner);
        
        m_Quartz_log = new G4LogicalVolume(Quartz_box,Pyrex,"Quartz_Log",0,0,0);
        
        m_Quartz_log->SetVisAttributes(m_VisQuartz);
        
        //??currently unused
        G4Tubs* QuartzCap_tube = new G4Tubs("QuartsCap_Tube",0,7.99*cm,0.3175*cm,0.0*deg,360.0*deg);
        m_QuartzCap_log = new G4LogicalVolume(QuartzCap_tube,Aluminum,"QuartzCap_Log",0,0,0);
        
        //Scintillator
        G4Box* Scintillator_box = new G4Box("Scintillator_Box",NeutronWall_NS::Scintillator_X*0.5,NeutronWall_NS::Scintillator_Y*0.5,NeutronWall_NS::Scintillator_Z*0.5);
        m_Scintillator_log = new G4LogicalVolume(Scintillator_box,ScintMaterial,"Scintillator_Log",0,0,0);
        
        m_Scintillator_log->SetVisAttributes(m_VisScintillator);
        m_Scintillator_log->SetSensitiveDetector(m_NeutronWallScorer);
        
        //Shadow bar construction??currently unused
        G4Trd* ShadowBar_trd = new G4Trd("ShadowBar_Trd",6.51/2*cm, 7.23/2*cm, 6.79/2*cm, 7.66/2*cm, 30./2*cm);
        m_ShadowBar_log = new G4LogicalVolume(ShadowBar_trd, Aluminum, "ShadowBar_Log");
        
        
        if (m_VetoWall[i] == 1){
/*
            //Initialize total height
            NeutronWall_NS::TotalHeightOfNeutronWall = m_Bars[i]*NeutronWall_NS::Py_Youter + (m_Bars[i]-1)*NeutronWall_NS::seperation_between_pyrex;
            NeutronWall_NS::TotalHeightOfVetoWall = (m_Bars[i]+1)/2*NeutronWall_NS::PlasticBar_Y+(m_Bars[i]-1)/2*(NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
            if (NeutronWall_NS::TotalHeightOfNeutronWall*NeutronWall_NS::ScaleDownFactor > NeutronWall_NS::TotalHeightOfVetoWall){
                G4cout << "\t*************************************************************" <<endl;
                G4cout << "\t* The shadow of VetoWall is not enough to cover NeutronWall *" <<endl;
                G4cout << "\t*            Re-input VWDistance or Overlap                 *" <<endl;
                G4cout << "\t*************************************************************" <<endl;
                exit(1);
            }*/


            //PlasticBar
            G4Material* Plastic = MaterialManager::getInstance()->GetMaterialFromLibrary(m_VWMaterial[i]);
            G4Box* PlasticBar_box = new G4Box("PlasticBar_Box", NeutronWall_NS::PlasticBar_X*0.5, NeutronWall_NS::PlasticBar_Y*0.5, NeutronWall_NS::PlasticBar_Z*0.5);
            m_PlasticBar_log = new G4LogicalVolume(PlasticBar_box, Plastic, "PlasticBar_Log");
            m_PlasticBar_log->SetSensitiveDetector(m_VetoWallScorer);
            m_PlasticBar_log->SetVisAttributes(m_VisPlasticBar);
        }
        
        
        //******************* Placement *******************//
        //----Aluminum Case----
        m_AlCase_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),m_AlCase_log,
                                          "AlCase_phys",m_NeutronWall_log,false,0);
        //----Scintillator and Quartz tube----
        for (int j = 0; j < 25; j++ ) {
	    double CenterOfScintillator_X = 0*mm;
	    double CenterOfScintillator_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5-j*(NeutronWall_NS::Py_Youter+NeutronWall_NS::seperation_between_pyrex));
   	    double CenterOfScintillator_Z = 0*mm;
	    G4ThreeVector ScintillatorDisplacement(CenterOfScintillator_X,CenterOfScintillator_Y,CenterOfScintillator_Z);
            m_ScintillatorTube_phys = new G4PVPlacement(0,ScintillatorDisplacement,m_Scintillator_log,"ScintillatorTube_phys",m_NeutronWall_log,false,j);
	    //Quartz center coincide with Scintillator's, therefore, they have the same displacement.
            m_Quartz_phys = new G4PVPlacement(0,ScintillatorDisplacement,m_Quartz_log, "Quartz_phys",m_NeutronWall_log,false,j);
	    if (m_VetoWall[i] == 1){
		
		//Even number is associated with 0th 2nd 4th ... plasticbar in vetowall which comprise of the backlayer.
		//Odd number is associated with 1st 3rd 5th ... plasticbar in vetowall which comprise of the frontlayer.
		double CenterOfVetoWall_Even_X = 0*mm;
		//double CenterOfVetoWall_Even_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-0.5*NeutronWall_NS::PlasticBar_Y-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		double CenterOfVetoWall_Even_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor - j*(NeutronWall_NS::PlasticBar_Y-ExtrudeY-m_Overlap[i]);
		//double CenterOfVetoWall_Even_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap)*NeutronWall_NS::ScaleDownFactor-0.5*NeutronWall_NS::PlasticBar_Y-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);

		double CenterOfVetoWall_Even_Z = -m_VWDistance[i];

		double CenterOfVetoWall_Odd_X = 0*mm;
		//double CenterOfVetoWall_Odd_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-1.5*NeutronWall_NS::PlasticBar_Y+m_Overlap[i]-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		double CenterOfVetoWall_Odd_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor - j*(NeutronWall_NS::PlasticBar_Y-ExtrudeY-m_Overlap[i]);
		//double CenterOfVetoWall_Odd_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap)*NeutronWall_NS::ScaleDownFactor-1.5*NeutronWall_NS::PlasticBar_Y+m_Overlap[i]-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		//double CenterOfVetoWall_Odd_Z = -m_VWDistance[i]-1*mm-NeutronWall_NS::PlasticBar_Z;
		double CenterOfVetoWall_Odd_Z = -m_VWDistance[i]-NeutronWall_NS::PlasticBar_Z-1*mm;


		G4ThreeVector VetoTransOfBackLayer(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z);
		G4ThreeVector VetoTransOfFrontLayer(CenterOfVetoWall_Odd_X,CenterOfVetoWall_Odd_Y,CenterOfVetoWall_Odd_Z);
                if (j%2 == 0){
                    m_PlasticBar_phys = new G4PVPlacement(0,G4ThreeVector(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z),m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,j,true);
                }
                else {
                    m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z-NeutronWall_NS::PlasticBar_Z-1*mm),m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,j,true);
		}
                    
            }
	    












		/*if (j <= 2){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor-13*mm - j*(NeutronWall_NS::PlasticBar_Y-4*mm-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}

		if ((j <= 5) && (j > 2)){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor-8*mm - j*(NeutronWall_NS::PlasticBar_Y-4*mm-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}		
		
		if ((j <= 7)&&(j > 5)){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor-3*mm - j*(NeutronWall_NS::PlasticBar_Y-4*mm-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}*/

		/*if ((j <= 9)&&(j > 7)){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor-3*mm - j*(NeutronWall_NS::PlasticBar_Y-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}*/

		/*if ((j < 17)&&(j >7)){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor - j*(NeutronWall_NS::PlasticBar_Y-4*mm-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}*/

		/*if ((j <19) && (j >= 15)){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor+3*mm - j*(NeutronWall_NS::PlasticBar_Y-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}*/


		/*if ((j < 19) && (j >=17 )){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor+3*mm - j*(NeutronWall_NS::PlasticBar_Y-4*mm-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}


		if ((j < 22) && (j >= 19)){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor+8*mm - j*(NeutronWall_NS::PlasticBar_Y-4*mm-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}
		
		if ((j <=24) && (j >= 22)){
		double MoveX = 0*mm;
		double MoveY = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor+13*mm - j*(NeutronWall_NS::PlasticBar_Y-4*mm-m_Overlap[i]);
		double MoveZ = -m_VWDistance[i];
		m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(MoveX,MoveY,MoveZ), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);		
		}

		







	    }
        }*/

	/*for (int j = 0; j < 25; j++){
	    double CenterOfVetoWall_X = 0*mm;
	    double CenterOfVetoWall_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap-NeutronWall_NS::Py_Youter*0.5)*NeutronWall_NS::ScaleDownFactor - j*(NeutronWall_NS::PlasticBar_Y-m_Overlap[i]);
	    double CenterOfVetoWall_Z = -m_VWDistance[i];
	    m_PlasticBar_phys = new G4PVPlacement(0, G4ThreeVector(CenterOfVetoWall_X,CenterOfVetoWall_Y,CenterOfVetoWall_Z), m_PlasticBar_log, "PlasticBar_phys", m_NeutronWall_log, false, j);
	}*/
        
        /*for(int j = 0; j<25; j++){
            if (m_VetoWall[i] == 1){
		//Even number is associated with 0th 2nd 4th ... plasticbar in vetowall which comprise of the backlayer.
		//Odd number is associated with 1st 3rd 5th ... plasticbar in vetowall which comprise of the frontlayer.
		double CenterOfVetoWall_Even_X = 0*mm;
		//double CenterOfVetoWall_Even_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-0.5*NeutronWall_NS::PlasticBar_Y-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		//double CenterOfVetoWall_Even_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap)*(1.0/NeutronWall_NS::ScaleDownFactor)*0.9-0.5*NeutronWall_NS::PlasticBar_Y-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		double CenterOfVetoWall_Even_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap)*NeutronWall_NS::ScaleDownFactor-0.5*NeutronWall_NS::PlasticBar_Y-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);

		double CenterOfVetoWall_Even_Z = -m_VWDistance[i];

		double CenterOfVetoWall_Odd_X = 0*mm;
		//double CenterOfVetoWall_Odd_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-1.5*NeutronWall_NS::PlasticBar_Y+m_Overlap[i]-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		double CenterOfVetoWall_Odd_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap)*(1.0/NeutronWall_NS::ScaleDownFactor)*0.9-1.5*NeutronWall_NS::PlasticBar_Y+m_Overlap[i]-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		//double CenterOfVetoWall_Odd_Y = (NeutronWall_NS::NS_Y*0.5-NeutronWall_NS::frame_thickness-NeutronWall_NS::upper_gap)*NeutronWall_NS::ScaleDownFactor-1.5*NeutronWall_NS::PlasticBar_Y+m_Overlap[i]-(j/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
		//double CenterOfVetoWall_Odd_Z = -m_VWDistance[i]-1*mm-NeutronWall_NS::PlasticBar_Z;
		double CenterOfVetoWall_Odd_Z = -m_VWDistance[i];


		G4ThreeVector VetoTransOfBackLayer(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z);
		G4ThreeVector VetoTransOfFrontLayer(CenterOfVetoWall_Odd_X,CenterOfVetoWall_Odd_Y,CenterOfVetoWall_Odd_Z);
                if (j%2 == 0){
                    m_PlasticBar_phys = new G4PVPlacement(0,VetoTransOfBackLayer,m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,j);
                }
                else {
                    m_PlasticBar_phys = new G4PVPlacement(0, VetoTransOfFrontLayer,m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,j);
                    
                }
            }*/
	
        }
	/*double CenterOfVetoWall_Even_X = 0*mm;
	double CenterOfVetoWall_Even_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-0.5*NeutronWall_NS::PlasticBar_Y-(0/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
	double CenterOfVetoWall_Even_Z = -m_VWDistance[i];
	G4ThreeVector VetoTransOfBackLayer1(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z);
	m_PlasticBar_phys = new G4PVPlacement(0,VetoTransOfBackLayer1,m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,0);

	CenterOfVetoWall_Even_X = 0*mm;
	CenterOfVetoWall_Even_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-0.5*NeutronWall_NS::PlasticBar_Y-(6/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
	CenterOfVetoWall_Even_Z = -m_VWDistance[i];
	G4ThreeVector VetoTransOfBackLayer4(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z);
	m_PlasticBar_phys = new G4PVPlacement(0,VetoTransOfBackLayer4,m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,6);

	CenterOfVetoWall_Even_X = 0*mm;
	CenterOfVetoWall_Even_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-0.5*NeutronWall_NS::PlasticBar_Y-(12/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
	CenterOfVetoWall_Even_Z = -m_VWDistance[i];
	G4ThreeVector VetoTransOfBackLayer3(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z);
	m_PlasticBar_phys = new G4PVPlacement(0,VetoTransOfBackLayer3,m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,12);

	CenterOfVetoWall_Even_X = 0*mm;
	CenterOfVetoWall_Even_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-0.5*NeutronWall_NS::PlasticBar_Y-(18/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
	CenterOfVetoWall_Even_Z = -m_VWDistance[i];
	G4ThreeVector VetoTransOfBackLayer5(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z);
	m_PlasticBar_phys = new G4PVPlacement(0,VetoTransOfBackLayer5,m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,18);
	
	CenterOfVetoWall_Even_X = 0*mm;
	CenterOfVetoWall_Even_Y = 0.5*NeutronWall_NS::TotalHeightOfVetoWall-0.5*NeutronWall_NS::PlasticBar_Y-(24/2)*(NeutronWall_NS::PlasticBar_Y+NeutronWall_NS::PlasticBar_Y-2*m_Overlap[i]);
	CenterOfVetoWall_Even_Z = -m_VWDistance[i];
	G4ThreeVector VetoTransOfBackLayer2(CenterOfVetoWall_Even_X,CenterOfVetoWall_Even_Y,CenterOfVetoWall_Even_Z);
	m_PlasticBar_phys = new G4PVPlacement(0,VetoTransOfBackLayer2,m_PlasticBar_log,"PlasticBar_phys",m_NeutronWall_log,false,24);*/
        /****************** Place the walls*************************/
        m_NeutronWall_phys = new G4PVPlacement(G4Transform3D(*Rot, Det_pos),
                                               m_NeutronWall_log,
                                               "NeutronWall_phys",world,false,i);
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void NeutronWall::InitializeRootOutput(){
    RootOutput *pAnalysis = RootOutput::getInstance();
    TTree *pTree = pAnalysis->GetTree();
    pTree->Branch("NeutronWall", "TNeutronWallData", &m_Event) ;
    pTree->SetBranchAddress("NeutronWall", &m_Event) ;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void NeutronWall::ReadSensitive(const G4Event* event){
    m_Event->Clear();
    
    ///////////
    // Calorimeter scorer
    G4THitsMap<G4double*>* CaloHitMap;
    std::map<G4int, G4double**>::iterator Calo_itr;
    
    G4int CaloCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("NeutronWallScorer/Calorimeter");
    CaloHitMap = (G4THitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(CaloCollectionID));
    
    // Loop on the Calo map
    for (Calo_itr = CaloHitMap->GetMap()->begin() ; Calo_itr != CaloHitMap->GetMap()->end() ; Calo_itr++){
        
        G4double* Info = *(Calo_itr->second);
        //(Info[0]/2.35)*((Info[0]*1.02)*pow((Info[0]*1.8),.5))
        // double Energy = RandGauss::shoot(Info[0],((Info[0]*1000*1.02/2.35)*pow((Info[0]*1000*1.8),.5)) );
        double Energy = RandGauss::shoot(Info[0],NeutronWall_NS::ResoEnergy);
        if(Energy>NeutronWall_NS::EnergyThreshold){
            double Time = RandGauss::shoot(Info[1],NeutronWall_NS::ResoTime);
            int DetectorNbr = (int) Info[3];
            int PadNbr = (int) Info[2];
            //cout << Info[2] << " " << Info[3] << endl;
            m_Event->SetEnergy(DetectorNbr,PadNbr,Energy);
            m_Event->SetTime(DetectorNbr,PadNbr,Time);
        }
    }
    // clear map for next event
    CaloHitMap->clear();
    
    ///////////
    // Veto wall scorer
    G4THitsMap<G4double*>* VetoHitMap;
    std::map<G4int, G4double**>::iterator Veto_itr;
    
    G4int VetoCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("VetoWallScorer/VetoCalorimeter");
    VetoHitMap = (G4THitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(VetoCollectionID));
    
    
    // Loop on the Calo map
    for (Veto_itr = VetoHitMap->GetMap()->begin() ; Veto_itr != VetoHitMap->GetMap()->end() ; Veto_itr++){
        
        G4double* Info = *(Veto_itr->second);
        //(Info[0]/2.35)*((Info[0]*1.02)*pow((Info[0]*1.8),.5))
        // double Energy = RandGauss::shoot(Info[0],((Info[0]*1000*1.02/2.35)*pow((Info[0]*1000*1.8),.5)) );
        double Energy = RandGauss::shoot(Info[0],NeutronWall_NS::ResoEnergy);
        if(Energy>NeutronWall_NS::EnergyThreshold){
            double Time = RandGauss::shoot(Info[1],NeutronWall_NS::ResoTime);
            int DetectorNbr = (int) Info[8];
            int PadNbr = (int) Info[7];
            
            m_Event->SetVetoEnergy(DetectorNbr,PadNbr,Energy);
            m_Event->SetVetoTime(DetectorNbr,PadNbr,Time);

	     // Interraction Coordinates
        ms_InterCoord->SetDetectedPositionX(Info[2]) ;
        ms_InterCoord->SetDetectedPositionY(Info[3]) ;
        ms_InterCoord->SetDetectedPositionZ(Info[4]) ;
        ms_InterCoord->SetDetectedAngleTheta(Info[5]/deg) ;
        ms_InterCoord->SetDetectedAnglePhi(Info[6]/deg) ;
        }
    }
    // clear map for next event
    VetoHitMap->clear();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void NeutronWall::InitializeScorers() {
    // Otherwise the scorer is initialised
    vector<int> level;
    level.push_back(0);
    level.push_back(1);
    
    // This check is necessary in case the geometry is reloaded
    bool already_exist = false;
    m_NeutronWallScorer = CheckScorer("NeutronWallScorer",already_exist) ;
    
    if(already_exist)
        return ;
    
    // Neutron Wall Scorer
    G4VPrimitiveScorer* Calorimeter= new CALORIMETERSCORERS::PS_Calorimeter("Calorimeter",level,1) ;
    //and register it to the multifunctional detector
    m_NeutronWallScorer->RegisterPrimitive(Calorimeter);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_NeutronWallScorer) ;
    
    
    // Veto Wall Scorer
    already_exist = false;
    m_VetoWallScorer = CheckScorer("VetoWallScorer",already_exist) ;
    if(already_exist)
        return;
    
    G4VPrimitiveScorer* VetoCalorimeter= new CALORIMETERSCORERS::PS_CalorimeterWithInteraction("VetoCalorimeter",level,1) ;
    //and register it to the multifunctional detector
    m_VetoWallScorer->RegisterPrimitive(VetoCalorimeter);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_VetoWallScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* NeutronWall::Construct(){
    return  (NPS::VDetector*) new NeutronWall();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
    class proxy_nps_plastic{
    public:
        proxy_nps_plastic(){
            NPS::DetectorFactory::getInstance()->AddToken("NeutronWall","NeutronWall");
            NPS::DetectorFactory::getInstance()->AddDetector("NeutronWall",NeutronWall::Construct);
        }
    };
    
    proxy_nps_plastic p_nps_plastic;
}