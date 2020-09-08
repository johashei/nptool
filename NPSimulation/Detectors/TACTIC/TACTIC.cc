/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Warren Lynch  contact address: Warren.Lynch@york.ac.uk                        *
 *                                                                           *
 * Creation Date  : June 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  TACTIC simulation                             *
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
#include "G4THitsMap.hh"
#include "G4SDParticleFilter.hh"

// NPTool header
#include "TACTIC.hh"
#include "TACTICScorer.hh"

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
namespace TACTIC_NS{
  // Energy and time Resolution

  //const double EnergyThreshold = 0.01*MeV;
  //const double ResoTime = 17.*ns ;
  //const double ResoEnergy = 1.0*MeV ;
  const double cathode_radius = 12.*mm;
  const double drift_radius = 50.*mm;
  const double anode_radius = 51.*mm;
  const double active_length = 251.9*mm;
  const double window_pos = 104.*mm; //from centre of TACTIC from https://elog.triumf.ca/Tactic/Documentation/18                                     
  const double window_radius = 12.*mm; //guess
  const double window_width = 1.5e-03*mm;
  const double vacuum_pos = (active_length/2. + window_pos + window_width/2.)/2.*mm;
  const double vacuum_width = active_length/2. - window_pos - window_width/2.*mm;
  const G4int  NumberOfStrips = 60;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// TACTIC Specific Method
TACTIC::TACTIC(){
  m_Event = new TTACTICData() ;
  m_Scorer = 0;
  m_CylindricalDetector = 0;

  m_ReactionRegion = NULL; 
  // RGB Color + Transparency

  m_VisChamber        = new G4VisAttributes(G4Colour(1., 1., 1., 1.0));
  m_VisWindows        = new G4VisAttributes(G4Colour(1, 1, 1, 1.0));
  m_VisGas            = new G4VisAttributes(G4Colour(1, 1, 1, 0.1));
  m_VisVacuum         = new G4VisAttributes(G4Colour(1,1,1,1.));

  m_VisChamber->SetForceWireframe(true);
  m_VisVacuum->SetForceWireframe(true);
}

TACTIC::~TACTIC(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TACTIC::AddDetector(G4ThreeVector POS, string  Shape){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TACTIC::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* TACTIC::BuildCylindricalDetector(){
  if(!m_CylindricalDetector){

    G4Material* Cu = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
    G4Material* Mylar = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");
    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    
    unsigned const int NumberOfGasMix = m_GasMaterial.size();

    double density=0;
    vector<G4Material*> GasComponent;
    
    for(unsigned int i=0; i<NumberOfGasMix; i++){
      if(m_GasMaterial[i] == "CO2") GasComponent.push_back(MaterialManager::getInstance()->GetGasFromLibrary(m_GasMaterial[i], 1.0/bar, m_Temperature));
      if(m_GasMaterial[i] == "P10_gas") {
	//G4Material *P10_gas = new G4Material("P10_gas", 0.00156 * g /cm3, 3, kStateGas, m_Temperature, 1.0/bar); //density for SRIM (1 atm);
	G4Material *P10_gas = new G4Material("P10_gas", 0.00156*g/cm3, 3);
	G4Element *elAr = new G4Element("Argon","Ar",18.,39.948*g/mole);
	G4Element *elC = new G4Element("Carbon","C",6.,12.0107*g/mole);
	G4Element *elH = new G4Element("Hydrogen","H",1.,1.00784*g/mole);
	P10_gas->AddElement(elAr,90);
	P10_gas->AddElement(elC,2);
	P10_gas->AddElement(elH,8);
	GasComponent.push_back(P10_gas);
      }
      else GasComponent.push_back(MaterialManager::getInstance()->GetMaterialFromLibrary(m_GasMaterial[i]));
    }
    
    for(unsigned int i=0; i<NumberOfGasMix; i++){
      density += ((double)m_GasFraction[i]/100)*GasComponent[i]->GetDensity()*(m_Pressure/bar)/(GasComponent[i]->GetPressure()/bar);
      // p2 = p1*(P2/P1) //e.g.  p2 = p1*(0.5/1) p scales with P, T=const 
    }

    G4Material* TACTIC_gas = new G4Material("TACTIC_gas", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);

    for(unsigned int i=0; i<NumberOfGasMix; i++) TACTIC_gas->AddMaterial(GasComponent[i], (double)m_GasFraction[i]/100);

    cout << TACTIC_gas << endl;

    G4Tubs* anode = new G4Tubs("anode",0,TACTIC_NS::anode_radius,TACTIC_NS::active_length*0.5,0,360*deg);
    G4Tubs* gas_volume = new G4Tubs("gas_volume",0,TACTIC_NS::drift_radius,TACTIC_NS::active_length*0.5,0,360*deg);
    G4Tubs* window = new G4Tubs("window",0,TACTIC_NS::window_radius,TACTIC_NS::window_width*0.5,0,360*deg);
    G4Tubs* vacuum = new G4Tubs("vacuum",0,TACTIC_NS::window_radius,TACTIC_NS::vacuum_width*0.5,0,360*deg);

    m_CylindricalDetector = new G4LogicalVolume(anode, Cu, "anode_log",0,0,0);
    gas_volume_log = new G4LogicalVolume(gas_volume, TACTIC_gas, "gas_volume_log",0,0,0);
    window_log = new G4LogicalVolume(window, Mylar, "window_log",0,0,0);
    //window_log = new G4LogicalVolume(window, TACTIC_gas, "window_log",0,0,0); //windows removed (effectively)
    vacuum_log = new G4LogicalVolume(vacuum, Vacuum, "vacuum_log",0,0,0);

    new G4PVPlacement(0,G4ThreeVector(0,0,0),gas_volume_log,"gas_volume_phys",m_CylindricalDetector,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,TACTIC_NS::window_pos),window_log,"window_phys",gas_volume_log,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,-TACTIC_NS::window_pos),window_log,"window_phys",gas_volume_log,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,TACTIC_NS::vacuum_pos),vacuum_log,"vacuum_phys",gas_volume_log,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,-TACTIC_NS::vacuum_pos),vacuum_log,"vacuum_phys",gas_volume_log,false,0);

    gas_volume_log->SetVisAttributes(m_VisGas);
    window_log->SetVisAttributes(m_VisWindows);
    vacuum_log->SetVisAttributes(m_VisVacuum);
    m_CylindricalDetector->SetVisAttributes(m_VisChamber);

    G4UserLimits *gas_volume_step = new G4UserLimits();
    G4double maxStep = 0.1*mm;
    gas_volume_step->SetMaxAllowedStep(maxStep);
    gas_volume_log->SetUserLimits(gas_volume_step);

    gas_volume_log->SetSensitiveDetector(m_Scorer);
  }
  return m_CylindricalDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void TACTIC::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("TACTIC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 
  
  vector<string> cart = {"POS","Shape","GasMaterial_1","GasMaterial_2","GasFraction_1","GasFraction_2","Temperature","Pressure"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TACTIC " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      string Shape = blocks[i]->GetString("Shape");
      m_GasMaterial.push_back(blocks[i]->GetString("GasMaterial_1"));
      m_GasMaterial.push_back(blocks[i]->GetString("GasMaterial_2"));
      m_GasFraction.push_back(blocks[i]->GetInt("GasFraction_1"));
      m_GasFraction.push_back(blocks[i]->GetInt("GasFraction_2"));
      m_Temperature = blocks[i]->GetDouble("Temperature","kelvin");
      m_Pressure = blocks[i]->GetDouble("Pressure","bar");
      AddDetector(Pos,Shape);
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
void TACTIC::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*TACTIC_NS::active_length*0.5;
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
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),BuildCylindricalDetector(),"TACTIC",world,false,i+1);

    if(!m_ReactionRegion){
      G4ProductionCuts* ecut = new G4ProductionCuts();
      ecut->SetProductionCut(1000,"e-"); //I think lowest is 900 eV for delta electron production, this is 1000 MeV to cut all electrons produced this way
      m_ReactionRegion= new G4Region("NPSimulationProcess");
      m_ReactionRegion->SetProductionCuts(ecut);
      m_ReactionRegion->AddRootLogicalVolume(gas_volume_log);
    }

    G4FastSimulationManager* mng = m_ReactionRegion->GetFastSimulationManager();
    unsigned int size = m_ReactionModel.size();
    for(unsigned int j = 0 ; j < size ; j++){
      mng->RemoveFastSimulationModel(m_ReactionModel[j]);
    }
    m_ReactionModel.clear();
    
    G4VFastSimulationModel* fsm;
    fsm = new NPS::BeamReaction("BeamReaction",m_ReactionRegion);
    m_ReactionModel.push_back(fsm);
        
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void TACTIC::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("TACTIC")){
    pTree->Branch("TACTIC", "TTACTICData", &m_Event) ;
  }
  pTree->SetBranchAddress("TACTIC", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void TACTIC::ReadSensitive(const G4Event* event ){
  m_Event->Clear();

  ofstream file;

  G4THitsMap<G4double*>* BeamHitMap;
  std::map<G4int, G4double**>::iterator Beam_itr;
  G4int BeamCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TACTICScorer/BeamScorer");
  BeamHitMap = (G4THitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(BeamCollectionID));

  file.open("signal.dat", std::ios::app);
  file << "Event" << endl;
  file.close();
  
  file.open("out.dat",std::ios::app);

  for (Beam_itr = BeamHitMap->GetMap()->begin(); Beam_itr != BeamHitMap->GetMap()->end(); Beam_itr++) {
    G4double* Info = *(Beam_itr->second);
    //file <<  floor(((Info[3]+TACTIC_NS::active_length*0.5)/(TACTIC_NS::active_length/TACTIC_NS::NumberOfStrips))) << "\t"; // To get PAD number
    file << event->GetEventID() << "\t";
    for(int s = 0; s<11; s++) {
      //if(s==12) file << Info[s] << endl;
      //else
      file << Info[s] << "\t";
    }
    file << Info[11] << endl;
  }

  BeamHitMap->clear();
  
  G4THitsMap<G4double*>* EjectHitMap;
  std::map<G4int, G4double**>::iterator Eject_itr;
  G4int EjectCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TACTICScorer/EjectScorer");
  EjectHitMap = (G4THitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(EjectCollectionID));

  for (Eject_itr = EjectHitMap->GetMap()->begin(); Eject_itr != EjectHitMap->GetMap()->end(); Eject_itr++) {
    G4double* Info = *(Eject_itr->second);
    //file <<  floor(((Info[3]+TACTIC_NS::active_length*0.5)/(TACTIC_NS::active_length/TACTIC_NS::NumberOfStrips))) << "\t"; // To get PAD number
    file << event->GetEventID() << "\t";
    for(int s = 0; s<11; s++) {
      //if(s==12) file << Info[s] << endl;
      //else
      file << Info[s] << "\t";
    }
    file << Info[11] << endl;
  }

  EjectHitMap->clear();

  file.close();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void TACTIC::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_Scorer = CheckScorer("TACTICScorer",already_exist) ;
  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  G4VPrimitiveScorer* EjectScorer = new TACTICScorer::Gas_Scorer("EjectScorer",1,TACTIC_NS::active_length,(int)TACTIC_NS::NumberOfStrips);
  G4VPrimitiveScorer* BeamScorer = new TACTICScorer::Gas_Scorer("BeamScorer",1,TACTIC_NS::active_length,(int)TACTIC_NS::NumberOfStrips);
  //G4SDParticleFilter* EjectFilter = new G4SDParticleFilter("EjectFilter","proton");
  G4SDParticleFilter* EjectFilter = new G4SDParticleFilter("EjectFilter","alpha"); //For studying alpha source data
  G4SDParticleFilter* BeamFilter = new G4SDParticleFilter("BeamFilter");
  BeamFilter->addIon(11,21);
  BeamFilter->addIon(10,18);
  EjectScorer->SetFilter(EjectFilter);
  BeamScorer->SetFilter(BeamFilter);

  m_Scorer->RegisterPrimitive(BeamScorer);
  m_Scorer->RegisterPrimitive(EjectScorer);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_Scorer);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* TACTIC::Construct(){
  return  (NPS::VDetector*) new TACTIC();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_TACTIC{
    public:
      proxy_nps_TACTIC(){
        NPS::DetectorFactory::getInstance()->AddToken("TACTIC","TACTIC");
        NPS::DetectorFactory::getInstance()->AddDetector("TACTIC",TACTIC::Construct);
      }
  };

  proxy_nps_TACTIC p_nps_TACTIC;
}
