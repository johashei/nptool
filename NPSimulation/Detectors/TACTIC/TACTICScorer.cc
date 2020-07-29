#include "TACTICScorer.hh"
#include "G4UnitsTable.hh"
using namespace TACTICScorer;

Gas_Scorer::Gas_Scorer(G4String name,G4int Level,G4double ScorerLength,G4int NumberOfSegments, G4int depth) //what do level and depth do?       
:G4VPrimitiveScorer(name, depth),HCID(-1){
  m_ScorerLength = ScorerLength;
  m_NumberOfSegments = NumberOfSegments;
  m_SegmentLength = m_ScorerLength / m_NumberOfSegments;
  m_Level = Level;

  m_Position = G4ThreeVector(-1000,-1000,-1000);
  m_SegmentNumber = -1;
  m_Index = -1;
}

Gas_Scorer::~Gas_Scorer(){}

G4bool Gas_Scorer::ProcessHits(G4Step* aStep, G4TouchableHistory*){

  G4double* Infos = new G4double[12];
  m_Position  = aStep->GetPreStepPoint()->GetPosition();
  
  Infos[0] = aStep->GetTrack()->GetTrackID();

  Infos[1] = -1.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "Ne18") Infos[1] = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "Na21") Infos[1] = 1.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "proton") Infos[1] = 2.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "alpha") Infos[1] = 3.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "Li8") Infos[1] = 4.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "B11") Infos[1] = 5;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron") Infos[1] = 6;

  Infos[2] = aStep->GetPreStepPoint()->GetGlobalTime();
  Infos[3] = aStep->GetPreStepPoint()->GetKineticEnergy();
  Infos[4] = aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit();

  m_SegmentNumber = (int)((m_Position.z() + m_ScorerLength / 2.) / m_SegmentLength ) + 1; //Pad number 
  Infos[5] = m_SegmentNumber;
  Infos[6] = m_Position.z();
  Infos[7] = pow(pow(m_Position.x(),2) + pow(m_Position.y(),2),0.5); //R
  Infos[8] = aStep->GetTrack()->GetVertexPosition()[2];
  Infos[9] = aStep->GetTrack()->GetVertexKineticEnergy();
  G4ThreeVector p_vec = aStep->GetTrack()->GetVertexMomentumDirection();
  Infos[10] = acos(p_vec[2]/pow(pow(p_vec[0],2)+pow(p_vec[1],2)+pow(p_vec[2],2),0.5))/deg; //angle relative to z axis (theta);   
  Infos[11] = aStep->GetTrack()->GetTrackLength();

  m_DetectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  m_Index = m_DetectorNumber * 1e3 + m_SegmentNumber * 1e6;

  map<G4int, G4double**>::iterator it;
  it= EvtMap->GetMap()->find(m_Index);
  if(it!=EvtMap->GetMap()->end()){ //accumulate ionisation energy deposit
    G4double* dummy = *(it->second);
    Infos[4]+=dummy[4];
    delete dummy;
  }
  
  EvtMap->set(m_Index, Infos);
  
  return TRUE;

}

void Gas_Scorer::Initialize(G4HCofThisEvent* HCE){
    EvtMap = new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
    if (HCID < 0) {
        HCID = GetCollectionID(0);
    }
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void Gas_Scorer::EndOfEvent(G4HCofThisEvent*){}

void Gas_Scorer::clear(){
    std::map<G4int, G4double**>::iterator    MapIterator;
    for (MapIterator = EvtMap->GetMap()->begin() ; MapIterator != EvtMap->GetMap()->end() ; MapIterator++){
        delete *(MapIterator->second);
    }

    EvtMap->clear();
}


void Gas_Scorer::DrawAll(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                              

void Gas_Scorer::PrintAll(){
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
    G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}

