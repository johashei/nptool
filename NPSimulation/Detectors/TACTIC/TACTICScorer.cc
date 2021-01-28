//#define USE_Garfield //only use if compiled with Garfield

#include "TACTICScorer.hh"
#include "G4UnitsTable.hh"

#ifdef USE_Garfield
#include "GARFDRIFT.h"
#endif

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
  //G4double w_value = 26.31*eV;
  m_Position  = aStep->GetPreStepPoint()->GetPosition();

  Infos[0] = aStep->GetTrack()->GetTrackID();

  /*
  Infos[1] = -1.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "Ne18") Infos[1] = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "Na21") Infos[1] = 1.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "proton") Infos[1] = 2.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "alpha") Infos[1] = 3.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "Li8") Infos[1] = 4.;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "B11") Infos[1] = 5;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron") Infos[1] = 6;
  */

  Infos[1] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();;
  
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

  //Infos[12] = 0.;
  //Infos[13] = 0.;
  
  //Infos[14] = 0.;

  //cout << "\n" << "PAD " << Infos[5] <<  " first_step = " << first_step << endl;
  
  //double TOA_min, TOA_max;
  //if(PAD==1000) PAD = Infos[5];
  //if(PAD!=Infos[5]) TOA_PAD.clear(); // if new pad clear the TOA vec
#ifdef USE_Garfield
  G4ThreeVector delta_Position = aStep->GetDeltaPosition();
  //  bool last_step;
  //if(aStep->GetTrack()->GetTrackStatus()!=fAlive) last_step = 1; else last_step = 0;
  //vector<double> TOA_vec = GARFDRIFT(Infos[4]/eV, Infos[2], m_Position/cm, delta_Position/cm, Infos[7]/cm, Infos[5], Infos[1], m_ScorerLength/cm, m_SegmentLength/cm, last_step); //Ionization Energy Deposit, Global Time, x, y, z,R, ParticleID //Garfield works in cm, both G4 and Garfield work in ns.

  GARFDRIFT(Infos[4]/eV, Infos[2], m_Position/cm, delta_Position/cm, Infos[7]/cm, Infos[5], Infos[1], m_ScorerLength/cm, m_SegmentLength/cm);
#endif
  /*
  for(int t=0; t<TOA_vec.size(); t++) TOA_PAD.push_back(TOA_vec[t]);
  sort(TOA_PAD.begin(), TOA_PAD.end());
  reverse(TOA_PAD.begin(), TOA_PAD.end());
  if(TOA_PAD.size() > 42) TOA_max = TOA_PAD[42], TOA_min = *min_element(TOA_PAD.begin(), TOA_PAD.end()); //threshold of 42 (420 e- ~10 mV threshold)
  else TOA_min = 1.e06, TOA_max = 0.;
  cout << "\n" << "PAD: " << Infos[5] << " vec size: " << TOA_PAD.size() << " TOA_min: " << TOA_min << " TOA_max " << TOA_max << endl; 

  
#else
  TOA_min =1.e06;
  TOA_max = 0.;
#endif
  */
  
  map<G4int, G4double**>::iterator it;
  it= EvtMap->GetMap()->find(m_Index);
  if(it!=EvtMap->GetMap()->end()){
    G4double* dummy = *(it->second);
    /*
    if(TOA_min < 1.e06) Infos[12] = TOA_max - TOA_min;
    if(Infos[12] < dummy[12]) Infos[12] = dummy[12]; //Ensures that max risetime is maintained
    PAD = Infos[5];
    */
    
    Infos[4]+=dummy[4]; //accumulate ionisation energy deposit to get total accross pad
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
