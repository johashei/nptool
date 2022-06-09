/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Freddy Flavigny  contact: flavigny@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : April 2022                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Custom scorer for Germanium detectors                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new type of scorer is for Ge detectors built using the geometry      *
 * conventions used in the AGATA/GRETINA/HiCARI simulation codes             *
 * file (euler, clust, solid, slice)                                         *
 *****************************************************************************/
#include "GeScorers.hh"
#include "G4UnitsTable.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"

using namespace GeScorers ;
vector<GeData>::iterator GeDataVector::find(const unsigned int& index){
  for(vector<GeData>::iterator it= m_Data.begin()  ; it !=m_Data.end() ; it++){
    if((*it).GetIndex()==index)
      return it;
  }
  return m_Data.end();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//PS_GeDetector::PS_GeDetector(G4String name, G4int level,G4int depth){
	//m_NestingLevel = NestingLevel;
//}
PS_GeDetector::PS_GeDetector(G4String name, G4int level,G4int depth)
	:G4VPrimitiveScorer(name, depth){
	m_Level = level;
	m_Depth = depth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_GeDetector::~PS_GeDetector(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_GeDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*){

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();

    G4int detCode = theTouchable->GetReplicaNumber(0);
    if(detCode == 0)
        detCode = theTouchable->GetReplicaNumber(m_Depth);

    G4int detNum = detCode%1000;

    G4String name = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

    //G4cout << name << ", " << detCode << ", " << detNum << G4endl;

    // Track everything (the primary gamma, all secondary electrons, 
    // and all secondary gammas).
    G4double edep = aStep->GetTotalEnergyDeposit();

    // Keep the gamma parent of pair-production tracks for hit processing.
    //const G4VProcess* test = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    //const G4String& testname = test->GetProcessName();
    //if(edep<0.001*eV && testname != "conv") 
    if(edep<0.001*eV && 
       aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() != "conv") 
        return false;

    G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();

    G4VPhysicalVolume* topVolume;
    topVolume = theTouchable->GetVolume(m_Depth);


    G4ThreeVector frameTrans = topVolume->GetFrameTranslation();

    const G4RotationMatrix* rot = topVolume->GetFrameRotation();

    G4RotationMatrix frameRot;
    if( rot )
        frameRot = *rot;

    G4ThreeVector posSol = frameRot( position );
    posSol += frameRot( frameTrans );


  G4String particlename
      = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    // cout << particlename << endl;
    // Contain Energy, Time + as many copy number as nested volume
    //unsigned int mysize = m_NestingLevel.size();
    t_Time  = aStep->GetPreStepPoint()->GetGlobalTime();
/*
    t_Level.clear();
    for (unsigned int i = 0; i < mysize; i++) {
      t_Level.push_back(
          aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(
              m_NestingLevel[i]));
    }

    // Check if the particle has interact before, if yes, add up the energies.
    
    vector<GeData>::iterator it;
    it = m_Hit.find(detCode);
    if (it != m_Hit.end()) {
      it->Add(edep);
    } else {
      m_Hit.Set(edep, t_Time, detCode, posSol);
    }
    */
    m_Hit.Set(edep, t_Time, detCode, posSol);

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GeDetector::Initialize(G4HCofThisEvent* HCE){
/*
	EvtMap = new NPS::HitsMap<StepInfo>(GetMultiFunctionalDetector()->GetName(), GetName());
	if (HCID < 0) {
		HCID = GetCollectionID(0);
	}
	HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
*/

    clear(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GeDetector::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GeDetector::clear(){   
	m_Hit.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GeDetector::DrawAll(){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GeDetector::PrintAll(){
/*
	G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
	G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
	G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
*/
}


/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Images::PS_Images(G4String name, string imageFront,string imageBack,double scalingFront, double scalingBack, double centerOffsetX,double centerOffsetY,unsigned int ignoreValue, G4int depth)  :G4VPrimitiveScorer(name, depth){
  m_ImageFront = new NPL::Image(imageFront,scalingFront,scalingFront);
  m_ImageBack  = new NPL::Image(imageBack,scalingBack,scalingBack);
  m_ScalingFront = scalingFront;
  m_ScalingBack  = scalingBack;
  m_CenterOffsetX = centerOffsetX;
  m_CenterOffsetY = centerOffsetY;
  m_IgnoreValue = ignoreValue;
  m_Level = depth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Images::ProcessHits(G4Step* aStep, G4TouchableHistory*){

  // contain Energy Time, DetNbr, PixelFront and PixelBack
  t_Energy = aStep->GetTotalEnergyDeposit();
  t_Time = aStep->GetPreStepPoint()->GetGlobalTime();

  t_DetectorNbr = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  t_Position  = aStep->GetPreStepPoint()->GetPosition();

  // transforming the position to the local volume
  t_Position = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(t_Position);
  t_PixelFront = m_ImageFront->GetPixelAtCoordinate(t_Position.x(),t_Position.y());
  t_PixelBack = m_ImageBack->GetPixelAtCoordinate(t_Position.x(),t_Position.y());

  // If front and back are in inactive part of the wafer,
  // nothing is added to the unordered_map
  if(t_PixelFront == m_IgnoreValue && t_PixelBack == m_IgnoreValue)
    return FALSE;


  // Check if the particle has interact before, if yes, add up the energies.
   vector<GeData>::iterator it;
  
  it= m_HitFront.find(GeData::CalculateIndex(t_PixelFront,t_DetectorNbr));
  if(it!=m_HitFront.end()){
    it->Add(t_Energy);
  }

  else{
    m_HitFront.Set(t_Energy,t_Time,t_PixelFront,t_DetectorNbr);
  }

  // Check if the particle has interact before, if yes, add up the energies.
  it= m_HitBack.find(GeData::CalculateIndex(t_PixelBack,t_DetectorNbr));
  if(it!=m_HitBack.end()){
    it->Add(t_Energy);
  }

  else{
    m_HitBack.Set(t_Energy,t_Time,t_PixelBack,t_DetectorNbr);
  }

  return TRUE;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Images::Initialize(G4HCofThisEvent*){
  clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Images::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Images::clear(){
  m_HitFront.clear();
  m_HitBack.clear();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Images::GetARGBFront(unsigned int& i,unsigned int& a,unsigned int& r,unsigned int& g,unsigned int& b){
  unsigned int Info = m_HitFront[i]->GetStrip();
  a = (Info>>24)&0xff;
  r = (Info>>16)&0xff;
  g = (Info>>8)&0xff;
  b = (Info>>0)&0xff;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Images::GetARGBBack(unsigned int& i,unsigned int& a,unsigned int& r,unsigned int& g,unsigned int& b){
  unsigned int Info = m_HitBack[i]->GetStrip();
  a = (Info>>24)&0xff;
  r = (Info>>16)&0xff;
  g = (Info>>8)&0xff;
  b = (Info>>0)&0xff;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Rectangle::PS_Rectangle(G4String name,G4int Level, G4double StripPlaneLength, G4double StripPlaneWidth, G4int NumberOfStripLength,G4int NumberOfStripWidth,G4int depth,G4String axis)
  :G4VPrimitiveScorer(name, depth){
    m_StripPlaneLength = StripPlaneLength;
    m_StripPlaneWidth = StripPlaneWidth;
    m_NumberOfStripLength = NumberOfStripLength;
    m_NumberOfStripWidth = NumberOfStripWidth;
    m_StripPitchLength = m_StripPlaneLength / m_NumberOfStripLength;
    m_StripPitchWidth = m_StripPlaneWidth / m_NumberOfStripWidth;
    m_Level = Level;
    if(axis=="xy")
      m_Axis=ps_xy;
    else if(axis=="yz")
      m_Axis=ps_yz;
    else if(axis=="xz")
      m_Axis=ps_xz;


  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Rectangle::~PS_Rectangle(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Rectangle::ProcessHits(G4Step* aStep, G4TouchableHistory*){

  // contain Energy Time, DetNbr, StripFront and StripBack
  t_Energy  = aStep->GetTotalEnergyDeposit();
  t_Time  = aStep->GetPreStepPoint()->GetGlobalTime();

  t_DetectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  t_Position  = aStep->GetPreStepPoint()->GetPosition();

  t_Position = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(t_Position);

  if(m_Axis==ps_xy){
    t_StripLengthNumber = (int)((t_Position.x() + m_StripPlaneLength / 2.) / m_StripPitchLength ) + 1  ;
    t_StripWidthNumber = (int)((t_Position.y() + m_StripPlaneWidth / 2.) / m_StripPitchWidth ) + 1  ;
  }
  else if(m_Axis==ps_yz){
    t_StripLengthNumber = (int)((t_Position.y() + m_StripPlaneLength / 2.) / m_StripPitchLength ) + 1  ;
    t_StripWidthNumber = (int)((t_Position.z() + m_StripPlaneWidth / 2.) / m_StripPitchWidth ) + 1  ;
  }
  else if(m_Axis==ps_xz){
    t_StripLengthNumber = (int)((t_Position.x() + m_StripPlaneLength / 2.) / m_StripPitchLength ) + 1  ;
    t_StripWidthNumber = (int)((t_Position.z() + m_StripPlaneWidth / 2.) / m_StripPitchWidth ) + 1  ;
  }

  //Rare case where particle is close to edge of silicon plan
  if (t_StripLengthNumber > m_NumberOfStripLength) t_StripLengthNumber = m_NumberOfStripLength;
  if (t_StripWidthNumber > m_NumberOfStripWidth) t_StripWidthNumber = m_NumberOfStripWidth;

  // Check if the particle has interact before, if yes, add up the energies.
  vector<GeData>::iterator it;
  // Length
  it = m_HitLength.find(GeData::CalculateIndex(t_StripLengthNumber,t_DetectorNumber));
  if(it!=m_HitLength.end()){
    it->Add(t_Energy);
  }
  else
    m_HitLength.Set(t_Energy,t_Time,t_StripLengthNumber,t_DetectorNumber);
  // Width
  it = m_HitWidth.find(GeData::CalculateIndex(t_StripWidthNumber,t_DetectorNumber));
  if(it!=m_HitWidth.end()){
    it->Add(t_Energy);
  }
  else
    m_HitWidth.Set(t_Energy,t_Time,t_StripWidthNumber,t_DetectorNumber);


return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Rectangle::Initialize(G4HCofThisEvent*){
  clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Rectangle::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Rectangle::clear(){
  m_HitLength.clear();
  m_HitWidth.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Rectangle::DrawAll(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Rectangle::PrintAll(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Annular::PS_Annular(G4String name,G4int Level, G4double StripPlaneInnerRadius, G4double StripPlaneOuterRadius, G4double PhiStart,G4double PhiStop, G4int NumberOfStripRing,G4int NumberOfStripSector,G4int NumberOfQuadrant,G4int depth)
  :G4VPrimitiveScorer(name, depth){

    m_StripPlaneInnerRadius = StripPlaneInnerRadius;
    m_StripPlaneOuterRadius = StripPlaneOuterRadius;
    m_PhiStart = PhiStart;
    m_PhiStop = PhiStop;
    m_NumberOfStripRing = NumberOfStripRing;
    m_NumberOfStripSector = NumberOfStripSector;
    m_NumberOfStripQuadrant  = NumberOfQuadrant;
    m_StripPitchRing =  (m_StripPlaneOuterRadius-m_StripPlaneInnerRadius)/m_NumberOfStripRing;
    m_StripPitchSector = (m_PhiStop-m_PhiStart)/m_NumberOfStripSector;
    m_StripPitchQuadrant = (m_PhiStop-m_PhiStart)/m_NumberOfStripQuadrant;
    m_Level = Level;

    m_uz = G4ThreeVector(0,0,1);
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Annular::~PS_Annular(){
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Annular::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  t_Energy = aStep->GetTotalEnergyDeposit();

  t_Time = aStep->GetPreStepPoint()->GetGlobalTime();

  t_DetectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  t_Position = aStep->GetPreStepPoint()->GetPosition();

  //Transform into local coordinates
  t_Position = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(t_Position);
  t_StripRingNumber = (int) ((t_Position.rho() - m_StripPlaneInnerRadius) / m_StripPitchRing ) + 1 ;

  // phi() from G4-CLHEP method return azimuth between [-pi;pi]
  // we need it in [0;2pi] to calculate sector nbr in [1,NSectors]
  // only add 360 if the value is negative
  double phi = (t_Position.phi()<0)?  t_Position.phi()+2*pi : t_Position.phi() ;

  // factor out the extra 360 degrees before strip/quad calculation
  t_StripSectorNumber   = (int) ( fmod((phi - m_PhiStart),2*pi)  / m_StripPitchSector  ) + 1 ;
  t_StripQuadrantNumber = (int) ( fmod((phi - m_PhiStart),2*pi)  / m_StripPitchQuadrant) + 1 ;

  //Rare case where particle is close to edge of silicon plan
  if (t_StripRingNumber > m_NumberOfStripRing) t_StripRingNumber = m_NumberOfStripRing;
  if (t_StripSectorNumber > m_NumberOfStripSector) t_StripSectorNumber = m_NumberOfStripSector;
  if (t_StripQuadrantNumber > m_NumberOfStripQuadrant) t_StripQuadrantNumber = m_NumberOfStripQuadrant;

  vector<GeData>::iterator it;
  // Ring
  it = m_HitRing.find(GeData::CalculateIndex(t_StripRingNumber,t_DetectorNumber));
  if(it!=m_HitRing.end()){
    it->Add(t_Energy);
  }
  else
    m_HitRing.Set(t_Energy,t_Time,t_StripRingNumber,t_DetectorNumber);

  //Sector
  it = m_HitSector.find(GeData::CalculateIndex(t_StripSectorNumber,t_DetectorNumber));
  if(it!=m_HitSector.end()){
    it->Add(t_Energy);
  }
  else
    m_HitSector.Set(t_Energy,t_Time,t_StripSectorNumber,t_DetectorNumber);

  //Quadrant
  it = m_HitQuadrant.find(GeData::CalculateIndex(t_StripQuadrantNumber,t_DetectorNumber));
  if(it!=m_HitQuadrant.end()){
    it->Add(t_Energy);
  }
  else
    m_HitQuadrant.Set(t_Energy,t_Time,t_StripQuadrantNumber,t_DetectorNumber);

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Annular::Initialize(G4HCofThisEvent*){
  clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Annular::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Annular::clear(){
  m_HitRing.clear();
  m_HitSector.clear();
  m_HitQuadrant.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Annular::DrawAll(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Annular::PrintAll(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Resistive::PS_Resistive(G4String name,G4int Level, G4double StripPlaneLength, G4double StripPlaneWidth, G4int NumberOfStripWidth,G4int depth)
  :G4VPrimitiveScorer(name, depth){
    m_StripPlaneLength = StripPlaneLength;
    m_StripPlaneWidth = StripPlaneWidth;
    m_NumberOfStripWidth = NumberOfStripWidth;
    m_StripPitchWidth = m_StripPlaneWidth / m_NumberOfStripWidth;
    m_Level = Level;

    t_Position = G4ThreeVector(-1000,-1000,-1000);
    t_DetectorNumber = -1;
    t_StripWidthNumber = -1;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Resistive::~PS_Resistive(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Resistive::ProcessHits(G4Step* aStep, G4TouchableHistory*){

  // contain Energy Total, E1, E2, Time, DetNbr,  and StripWidth
  t_Energy = aStep->GetTotalEnergyDeposit(); 
  t_Time = aStep->GetPreStepPoint()->GetGlobalTime();
  
  t_DetectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  t_Position = aStep->GetPreStepPoint()->GetPosition();
  t_Position = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(t_Position);

  t_StripWidthNumber = (int)((t_Position.x() + m_StripPlaneWidth / 2.) / m_StripPitchWidth ) + 1  ;

  // The energy is divided in two depending on the position
  // position along the resistive strip
  double P = (t_Position.z())/(0.5*m_StripPlaneLength);
  // Energy
  t_EnergyUp = aStep->GetTotalEnergyDeposit()*(1+P)*0.5; 
  t_EnergyDown = t_Energy-t_EnergyUp; 

  //Rare case where particle is close to edge of silicon plan
  if (t_StripWidthNumber > m_NumberOfStripWidth) t_StripWidthNumber = m_NumberOfStripWidth;


  // Up
  vector<GeData>::iterator it;
  it = m_HitUp.find(GeData::CalculateIndex(t_DetectorNumber,t_StripWidthNumber));
  if(it!=m_HitUp.end())
    it->Add(t_EnergyUp);
  else
    m_HitUp.Set(t_EnergyUp,t_Time,t_StripWidthNumber,t_DetectorNumber);
    
  // Down
  it = m_HitDown.find(GeData::CalculateIndex(t_DetectorNumber,t_StripWidthNumber));
  if(it!=m_HitDown.end())
    it->Add(t_EnergyDown);
  else
    m_HitDown.Set(t_EnergyDown,t_Time,t_StripWidthNumber,t_DetectorNumber);
  
   // Back
  it = m_HitBack.find(GeData::CalculateIndex(t_DetectorNumber,t_StripWidthNumber));
  if(it!=m_HitBack.end())
    it->Add(t_Energy);
  else
    m_HitBack.Set(t_Energy,t_Time,1,t_DetectorNumber);
  
  
  
  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Resistive::Initialize(G4HCofThisEvent* ){
  clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Resistive::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Resistive::clear(){
  m_HitUp.clear();
  m_HitDown.clear();
  m_HitBack.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Resistive::DrawAll(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Resistive::PrintAll(){
}
*/
