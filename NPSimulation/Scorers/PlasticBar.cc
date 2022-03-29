/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : February 2013                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  File old the scorer specific to the Sharc Detector                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new type of scorer is aim to become the standard for DSSD,SSSD and   *
 * PAD detector (any Silicon Detector)                                       *
 *****************************************************************************/
#include "PlasticBar.hh"
#include "G4UnitsTable.hh"
using namespace PlasticBar;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
unsigned int
PlasticBarData::CalculateIndex(const vector<unsigned int>& level) {

  unsigned int size       = level.size();
  unsigned int result     = 0;
  unsigned int multiplier = 1;
  for (unsigned int i = 0; i < size; i++) {
    result += level[i] * multiplier;
    multiplier *= 1000;
  }
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<PlasticBarData>::iterator
PlasticBarDataVector::find(const unsigned int& index) {
  for (vector<PlasticBarData>::iterator it = m_Data.begin();
       it != m_Data.end(); it++) {
    if ((*it).GetIndex() == index)
      return it;
  }
  return m_Data.end();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_PlasticBar::PS_PlasticBar(G4String name, vector<G4int> NestingLevel,
                               G4int depth)
    : G4VPrimitiveScorer(name, depth) {
  m_NestingLevel = NestingLevel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_PlasticBar::~PS_PlasticBar() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_PlasticBar::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  G4String particlename
      = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    cout << particlename << " has left ";
    // Contain Energy, Time + as many copy number as nested volume
    unsigned int mysize = m_NestingLevel.size();
    t_Energy            = aStep->GetTotalEnergyDeposit();
    cout << t_Energy << " MeV";
    t_Time              = aStep->GetPreStepPoint()->GetGlobalTime();
    cout << " after " << t_Time << " ms" << endl;
    t_Level.clear();
    for (unsigned int i = 0; i < mysize; i++) {
      t_Level.push_back(
          aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(
              m_NestingLevel[i]));
    }
    // Check if the particle has interact before, if yes, add up the energies.
    vector<PlasticBarData>::iterator it;
    it = m_Data.find(PlasticBarData::CalculateIndex(t_Level));
    if (it != m_Data.end()) {
      it->Add(t_Energy);
    } else {
      m_Data.Set(t_Energy, t_Time, t_Level);
    }
    return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::Initialize(G4HCofThisEvent*) { clear(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::EndOfEvent(G4HCofThisEvent*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::clear() {
  m_Data.clear();
  t_Level.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::PrintAll() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
