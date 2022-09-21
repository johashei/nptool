/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2012                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
using namespace std;
#include "Analysis.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPFunction.h"
#include "NPOptionManager.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis() {}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis() {}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
  InitOutputBranch();
  InitInputBranch();

  QQQ3 = (TQQQ3Physics*)m_DetectorManager->GetDetector("QQQ3");
  myReaction = new NPL::Reaction();
  myReaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();

  // energy losses
  string light = NPL::ChangeNameToG4Standard(myReaction->GetNucleus3().GetName());
  string beam = NPL::ChangeNameToG4Standard(myReaction->GetNucleus1().GetName());

  LightCD2 = NPL::EnergyLoss(light + "_" + TargetMaterial + ".G4table", "G4Table", 100);
  LightAl = NPL::EnergyLoss(light + "_Al.G4table", "G4Table", 100);
  LightSi = NPL::EnergyLoss(light + "_Si.G4table", "G4Table", 100);
  BeamCD2 = NPL::EnergyLoss(beam + "_" + TargetMaterial + ".G4table", "G4Table", 100);

  OriginalBeamEnergy = myReaction->GetBeamEnergy();
  Rand = TRandom3();
  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
  Si_E_M2 = 0;
  CsI_E_M2 = 0;
  Energy = 0;
  E_M2 = 0;

  ThetaQQQ3Surface = 0;
  X_QQQ3 = 0;
  Y_QQQ3 = 0;
  Z_QQQ3 = 0;
  Si_E_QQQ3 = 0;
  E_QQQ3 = 0;
  Si_X_QQQ3 = 0;
  Si_Y_QQQ3 = 0;

  double BeamEnergy = BeamCD2.Slow(OriginalBeamEnergy, TargetThickness * 0.5, 0);
  myReaction->SetBeamEnergy(BeamEnergy);
  cout << "Beam energy set at " << BeamEnergy << " MeV" << endl;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
  // Reinitiate calculated variable
  ReInitValue();
  double XTarget = 0;
  double YTarget = 0;
  TVector3 BeamDirection = TVector3(0, 0, 1);
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////// LOOP on QQQ3//////////////////
  if (QQQ3->Strip_E.size() == 1) {
    /************************************************/
    // Part 1 : Impact Angle
    ThetaQQQ3Surface = 0;
    ThetaNormalTarget = 0;
    if (XTarget > -1000 && YTarget > -1000) {
      TVector3 HitDirection = QQQ3->GetPositionOfInteraction(0);
      ThetaLab = HitDirection.Angle(BeamDirection);

      ThetaQQQ3Surface = HitDirection.Angle(TVector3(0, 0, 1));
      ThetaNormalTarget = HitDirection.Angle(TVector3(0, 0, 1));
    }

    else {
      BeamDirection = TVector3(-1000, -1000, -1000);
      ThetaQQQ3Surface = -1000;
      ThetaNormalTarget = -1000;
    }

    /************************************************/

    /************************************************/
    // Part 2 : Impact Energy

    Energy = 0;
    if (QQQ3->PAD_E[0] > 0) {
      Energy = QQQ3->PAD_E[0];
    }

    Energy += QQQ3->Strip_E[0];
    // Target Correction

    ELab = LightCD2.EvaluateInitialEnergy(Energy, TargetThickness * 0.5, ThetaNormalTarget);
    /************************************************/

    /************************************************/
    // Part 3 : Excitation Energy Calculation
    Ex = myReaction->ReconstructRelativistic(ELab, ThetaLab);

    /************************************************/

    /************************************************/
    // Part 4 : Theta CM Calculation
    ThetaCM = myReaction->EnergyLabToThetaCM(ELab, ThetaLab) / deg;
    ThetaLab = ThetaLab / deg;
    /************************************************/
  } // end loop QQQ3
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Ex", &Ex, "Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab", &ELab, "ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab", &ThetaLab, "ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM", &ThetaCM, "ThetaCM/D");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch() {}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue() {
  Ex = -1000;
  ELab = -1000;
  ThetaLab = -1000;
  ThetaCM = -1000;
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the AnalysisFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct() { return (NPL::VAnalysis*)new Analysis(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy {
 public:
  proxy() { NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct); }
};

proxy p;
}

