#ifndef Analysis_h
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2025                                               *
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
#include "NPEnergyLoss.h"
#include "NPVAnalysis.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TReactionConditions.h"
#include "TSuperX3Physics.h"
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>

class Analysis : public NPL::VAnalysis {
 public:
  Analysis();
  ~Analysis();

 public:
  void Init();
  void TreatEvent();
  void End();
  void InitOutputBranch();
  void InitInputBranch();
  void ReInitValue();
  static NPL::VAnalysis* Construct();

 private:
  TSuperX3Physics* SuperX3;
};
#endif
