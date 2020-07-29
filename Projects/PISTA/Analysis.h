#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  PISTA analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "NPVAnalysis.h"
#include "TPISTAPhysics.h"
#include "TInitialConditions.h"
#include "TReactionConditions.h"
#include "TInteractionCoordinates.h"
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "TRandom3.h"
class Analysis: public NPL::VAnalysis{
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
    double OriginalBeamEnergy;
    double BeamEnergy;
    double R;
    double XTarget;
    double YTarget;
    double ZTarget;
    double OriginalElab;
    double Elab;
    double ThetaLab;
    double ThetaCM;
    double OptimumEx;
    double Ex;
    double PID;
    double OriginalThetaLab;
    NPL::Reaction* Transfer;

    TRandom3 Rand;
    double ThetaNormalTarget;
    double ThetaDetectorSurface;
    double TargetThickness;

    NPL::EnergyLoss Be10C;
    NPL::EnergyLoss U238C;

  private:
    TPISTAPhysics* PISTA;
    TInteractionCoordinates* InteractionCoordinates;
    TInitialConditions* InitialConditions;
    TReactionConditions* ReactionConditions;

};
#endif
