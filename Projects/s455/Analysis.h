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
 *  This class describe  Sofia analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
// Root
#include "TCutG.h"
#include "TRandom3.h"

// NPTool
#include"NPVAnalysis.h"
#include"TSofTofWPhysics.h"
#include"TSofTrimPhysics.h"
#include"TSofTwimPhysics.h"
#include"TSofSciPhysics.h"
#include"TSofBeamID.h"

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();
    void InitOutputBranch();
    void InitParameter();
    void ReInitValue();
    void BeamAnalysis();

    static NPL::VAnalysis* Construct();

  public:
    void LoadCut();
    int DetermineQmax();

  private:
    TSofSciPhysics* SofSci;
    TSofTrimPhysics* SofTrim;
    TSofTwimPhysics* SofTwim;
    TSofTofWPhysics* SofTofW;
    TSofBeamID* SofBeamID;

  private:
    double fLS2_0; 
    double fBrho0;
    double fDS2;
    double fDCC;
    double fK_LS2;

    TRandom3 rand;
    TCutG* cutQ78[3];
    TCutG* cutQ79[3];
    TCutG* cutQ80[3];
    TCutG* cutQ81[3];
};
#endif