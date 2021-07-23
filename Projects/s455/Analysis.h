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
#include "TSpline.h"

// NPTool
#include"NPVAnalysis.h"
#include"TSofTofWPhysics.h"
#include"TSofTrimPhysics.h"
#include"TSofAtPhysics.h"
#include"TSofTwimPhysics.h"
#include"TSofSciPhysics.h"
#include"TSofMwpcPhysics.h"
#include"TSofBeamID.h"
#include"TSofFissionFragment.h"

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
    void FissionFragmentAnalysis();

    static NPL::VAnalysis* Construct();

  public:
    void LoadCut();
    void LoadSpline();
    int DetermineQmax();

  private:
    TSofSciPhysics* SofSci;
    TSofMwpcPhysics* SofMwpc;
    TSofTrimPhysics* SofTrim;
    TSofAtPhysics* SofAt;
    TSofTwimPhysics* SofTwim;
    TSofTofWPhysics* SofTofW;
    TSofBeamID* SofBeamID;
    TSofFissionFragment* SofFF;

  private:
    double fLS2_0; 
    double fBrho0;
    double fDS2;
    double fDCC;
    double fK_LS2;

    TRandom3 rand;
    TCutG* cutQ77[3];
    TCutG* cutQ78[3];
    TCutG* cutQ79[3];
    TCutG* cutQ80[3];
    TCutG* cutQ81[3];
    
    TSpline3* fcorr_z_beta[4];

};
#endif
