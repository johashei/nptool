#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr *
 *                                                                             *
 * Creation Date  : 06/2021                                             *
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
    void BeamFragmentAnalysis();

    static NPL::VAnalysis* Construct();

  public:
    void LoadSpline();

  private:
    TSofSciPhysics* SofSci;
    TSofMwpcPhysics* SofMwpc;
    TSofTrimPhysics* SofTrim;
    TSofAtPhysics* SofAt;
    TSofTwimPhysics* SofTwim;
    TSofTofWPhysics* SofTofW;
    TSofBeamID* SofBeamID;
    TSofFissionFragment* SofFF;
    int RunID;

  private:
    int fRunID;
    double fLS2_0; 
    double fBrho0;
    double fDS2;
    double fDCC;
    double fK_LS2;
    double fZbeam_p0; 
    double fZbeam_p1; 
    double fZbeam_p2; 
    double fZff_p0;
    double fZff_p1;
    double fZff_p2;
    double fZBeta_p0;
    double fZBeta_p1;

    TRandom3 rand;
    
    TSpline3* fcorr_z_beta[4];
    TSpline3* fcorr_z_dt[4];

};
#endif
