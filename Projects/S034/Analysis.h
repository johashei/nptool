#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2021    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta contact address: matta@lpccaen.in2p3.fr         *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  S034 analysis project                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include"NPVAnalysis.h"
#include"TMinosPhysics.h"
#include"TNebulaPhysics.h"
#include"TSamuraiBDCPhysics.h"
#include"TSamuraiFDC0Physics.h"
#include"TSamuraiFDC2Physics.h"
#include"TSamuraiHodoscopePhysics.h"

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();

    static NPL::VAnalysis* Construct();

  private:
    TMinosPhysics* Minos;
    TNebulaPhysics* Nebula;
    TSamuraiBDCPhysics* BDC;
    TSamuraiFDC0Physics* FDC0;
    TSamuraiFDC2Physics* FDC2;
    TSamuraiHodoscopePhysics* Hodo;
 
  private: // output variable
    double Brho;
    double Beta_f;
    int    Trigger;
  public:
    void  Clear();
    void  InitOutputBranch();
    void  InitInputBranch();
};

// use for broh calculation
double r_fit(double *x) ;
#endif
