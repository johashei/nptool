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
 *  This class describe  Scone analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include"NPVAnalysis.h"
#include"TSconePhysics.h"
#include"TInitialConditions.h"
#include"TInteractionCoordinates.h"

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();
    void InitInputBranch();
    void InitOutputBranch();
    void ReInitValue();

    static NPL::VAnalysis* Construct();

  private:
    double Time_max;
    double E_init;
    double E_max;
    double E_sum;
    double E_sum_gamma;
    double E_sum_proton;
    double E_mean_proton;
  
  private:
    TSconePhysics* Scone;
    TInitialConditions* InitialConditions;

  private:
    double m_DetectedNeutron;
    int m_entries;
    vector<double> vDetectedNeutron;

  private:
    bool new_energy;
    double E_new;
    vector<double> vE_init;

};

#endif
