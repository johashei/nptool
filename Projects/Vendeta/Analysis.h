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
 *  This class describe  Vendeta analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include"NPVAnalysis.h"
#include"TVendetaPhysics.h"
#include"TFissionChamberPhysics.h"
#include "NPParticle.h"
#include "TRandom3.h"

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
    vector<double> inToF;
    vector<double> inEnergy;

				vector<double> LG_DT;
				vector<double> LG_ID;
    vector<double> LG_Anode_ID;
    vector<double> LG_ThetaLab;
    vector<double> LG_ELab;
    vector<double> LG_Tof;
    vector<double> LG_Q1;
    vector<double> LG_Q2;
    vector<double> LG_Qmax;
    vector<bool>   LG_FakeFission;

    vector<double> HG_DT;
    vector<double> HG_ID;
    vector<double> HG_Anode_ID;
    vector<double> HG_ThetaLab;
    vector<double> HG_ELab;
    vector<double> HG_Tof;
    vector<double> HG_Q1;
    vector<double> HG_Q2;
    vector<double> HG_Qmax;
    vector<bool>   HG_FakeFission;

		vector<double> FC_Q1;    
		vector<double> FC_Q2;    

  private:
    TVendetaPhysics* Vendeta;
    TFissionChamberPhysics* FC;

    NPL::Particle* neutron;
    TRandom3 ra;
};
#endif
