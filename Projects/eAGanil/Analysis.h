#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adriment Matta contact address: matta@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : Octover 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  eAGanil analysis project                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include"NPVAnalysis.h"
#include"TeAGanilPhysics.h"
#include"TInteractionCoordinates.h"
#include"TReactionConditions.h"
#include"TRandom3.h"
#include"NPReaction.h"
#include<vector>
using namespace std;
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
   TRandom3 Rand; 
   NPL::Reaction m_reaction;
   vector<double> Ex,ELab,ThetaLab,Resolution;
   vector<int> Detected;
  private:
   TeAGanilPhysics* eAGanil;
   TInteractionCoordinates* Inter;   
   TReactionConditions* Initial;   
};
#endif
