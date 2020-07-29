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
#include "NPVAnalysis.h"
#include "TChiNuPhysics.h"
#include "TInitialConditions.h"
#include "TInteractionCoordinates.h"
#include "TReactionConditions.h"
#include "NPNucleus.h"
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
	double Einit;
	vector<double> ThetaLab;
        vector<double> Ex;
        vector<double> Elab;
	vector<double> E;
	vector<double> T;
	vector<double> Distance;
	vector<int>    Det;

  private:
	TChiNuPhysics* m_ChiNu;
	TInitialConditions* InitialConditions;
	TInteractionCoordinates* InteractionCoordinates;
        TReactionConditions* ReactionConditions;

        NPL::Reaction* my_Reaction;
   	NPL::Nucleus* neutron;

        TRandom3 ra;
};
#endif
