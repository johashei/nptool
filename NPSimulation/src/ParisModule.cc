/*****************************************************************************
 * Copyright (C) 2009   this file is part of the NPTool Project              *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 04/12/09                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: This class is an Abstract Base Class (ABC) from which should  *
 *             derive all different modules from the Paris detector.         *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "ParisModule.hh"
#include "RootOutput.h"


TParisData *ParisModule::ms_Event = 0;



ParisModule::ParisModule()
{
   if (ms_Event == 0) ms_Event = new TParisData();

   InitializeRootOutput();
   InitializeIndex();
}



ParisModule::~ParisModule()
{
}



void ParisModule::InitializeRootOutput()
{
   RootOutput *pAnalysis = RootOutput::getInstance();
   TTree *pTree = pAnalysis->GetTree();
   // if the branch does not exist yet, create it
   if (!pTree->GetBranch("PARIS"))
      pTree->Branch("PARIS", "TParisData", &ms_Event);
}



void ParisModule::InitializeIndex()
{
   m_index["Phoswitch"]     =    0;
}
