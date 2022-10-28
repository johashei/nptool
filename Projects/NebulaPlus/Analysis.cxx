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
 *  This class describe  Nebula analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
using namespace std;
#include"Analysis.h"
#include"TGraph.h"
#include"TCanvas.h"
#include"TPad.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
   NebulaPlus= (TNebulaPlusPhysics*) m_DetectorManager->GetDetector("NEBULAPLUS");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
/*
    TGraph* g = new TGraph();
    for(int i=0; i<NebulaPlus->DetectorNumber.size(); i++){
        //if(i==0)std::cout << "------------------" << std::endl;
        //std::cout << NebulaPlus->GetPos(i).X() << std::endl;
        g->SetPoint(i,NebulaPlus->GetPos(i).X(),NebulaPlus->GetPos(i).Y())
    }
*/
//    TGraph* g = new TGraph(NebulaPlus->DetectorNumber.size(),&(NebulaPlus->PosX[0]),&(NebulaPlus->PosY[0]));
    //    g->Fit("pol1");
    
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy{
  public:
    proxy(){
      NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
    }
};

proxy p;
}

