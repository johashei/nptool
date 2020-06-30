/*****************************************************************************
 * Copyright (C) 2009   this file is part of the NPTool Project              *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 22/07/09                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *    + This macro displays everything concerning the incident beam and the  *
 *      emitted particle from NPSimulation                                   *
 *                                                                           *
 *    + Use in a ROOT session:                                               *
 *      .x ControlSimu.C("FileToAnalyse")                                    *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <iostream>
#include <ctime>
#include <cstdlib>
using namespace std;

// ROOT headers
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

// nptool headers
#include "TInitialConditions.h"
#include "TInteractionCoordinates.h"
//#include "NPDetectorManager.h"
#include "NPReaction.h"
using namespace NPL;

void ShowResults(){
  // get tree   
  TFile* f = new TFile("../../Outputs/Analysis/PhysicsTree.root");
  TTree* t = (TTree*) f->Get("PhysicsTree");

  // draw kinematic information
  // canvas
  TCanvas *c1 = new TCanvas("c1", "Results", 1000, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  // kinematic line
  TH2F* hk = new TH2F("hk", "hk", 180*3, 0, 180, 1000, 0, 60);
  t->Draw("ELab:ThetaLab>>hk","","col");
  hk->GetXaxis()->SetTitle("#Theta_{lab} (deg)");
  hk->GetYaxis()->SetTitle("E_{p} (MeV)");
  NPL::Reaction* reaction = new NPL::Reaction();
  reaction->ReadConfigurationFile("16Odp17O_870keV_12.reaction");
  reaction->GetKinematicLine3()->Draw("c");

  c1->cd(2);
  TH1F* hEx = new TH1F("hEx", "hEx",1000, -1, 5);
  t->Draw("Ex>>hEx","ThetaLab>100 && ThetaLab<156","col");
  hEx->GetXaxis()->SetTitle("Ex (MeV)");
  hEx->GetYaxis()->SetTitle("counts/100 keV");
 
  c1->cd(3);
  TH1F *hCM = new TH1F("hCM", "hCM", 36, 0, 180); 
  t->Draw("ThetaCM>>hCM","Ex>0&&Ex<6","");
  for(int i = 0 ; i < hCM->GetNbinsX();i++){
    if(hCM->GetBinCenter(i)==37.5 || hCM->GetBinCenter(i)==97.5|| hCM->GetBinCenter(i)==167.5|| hCM->GetBinCenter(i)==42.5){
      hCM->SetBinContent(i,0);

    }
  }

  
  TCanvas *c2 = new TCanvas("c2", "Control", 1000, 1000);
  c2->Divide(2,2);
  c2->cd(1);
  TH1F* hcT = new TH1F("hcT", "hcT", 180*3, -1,1);
  t->Draw("ThetaLab-OriginalThetaLab>>hcT","","col");
  TLine* lT = new TLine(0,0,180,180);
  //lT->Draw();
  c2->cd(2);
  TH1F* hcE = new TH1F("hcE", "hcE", 1000, -1, 1);
  t->Draw("ELab-OriginalELab>>hcE","","col");
  TLine* lE = new TLine(0,0,60,60);
  //lE->Draw();


}

