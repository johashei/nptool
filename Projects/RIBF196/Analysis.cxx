/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Freddy Flavigny  contact: flavigny@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : may 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Basic AoQ,Z reconstruction for RIBF196 experiment                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "Analysis.h"
#include "TMatrixD.h"
#include "NPFunction.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPOptionManager.h"

#include<iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {

  InitOutputBranch();
  InitInputBranch();
  // get PPAC, PL, and IC  objects
  PPAC = (TBigRIPSPPACPhysics*)  m_DetectorManager -> GetDetector("BigRIPSPPAC");
  PL = (TBigRIPSPlasticPhysics*)  m_DetectorManager -> GetDetector("BigRIPSPlastic");
  IC = (TBigRIPSICPhysics*)  m_DetectorManager -> GetDetector("BigRIPSIC");

  ReadXmls();

  // initialize various parameters
  //Rand = TRandom3();
  //DetectorNumber = 0;

}

////////////////////////////////////////////////////////////////////////////////
// Temporary Function to read some parameter necessary for reconstruction
// Could be improved by using directly xml files defined in detector file (same)
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReadXmls() {

  NPL::XmlParser parser;
  parser.LoadFile("FocalPlane.xml");
  vector<NPL::XML::block *> b = parser.GetAllBlocksWithName("FocalPlane");

  unsigned int size = b.size();
  for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = b[i]->AsInt("ID"); 
      FPL_IDtoFP[ID] = b[i]->AsInt("FPL"); 
      FPL_Z[ID] = b[i]->AsDouble("zpos"); 
      FPL_Zoffset[ID] = b[i]->AsDouble("zshift"); 
      if(FPL_Zoffset[ID]==-1000) FPL_Zoffset[ID] = 0; 
  }

  parser.LoadFile("BigRIPSPPAC.xml");
  vector<NPL::XML::block *> bb = parser.GetAllBlocksWithName("BigRIPSPPAC");

  size = bb.size();
  for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = bb[i]->AsInt("ID"); 
      PPAC_IDtoFP[ID] = bb[i]->AsInt("FPL"); 
      PPAC_XZoffset[ID] = bb[i]->AsDouble("xzpos"); 
      PPAC_YZoffset[ID] = bb[i]->AsDouble("yzpos"); 
  }
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
    // Reinitiate calculated variable
    ReInitValue();

    FP_PPACHitList.clear();

    for(int i=0;i<PPAC->ID.size();i++){ 
        FP_PPACHitList[PPAC->FP[i]].AddPPACHit(PPAC->X[i],PPAC->Y[i],PPAC->ID[i]);
    }

    for(auto it = FP_PPACHitList.begin();it!=FP_PPACHitList.end();++it){
        std::vector<double> Rec;
        Rec = RecFPposition(it->second.PPACHit,it->second.PPACID);
        //it->second.Print();

        if(it->first==8){
            fX=Rec[0];
            fA=Rec[1];
            fY=Rec[2];
            fB=Rec[3];
        }
    } //end for loop on FP
}

/*
   std::cout << "______________________" << std::endl;
   std::cout << "ID:" << PPAC->ID[i]  << std::endl;
   std::cout << "FP:" << PPAC->FP[i]  << std::endl;
   std::cout << "X:" << PPAC->X[i]  << std::endl;
   std::cout << "Y:" << PPAC->Y[i]  << std::endl;
   std::cout << "------------------:" << std::endl;
   std::cout << "FP:" << it->first << std::endl;
   it->second.Print();
   */

////////////////////////////////////////////////////////////////////////////////
// Function transforming a list of PPAC Hit positions and IDs for a given FP
// into  positions and angles at the FP position (FX,FY,FA,FB)
// requires PPAC_XZoffset, PPAC_YZoffset and FPL_Zoffset to be defined for all id's;
////////////////////////////////////////////////////////////////////////////////
std::vector<double> Analysis::RecFPposition(std::vector<TVector2> HitList,std::vector<int> IdList){
    TMatrixD xvec(2,1); xvec.Zero();
    TMatrixD yvec(2,1); yvec.Zero();
    TMatrixD xmat(2,2); xmat.Zero();
    TMatrixD ymat(2,2); ymat.Zero();

    Double_t sum_zx[2] = {0,0};
    Double_t sum_x[2] = {0,0};
    Double_t sum_z[2] = {0,0};
    Double_t sum_zz[2] = {0,0};
    Double_t x, y, z_x, z_y;
    Double_t id;

    Int_t nfired_ppacx = 0;
    Int_t nfired_ppacy = 0;
    Int_t nfired_ppacx_u = 0; // up stream side
    Int_t nfired_ppacy_u = 0;
    Int_t nfired_ppacx_d = 0; // down stream side
    Int_t nfired_ppacy_d = 0;

    int nppac = HitList.size(); 

    //XY focal plane rec
    std::vector<double> Rec;
    for(int i=0; i<4; i++) Rec.push_back(-99999);

    if(nppac==1){
        Rec[0]= HitList[0].X();
        Rec[1]= -99999;
        Rec[2]= HitList[0].Y();
        Rec[3]= -99999;
    }        
    else if(nppac>1){        
        for(int i=0;i<nppac;i++){
            x = HitList[i].X();
            y = HitList[i].Y();
            id = IdList[i];
            z_x = PPAC_XZoffset[id] - FPL_Zoffset[id];
            z_y = PPAC_YZoffset[id] - FPL_Zoffset[id];
            if(HitList[i].X()!=-99999){
                xvec(0,0) += z_x * x; // b(1) in rayfit
                xvec(1,0) += x; // b(2) in rayfit
                xmat(0,1) += z_x; // a(1,2) in rayfit
                xmat(1,0) += z_x; // a(1,2) in rayfit
                xmat(0,0) += z_x * z_x;   // a(1,1) in rayfit
                xmat(1,1) ++;  // a(2,2) in rayfit
                nfired_ppacx ++;  // a(2,2) in rayfit
                if(i<2) nfired_ppacx_u ++; else nfired_ppacx_d ++; 
            }
            if(HitList[i].Y()!=-99999){
                yvec(0,0) += z_y * y; // b(1) in rayfit
                yvec(1,0) += y; // b(2) in rayfit
                ymat(0,1) += z_y; // a(1,2) in rayfit
                ymat(1,0) += z_y; // a(1,2) in rayfit
                ymat(0,0) += z_y * z_y;   // a(1,1) in rayfit
                ymat(1,1) ++;  // a(2,2) in rayfit
                nfired_ppacy ++;  // a(2,2) in rayfit
                if(i<2) nfired_ppacy_u ++; else nfired_ppacy_d ++; 
            }
        }

        // determine a and b in x = a + b * z (i.e.); 
        if((nppac==2&&nfired_ppacx>1)||
                (nppac==4&&nfired_ppacx_u>0&&nfired_ppacx_d>0) // 20141023 by TI 
          ){
            TMatrixD rxvec = xmat.Invert() * xvec;
            Rec[0] = rxvec(1,0);
            Rec[1] = TMath::ATan(rxvec(0,0)) * 1000;
        }
        else{
            Rec[0] = -99999;
            Rec[1] = -99999;
        }
        if((nppac==2&&nfired_ppacy>1)||
                (nppac==4&&nfired_ppacy_u>0&&nfired_ppacy_d>0) // 20141023 by TI 
          ){
            TMatrixD ryvec = ymat.Invert() * yvec;
            Rec[2] = ryvec(1,0);
            Rec[3]= TMath::ATan(ryvec(0,0)) * 1000; // the value is given in mrad.;
        }
        else{
            Rec[2] = -99999;
            Rec[3] = -99999;
        }

    }//end if(nppac>1)

    return Rec;

}
////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {

  RootOutput::getInstance()->GetTree()->Branch("fX",&fX,"fX/D");
  RootOutput::getInstance()->GetTree()->Branch("fY",&fY,"fY/D");
/*
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("EDC",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab,"ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Run",&Run,"Run/I");
  RootOutput::getInstance()->GetTree()->Branch("X",&X,"X/D");
  RootOutput::getInstance()->GetTree()->Branch("Y",&Y,"Y/D");
  RootOutput::getInstance()->GetTree()->Branch("Z",&Z,"Z/D");
  RootOutput::getInstance()->GetTree()->Branch("dE",&dE,"dE/D");
  if(!simulation){
  // Vamos 
  RootOutput::getInstance()->GetTree()->Branch("LTS",&LTS,"LTS/l");

  // Agata
  // Time stamp of the agata trigger
  RootOutput::getInstance()->GetTree()->Branch("TStrack",&TStrack,"TStrack/l");

  // Array of reconstructed tracks
  RootOutput::getInstance()->GetTree()->Branch("nbTrack",&nbTrack,"nbTrack/I");
  RootOutput::getInstance()->GetTree()->Branch("trackE",trackE,"trackE[nbTrack]/F");
  RootOutput::getInstance()->GetTree()->Branch("trackX1",trackX1,"trackX1[nbTrack]/F");
  RootOutput::getInstance()->GetTree()->Branch("trackY1",trackY1,"trackY1[nbTrack]/F");
  RootOutput::getInstance()->GetTree()->Branch("trackZ1",trackZ1,"trackZ1[nbTrack]/F");
  RootOutput::getInstance()->GetTree()->Branch("trackT",trackT,"trackT[nbTrack]/F");
  RootOutput::getInstance()->GetTree()->Branch("trackCrystalID",trackCrystalID,"trackCrystalID[nbTrack]/I");

  // Array of reconstructed core
  RootOutput::getInstance()->GetTree()->Branch("nbCores",&nbCores,"nbCores/I");
  RootOutput::getInstance()->GetTree()->Branch("coreId",coreId,"coreId[nbCores]/I");
  RootOutput::getInstance()->GetTree()->Branch("coreTS",coreTS,"coreTS[nbCores]/l");
  RootOutput::getInstance()->GetTree()->Branch("coreE0",coreE0,"coreE0[nbCores]/F");
  //
  }
  else{
    RootOutput::getInstance()->GetTree()->Branch("OriginalELab",&OriginalELab,"OriginalELab/D");
    RootOutput::getInstance()->GetTree()->Branch("OriginalThetaLab",&OriginalThetaLab,"OriginalThetaLab/D");
    RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",&BeamEnergy,"BeamEnergy/D");
  }
*/
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  // RootInput:: getInstance()->GetChain()->SetBranchAddress("GATCONF",&vGATCONF);
/*
 *
  if(!simulation){
    // Vamos
    RootInput::getInstance()->GetChain()->SetBranchAddress("LTS",&LTS);
    // Agata
    RootInput::getInstance()->GetChain()->SetBranchAddress("TStrack",&TStrack);
    RootInput::getInstance()->GetChain()->SetBranchAddress("nbTrack",&nbTrack);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackE",trackE);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackX1",trackX1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackY1",trackY1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackZ1",trackZ1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackT",trackT);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackCrystalID",trackCrystalID);
    RootInput::getInstance()->GetChain()->SetBranchAddress("nbCores",&nbCores);
    RootInput::getInstance()->GetChain()->SetBranchAddress("coreId",coreId);
    RootInput::getInstance()->GetChain()->SetBranchAddress("coreTS",coreTS);
    RootInput::getInstance()->GetChain()->SetBranchAddress("coreE0",coreE0);

  }
  else{

    RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
    RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&Initial);
    RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
    RootInput:: getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&ReactionConditions);

  }
*/
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  fX = -1000 ;
  fY = -1000 ;
  fA = -1000 ;
  fB = -1000 ;
/*
  Ex = -1000 ;
  EDC= -1000;
  ELab = -1000;
  BeamEnergy = -1000;
  ThetaLab = -1000;
  ThetaCM = -1000;
  X = -1000;
  Y = -1000;
  Z = -1000;
  dE= -1000;
*/
}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the AnalysisFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_analysis{
    public:
      proxy_analysis(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
      }
  };

  proxy_analysis p_analysis;
}

