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

  NmaxFP=11;
  //FP = new TBigRIPSFocal();
  //FP->Init(NmaxFP+1);
  Rec35 = new TBigRIPSReco();
  Rec57 = new TBigRIPSReco();
  Rec89 = new TBigRIPSReco();
  Rec911 = new TBigRIPSReco();

  InitOutputBranch();
  InitInputBranch();
  // get PPAC, PL, and IC  objects
  PPAC = (TBigRIPSPPACPhysics*)  m_DetectorManager -> GetDetector("BigRIPSPPAC");
  PL = (TBigRIPSPlasticPhysics*)  m_DetectorManager -> GetDetector("BigRIPSPlastic");
  IC = (TBigRIPSICPhysics*)  m_DetectorManager -> GetDetector("BigRIPSIC");

  ReadXmls();
/*
  matF35.ResizeTo(6,5); matF35 = RecReadTransferMatrix("mat1.mat");
  matF57.ResizeTo(6,5); matF57 = RecReadTransferMatrix("mat2.mat");
  matF89.ResizeTo(6,5); matF89 = RecReadTransferMatrix("F8F9_LargeAccAchr.mat");
  matF911.ResizeTo(6,5); matF911 = RecReadTransferMatrix("F9F11_LargeAccAchr.mat");
*/
  matrixF35  = RecReadTransferMatrix2("mat1.mat");
  matrixF57  = RecReadTransferMatrix2("mat2.mat");
  matrixF89  = RecReadTransferMatrix2("F8F9_LargeAccAchr.mat");
  matrixF911 = RecReadTransferMatrix2("F9F11_LargeAccAchr.mat");

  length37  = (FP_Z[FPtoID[7]]+PL_NametoZ["F7pl"]) 
             -(FP_Z[FPtoID[3]]+PL_NametoZ["F3pl"]); 
  length35  =  FP_Z[FPtoID[5]] - (FP_Z[FPtoID[3]]+PL_NametoZ["F3pl"]);
  length57  = (FP_Z[FPtoID[7]]+PL_NametoZ["F7pl"]) - FP_Z[FPtoID[5]];
  length811 = (FP_Z[FPtoID[11]]+PL_NametoZ["F11pl-1"]) 
             -(FP_Z[FPtoID[8]]+PL_NametoZ["F8pl"]); 
  length89  =  FP_Z[FPtoID[9]] - (FP_Z[FPtoID[8]]+PL_NametoZ["F8pl"]);
  length911 = (FP_Z[FPtoID[11]]+PL_NametoZ["F11pl-1"]) - FP_Z[FPtoID[9]];

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
      FP_IDtoFP[ID] = b[i]->AsInt("FPL"); 
      int fp = b[i]->AsInt("FPL");
      FPtoID[fp] = ID; 
      FP_Z[ID] = b[i]->AsDouble("zpos"); 
      FP_Zoffset[ID] = b[i]->AsDouble("zshift"); 
      if(FP_Zoffset[ID]==-1000) FP_Zoffset[ID] = 0; 
  }

  parser.LoadFile("BigRIPSPPAC.xml");
  vector<NPL::XML::block *> bppac = parser.GetAllBlocksWithName("BigRIPSPPAC");

  size = bppac.size();
  for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = bppac[i]->AsInt("ID"); 
      PPAC_IDtoFP[ID] = bppac[i]->AsInt("FPL"); 
      PPAC_XZoffset[ID] = bppac[i]->AsDouble("xzpos"); 
      PPAC_YZoffset[ID] = bppac[i]->AsDouble("yzpos"); 
  }

  parser.LoadFile("BigRIPSPlastic.xml");
  vector<NPL::XML::block *> bpl = parser.GetAllBlocksWithName("BigRIPSPlastic");
  size = bpl.size();
  for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = bpl[i]->AsInt("ID"); 
      PL_IDtoFP[ID] = bpl[i]->AsInt("FPL"); 
      PL_IDtoName[ID] = bpl[i]->AsString("NAME"); 
      string name = bpl[i]->AsString("NAME");
      PL_NametoID[name] = ID; 
      PL_NametoZ[name] = bpl[i]->AsDouble("zpos"); 
  }

  parser.LoadFile("BigRIPSIC.xml");
  vector<NPL::XML::block *> bic = parser.GetAllBlocksWithName("BigRIPSIC");
  size = bic.size();
  for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = bic[i]->AsInt("ID"); 
      IC_IDtoFP[ID] = bic[i]->AsInt("FPL"); 
      IC_IDtoName[ID] = bic[i]->AsString("NAME"); 
      string name = bic[i]->AsString("NAME");
      IC_NametoID[name] = ID; 
      IC_Zcoef[name].push_back(bic[i]->AsDouble("zcoef_0")); 
      IC_Zcoef[name].push_back(bic[i]->AsDouble("zcoef_1")); 
      IC_Zcoef[name].push_back(bic[i]->AsDouble("zcoef_2")); 
      IC_Ionpair[name] = bic[i]->AsDouble("ionpair"); 
  }
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
    // Reinitiate and clear variables
    ReInitValue();
    FP_PPACHitList.clear();

    /////////////////////////////////////////////
    /// STEP1: Trajectory/Brho reconstruction ///
    /////////////////////////////////////////////

    // Get all PPAC X,Y hits and store them in a HitList indexed by Focalplane
    // For a given Focal plane, the HitList size depends on how many PPAC 
    // are installed at this FP
    for(int i=0;i<PPAC->ID.size();i++){ 
        FP_PPACHitList[PPAC->FP[i]].AddPPACHit(PPAC->X[i],PPAC->Y[i],PPAC->ID[i]);
    }

    // For each FP, perform reco. and store in vector<vector>
    // RecFP[FPNbr][0,1,2,3 for x,a,y,b]
    for(auto it = FP_PPACHitList.begin();it!=FP_PPACHitList.end();++it){
        if(it->first < NmaxFP){
           RecFP[it->first] = RecFPposition(it->second.PPACHit,it->second.PPACID);
        }
    }
    // Same as above but using an object from class TBigRIPSFocal
    /*
    for(auto it = FP_PPACHitList.begin();it!=FP_PPACHitList.end();++it){
        if(it->first < NmaxFP){
           FP->SetFPTrack(it->first, RecFPposition(it->second.PPACHit,it->second.PPACID) );
           //FP->Print(it->first);
        }
    }
*/

    // Brho/delta reconstruction
    Rec35->RecBrho(RecFP[3],RecFP[5],matrixF35,BrhoD3);
    Rec57->RecBrho(RecFP[5],RecFP[7],matrixF57,BrhoD5);
    Rec89->RecBrho(RecFP[8],RecFP[9],matrixF89,BrhoD7);
    Rec911->RecBrho(RecFP[9],RecFP[11],matrixF911,BrhoD8);


    /////////////////////////////////
    /// STEP2: TOF, flight length ///
    /////////////////////////////////
    
    std::map<string,double> PL_T; //index is plastic name in xml 
    for(int i=0;i<PL->ID.size();i++) {PL_T[PL_IDtoName[PL->ID[i]]] = PL->T[i];}
    tf3 = PL_T["F3pl"];
    tf7 = PL_T["F7pl"];
    if(PL_T["F7pl"]!=-9999 && PL_T["F3pl"]!=-9999) tof37 = PL_T["F7pl"] - PL_T["F3pl"];
    tof37 += tof_offset37;

    tf8 = PL_T["F8pl"];
    tf11 = PL_T["F11pl-1"];
    if(PL_T["F11pl-1"]!=-9999 && PL_T["F8pl"]!=-9999) tof811 = PL_T["F11pl-1"] - PL_T["F8pl"];
    tof811 += tof_offset811;

    
    /////////////////////////////////
    /// STEP3: Beta, Gamma, AoQ   ///
    /////////////////////////////////

    // Calculate Beta from length/tof and combine it with Brho to get AoQ
    Rec35->RecAoqOne(tof37,length37);
    Rec57->RecAoqOne(tof37,length37);
    Rec89->RecAoqOne(tof811,length811);
    Rec911->RecAoqOne(tof811,length811);
    
    ///////////////////////////////////////////////
    /// STEP4: Z reco. from dE in IC and Beta   ///
    ///////////////////////////////////////////////
    
    std::map<string,double> IC_dE; //index is IC name in xml 
    for(int i=0;i<IC->ID.size();i++) {
        IC_dE[IC_IDtoName[IC->ID[i]]] = IC->CalSqSum[i];
        //if(IC->ID[i]==1) std::cout << "IC ID1 CalSq:" << IC->CalSqSum[i]  << std::endl;
    }
    dE_ICF7 = IC_dE["F7IC"];
    Rec35->RecZet(dE_ICF7,IC_Ionpair["F7IC"],IC_Zcoef["F7IC"]);
    Rec57->RecZet(dE_ICF7,IC_Ionpair["F7IC"],IC_Zcoef["F7IC"]);

    dE_ICF11 = IC_dE["F11IC"];
    Rec89->RecZet(dE_ICF11,IC_Ionpair["F11IC"],IC_Zcoef["F11IC"]);
    Rec911->RecZet(dE_ICF11,IC_Ionpair["F11IC"],IC_Zcoef["F11IC"]);

}

////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<double>> Analysis::RecReadTransferMatrix2(string matfile){
  char buffer[256];
  std::ifstream matfin(matfile); 
  if(matfin.fail()) cout<<"Failed to open optical matrix file: "<<matfile<<endl;
  matfin.getline(buffer,256);

  std::vector<std::vector<double>> matrix;
  //matrix.resize(5, vector<double>(6));

  double val[6];  
  std::vector<double> row;

  for(Int_t i=0;i<6;i++){
    row.clear();
    matfin>>val[0]>>val[1]>>val[2]>>val[3]>>val[4]>>val[5];
    //cout<<val[0]<<val[1]<<val[2]<<val[3]<<val[4]<<val[5]<<endl;
    for(Int_t j=0;j<5;j++) row.push_back(val[j]);
    matrix.push_back(row);
  }
  matfin.close();
  return matrix;
}

////////////////////////////////////////////////////////////////////////////////
// Function transforming a list of PPAC Hit positions and IDs for a given FP
// into  positions and angles at the FP position (FX,FY,FA,FB)
// requires PPAC_XZoffset, PPAC_YZoffset and FP_Zoffset to be defined for all id's;
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
            z_x = PPAC_XZoffset[id] - FP_Zoffset[id];
            z_y = PPAC_YZoffset[id] - FP_Zoffset[id];
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

  RootOutput::getInstance()->GetTree()->Branch("RunNumber",&RunNumber,"RunNumber/I");
  RootOutput::getInstance()->GetTree()->Branch("EventNumber",&EventNumber,"EventNumber/I");
  RootOutput::getInstance()->GetTree()->Branch("Trigger",&Trigger,"Trigger/I");
  RootOutput::getInstance()->GetTree()->Branch("TimeStamp",&TimeStamp,"TimeStamp/I");
  RootOutput::getInstance()->GetTree()->Branch("RecFP",&RecFP);
  //RootOutput::getInstance()->GetTree()->Branch("FP","TBigRIPSFocal",&FP);
  RootOutput::getInstance()->GetTree()->Branch("Rec35","TBigRIPSReco",&Rec35);
  RootOutput::getInstance()->GetTree()->Branch("Rec57","TBigRIPSReco",&Rec57);
  RootOutput::getInstance()->GetTree()->Branch("Rec89","TBigRIPSReco",&Rec89);
  RootOutput::getInstance()->GetTree()->Branch("Rec911","TBigRIPSReco",&Rec911);
  RootOutput::getInstance()->GetTree()->Branch("tof37",&tof37,"tof37/D");
  RootOutput::getInstance()->GetTree()->Branch("tof811",&tof811,"tof811/D");
  RootOutput::getInstance()->GetTree()->Branch("tf7",&tf7,"tf7/D");
  RootOutput::getInstance()->GetTree()->Branch("tf8",&tf8,"tf8/D");
  RootOutput::getInstance()->GetTree()->Branch("tf11",&tf11,"tf11/D");
  RootOutput::getInstance()->GetTree()->Branch("dE_ICF7",&dE_ICF7,"dE_ICF7/D");
  RootOutput::getInstance()->GetTree()->Branch("dE_ICF11",&dE_ICF11,"dE_ICF11/D");
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput:: getInstance()->GetChain()->SetBranchAddress("RunNumber",&RunNumber);
  RootInput:: getInstance()->GetChain()->SetBranchAddress("EventNumber",&EventNumber);
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){

  RecFP.clear();
  std::vector<double> reset;
  for(int j=0;j<4;j++) reset.push_back(-9999);
  for(int i=0;i<NmaxFP+1;i++) RecFP.push_back(reset);

  //FP->Clear();
  //FP->Init(NmaxFP+1);
  Rec35->Init();
  Rec57->Init();
  Rec89->Init();
  Rec911->Init();
  tf3 =-9999;
  tf7 =-9999;
  tof37 = -9999 ;
  dE_ICF7 = -9999;
  tf8 =-9999;
  tf11 =-9999;
  tof811 = -9999 ;
  dE_ICF11 = -9999;

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

