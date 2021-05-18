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

  matF35.ResizeTo(6,5); matF35 = RecReadTransferMatrix("mat1.mat");
  matF57.ResizeTo(6,5); matF57 = RecReadTransferMatrix("mat2.mat");

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
  }
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
    // Reinitiate calculated variable
    ReInitValue();
    FP_PPACHitList.clear();

    /////////////////////////////////////////////
    /// STEP1: Trajectory/Brho reconstruction ///
    /////////////////////////////////////////////

    // Get all PPAC X,Y hits and store them in a HitList indexed by Focalplane
    for(int i=0;i<PPAC->ID.size();i++){ 
        FP_PPACHitList[PPAC->FP[i]].AddPPACHit(PPAC->X[i],PPAC->Y[i],PPAC->ID[i]);
    }

    // Perform FP reco. and store in map: index FP, content vector(x,a,y,b)
    std::map<int,vector<double>> RecFP;
    double NmaxFP=12;
    for(int i=0;i<NmaxFP;i++){
        for(int j=0;j<4;j++){
            RecFP[i].push_back(-9999);}}

    for(auto it = FP_PPACHitList.begin();it!=FP_PPACHitList.end();++it){
        RecFP[it->first] = RecFPposition(it->second.PPACHit,it->second.PPACID);
    }

    fX=RecFP[8][0];
    fA=RecFP[8][1];
    fY=RecFP[8][2];
    fB=RecFP[8][3];

    delta35 = RecDeltaBrho(RecFP[3],RecFP[5],matF35);
    Brho35 = BrhoD3 * (1.0 + delta35*0.01);

    delta57 = RecDeltaBrho(RecFP[5],RecFP[7],matF57);
    Brho57 = BrhoD5 * (1.0 + delta57*0.01);

    /////////////////////////////////
    /// STEP2: TOF, flight length ///
    /////////////////////////////////
    
    std::map<string,double> PL_T; //index is plastic name in xml 
    for(int i=0;i<PL->ID.size();i++) {PL_T[PL_IDtoName[PL->ID[i]]] = PL->T[i];
        //if(PL->ID[i]==2) std::cout << "PL ID2 T::" << PL->T[i]  << std::endl;
    }
    tf3 = PL_T["F3pl"];
    tf7 = PL_T["F7pl"];
    if(PL_T["F7pl"]!=-9999 && PL_T["F3pl"]!=-9999)  tof37 = PL_T["F7pl"] - PL_T["F3pl"];
    tof37 += tof_offset37;
    
    //double length35 = (FP_Z[5]+PL_NametoZ["F5pl"]) - (FP_Z[3]+PL_NametoZ["F3pl"]) ; 
    //double length57 = (FP_Z[7]+PL_NametoZ["F7pl"]) - (FP_Z[5]+PL_NametoZ["F5pl"]) ; 
    //double length37 = length35 + length57 ; 

    double length37 =  (FP_Z[FPtoID[7]]+PL_NametoZ["F7pl"]) - (FP_Z[FPtoID[3]]+PL_NametoZ["F3pl"]); 
    double length35 =   FP_Z[FPtoID[5]] - (FP_Z[FPtoID[3]]+PL_NametoZ["F3pl"]);
    double length57 =  (FP_Z[FPtoID[7]]+PL_NametoZ["F7pl"]) - FP_Z[FPtoID[5]];

    //std::cout << "length37:" << length37  << std::endl;
    //std::cout << "Z FP7:" << FP_Z[FPtoID[7]] << std::endl;
    //std::cout << "Zoff PlasticF7:" << PL_NametoZ["F7pl"] << std::endl;
    //std::cout << "Z FP3:" << FP_Z[FPtoID[3]] << std::endl;
    //std::cout << "Zoff PlasticF3:" << PL_NametoZ["F3pl"] << std::endl;

    beta35 = length35 /(tof37/2. * clight);
    beta57 = length57 /(tof37/2. * clight);
    beta37 = length37 /(tof37 * clight);

    /////////////////////////////////
    /// STEP3: Beta, Gamma, AoQ   ///
    /////////////////////////////////

    //   std::cout << "______________________" << std::endl;
    //   std::cout << "Brho35:" << Brho35  << std::endl;
    //   std::cout << "Tof37:" << tof37  << std::endl;
    //   std::cout << "length35:" << length35  << std::endl;
    aoq35 = RecAoqOneFold(Brho35, tof37, length37);
    //   std::cout << "aoq35:" << aoq35  << std::endl;

    //   std::cout << "______________________" << std::endl;
    //   std::cout << "Brho57:" << Brho57  << std::endl;
    //   std::cout << "Tof37:" << tof37  << std::endl;
    //   std::cout << "length57:" << length57  << std::endl;
    aoq57 = RecAoqOneFold(Brho57, tof37, length37);
    //   std::cout << "aoq57:" << aoq57  << std::endl;
    
    //////////////////
    /// STEP4: Z   ///
    //////////////////
    
    std::map<string,double> IC_dE; //index is IC name in xml 
    for(int i=0;i<IC->ID.size();i++) {IC_dE[IC_IDtoName[IC->ID[i]]] = IC->CalSqSum[i];
        //if(IC->ID[i]==1) std::cout << "IC ID1 CalSq:" << IC->CalSqSum[i]  << std::endl;
    }
    dE_ICF7 = IC_dE["F7IC"];
    z_BR = RecZ(dE_ICF7, tof37, length37);

}

////////////////////////////////////////////////////////////////////////////////
TMatrixD Analysis::RecReadTransferMatrix(string matfile){
  char buffer[256];
  std::ifstream matfin(matfile); 
  if(matfin.fail()) cout<<"Failed to open optical matrix file: "<<matfile<<endl;
  matfin.getline(buffer,256);

  TMatrixD matrix(6,5);
  double val[6];  

  for(Int_t i=0;i<6;i++){
    matfin>>val[0]>>val[1]>>val[2]>>val[3]>>val[4]>>val[5];
    //cout<<val[0]<<val[1]<<val[2]<<val[3]<<val[4]<<val[5]<<endl;
    for(Int_t j=0;j<5;j++)
      matrix(i,j) = val[j];
  }
  matfin.close();

  matrix.Print();
  return matrix;
}
////////////////////////////////////////////////////////////////////////////////
double Analysis::RecDeltaBrho(std::vector<double> RecFPUpstream,std::vector<double> RecFPDownstream, TMatrixD matrix){

    if(matrix.GetNrows()!=6 && matrix.GetNcols()!=5){
        std::cout << "MATRIX used for Brho Rec. is NULL" << std::endl;
        return -9999;
    }

    double x_a = matrix(0,0);

    TMatrixD mat1(2,1);
    mat1(0,0) = matrix(0,0); mat1(1,0) = matrix(0,1);
    //mat1.Print();

    TMatrixD mat2(2,2);
    mat2(0,0) = matrix(1,0); mat2(0,1) = matrix(5,0);
    mat2(1,0) = matrix(1,1); mat2(1,1) = matrix(5,1);
    //mat2.Print();
    mat2.Invert();

    //TVectorD * fvec = (TVectorD *)fUpstreamFplArrayBuffer[i]->GetOptVector();
    //TVectorD * bvec = (TVectorD *)fDownstreamFplArrayBuffer[i]->GetOptVector();

    if(RecFPUpstream[0]==-9999 || 
       RecFPUpstream[1]==-9999 || 
       RecFPDownstream[0]==-9999)
        return -9999;

    TMatrixD xvec(2,1);
    //xvec(0,0) = (*bvec)(0);
    //xvec(1,0) = (*bvec)(1);
    xvec(0,0) = RecFPDownstream[0]; // downst. X
    xvec(1,0) = RecFPDownstream[1]; // downst. A

    TMatrixD rvec = mat2 * (xvec - RecFPUpstream[0] * mat1);
    Double_t angle = rvec(0,0); // change to mrad
    Double_t delta = rvec(1,0);
    // simple_delta = (*bvec)(0) / 33.; // only for fpl-5

    //double BrhoRec = BrhoCentral*(1.0+delta*0.01);
    //rips->SetDelta(delta);
    //rips->SetAngle(angle);


    //fReconstructed = true;
    return delta;
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
double Analysis::RecAoqOneFold(double Brho, double tof, double length){
    double beta = length /(tof * clight);
    double gamma = 1./sqrt(1 - beta*beta);
    double aoq = (Brho * clight) / (mnucleon * beta * gamma);
    return aoq;
}

////////////////////////////////////////////////////////////////////////////////
double Analysis::RecZ(double dE, double tof, double length){
 double ionpair = 4866.;
 double beta = length /(tof * clight);
 Double_t de_v = log(ionpair*beta*beta) - log((1-beta*beta)) - beta*beta;
 double z = 17.9143234533*sqrt(dE/de_v)*beta -13.0990490522;   
 return z;
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {

  RootOutput::getInstance()->GetTree()->Branch("EventNumber",&EventNumber,"EventNumber/I");
  RootOutput::getInstance()->GetTree()->Branch("RunNumber",&RunNumber,"RunNumber/I");
  RootOutput::getInstance()->GetTree()->Branch("fX",&fX,"fX/D");
  RootOutput::getInstance()->GetTree()->Branch("fY",&fY,"fY/D");
  RootOutput::getInstance()->GetTree()->Branch("delta35",&delta35,"delta35/D");
  RootOutput::getInstance()->GetTree()->Branch("delta57",&delta57,"delta57/D");
  RootOutput::getInstance()->GetTree()->Branch("Brho35",&Brho35,"Brho35/D");
  RootOutput::getInstance()->GetTree()->Branch("Brho57",&Brho57,"Brho57/D");
  RootOutput::getInstance()->GetTree()->Branch("tof37",&tof37,"tof37/D");
  RootOutput::getInstance()->GetTree()->Branch("aoq35",&aoq35,"aoq35/D");
  RootOutput::getInstance()->GetTree()->Branch("aoq57",&aoq57,"aoq57/D");
  RootOutput::getInstance()->GetTree()->Branch("aoq37",&aoq57,"aoq37/D");
  RootOutput::getInstance()->GetTree()->Branch("tf3",&tf3,"tf3/D");
  RootOutput::getInstance()->GetTree()->Branch("tf7",&tf7,"tf7/D");
  RootOutput::getInstance()->GetTree()->Branch("beta35",&beta35,"beta35/D");
  RootOutput::getInstance()->GetTree()->Branch("beta57",&beta57,"beta57/D");
  RootOutput::getInstance()->GetTree()->Branch("beta37",&beta57,"beta37/D");
  RootOutput::getInstance()->GetTree()->Branch("z_BR",&z_BR,"z_BR/D");
  RootOutput::getInstance()->GetTree()->Branch("dE_ICF7",&dE_ICF7,"dE_ICF7/D");
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput:: getInstance()->GetChain()->SetBranchAddress("RunNumber",&RunNumber);
  RootInput:: getInstance()->GetChain()->SetBranchAddress("EventNumber",&EventNumber);
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
  fX = -9999 ;
  fY = -9999 ;
  fA = -9999 ;
  fB = -9999 ;
  delta35 = -9999 ;
  delta57 = -9999 ;
  Brho35 = -9999 ;
  Brho57 = -9999 ;
  tf3 =-9999;
  tf7 =-9999;
  tof37 = -9999 ;
  aoq35 = -9999 ;
  aoq57 = -9999 ;
  aoq37 = -9999 ;

  beta35 = -9999 ;
  beta57 = -9999 ;
  beta37 = -9999 ;
  dE_ICF7 = -9999;
  z_BR = -9999 ;

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

