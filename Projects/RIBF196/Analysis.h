#ifndef Analysis_h 
#define Analysis_h
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
#include"NPVAnalysis.h"
//#include"NPEnergyLoss.h"
//#include"NPReaction.h"
#include"RootOutput.h"
#include"RootInput.h"
//#include "TInitialConditions.h"
//#include "TReactionConditions.h"
#include "TBigRIPSPPACPhysics.h"
#include "TBigRIPSPlasticPhysics.h"
#include "TBigRIPSICPhysics.h"
#include <TRandom3.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TMath.h>

//Little class to gather all PPACHits for a given Focal plane
class BigRIPSPPACHitList{
  public:
   BigRIPSPPACHitList(){Clear();};  
   ~BigRIPSPPACHitList(){};  

  public:
    std::vector<TVector2> PPACHit;
    std::vector<int> PPACID;

    void AddPPACHit(double x, double y, double ID){
            TVector2 tmp(x,y);
            PPACHit.push_back(tmp);
            PPACID.push_back(ID);
    };
    TVector2 GetPPACHit(int i){return PPACHit[i];}
    int GetPPACID(int i){return PPACID[i];}
    void Clear(){PPACHit.clear();PPACID.clear();};

    void Print(){
        cout << "----PPACHitList----size:" << PPACHit.size()<<endl;
        for (UShort_t i = 0; i < PPACHit.size(); i++){
            cout << "ID: " << PPACID[i] << endl;
            cout << "X: " << PPACHit[i].X() << endl;
            cout << "Y: " << PPACHit[i].Y() << endl;
        }
    }
};

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();
    void ReadXmls();
    std::vector<double> RecFPposition(std::vector<TVector2>,std::vector<int>);

    void InitOutputBranch();
    void InitInputBranch();
    void ReInitValue();
    static NPL::VAnalysis* Construct();

  private:
    double aoq;
    double zet;
    double beta;
    double delta;
    double fX;
    double fY;
    double fA;
    double fB;

    // intermediate variable


    // Branches and detectors
    TBigRIPSPPACPhysics* PPAC;
    TBigRIPSPlasticPhysics* PL;
    TBigRIPSICPhysics* IC;
    map<int,BigRIPSPPACHitList> FP_PPACHitList  ;

    //maps containings infos from input xml config files
    // for PPACs
    map<int,int>  PPAC_IDtoFP;//! Focal plane where the PPAC is located
    map<int,double> PPAC_XZoffset;
    map<int,double> PPAC_YZoffset;
    // for FocalPlanes
    map<int,int>  FPL_IDtoFP;//! Focal plane where the PPAC is located
    map<int,double> FPL_Z;
    map<int,double> FPL_Zoffset;

};




#endif
