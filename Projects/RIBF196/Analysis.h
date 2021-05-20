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
#include "TBigRIPSFocal.h"
#include "TBigRIPSReco.h"
#include <TMatrixD.h>
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
    //TMatrixD* RecReadTransferMatrix(char*);
    //TMatrixD RecReadTransferMatrix(string);
    std::vector<std::vector<double>> RecReadTransferMatrix2(string);
    std::vector<double> RecFPposition(std::vector<TVector2>,std::vector<int>);
    //double RecDeltaBrho(std::vector<double>,std::vector<double>, TMatrixD matrix);
    //double RecAoqOneFold(double, double, double);
    double RecZ(double,double,double, bool);

    void InitOutputBranch();
    void InitInputBranch();
    void ReInitValue();
    static NPL::VAnalysis* Construct();

  private:
    //calculated variables
    double tof37;
    double tof811;
    double dE_ICF7;
    double dE_ICF11;
    double z_BR;
    double z_ZD;
    double tf3;
    double tf7;
    double tf8;
    double tf11;

    // List of Dipole Brho, manually set for the moment 
    // But method to read from root file should be implemented
    double BrhoD1 = 7.9677;
    double BrhoD2 = 7.0945;
    double BrhoD3 = 7.0658;
    double BrhoD4 = 7.0658;
    double BrhoD5 = 6.8037;
    double BrhoD6 = 6.8037;
    double BrhoD7 = 5.9429;
    double BrhoD8 = 5.9381;

    double length35;
    double length57;
    double length37;
    double length89;
    double length911;
    double length811;

    double tof_offset37 = 300.06;
    double tof_offset811 = -134.503;

    const double clight = 299.7792458; // in mm/ns
    const double mnucleon = 931.49432;  // MeV

    //5 columns/ 6 rows transfer matrices
    vector<vector<double>> matrixF35;
    vector<vector<double>> matrixF57;
    vector<vector<double>> matrixF89;
    vector<vector<double>> matrixF911;


    // Branches and detectors
    TBigRIPSPPACPhysics* PPAC;
    TBigRIPSPlasticPhysics* PL;
    TBigRIPSICPhysics* IC;
    map<int,BigRIPSPPACHitList> FP_PPACHitList  ;
    int EventNumber;
    int RunNumber;
    int Trigger;
    unsigned long long int TimeStamp;
    //Focal* FP;
    //TBigRIPSFocal* FP;
    TBigRIPSReco* Rec35;
    TBigRIPSReco* Rec57;
    TBigRIPSReco* Rec89;
    TBigRIPSReco* Rec911;
    std::vector<std::vector<double>> RecFP;

    //maps containings infos from input xml config files
    // for PPACs
    map<int,int>  PPAC_IDtoFP;//! Focal plane where the PPAC is located
    map<int,double> PPAC_XZoffset;
    map<int,double> PPAC_YZoffset;
    // for Plastics
    map<int,int>  PL_IDtoFP;
    map<int, string>  PL_IDtoName; 
    map<string,int>  PL_NametoID; 
    map<string,double>  PL_NametoZ; 
    // for IC
    map<int,int>  IC_IDtoFP;
    map<int, string>  IC_IDtoName; 
    map<string,int>  IC_NametoID; 
    map<string,vector<double>>  IC_Zcoef;
    map<string,double>  IC_Ionpair;
    // for FocalPlanes
    map<int,int>  FP_IDtoFP;
    map<int,int>  FPtoID;
    map<int,double> FP_Z;
    map<int,double> FP_Zoffset;
    int NmaxFP;

};




#endif
