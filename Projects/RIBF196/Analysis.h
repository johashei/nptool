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
    TMatrixD RecReadTransferMatrix(string);
    std::vector<double> RecFPposition(std::vector<TVector2>,std::vector<int>);
    double RecDeltaBrho(std::vector<double>,std::vector<double>, TMatrixD matrix);
    double RecAoqOneFold(double, double, double);
    double RecZ(double,double,double);

    void InitOutputBranch();
    void InitInputBranch();
    void ReInitValue();
    static NPL::VAnalysis* Construct();

  private:
    //calculated variables
    double aoq;
    double zet;
    double beta;
    double delta35;
    double delta57;
    double Brho35;
    double Brho57;
    double tof37;
    double aoq35;
    double aoq57;
    double aoq37;
    double beta35;
    double beta57;
    double beta37;
    double dE_ICF7;
    double z_BR;
    double tf3;
    double tf7;
    double fX;
    double fY;
    double fA;
    double fB;

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

    //double tof_offset37 = 0.;
    //double tof_offset37 = 300.06;
    double tof_offset37 = 299.968;

    const double clight = 299.7792458; // in mm/ns
    const double mnucleon = 931.49432;  // MeV


    // intermediate variables
    TMatrixD matF35;
    TMatrixD matF57;


    // Branches and detectors
    TBigRIPSPPACPhysics* PPAC;
    TBigRIPSPlasticPhysics* PL;
    TBigRIPSICPhysics* IC;
    map<int,BigRIPSPPACHitList> FP_PPACHitList  ;
    int EventNumber;
    int RunNumber;

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
    // for FocalPlanes
    map<int,int>  FP_IDtoFP;
    map<int,int>  FPtoID;
    map<int,double> FP_Z;
    map<int,double> FP_Zoffset;

};




#endif
