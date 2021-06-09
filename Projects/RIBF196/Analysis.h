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
 * AoQ,Z reconstruction for RIBF196 experiment                               *
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
    double aoqc37;
    double aoqc811;
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
    TBigRIPSReco* Rec37;
    TBigRIPSReco* Rec89;
    TBigRIPSReco* Rec911;
    TBigRIPSReco* Rec811;
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

    // Higher order Optical corrections
    
    double ZD_OptCorr_Ga(double *);
    constexpr static int    gGaZDNVariables    = 8;
    constexpr static int    gGaZDNCoefficients = 10;
    constexpr static double gGaZDDMean         = 2.67982;
    // Assignment to mean vector.
    constexpr static double gGaZDXMean[] = {
        -0.151618, -0.97249, -0.226697, -0.0028457, -0.0250826, 0.00288959, -0.0263575, 0.180498 };
    // Assignment to minimum vector.
    constexpr static double gGaZDXMin[] = {
        -112.762, -67.3691, -21.6583, -1461.26, -19.009, -25.9623, -63.8841, -82.8509 };

    // Assignment to maximum vector.
    constexpr static double gGaZDXMax[] = {
        110.757, 6247.9, 37.3008, 1461.03, 19.6172, 25.5737, 55.9905, 60.4835 };

    // Assignment to coefficients vector.
    constexpr static double gGaZDCoefficient[] = {
        0.0045763,
        0.0260008,
        0.0202036,
        -0.00373027,
        0.0117822,
        -0.00717916,
        -0.010714,
        -6.11157e-06,
        0.0045109,
        -0.00397418
    };

    // Assignment to error coefficients vector.
    constexpr static double gGaZDCoefficientRMS[] = {
        0.0263103,
        0.0627114,
        0.0645658,
        0.0289265,
        0.0959883,
        0.0450067,
        0.125588,
        0.0482148,
        0.0962523,
        0.0798097
    };

    // Assignment to powers vector.
    // The powers are stored row-wise, that is
    //  p_ij = gPower[i * NVariables + j];
    constexpr static int    gGaZDPower[] = {
        1,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  1,  2,  1,
        1,  1,  1,  1,  2,  1,  1,  1,
        1,  1,  1,  1,  1,  2,  1,  1,
        1,  1,  2,  1,  1,  1,  1,  1,
        3,  1,  1,  1,  1,  1,  1,  1,
        2,  1,  1,  1,  1,  1,  3,  1,
        1,  1,  1,  1,  1,  1,  1,  2,
        1,  1,  1,  1,  1,  1,  1,  3,
        1,  1,  1,  1,  1,  3,  1,  1
    };
    

   double BR_OptCorr_Ga(double *);
   constexpr static int    gGaBRNVariables    = 8;
   constexpr static int    gGaBRNCoefficients = 25;
   constexpr static double gGaBRDMean         = 2.67489;
    // Assignment to mean vector.
   constexpr static double gGaBRXMean[] = {
        -0.0642977, -0.00307703, -0.00100196, -0.0107261, -0.0131277, -0.00317987, -0.158857, -0.0435761 };

    // Assignment to minimum vector.
   constexpr static double gGaBRXMin[] = {
        -107.29, -25.9837, -31.9689, -42.8947, -29.3661, -31.5436, -20.0646, -37.375 };

    // Assignment to maximum vector.
   constexpr static double gGaBRXMax[] = {
        64.5168, 21.8811, 30.0843, 43.7818, 35.0687, 23.5787, 29.204, 35.7001 };

    // Assignment to coefficients vector.
   constexpr static double gGaBRCoefficient[] = {
        0.000721757,
        0.00503421,
        -0.00479955,
        0.00219026,
        0.00552134,
        0.00247818,
        0.0026914,
        -0.000758646,
        0.00743853,
        0.0162762,
        0.0157268,
        0.00235257,
        -0.0029327,
        -0.0115376,
        -0.0139695,
        -0.0526337,
        0.000312724,
        -0.00270404,
        0.00261865,
        -0.00281391,
        0.00284842,
        -0.00417684,
        0.029841,
        -0.0494852,
        0.00822339
    };

    // Assignment to error coefficients vector.
   constexpr static double gGaBRCoefficientRMS[] = {
        0.0221415,
        0.0946943,
        0.195954,
        0.0827622,
        0.10566,
        0.0550079,
        0.091793,
        0.0767151,
        0.172153,
        0.94549,
        0.578329,
        0.49222,
        0.16515,
        0.50208,
        0.489808,
        1.89549,
        0.0284801,
        0.100958,
        0.121855,
        0.161108,
        0.277453,
        0.292148,
        4.02157,
        3.97904,
        0.647432
    };

    // Assignment to powers vector.
    // The powers are stored row-wise, that is
    //  p_ij = gPower[i * NVariables + j];
   constexpr static int    gGaBRPower[] = {
        1,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  1,  2,  1,
        1,  1,  1,  1,  1,  2,  1,  2,
        1,  1,  1,  1,  1,  1,  3,  1,
        2,  1,  1,  1,  1,  1,  2,  1,
        2,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  2,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  2,  2,  1,
        1,  1,  2,  1,  2,  1,  1,  1,
        1,  1,  1,  1,  2,  1,  1,  3,
        2,  1,  1,  1,  2,  1,  2,  1,
        1,  2,  1,  2,  1,  1,  2,  1,
        4,  1,  1,  1,  1,  1,  1,  1,
        2,  1,  1,  1,  3,  1,  1,  1,
        1,  1,  1,  1,  3,  1,  3,  1,
        1,  2,  4,  1,  3,  2,  1,  1,
        1,  1,  1,  1,  1,  2,  1,  1,
        3,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  3,  1,  1,
        2,  1,  1,  1,  2,  1,  1,  1,
        1,  2,  1,  1,  1,  2,  2,  1,
        3,  1,  1,  1,  2,  1,  1,  1,
        1,  1,  1,  1,  2,  1,  2,  3,
        1,  1,  2,  1,  2,  1,  1,  3,
        3,  1,  1,  1,  1,  1,  3,  1
    };

};




#endif
