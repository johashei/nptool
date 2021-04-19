#ifndef BIGRIPSPPACVariables_H
#define BIGRIPSPPACVariables_H

/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiFDC0 treated data                                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *  little class to index each of the DC wire                                *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
using namespace std;

class BigRIPSPPACVariables{
  public:
   BigRIPSPPACVariables(){Clear();};  
   ~BigRIPSPPACVariables(){};  

  public:
    std::vector<double> FTX1;
    std::vector<double> FTX2;
    std::vector<double> FTY1;
    std::vector<double> FTY2;
    std::vector<double> FTA;
    int FmultiHit[5];

    void Clear(){
        FTX1.clear();
        FTX2.clear();
        FTY1.clear();
        FTY2.clear();
        FTA.clear();
        for(int i=0; i<5; i++) FmultiHit[i]=0;
    };

    void Print(){
        //cout << "XXXXXXXXXXXXXXXXXXXXXXXX PPAC Event XXXXXXXXXXXXXXXXX" << endl;
        cout << "FTX1_Mult = " << FTX1.size();
        for (UShort_t i = 0; i < FTX1.size(); i++){cout << "\tFTX1: " << FTX1[i] << endl;}
        cout << "FTX2_Mult = " << FTX2.size();
        for (UShort_t i = 0; i < FTX2.size(); i++){cout << "\tFTX2: " << FTX2[i] << endl;}
        cout << "FTY1_Mult = " << FTY1.size();
        for (UShort_t i = 0; i < FTY1.size(); i++){cout << "\tFTY1: " << FTY1[i] << endl;}
        cout << "FTY2_Mult = " << FTY2.size();
        for (UShort_t i = 0; i < FTY2.size(); i++){cout << "\tFTY2: " << FTY2[i] << endl;}
        cout << "FTA_Mult = " << FTA.size();
        for (UShort_t i = 0; i < FTA.size(); i++){cout << "\tFTA: " << FTA[i] << endl;}
        cout << "MultHit = " <<endl;
        for (UShort_t i = 0; i <5; i++){cout << FmultiHit[i] << endl;}
    }

    bool HasTXs(){
        if(FTX1.size()==1 && FTX2.size()==1){return true;}
        else{return false;}
    }
    bool HasTYs(){
        if(FTY1.size()==1 && FTY2.size()==1){return true;}
        else{return false;}
    }
    bool HasTA(){
        if(FTA.size()==1){return true;}
        else{return false;}
    }
    bool HasEverything(){
        if(FTX1.size()==1 && FTX2.size()==1 &&
           FTY1.size()==1 && FTY2.size()==1 &&
           FTA.size()==1){
            return true;
        }else{return false;}
    }
/*
  public:
    void SetTX1(double value){TX1=value;}
    void SetTX2(double value){TX2=value;}
    void SetTY1(double value){TY1=value;}
    void SetTY2(double value){TY2=value;}
    void SetTA(double value){TA=value;}


    int GetTX1Mult() const {return TX1.size();}
    int const GetTX1(const unsigned int& i){return fPPAC_TX1[i];};

    int GetTX2Mult() const {return fPPAC_TX2_ID.size();}
    int const GetTX2(const unsigned int& i){return fPPAC_TX2[i];};

    int GetTY1Mult() const {return fPPAC_TY1_ID.size();}
    int const GetTY1(const unsigned int& i){return fPPAC_TY1[i];};

    int GetTY2Mult() const {return fPPAC_TY2_ID.size();}
    int const GetTY2(const unsigned int& i){return fPPAC_TY2[i];};

    int GetTAMult() const {return fPPAC_TA_ID.size();}
    int const GetTA(const unsigned int& i){return fPPAC_TA[i];};


    double GetTX1(){return TX1;}
    double GetTX2(){return TX2;}
    double GetTY1(){return TY1;}
    double GetTY2(){return TY2;}
    double GetTA(){return TA;}
 */  
  };

#endif
