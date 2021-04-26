#ifndef BIGRIPSPlasticVariables_H
#define BIGRIPSPlasticVariables_H

/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Freddy Flavigny  contact: flavigny@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class aims at holding Plastic data after pretreatment               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *  Intermediate class necessary to hold all variables per detector per event*
 *  Different from TPlasticData whose variable (vectors) are independent     *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
using namespace std;

class BigRIPSPlasticVariables{
  public:
   BigRIPSPlasticVariables(){Clear();};  
   ~BigRIPSPlasticVariables(){};  

  public:
    std::vector<double> FTL;
    std::vector<double> FTR;
    std::vector<double> FQL;
    std::vector<double> FQR;
    int FmultiHit[4];

    void Clear(){
        FTL.clear();
        FTR.clear();
        FQL.clear();
        FQR.clear();
        for(int i=0; i<4; i++) FmultiHit[i]=0;
    };

    void Print(){
        //cout << "XXXXXXXXXXXXXXXXXXXXXXXX Plastic Event XXXXXXXXXXXXXXXXX" << endl;
        cout << "FTL_Mult = " << FTL.size();
        for (UShort_t i = 0; i < FTL.size(); i++){cout << "\tFTL: " << FTL[i] << endl;}
        cout << "FTR_Mult = " << FTR.size();
        for (UShort_t i = 0; i < FTR.size(); i++){cout << "\tFTR: " << FTR[i] << endl;}
        cout << "FQL_Mult = " << FQL.size();
        for (UShort_t i = 0; i < FQL.size(); i++){cout << "\tFQL: " << FQL[i] << endl;}
        cout << "FQR_Mult = " << FQR.size();
        for (UShort_t i = 0; i < FQR.size(); i++){cout << "\tFQR: " << FQR[i] << endl;}
        cout << "MultHit = " <<endl;
        for (UShort_t i = 0; i <4; i++){cout << FmultiHit[i] << endl;}
    }

    bool HasTLandQL(){
        if(FTL.size()==1 && FQL.size()==1){return true;}
        else{return false;}
    }
    bool HasTRandQR(){
        if(FTR.size()==1 && FQR.size()==1){return true;}
        else{return false;}
    }
    bool HasTLandTR(){
        if(FTL.size()==1 && FTR.size()==1){return true;}
        else{return false;}
    }
    bool HasQLandQR(){
        if(FQL.size()==1 && FQR.size()==1){return true;}
        else{return false;}
    }
    bool HasEverything(){
        if(FTL.size()==1 && FTR.size()==1 &&
           FQL.size()==1 && FQR.size()==1){
            return true;
        }else{return false;}
    }
  };

#endif
