#ifndef __SamuraiHodoscopeDATA__
#define __SamuraiHodoscopeDATA__
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : April 2021                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiHodoscope Raw data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TSamuraiHodoscopeData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // UP // 
    // Charge 
    vector<UShort_t>   fSamuraiHodoscope_Qu_ID;
    vector<Double_t>   fSamuraiHodoscope_Qu_Charge;
    
    // Time
    vector<UShort_t>   fSamuraiHodoscope_Tu_ID;
    vector<Double_t>   fSamuraiHodoscope_Tu_Time;
    
    // DOWN // 
    // Charge 
    vector<UShort_t>   fSamuraiHodoscope_Qd_ID;
    vector<Double_t>   fSamuraiHodoscope_Qd_Charge;
    
    // Time
    vector<UShort_t>   fSamuraiHodoscope_Td_ID;
    vector<Double_t>   fSamuraiHodoscope_Td_Time;
 

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSamuraiHodoscopeData();
    ~TSamuraiHodoscopeData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};
    void Dump() const;


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    // UP // 
    // Charge
    inline void SetChargeUp(const Double_t& ID, const Double_t& Charge){
    fSamuraiHodoscope_Qu_ID.push_back(ID);
    fSamuraiHodoscope_Qu_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeUp(const Double_t& ID, const Double_t& Time){
    fSamuraiHodoscope_Tu_ID.push_back(ID);
    fSamuraiHodoscope_Tu_Time.push_back(Time);
    };//!

    // DOWN // 
    // Charge
    inline void SetChargeDown(const Double_t& ID, const Double_t& Charge){
    fSamuraiHodoscope_Qd_ID.push_back(ID);
    fSamuraiHodoscope_Qd_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeDown(const Double_t& ID, const Double_t& Time){
    fSamuraiHodoscope_Td_ID.push_back(ID);
    fSamuraiHodoscope_Td_Time.push_back(Time);
    };//!
/*
    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy() const
      {return fSamuraiHodoscope_E_DetectorNbr.size();}
    inline UShort_t GetE_DetectorNbr(const unsigned int &i) const 
      {return fSamuraiHodoscope_E_DetectorNbr[i];}//!
    inline Double_t Get_Energy(const unsigned int &i) const 
      {return fSamuraiHodoscope_Energy[i];}//!

    // Time
    inline UShort_t GetMultTime() const
      {return fSamuraiHodoscope_T_DetectorNbr.size();}
    inline UShort_t GetT_DetectorNbr(const unsigned int &i) const 
      {return fSamuraiHodoscope_T_DetectorNbr[i];}//!
    inline Double_t Get_Time(const unsigned int &i) const 
      {return fSamuraiHodoscope_Time[i];}//!
*/

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSamuraiHodoscopeData,1)  // SamuraiHodoscopeData structure
};

#endif
