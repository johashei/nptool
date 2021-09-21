#ifndef __SAMURAIDATA__
#define __SAMURAIDATA__
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto  contact address: pilottoelia@gmail.com                        *
 *                                                                           *
 * Creation Date  : septembre 2021                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SAMURAI Raw data                                    *
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

class TSAMURAIData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // Energy
    vector<UShort_t>   fSAMURAI_E_DetectorNbr;
    vector<Double_t>   fSAMURAI_Energy;

    // Time
    vector<UShort_t>   fSAMURAI_T_DetectorNbr;
    vector<Double_t>   fSAMURAI_Time;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSAMURAIData();
    ~TSAMURAIData();
    

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
    // Energy
    inline void SetEnergy(const UShort_t& DetNbr,const Double_t& Energy){
      fSAMURAI_E_DetectorNbr.push_back(DetNbr);
      fSAMURAI_Energy.push_back(Energy);
    };//!

    // Time
    inline void SetTime(const UShort_t& DetNbr,const Double_t& Time)	{
      fSAMURAI_T_DetectorNbr.push_back(DetNbr);     
      fSAMURAI_Time.push_back(Time);
    };//!


    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy() const
      {return fSAMURAI_E_DetectorNbr.size();}
    inline UShort_t GetE_DetectorNbr(const unsigned int &i) const 
      {return fSAMURAI_E_DetectorNbr[i];}//!
    inline Double_t Get_Energy(const unsigned int &i) const 
      {return fSAMURAI_Energy[i];}//!

    // Time
    inline UShort_t GetMultTime() const
      {return fSAMURAI_T_DetectorNbr.size();}
    inline UShort_t GetT_DetectorNbr(const unsigned int &i) const 
      {return fSAMURAI_T_DetectorNbr[i];}//!
    inline Double_t Get_Time(const unsigned int &i) const 
      {return fSAMURAI_Time[i];}//!


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSAMURAIData,1)  // SAMURAIData structure
};

#endif
