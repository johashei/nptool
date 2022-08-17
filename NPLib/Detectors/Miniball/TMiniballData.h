#ifndef __MiniballDATA__
#define __MiniballDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Modified by: Johannes Sorby Heines                                        * 
 *         contact address: j.s.heines@protonmailch                          *
 *                                                                           *
 * Creation Date  : January 2016                                             *
 * Last update    : August 2022                                              *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Miniball Raw data                                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
#include <iostream>
using namespace std;

// ROOT
#include "TObject.h"

class TMiniballData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // Energy
    vector<UShort_t>   fMiniball_E_DetectorNbr;
    vector<Double_t>   fMiniball_Energy;	

    // Angle
    vector<UShort_t>   fMiniball_A_DetectorNbr;
    vector<Double_t>   fMiniball_Angle;

    // Time
    vector<UShort_t>   fMiniball_T_DetectorNbr;
    vector<Double_t>   fMiniball_Time;

    vector<UShort_t>   fMiniball_E_CrystalNbr;
    vector<UShort_t>   fMiniball_A_CrystalNbr;
    vector<UShort_t>   fMiniball_T_CrystalNbr;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TMiniballData();
    virtual ~TMiniballData();
    

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
    inline void SetE_DetectorNbr(const UShort_t& DetNbr)
      {fMiniball_E_DetectorNbr.push_back(DetNbr);} //!
    inline void SetE_CrystalNbr(const UShort_t& CrystalNbr)
      {fMiniball_E_CrystalNbr.push_back(CrystalNbr);} //!
    inline void Set_Energy(const Double_t& Energy)
      {fMiniball_Energy.push_back(Energy);}//!
    // Prefer global setter so that all vectors have the same size
    inline void SetEnergy(const UShort_t& DetNbr,const UShort_t& CrystalNbr,const Double_t& Energy) {
      SetE_DetectorNbr(DetNbr);
      SetE_CrystalNbr(CrystalNbr);
      Set_Energy(Energy);
//      std::cout << "\nValues added to Miniball:\nfMiniball_E_DetectorNbr : " << fMiniball_E_DetectorNbr.back() << "\nfMiniball_E_CrystalNbr : " << fMiniball_E_CrystalNbr.back() << "\nfMiniball_Energy : " << fMiniball_Energy.back() << "\n" << std::endl;
    };//!

    // Angle
    inline void SetA_DetectorNbr(const UShort_t& DetNbr)
      {fMiniball_A_DetectorNbr.push_back(DetNbr);} //!
    inline void SetA_CrystalNbr(const UShort_t& CrystalNbr)
      {fMiniball_A_CrystalNbr.push_back(CrystalNbr);} //!
    inline void Set_Angle(const Double_t& Angle)
      {fMiniball_Angle.push_back(Angle);}//!
    // Prefer global setter so that all vectors have the same size
    inline void SetAngle(const UShort_t& DetNbr,const UShort_t& CrystalNbr,const Double_t& Angle) {
      SetA_DetectorNbr(DetNbr);
      SetA_CrystalNbr(CrystalNbr);
      Set_Angle(Angle);
//      std::cout << "\nValues added to Miniball:\nfMiniball_A_DetectorNbr : " << fMiniball_A_DetectorNbr.back() << "\nfMiniball_A_CrystalNbr : " << fMiniball_A_CrystalNbr.back() << "\nfMiniball_Angle : " << fMiniball_Angle.back() << "\n" << std::endl;
    };//!


    // Time
    inline void SetT_DetectorNbr(const UShort_t& DetNbr)
      {fMiniball_T_DetectorNbr.push_back(DetNbr);} //!
    inline void SetT_CrystalNbr(const UShort_t& CrystalNbr)
      {fMiniball_T_CrystalNbr.push_back(CrystalNbr);} //!
    inline void Set_Time(const Double_t& Time)
      {fMiniball_Time.push_back(Time);}//!
    // Prefer global setter so that all vectors have the same size
    inline void SetTime(const UShort_t& DetNbr,const UShort_t& CrystalNbr,const Double_t& Time)	{
      SetT_DetectorNbr(DetNbr);
      SetT_CrystalNbr(CrystalNbr);
      Set_Time(Time);
//      std::cout << "\nValues added to Miniball:\nfMiniball_T_DetectorNbr : " << fMiniball_T_DetectorNbr.back() << "\nfMiniball_T_CrystalNbr : " << fMiniball_T_CrystalNbr.back() << "\nfMiniball_Time : " << fMiniball_Time.back() << "\n" << std::endl;
    };//!


    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy() const
      {return fMiniball_E_DetectorNbr.size();}
    inline UShort_t GetE_DetectorNbr(const unsigned int &i) const 
      {return fMiniball_E_DetectorNbr[i];}//!
    inline UShort_t GetE_CrystalNbr(const unsigned int &i) const 
      {return fMiniball_E_CrystalNbr[i];}//!
    inline Double_t Get_Energy(const unsigned int &i) const 
      {return fMiniball_Energy[i];}//!

    // Angle
    inline UShort_t GetMultAngle() const
      {return fMiniball_A_DetectorNbr.size();}
    inline UShort_t GetA_DetectorNbr(const unsigned int &i) const 
      {return fMiniball_A_DetectorNbr[i];}//!
    inline UShort_t GetA_CrystalNbr(const unsigned int &i) const 
      {return fMiniball_A_DetectorNbr[i];}//!
    inline Double_t Get_Angle(const unsigned int &i) const 
      {return fMiniball_Angle[i];}//!


    // Time
    inline UShort_t GetMultTime() const
      {return fMiniball_T_DetectorNbr.size();}
    inline UShort_t GetT_DetectorNbr(const unsigned int &i) const 
      {return fMiniball_T_DetectorNbr[i];}//!
    inline UShort_t GetT_CrystalNbr(const unsigned int &i) const 
      {return fMiniball_T_DetectorNbr[i];}//!
    inline Double_t Get_Time(const unsigned int &i) const 
      {return fMiniball_Time[i];}//!


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TMiniballData,1)  // MiniballData structure
};

#endif
