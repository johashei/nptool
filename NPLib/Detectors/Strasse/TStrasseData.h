#ifndef __StrasseDATA__
#define __StrasseDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny    contact : flavigny@lpccaen.in2p3.fr       *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Strasse Raw data                                         *
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

class TStrasseData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // First Stage Front Energy
    vector<unsigned short>  fInner_XE_DetectorNbr;
    vector<unsigned short>  fInner_XE_StripNbr;
    vector<double>          fInner_XE_Energy;
    // First Stage Back Energy
    vector<unsigned short>  fInner_YE_DetectorNbr;
    vector<unsigned short>  fInner_YE_StripNbr;
    vector<double>          fInner_YE_Energy;

    // Second Stage Front Energy
    vector<unsigned short>  fOutter_XE_DetectorNbr;
    vector<unsigned short>  fOutter_XE_StripNbr;
    vector<double>          fOutter_XE_Energy;
    // Second Stage Back Energy
    vector<unsigned short>  fOutter_YE_DetectorNbr;
    vector<unsigned short>  fOutter_YE_StripNbr;
    vector<double>          fOutter_YE_Energy;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TStrasseData();
    ~TStrasseData();
    

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
    // First Stage Energy Front
    inline void SetInnerXE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fInner_XE_DetectorNbr.push_back(DetNbr);
      fInner_XE_StripNbr.push_back(StripNbr);
      fInner_XE_Energy.push_back(Energy);
    };//!
    // First Stage Energy Back
    inline void SetInnerYE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fInner_YE_DetectorNbr.push_back(DetNbr);
      fInner_YE_StripNbr.push_back(StripNbr);
      fInner_YE_Energy.push_back(Energy);
    };//!
   
    //////
    // Second Stage Energy Front
    inline void SetOutterXE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fOutter_XE_DetectorNbr.push_back(DetNbr);
      fOutter_XE_StripNbr.push_back(StripNbr);
      fOutter_XE_Energy.push_back(Energy);
    };//!
    // Second Stage Energy Back
    inline void SetOutterYE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fOutter_YE_DetectorNbr.push_back(DetNbr);
      fOutter_YE_StripNbr.push_back(StripNbr);
      fOutter_YE_Energy.push_back(Energy);
    };//!

    //////////////////////    GETTERS    ////////////////////////
    // First Stage Energy X
    inline UShort_t GetInnerMultXEnergy() const
      {return fInner_XE_DetectorNbr.size();}
    inline UShort_t GetInner_XE_DetectorNbr(const unsigned int &i) const 
      {return fInner_XE_DetectorNbr[i];}//!
    inline UShort_t GetInner_XE_StripNbr(const unsigned int &i) const 
      {return fInner_XE_StripNbr[i];}//!
    inline Double_t GetInner_XE_Energy(const unsigned int &i) const 
      {return fInner_XE_Energy[i];}//!
    // First Stage Energy Y
    inline UShort_t GetInnerMultYEnergy() const
      {return fInner_YE_DetectorNbr.size();}
    inline UShort_t GetInner_YE_DetectorNbr(const unsigned int &i) const 
      {return fInner_YE_DetectorNbr[i];}//!
    inline UShort_t GetInner_YE_StripNbr(const unsigned int &i) const 
      {return fInner_YE_StripNbr[i];}//!
    inline Double_t GetInner_YE_Energy(const unsigned int &i) const 
      {return fInner_YE_Energy[i];}//!
   
    //////
    // Second Stage Energy X
    inline UShort_t GetOutterMultXEnergy() const
      {return fOutter_XE_DetectorNbr.size();}
    inline UShort_t GetOutter_XE_DetectorNbr(const unsigned int &i) const 
      {return fOutter_XE_DetectorNbr[i];}//!
    inline UShort_t GetOutter_XE_StripNbr(const unsigned int &i) const 
      {return fOutter_XE_StripNbr[i];}//!
    inline Double_t GetOutter_XE_Energy(const unsigned int &i) const 
      {return fOutter_XE_Energy[i];}//!
    // Second Stage Energy Y
    inline UShort_t GetOutterMultYEnergy() const
      {return fOutter_YE_DetectorNbr.size();}
    inline UShort_t GetOutter_YE_DetectorNbr(const unsigned int &i) const 
      {return fOutter_YE_DetectorNbr[i];}//!
    inline UShort_t GetOutter_YE_StripNbr(const unsigned int &i) const 
      {return fOutter_YE_StripNbr[i];}//!
    inline Double_t GetOutter_YE_Energy(const unsigned int &i) const 
      {return fOutter_YE_Energy[i];}//!

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TStrasseData,1)  // StrasseData structure
};

#endif
