#ifndef __VendetaDATA__
#define __VendetaDATA__
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Vendeta Raw data                                    *
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

class TVendetaData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<UShort_t>   fVendeta_DetectorNbr;
    vector<Double_t>   fVendeta_Q1;
    vector<Double_t>   fVendeta_Q2;
    vector<Double_t>   fVendeta_Time;
    vector<Bool_t>     fVendeta_isHG; 


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TVendetaData();
    ~TVendetaData();
    

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
    inline void SetDetectorNbr(const UShort_t& DetNbr) {fVendeta_DetectorNbr.push_back(DetNbr);};//!
    inline void SetQ1(const Double_t& Q1) {fVendeta_Q1.push_back(Q1);};//!
    inline void SetQ2(const Double_t& Q2) {fVendeta_Q2.push_back(Q2);};//!
    inline void SetTime(const Double_t& Time) {fVendeta_Time.push_back(Time);};//!
    inline void SetHighGainStatus(const Bool_t& isHG) {fVendeta_isHG.push_back(isHG);};//

    //////////////////////    GETTERS    ////////////////////////
    inline UShort_t GetMultEnergy() const {return fVendeta_DetectorNbr.size();}
    inline UShort_t GetDetectorNbr(const unsigned int &i) const {return fVendeta_DetectorNbr[i];}//!
    inline Double_t GetQ1(const unsigned int &i) const {return fVendeta_Q1[i];}//!
    inline Double_t GetQ2(const unsigned int &i) const {return fVendeta_Q2[i];}//!
    inline Double_t GetTime(const unsigned int &i) const {return fVendeta_Time[i];}//!
    inline Bool_t GetHighGainStatus(const unsigned int &i) const {return fVendeta_isHG[i];}//!

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TVendetaData,1)  // VendetaData structure
};

#endif
