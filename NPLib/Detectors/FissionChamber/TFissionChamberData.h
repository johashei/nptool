#ifndef __FissionChamberDATA__
#define __FissionChamberDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : September 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FissionChamber Raw data                                    *
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

class TFissionChamberData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<UShort_t>   fFC_AnodeNbr;
    vector<Double_t>   fFC_Energy;
    vector<Double_t>   fFC_Time;
    vector<Bool_t>     fFC_isFakeFission;
    vector<Double_t>   fFC_Time_HF;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TFissionChamberData();
    ~TFissionChamberData();
    

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
    inline void SetAnodeNbr(const UShort_t& AnodeNbr){fFC_AnodeNbr.push_back(AnodeNbr);}//!
    inline void SetEnergy(const Double_t& Energy){fFC_Energy.push_back(Energy);}//!
    inline void SetTime(const Double_t& Time){fFC_Time.push_back(Time);}//!
    inline void SetTimeHF(const Double_t& Time){fFC_Time_HF.push_back(Time);}//!
    inline void SetFakeFissionStatus(const Bool_t& isFF){fFC_isFakeFission.push_back(isFF);}//!

    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy() const
      {return fFC_Energy.size();}
    inline UShort_t GetAnodeNbr(const unsigned int &i) const 
      {return fFC_AnodeNbr[i];}//!
    inline Double_t GetEnergy(const unsigned int &i) const 
      {return fFC_Energy[i];}//!
    inline Double_t GetTime(const unsigned int &i) const 
      {return fFC_Time[i];}//!
    inline Double_t GetTimeHF(const unsigned int &i) const 
      {return fFC_Time_HF[i];}//!
    inline Bool_t GetFakeFissionStatus(const unsigned int &i) const 
      {return fFC_isFakeFission[i];}//!




  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TFissionChamberData,1)  // FissionChamberData structure
};

#endif
