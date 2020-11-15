#ifndef __SofiaDATA__
#define __SofiaDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : November 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Sofia Raw data                                    *
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

class TSofiaData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int>      fTOF_DetectorNbr;
    vector<int>      fTOF_PlasticNbr;
    vector<double>   fTOF_Energy;
    vector<double>   fTOF_Time;



  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofiaData();
    ~TSofiaData();
    

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
    inline void SetDetectorNbr(int det){fTOF_DetectorNbr.push_back(det);};//!
    inline void SetPlasticNbr(int plastic){fTOF_PlasticNbr.push_back(plastic);};//!
    inline void SetEnergy(double Energy){fTOF_Energy.push_back(Energy);};//!
    inline void SetTime(double Time){fTOF_Time.push_back(Time);};//!



    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline int GetMultiplicity() const {return fTOF_PlasticNbr.size();}
    inline int GetDetectorNbr(const unsigned int &i) const {return fTOF_DetectorNbr[i];}//! 
    inline int GetPlasticNbr(const unsigned int &i) const {return fTOF_PlasticNbr[i];}//!     
    inline double GetEnergy(const unsigned int &i) const {return fTOF_Energy[i];}//!     
    inline double GetTime(const unsigned int &i) const {return fTOF_Time[i];}//!     


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofiaData,1)  // SofiaData structure
};

#endif
