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
    // TOF //
    vector<int>      fTOF_DetectorNbr;
    vector<int>      fTOF_PlasticNbr;
    vector<double>   fTOF_Energy;
    vector<double>   fTOF_Time;
  
    // TWIN MUSIC //
    vector<int>      fTWIN_SectorNbr;
    vector<int>      fTWIN_AnodeNbr;
    vector<double>   fTWIN_AnodeEnergy;
    vector<double>   fTWIN_AnodeTime;
    double           fTWIN_Esum1;
    double           fTWIN_Esum2;
    double           fTWIN_Esum3;
    double           fTWIN_Esum4;

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
    // TOF
    inline void SetDetectorNbr(int det){fTOF_DetectorNbr.push_back(det);};//!
    inline void SetPlasticNbr(int plastic){fTOF_PlasticNbr.push_back(plastic);};//!
    inline void SetEnergy(double Energy){fTOF_Energy.push_back(Energy);};//!
    inline void SetTime(double Time){fTOF_Time.push_back(Time);};//!

    // TWIN
    inline void SetTwinSectorNbr(int Sector){fTWIN_SectorNbr.push_back(Sector);};//!
    inline void SetTwinAnodeNbr(int Anode){fTWIN_AnodeNbr.push_back(Anode);};//!
    inline void SetTwinAnodeEnergy(double Energy){fTWIN_AnodeEnergy.push_back(Energy);};//!
    inline void SetTwinAnodeTime(double Time){fTWIN_AnodeTime.push_back(Time);};//!
    inline void SetTwinEsum1(double E){fTWIN_Esum1=E;};//!
    inline void SetTwinEsum2(double E){fTWIN_Esum2=E;};//!
    inline void SetTwinEsum3(double E){fTWIN_Esum3=E;};//!
    inline void SetTwinEsum4(double E){fTWIN_Esum4=E;};//!

    //////////////////////    GETTERS    ////////////////////////
    // TOF
    inline int GetMultiplicity() const {return fTOF_PlasticNbr.size();}//!
    inline int GetDetectorNbr(const unsigned int &i) const {return fTOF_DetectorNbr[i];}//! 
    inline int GetPlasticNbr(const unsigned int &i) const {return fTOF_PlasticNbr[i];}//!     
    inline double GetEnergy(const unsigned int &i) const {return fTOF_Energy[i];}//!     
    inline double GetTime(const unsigned int &i) const {return fTOF_Time[i];}//!     

    // TWIN
    inline int GetTwinMult() const {return fTWIN_AnodeNbr.size();}//!
    inline int GetTwinSectorNbr(const unsigned int &i) const {return fTWIN_SectorNbr[i];}//!
    inline int GetTwinAnodeNbr(const unsigned int &i) const {return fTWIN_AnodeNbr[i];}//!
    inline double GetTwinAnodeEnergy(const unsigned int &i) const {return fTWIN_AnodeEnergy[i];}//!
    inline double GetTwinAnodeTime(const unsigned int &i) const {return fTWIN_AnodeTime[i];}//!
    inline double GetTwinEsum1() const {return fTWIN_Esum1;}//!
    inline double GetTwinEsum2() const {return fTWIN_Esum2;}//!
    inline double GetTwinEsum3() const {return fTWIN_Esum3;}//!
    inline double GetTwinEsum4() const {return fTWIN_Esum4;}//!

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofiaData,1)  // SofiaData structure
};

#endif
