#ifndef __SofBeamIDDATA__
#define __SofBeamIDDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SofBeamID Raw data                                    *
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

class TSofBeamID : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private:
    double Zbeam; 
    double Qmax;
    double AoQ;
    double Abeam;
    double Beta;
    double Gamma;
    double Brho;
    double XS2;
    double XCC;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofBeamID();
    ~TSofBeamID();
    

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
    inline void SetZbeam(double val){Zbeam = val;};//!
    inline void SetQmax(double val){Qmax = val;};//!
    inline void SetAoQ(double val){AoQ = val;};//!
    inline void SetAbeam(double val){Abeam = val;};//!
    inline void SetBeta(double val){Beta = val;};//!
    inline void SetGamma(double val){Gamma = val;};//!
    inline void SetBrho(double val){Brho = val;};//!
    inline void SetXS2(double val){XS2 = val;};//!
    inline void SetXCC(double val){XCC = val;};//!

    //////////////////////    GETTERS    ////////////////////////
    inline double GetZbeam() const {return Zbeam;}//! 
    inline double GetQmax() const {return Qmax;}//! 
    inline double GetAoQ() const {return AoQ;}//! 
    inline double GetAbeam() const {return Abeam;}//! 
    inline double GetBeta() const {return Beta;}//! 
    inline double GetGamma() const {return Gamma;}//! 
    inline double GetBrho() const {return Brho;}//! 
    inline double GetXS2() const {return XS2;}//! 
    inline double GetXCC() const {return XCC;}//! 

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofBeamID,1)  // SofBeamID structure
};

#endif
