#ifndef __FissionFragmentDATA__
#define __FissionFragmentDATA__
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
 *  This class hold FissionFragment Raw data                                    *
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

class TSofFissionFragment : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private:
    vector<double> fFF_Z; 
    vector<int> fFF_iZ;
    vector<double> fFF_AoQ;
    vector<double> fFF_A;
    vector<double> fFF_Beta;
    vector<double> fFF_TOF;
    vector<double> fFF_Gamma;
    vector<double> fFF_Brho;
    vector<double> fFF_DT;
    vector<double> fFF_ThetaIn;
    vector<double> fFF_TofPosX;
    vector<double> fFF_TofPosY;
    vector<double> fFF_PosX1;
    vector<double> fFF_PosX2;
    vector<double> fFF_PosX3;
    vector<double> fFF_PosY1;
    vector<double> fFF_PosY2;
    vector<double> fFF_PosY3;
    double fFF_Zsum;
    int fFF_IntZsum;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofFissionFragment();
    ~TSofFissionFragment();
    

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
    inline void SetZsum(double val){fFF_Zsum = val;};//!
    inline void SetIntZsum(int val){fFF_IntZsum = val;};//!
    inline void SetZ(double val){fFF_Z.push_back(val);};//!
    inline void SetiZ(int val){fFF_iZ.push_back(val);};//!
    inline void SetAoQ(double val){fFF_AoQ.push_back(val);};//!
    inline void SetA(double val){fFF_A.push_back(val);};//!
    inline void SetBeta(double val){fFF_Beta.push_back(val);};//!
    inline void SetTOF(double val){fFF_TOF.push_back(val);};//!
    inline void SetGamma(double val){fFF_Gamma.push_back(val);};//!
    inline void SetBrho(double val){fFF_Brho.push_back(val);};//!
    inline void SetDT(double val){fFF_DT.push_back(val);};//!
    inline void SetThetaIn(double val){fFF_ThetaIn.push_back(val);};//!
    inline void SetTofPosX(double val){fFF_TofPosX.push_back(val);};//!
    inline void SetTofPosY(double val){fFF_TofPosY.push_back(val);};//!
    inline void SetPosX1(double val){fFF_PosX1.push_back(val);};//!
    inline void SetPosX2(double val){fFF_PosX2.push_back(val);};//!
    inline void SetPosX3(double val){fFF_PosX3.push_back(val);};//!
    inline void SetPosY1(double val){fFF_PosY1.push_back(val);};//!
    inline void SetPosY2(double val){fFF_PosY2.push_back(val);};//!
    inline void SetPosY3(double val){fFF_PosY3.push_back(val);};//!


    //////////////////////    GETTERS    ////////////////////////
    int GetMult() {return fFF_Z.size();}//!
    int GetMultTofPos() {return fFF_TofPosY.size();}//!
    int GetMultMwpc1() {return fFF_PosY1.size();}//!
    int GetMultMwpc2() {return fFF_PosY2.size();}//!
    int GetMultMwpc3() {return fFF_PosY3.size();}//!
    inline double GetZsum() const {return fFF_Zsum;}//! 
    inline int GetIntZsum() const {return fFF_IntZsum;}//! 
    inline double GetZ(int i) const {return fFF_Z[i];}//! 
    inline int GetiZ(int i) const {return fFF_iZ[i];}//! 
    inline double GetAoQ(int i) const {return fFF_AoQ[i];}//! 
    inline double GetA(int i) const {return fFF_A[i];}//! 
    inline double GetBeta(int i) const {return fFF_Beta[i];}//! 
    inline double GetTOF(int i) const {return fFF_TOF[i];}//! 
    inline double GetGamma(int i) const {return fFF_Gamma[i];}//! 
    inline double GetBrho(int i) const {return fFF_Brho[i];}//! 
    inline double GetDT(int i) const {return fFF_DT[i];}//! 
    inline double GetThetaIn(int i) const {return fFF_ThetaIn[i];}//! 
    inline double GetTofPosX(int i) const {return fFF_TofPosX[i];}//! 
    inline double GetTofPosY(int i) const {return fFF_TofPosY[i];}//! 
    inline double GetPosX1(int i) const {return fFF_PosX1[i];}//! 
    inline double GetPosX2(int i) const {return fFF_PosX2[i];}//! 
    inline double GetPosX3(int i) const {return fFF_PosX3[i];}//! 
    inline double GetPosY1(int i) const {return fFF_PosY1[i];}//! 
    inline double GetPosY2(int i) const {return fFF_PosY2[i];}//! 
    inline double GetPosY3(int i) const {return fFF_PosY3[i];}//! 


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofFissionFragment,1)  // FissionFragment structure
};

#endif
