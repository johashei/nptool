#ifndef NPDCRECONSTRUCTION_H
#define NPDCRECONSTRUCTION_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta   contact address: matta@lpcaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : October 2020                                            *
 * Last update     : October 2020                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class have all the method needed to analyse Drift Chambers          *
 *****************************************************************************/
#include<vector>
#include"TVector3.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
namespace NPL{
  
  class DCReconstruction{
    public:
      DCReconstruction();
      ~DCReconstruction();
    
    public:
    // Build a track in 2D based on drift circle of Radius R and position X,Z
    // return X0(X100) the X position at Z=0 (Z=100)
    // return a and b the coeff of the 2D line
    void BuildTrack2D(const std::vector<double>& X,const std::vector<double>& Z,const std::vector<double>& R,double& X0,double& X100,double& a, double& b );
    
    // Compute X and Y crossing coordinate of 2 plane of Wire
    void ResolvePlane(const TVector3& L,const double& ThetaU ,const TVector3& H, const double& ThetaV, TVector3& PosXY);
    
    // Function used by the minimizer in BuildTrack2D
    double SumD(const double* parameter );

    private: // private member used by SumD
      ROOT::Math::Minimizer* m_min;
      ROOT::Math::Functor    m_func;
      const std::vector<double>* fitX;
      const std::vector<double>* fitZ;
      const std::vector<double>* fitR;
    };
  }

#endif
