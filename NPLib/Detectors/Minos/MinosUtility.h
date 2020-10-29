#ifndef MINOSUTILITY_H
#define MINOSUTILITY_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta   contact address: matta@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : october 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Minos class utility. Hold various usefull fonction to analyse MINOS      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include<vector>
#include "Math/Minimizer.h"
#include "Math/Functor.h"

namespace NPL{

  class MinosUtility{
    public:
      MinosUtility();
      ~MinosUtility();

    public:
      // take the vector describing the charge, perform a fit and return the charge and time of the trace
      double Calibrate(const std::vector<double>& T, const std::vector<double>& Q, double& time, double& charge); 
      // This function describe the shape of signal 
      double signal_shape(const double& x, const double* p);
      double signal_chi2(const double* p);

    private:
      double m_signal_offset; // offset apply to the signal before fitting
      ROOT::Math::Minimizer* m_signal_min; // minimiser used for the signal fitting
      ROOT::Math::Functor    m_signal_func; // functor passed to the miniser
      const std::vector<double>* m_fitSignalT;
      const std::vector<double>* m_fitSignalQ;
      unsigned int m_signal_size;
      double t0;
      double m_signal_parameter[4];
      double m_Dp1; // mean distance between t0 and T[0]
      double m_p2;// mean tau
      double m_Ap0;// mean ratio between q max and A
      unsigned int m_counter;
  };



}  

#endif
