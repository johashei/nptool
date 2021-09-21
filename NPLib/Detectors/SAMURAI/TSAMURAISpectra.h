#ifndef TSAMURAISPECTRA_H
#define TSAMURAISPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto  contact address: pilottoelia@gmail.com                        *
 *                                                                           *
 * Creation Date  : septembre 2021                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SAMURAI Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TSAMURAIData.h"
#include "TSAMURAIPhysics.h"

// Forward Declaration
class TSAMURAIPhysics;


class TSAMURAISpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TSAMURAISpectra();
    TSAMURAISpectra(unsigned int NumberOfDetectors);
    ~TSAMURAISpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TSAMURAIData*);
    void FillPreTreatedSpectra(TSAMURAIData*);
    void FillPhysicsSpectra(TSAMURAIPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
