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
#include "TStrasseData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TStrasseData)


//////////////////////////////////////////////////////////////////////
TStrasseData::TStrasseData() {
}



//////////////////////////////////////////////////////////////////////
TStrasseData::~TStrasseData() {
}



//////////////////////////////////////////////////////////////////////
void TStrasseData::Clear() {
  // Energy X
  fInner_XE_DetectorNbr.clear();
  fInner_XE_StripNbr.clear();
  fInner_XE_Energy.clear();
  // Energy Y
  fInner_YE_DetectorNbr.clear();
  fInner_YE_StripNbr.clear();
  fInner_YE_Energy.clear();

  // Energy X
  fOuter_XE_DetectorNbr.clear();
  fOuter_XE_StripNbr.clear();
  fOuter_XE_Energy.clear();
  // Energy Y
  fOuter_YE_DetectorNbr.clear();
  fOuter_YE_StripNbr.clear();
  fOuter_YE_Energy.clear();

}



//////////////////////////////////////////////////////////////////////
void TStrasseData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TStrasseData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fInner_XE_DetectorNbr.size();
  cout << "Inner Strasse_XE_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "X-DetNbr: " << fInner_XE_DetectorNbr[i]
         << " X-Energy: " << fInner_XE_Energy[i];
  }
  
}
