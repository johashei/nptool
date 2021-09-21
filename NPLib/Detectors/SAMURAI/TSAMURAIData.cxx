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
 *  This class hold SAMURAI Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSAMURAIData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSAMURAIData)


//////////////////////////////////////////////////////////////////////
TSAMURAIData::TSAMURAIData() {
}



//////////////////////////////////////////////////////////////////////
TSAMURAIData::~TSAMURAIData() {
}



//////////////////////////////////////////////////////////////////////
void TSAMURAIData::Clear() {
  // Energy
  fSAMURAI_E_DetectorNbr.clear();
  fSAMURAI_Energy.clear();
  // Time
  fSAMURAI_T_DetectorNbr.clear();
  fSAMURAI_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TSAMURAIData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSAMURAIData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fSAMURAI_E_DetectorNbr.size();
  cout << "SAMURAI_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fSAMURAI_E_DetectorNbr[i]
         << " Energy: " << fSAMURAI_Energy[i];
  }
  
  // Time
  mysize = fSAMURAI_T_DetectorNbr.size();
  cout << "SAMURAI_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fSAMURAI_T_DetectorNbr[i]
         << " Time: " << fSAMURAI_Time[i];
  }
}
