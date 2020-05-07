/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TPISTAData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TPISTAData)


//////////////////////////////////////////////////////////////////////
TPISTAData::TPISTAData() {
}



//////////////////////////////////////////////////////////////////////
TPISTAData::~TPISTAData() {
}



//////////////////////////////////////////////////////////////////////
void TPISTAData::Clear() {
  // Energy
  fPISTA_E_DetectorNbr.clear();
  fPISTA_Energy.clear();
  // Time
  fPISTA_T_DetectorNbr.clear();
  fPISTA_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TPISTAData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TPISTAData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fPISTA_E_DetectorNbr.size();
  cout << "PISTA_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fPISTA_E_DetectorNbr[i]
         << " Energy: " << fPISTA_Energy[i];
  }
  
  // Time
  mysize = fPISTA_T_DetectorNbr.size();
  cout << "PISTA_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fPISTA_T_DetectorNbr[i]
         << " Time: " << fPISTA_Time[i];
  }
}
