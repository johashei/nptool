/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Vendeta Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TVendetaData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TVendetaData)


//////////////////////////////////////////////////////////////////////
TVendetaData::TVendetaData() {
}



//////////////////////////////////////////////////////////////////////
TVendetaData::~TVendetaData() {
}



//////////////////////////////////////////////////////////////////////
void TVendetaData::Clear() {
  // Energy
  fVendeta_E_DetectorNbr.clear();
  fVendeta_Energy.clear();
  // Time
  fVendeta_T_DetectorNbr.clear();
  fVendeta_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TVendetaData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TVendetaData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fVendeta_E_DetectorNbr.size();
  cout << "Vendeta_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fVendeta_E_DetectorNbr[i]
         << " Energy: " << fVendeta_Energy[i];
  }
  
  // Time
  mysize = fVendeta_T_DetectorNbr.size();
  cout << "Vendeta_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fVendeta_T_DetectorNbr[i]
         << " Time: " << fVendeta_Time[i];
  }
}
