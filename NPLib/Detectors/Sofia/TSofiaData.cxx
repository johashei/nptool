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
#include "TSofiaData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofiaData)


//////////////////////////////////////////////////////////////////////
TSofiaData::TSofiaData() {
}



//////////////////////////////////////////////////////////////////////
TSofiaData::~TSofiaData() {
}



//////////////////////////////////////////////////////////////////////
void TSofiaData::Clear() {
  // TOF
  fTOF_DetectorNbr.clear();
  fTOF_PlasticNbr.clear();
  fTOF_Energy.clear();
  fTOF_Time.clear();
  
  // TWIN
  fTWIN_SectorNbr.clear();
  fTWIN_AnodeNbr.clear();
  fTWIN_AnodeEnergy.clear();
  fTWIN_AnodeTime.clear();
  
  fTWIN_Esum1 = -10;
  fTWIN_Esum2 = -10;
  fTWIN_Esum3 = -10;
  fTWIN_Esum4 = -10;
}



//////////////////////////////////////////////////////////////////////
void TSofiaData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofiaData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fTOF_DetectorNbr.size();
  cout << "TWIN_Mult: " << GetTwinMult() << endl;
  cout << "TOF_Mult: " << mysize << endl;
 
}
