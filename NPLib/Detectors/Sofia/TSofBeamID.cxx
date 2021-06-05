/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 * Creation Date  : May 2021                                           *
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
#include "TSofBeamID.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofBeamID)


//////////////////////////////////////////////////////////////////////
TSofBeamID::TSofBeamID() {
}



//////////////////////////////////////////////////////////////////////
TSofBeamID::~TSofBeamID() {
}



//////////////////////////////////////////////////////////////////////
void TSofBeamID::Clear() {
  Zbeam = -1;
  Qmax  = -1;
  AoQ   = -1;
  Abeam = -1;
  Beta  = -1;
  Gamma = -1;
  Brho  = -1;
  XS2   = -1000;
  XCC   = -1000;

}



//////////////////////////////////////////////////////////////////////
void TSofBeamID::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofBeamID::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  cout << "Zbeam: " << Zbeam << endl;
  cout << "AoQ: " << AoQ << endl;
  cout << "Abeam: " << Abeam << endl;
  cout << "Beta: " << Beta << endl;
  cout << "Gamma: " << Gamma << endl;
  cout << "Brho: " << Brho << endl;
  cout << "XS2: " << XS2 << endl;
  cout << "XCC: " << XCC << endl;
 
}
