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
 *  This class hold Sofia Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TSofiaSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TSofiaSpectra::TSofiaSpectra() 
   : fNumberOfDetectors(0) {
  SetName("Sofia");
}



////////////////////////////////////////////////////////////////////////////////
TSofiaSpectra::TSofiaSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TSofiaSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("Sofia");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TSofiaSpectra::~TSofiaSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TSofiaSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Sofia"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Sofia/RAW");
    // Time 
    name = "Sofia"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Sofia/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TSofiaSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Sofia"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Sofia/CAL");
    // Time
    name = "Sofia"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Sofia/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TSofiaSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "Sofia_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "Sofia/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TSofiaSpectra::FillRawSpectra(TSofiaData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultiplicity();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Sofia"+NPL::itoa(RawData->GetPlasticNbr(i))+"_ENERGY_RAW";
    family = "Sofia/RAW";

    FillSpectra(family,name,RawData->GetEnergy(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TSofiaSpectra::FillPreTreatedSpectra(TSofiaData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultiplicity();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Sofia"+NPL::itoa(PreTreatedData->GetPlasticNbr(i))+"_ENERGY_CAL";
    family = "Sofia/CAL";

    FillSpectra(family,name,PreTreatedData->GetEnergy(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TSofiaSpectra::FillPhysicsSpectra(TSofiaPhysics* Physics) {
  static string name;
  static string family;
  family= "Sofia/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "Sofia_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

