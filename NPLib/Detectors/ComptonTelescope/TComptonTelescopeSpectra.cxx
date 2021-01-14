/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                  B. Le Crom                        lecrom@ipno.in2p3.fr   *
 * Creation Date  : April 2014                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for ComptonTelescope      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + first version (not complete yet)                                     *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// NPL
#include "TComptonTelescopeSpectra.h"
#include "NPOptionManager.h"
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

// STL
#include <stdexcept>
#include <iostream>  
#include <cstdlib>
#include <string>
using namespace std;


////////////////////////////////////////////////////////////////////////////////
TComptonTelescopeSpectra::TComptonTelescopeSpectra(){
  SetName("ComptonTelescope");
  fNumberOfTelescope = 0;
  fNumberOfDetectors = 1;
  fNumberOfStripsFront=32;
  fNumberOfStripsBack=32;
  fNumberOfCounters=50;
  fCalorimeterNPixels=64;
}



////////////////////////////////////////////////////////////////////////////////
TComptonTelescopeSpectra::TComptonTelescopeSpectra(unsigned int NumberOfTelescope)
{
   if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
      cout << "************************************************" << endl
         << "TComptonTelescopeSpectra : Initalising control spectra for " 
         << NumberOfTelescope << " Telescopes" << endl
         << "************************************************" << endl ;

   SetName("ComptonTelescope");
   fNumberOfTelescope = NumberOfTelescope;
   fNumberOfDetectors = 1;
   fNumberOfStripsFront=32;
   fNumberOfStripsBack=32;
   fNumberOfCounters=50;
   fCalorimeterNPixels=64;

   InitRawSpectra();
   InitPreTreatedSpectra();
   InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TComptonTelescopeSpectra::~TComptonTelescopeSpectra()
{
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::InitRawSpectra()
{
  string name;
  int ntot = (fNumberOfStripsFront + fNumberOfStripsBack);

  for (unsigned int i = 0; i < fNumberOfTelescope; i++) { // loop on number of telescope

    // DSSSD
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {

     // FRONT_E_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_E_RAW";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 512, 0, 1024, "COMPTONTELESCOPE/RAW/FRONTE");

      // BACK_E_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_E_RAW";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 512, 0, 1024, "COMPTONTELESCOPE/RAW/BACKE");

      // E Front and Back vs Strip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_BACK_E_RAW";
      AddHisto2D(name, name, ntot, 0, ntot, 512, 0, 1024, "COMPTONTELESCOPE/RAW/FRONT_BACK_E");
 
      // FRONT_T_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_T_RAW";
      AddHisto1D(name, name, 10000, 0, 1e10, "COMPTONTELESCOPE/RAW/FRONTT");

      // BACK_T_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_T_RAW";
      AddHisto1D(name, name, 10000, 0, 1e10, "COMPTONTELESCOPE/RAW/BACKT");

      // FRONT_RAW_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/RAW/MULT");

      // BACK_RAW_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/RAW/MULT");
    }

    // CALORIMETER
    name = "CT"+NPL::itoa(i+1)+"_CALOR_RAW_TRIGGER";
    AddHisto1D(name, name, fCalorimeterNPixels, 1, fCalorimeterNPixels+1, "COMPTONTELESCOPE/RAW/CALORTRIGGER");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::InitPreTreatedSpectra()
{
  string name;
  int ntot = (fNumberOfStripsFront + fNumberOfStripsBack);

  for (unsigned int i = 0; i < fNumberOfTelescope; i++) { // loop on number of telescope

    // DSSSD
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {

      // FRONT_E_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_E_CAL";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 1400, 0, 1.4, "COMPTONTELESCOPE/CAL/FRONTE");

      // BACK_E_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_E_CAL";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 1400, 0, 1.4, "COMPTONTELESCOPE/CAL/BACKE");

      // E Front and Back vs Strip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_BACK_E_CAL";
      AddHisto2D(name, name, ntot, 0, ntot, 1400, 0, 1.4, "COMPTONTELESCOPE/CAL/FRONT_BACK_E");
      
      // FRONT_T_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_T_CAL";
      AddHisto1D(name, name, 10000, 0, 1e10, "COMPTONTELESCOPE/CAL/FRONTT");

      // BACK_T_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_T_CAL";
      AddHisto1D(name, name, 10000, 0, 1e10, "COMPTONTELESCOPE/CAL/BACKT");

      // FRONT_CAL_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/CAL/MULT");

      // BACK_CAL_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/CAL/MULT");

      // Front-Back Energy Correlation
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FB_COR_CAL";
      AddHisto2D(name, name, 1400,0,1.4, 1400,0,1.4, "COMPTONTELESCOPE/CAL/FB");

    }  // end loop on number of detectors
  }// end loop on number of telescopes
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::InitPhysicsSpectra()
{
  string name;
  int ntot = (fNumberOfStripsFront+fNumberOfStripsBack);

  // counters
  name = "CT_DSSSD_COUNTERS_EVTS";
  AddHisto1D(name, name, fNumberOfCounters, 0, fNumberOfCounters, "COMPTONTELESCOPE/PHY");

  name = "CT_DSSSD_COUNTERS_HITS";
  AddHisto1D(name, name, fNumberOfCounters, 0, fNumberOfCounters, "COMPTONTELESCOPE/PHY");

  // loop on number of telescopes
  for (unsigned int i = 0 ; i < fNumberOfTelescope ; i++) {

    //// DSSSD
    // loop on number of detectors
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {

      // E Front + Back vs Strip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_BACK_E_PHY";
      AddHisto2D(name, name, ntot, 0, ntot, 1400, 0, 1.4, "COMPTONTELESCOPE/PHY");

      // front back Energy Correlation
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FB_COR_PHY";
      AddHisto2D(name, name, 1400,0,1.4, 1400,0,1.4, "COMPTONTELESCOPE/PHY");
 
      // Front E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONTE_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1.4, "COMPTONTELESCOPE/PHY");

      // Back E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACKE_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1.4, "COMPTONTELESCOPE/PHY");

      // Half E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_HALFE_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1.4, "COMPTONTELESCOPE/PHY");

      // Multiplicity
      int MultMax = fNumberOfStripsFront*fNumberOfStripsBack;
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_MULT_PHYS";
      AddHisto1D(name, name, MultMax, 0, MultMax, "COMPTONTELESCOPE/PHY");
 
      // X-Y Impact Matrix
      //  name = "CT_IMPACT_MATRIX";
      //  AddHisto2D(name, name,500,-150,150,500,-150,150, "COMPTONTELESCOPE/PHY");

    }

    //// Calorimeter

    // Calorimeter energy spectrum
    name = "CT"+NPL::itoa(i+1)+"_CALOR_SPECTRUM";
    AddHisto1D(name, name, 1000, 1, 2000, "COMPTONTELESCOPE/PHY/CALOR");

    // Position on calorimeter
    name = "CT"+NPL::itoa(i+1)+"_CALOR_POS";
    AddHisto2D(name, name, 8, -24, 24, 8, -24, 24, "COMPTONTELESCOPE/PHY/CALOR_POS");

    // Sum spectrum
    name = "CT"+NPL::itoa(i+1)+"_SUM_SPECTRUM";
    AddHisto1D(name, name, 1000, 1, 2000, "COMPTONTELESCOPE/PHY/CALOR");

  }
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::FillRawSpectra(TComptonTelescopeData* RawData)
{
  string name;
  string family;

  // FRONT_E 
  for (unsigned int i = 0; i < RawData->GetCTTrackerFrontEMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_E_RAW";
    family = "COMPTONTELESCOPE/RAW/FRONTE";

    FillSpectra(family,name,
          RawData->GetCTTrackerFrontEStripNbr(i),
          RawData->GetCTTrackerFrontEEnergy(i));
  }

  // BACK_E
  for (unsigned int i = 0; i < RawData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerBackEDetectorNbr(i))+"_BACK_E_RAW";
    family = "COMPTONTELESCOPE/RAW/BACKE";

    FillSpectra(family,name,
          RawData->GetCTTrackerBackEStripNbr(i),
          RawData->GetCTTrackerBackEEnergy(i));
  }

  // E Front Back vs Strip
  family = "COMPTONTELESCOPE/RAW/FRONT_BACK_E";
  for (unsigned int i = 0; i < RawData->GetCTTrackerFrontEMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_BACK_E_RAW";
    FillSpectra(family, name,
          RawData->GetCTTrackerFrontEStripNbr(i),
          RawData->GetCTTrackerFrontEEnergy(i));
  }
  for (unsigned int i = 0; i < RawData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerBackEDetectorNbr(i))+"_FRONT_BACK_E_RAW";
    FillSpectra(family, name,
        fNumberOfStripsFront+RawData->GetCTTrackerBackEStripNbr(i),
        RawData->GetCTTrackerBackEEnergy(i));
  }

  // FRONT_T
  for (unsigned int i = 0; i < RawData->GetCTTrackerFrontTMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerFrontTDetectorNbr(i))+"_FRONT_T_RAW";
    family = "COMPTONTELESCOPE/RAW/FRONTT";

    FillSpectra(family,name,RawData->GetCTTrackerFrontTTime(i));
    //FillSpectra(family,name,RawData->GetCTTrackerFrontTStripNbr(i),RawData->GetCTTrackerFrontTTime(i));
  }
  // BACK_T
  for (unsigned int i = 0; i < RawData->GetCTTrackerBackTMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerBackTDetectorNbr(i))+"_BACK_T_RAW";
    family = "COMPTONTELESCOPE/RAW/BACKT";

    FillSpectra(family,name,RawData->GetCTTrackerBackTTime(i));
    //FillSpectra(family,name,RawData->GetCTTrackerBackTStripNbr(i),RawData->GetCTTrackerBackTTime(i));
  }

  // FRONT MULT
  int myMULT[fNumberOfTelescope][fNumberOfDetectors];
  //int myMULT[fNumberOfTelescope];
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }

  for (unsigned int i = 0 ; i < RawData->GetCTTrackerFrontEMult(); i++) { 
     myMULT[RawData->GetCTTrackerFrontETowerNbr(i)-1][RawData->GetCTTrackerFrontEDetectorNbr(i)-1] += 1;  
     //myMULT[RawData->GetCTTrackerFrontEDetectorNbr(i)-1] += 1;  
  }

  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_RAW_MULT";
      family= "COMPTONTELESCOPE/RAW/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // BACK MULT
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }

  for (unsigned int i = 0 ; i < RawData->GetCTTrackerBackEMult(); i++) {
     myMULT[RawData->GetCTTrackerBackETowerNbr(i)-1][RawData->GetCTTrackerBackEDetectorNbr(i)-1] += 1;  
  }

  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_RAW_MULT";
      family= "COMPTONTELESCOPE/RAW/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // CALORIMETER TRIGGERS
  for (unsigned int i = 0; i < RawData->GetCTCalorimeterTMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTCalorimeterEDetectorNbr(i))+"_CALOR_RAW_TRIGGER";
    family = "COMPTONTELESCOPE/RAW/CALORTRIGGER";
    FillSpectra(family, name, RawData->GetCTCalorimeterTChannelNbr(i)+1);
  }
}

////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::FillPreTreatedSpectra(TComptonTelescopeData* PreTreatedData)
{
  string name;
  string family;
  
  // FRONT_E
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_E_CAL";
    family = "COMPTONTELESCOPE/CAL/FRONTE";

    FillSpectra(family,name
      ,PreTreatedData->GetCTTrackerFrontEStripNbr(i), 
          PreTreatedData->GetCTTrackerFrontEEnergy(i));
  }

  // BACK_E
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackEDetectorNbr(i))+"_BACK_E_CAL";
    family = "COMPTONTELESCOPE/CAL/BACKE";

    FillSpectra(family,name
      ,PreTreatedData->GetCTTrackerBackEStripNbr(i), 
          PreTreatedData->GetCTTrackerBackEEnergy(i));
  }

  // E Front Back vs Strip
  family = "COMPTONTELESCOPE/CAL/FRONT_BACK_E";
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_BACK_E_CAL";
    FillSpectra(family, name,
        PreTreatedData->GetCTTrackerFrontEStripNbr(i),
        PreTreatedData->GetCTTrackerFrontEEnergy(i));
  }
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackEDetectorNbr(i))+"_FRONT_BACK_E_CAL";
    FillSpectra(family, name,
        fNumberOfStripsFront+PreTreatedData->GetCTTrackerBackEStripNbr(i),
        PreTreatedData->GetCTTrackerBackEEnergy(i));
  }        

  // E Front Back correlation
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
    for (unsigned int j = 0; j < PreTreatedData->GetCTTrackerBackEMult(); j++) {
      if(PreTreatedData->GetCTTrackerFrontETowerNbr(i) == PreTreatedData->GetCTTrackerBackETowerNbr(i) &&
         PreTreatedData->GetCTTrackerFrontEDetectorNbr(i) == PreTreatedData->GetCTTrackerBackEDetectorNbr(j)) {
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(i))+"_FB_COR_CAL";
        family = "COMPTONTELESCOPE/CAL/FB";
        FillSpectra(family,name,PreTreatedData->GetCTTrackerFrontEEnergy(i), PreTreatedData->GetCTTrackerBackEEnergy(j));
      }
    }
  }

  // FRONT_T
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontTMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontTDetectorNbr(i))+"_FRONT_T_CAL";
    family = "COMPTONTELESCOPE/CAL/FRONTT";
    FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontTTime(i));
  }
  // BACK_T
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerBackTMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackTDetectorNbr(i))+"_BACK_T_CAL";
    family = "COMPTONTELESCOPE/CAL/BACKT";
    FillSpectra(family,name, PreTreatedData->GetCTTrackerBackTTime(i));
  }

  // FRONT MULT
  int myMULT[fNumberOfTelescope][fNumberOfDetectors];
  for( unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }

  for(unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerFrontEMult();i++){
    myMULT[PreTreatedData->GetCTTrackerFrontETowerNbr(i)-1][PreTreatedData->GetCTTrackerFrontEDetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfTelescope; i++){
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // BACK MULT
  for(unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }

  for(unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerBackEMult();i++){
    myMULT[PreTreatedData->GetCTTrackerBackETowerNbr(i)-1][PreTreatedData->GetCTTrackerBackEDetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfTelescope; i++){
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  
}

////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::FillPhysicsSpectra(TComptonTelescopePhysics* Physics){
  string name;
  string family = "COMPTONTELESCOPE/PHY";

  //// DSSSD

  // counters
  // for events
  name = "CT_DSSSD_COUNTERS_EVTS";
  for (unsigned int i = 0; i < fNumberOfCounters; i++) {
    FillSpectra(family, name, i, Physics->m_CounterEvt[i]);
  }
  // for hits
  name = "CT_DSSSD_COUNTERS_HITS";
  for (unsigned int i = 0; i < fNumberOfCounters; i++) {
    FillSpectra(family, name, i, Physics->m_CounterHit[i]);
  }

  // Front + Back E per strip
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_FRONT_BACK_E_PHY";
    FillSpectra(family, name, Physics->GetFrontStrip(i), Physics->GetFrontEnergy(i));
    FillSpectra(family, name, fNumberOfStripsFront+Physics->GetBackStrip(i), Physics->GetBackEnergy(i));
  }

  // front back Energy correlation
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_FB_COR_PHY";
    FillSpectra(family, name, Physics->GetFrontEnergy(i), Physics->GetBackEnergy(i)); 
  }

  // Energy spectrum
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_FRONTE_SPECTRUM";
    FillSpectra(family, name, Physics->GetFrontEnergy(i));

    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_BACKE_SPECTRUM";
    FillSpectra(family, name, Physics->GetBackEnergy(i));

    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_HALFE_SPECTRUM";
    FillSpectra(family, name, Physics->GetHalfEnergy(i));
  }

  // Multiplicity
  int myMULT[fNumberOfTelescope][fNumberOfDetectors];
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ;
    }
  }
  for (unsigned int i = 0 ; i < Physics->GetEventMultiplicity(); i++) {
    myMULT[Physics->GetTowerNumber(i)-1][Physics->GetDetectorNumber(i)-1] += 1;
  }
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_MULT_PHYS";
      FillSpectra(family, name, myMULT[i][j]);
    }
  }


  // X-Y Impact Matrix
//  name = "CT_IMPACT_MATRIX";

  // T

  //// Calorimeter

  // Position on calorimeter
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    name = "CT"+NPL::itoa(i+1)+"_CALOR_POS";
    if (Physics->CalorPosX.size() == Physics->CalorPosY.size()) {
      for (int j = 0; j < Physics->CalorPosX.size(); j++) {
        FillSpectra(family, name, Physics->CalorPosX[j], Physics->CalorPosY[j]);
      }
    } else {
      cout << "Position not treated because size of x and y position vectors differs." << endl;
    }
  }

  // Calorimeters spectra
  double energy = 0;
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    name = "CT"+NPL::itoa(i+1)+"_CALOR_SPECTRUM";
    energy = 0;
    for (unsigned int j = 0; j < Physics->Calor_E.size(); j++) {
      energy += Physics->Calor_E[j];
    }
    FillSpectra(family, name, energy);
  }

  // Sum spectrum
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    name = "CT"+NPL::itoa(i+1)+"_SUM_SPECTRUM";
    energy = 0;
    for (unsigned int j = 0; j < Physics->Strip_E.size();j++) {
      energy += Physics->Strip_E[j];
    }
    for (unsigned int j = 0; j < Physics->Calor_E.size(); j++) {
      energy += Physics->Calor_E[j];
    }
    FillSpectra(family, name, energy);
  }


/*  string name;
  string family= "COMPTONTELESCOPE/PHY";
  // X-Y Impact Matrix

  for(unsigned int i = 0 ; i < Physics->Si_E.size(); i++){
    name = "CT_IMPACT_MATRIX";
    double x = Physics->GetPositionOfInteraction(i).x();
    double y = Physics->GetPositionOfInteraction(i).y();
    FillSpectra(family,name),x,y);

    name = "CT_THETA_E";
    double Theta = Physics->GetPositionOfInteraction(i).Angle(TVector3(0,0,1));
    Theta = Theta/deg;
    FillSpectra(family,name),Theta,Physics->Si_E[i]);

    // FRONT_E_CAL
    name = "CT"+NPL::itoa( Physics->TelescopeNumber[i])+"_XY_COR";
    FillSpectra(family,name),Physics->Si_EX[i],Physics->Si_EY[i]);


    // Fill only for particle stopped in the first stage
    if(Physics->SiLi_E[i]<0 && Physics->CsI_E[i]<0){
      // E-TOF:
      name = "CT_E_TOF";
      FillSpectra(family,name)->Fill(Physics->Si_E[i],Physics->Si_T[i]);

      name = "CT"+NPL::itoa(Physics->TelescopeNumber[i])+"_E_TOF";
      FillSpectra(family,name)->Fill(Physics->Si_E[i],Physics->Si_T[i]);
    }

    double Etot=0;
    if(Physics->SiLi_E[i]>0){
      name = "CT_SILIE_E";
      Etot = Physics->SiLi_E[i];
      FillSpectra(family,name)->Fill(Physics->SiLi_E[i],Physics->Si_E[i]);

      name = "CT"+NPL::itoa(Physics->TelescopeNumber[i])+"_SILIE_E";
      FillSpectra(family,name)->Fill(Physics->SiLi_E[i],Physics->Si_E[i]);
    }

    if(Physics->CsI_E[i]>0){
      name = "CT_CSIE_E";
      Etot += Physics->CsI_E[i];
      FillSpectra(family,name)->Fill(Physics->CsI_E[i],Physics->Si_E[i]);
      name = "CT"+NPL::itoa(Physics->TelescopeNumber[i])+"_CSIE_E";
      FillSpectra(family,name)->Fill(Physics->CsI_E[i],Physics->Si_E[i]);

    }

    if(Etot>0){
      name = "CT_Etot_E";
      FillSpectra(family,name)->Fill(Etot,Physics->Si_E[i]);
      name = "CT"+NPL::itoa(Physics->TelescopeNumber[i])+"_Etot_E";
      FillSpectra(family,name)->Fill(Etot,Physics->Si_E[i]);
    }

  }
*/
}

