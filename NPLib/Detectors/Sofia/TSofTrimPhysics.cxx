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
 *  This class hold SofTrim Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSofTrimPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

//   ROOT
#include "TChain.h"

ClassImp(TSofTrimPhysics)


  ///////////////////////////////////////////////////////////////////////////
TSofTrimPhysics::TSofTrimPhysics()
  : m_EventData(new TSofTrimData),
  m_PreTreatedData(new TSofTrimData),
  m_EventPhysics(this),
  m_NumberOfDetectors(0), 
  m_NumberOfSections(3), 
  m_NumberOfAnodesPerSection(6) {
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSofTrimPhysics::AddDetector(TVector3 ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultiplicity();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    //to do 
  }
}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysize = m_EventData->GetMultiplicity();
  for (unsigned int i = 0; i < mysize ; ++i) {
    Double_t Energy = Cal->ApplyCalibration("SofTrim/ENERGY_SEC"+NPL::itoa(m_EventData->GetSectionNbr(i))+"_ANODE"+NPL::itoa(m_EventData->GetAnodeNbr(i))+"_ENERGY",m_EventData->GetEnergy(i));
    Double_t Time = Cal->ApplyCalibration("SofTrim/TIME_SEC"+NPL::itoa(m_EventData->GetSectionNbr(i))+"_ANODE"+NPL::itoa(m_EventData->GetAnodeNbr(i))+"_TIME",m_EventData->GetDriftTime(i));

    m_PreTreatedData->SetSectionNbr(m_EventData->GetSectionNbr(i));
    m_PreTreatedData->SetAnodeNbr(m_EventData->GetAnodeNbr(i));
    m_PreTreatedData->SetEnergy(Energy);
    m_PreTreatedData->SetDriftTime(Time);
    m_PreTreatedData->SetPileUp(m_EventData->GetPileUp(i));
    m_PreTreatedData->SetOverflow(m_EventData->GetOverflow(i));
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSofTrim.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSofTrim.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSofTrim.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSofTrim.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSofTrim";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      /*else if (whatToDo=="E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
        }

        else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
        }*/

      else {
        ReadingStatus = false;
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::Clear() {
  SectionNbr.clear();
  EnergyPair1.clear();
  EnergyPair2.clear();
  EnergyPair2.clear();
  DriftTimePair1.clear();
  DriftTimePair2.clear();
  DriftTimePair3.clear();
  EnergySum.clear();
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofTrim");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTrim " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTrim " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetector(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}


///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for(int sec = 0; sec < m_NumberOfSections; sec++){
    for(int anode = 0; anode < m_NumberOfAnodesPerSection; anode++){
      Cal->AddParameter("SofTrim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ENERGY","SofTrim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ENERGY");
      Cal->AddParameter("SofTrim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_TIME","SofTrim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_TIME");

    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SofTrim",  true );
  inputChain->SetBranchAddress("SofTrim", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SofTrim", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("SofTrim", "TSofTrimPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSofTrimPhysics::Construct() {
  return (NPL::VDetector*) new TSofTrimPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_SofTrim{
    public:
      proxy_SofTrim(){
        NPL::DetectorFactory::getInstance()->AddToken("SofTrim","SofTrim");
        NPL::DetectorFactory::getInstance()->AddDetector("SofTrim",TSofTrimPhysics::Construct);
      }
  };

  proxy_SofTrim p_SofTrim;
}

