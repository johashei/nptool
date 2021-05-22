/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Nebula Treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TNebulaPhysics.h"

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

ClassImp(TNebulaPhysics)


///////////////////////////////////////////////////////////////////////////
TNebulaPhysics::TNebulaPhysics()
   : m_EventData(new TNebulaData),
     m_PreTreatedData(new TNebulaData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfBars(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TNebulaPhysics::ReadXML(NPL::XmlParser xml){ 
  m_NumberOfBars++;
} 

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();
/*
  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  unsigned int mysizeT = m_PreTreatedData->GetMultTime();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    for (UShort_t t = 0; t < mysizeT ; t++) {
      if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t)) {
        DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
        Energy.push_back(m_PreTreatedData->Get_Energy(e));
        Time.push_back(m_PreTreatedData->Get_Time(t));
      }
    }
  }*/
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();
/*
  // Energy
  unsigned int mysize = m_EventData->GetMultEnergy();
  for (UShort_t i = 0; i < mysize ; ++i) {
    if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = Cal->ApplyCalibration("Nebula/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("Nebula/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }
  */
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigNebula.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigNebula.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigNebula.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigNebula.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigNebula";
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

      else if (whatToDo=="E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("NEBULA");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML","Offset","InvertX","InvertY"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      cout << endl << "////  Nebula (" << i+1 << ")" << endl;
      unsigned int det = std::atoi(blocks[i]->GetMainValue().c_str());
      string xmlpath = blocks[i]->GetString("XML");
      NPL::XmlParser xml;
      xml.LoadFile(xmlpath);
      ReadXML(xml);
      TVector3 offset = blocks[i]->GetTVector3("Offset","mm"); 
      bool invertX = blocks[i]->GetInt("InvertX"); 
      bool invertY = blocks[i]->GetInt("InvertY"); 
      m_offset[det] = offset;
      m_invertX[det] = invertX;
      m_invertY[det] = invertY;
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitSpectra() {
  m_Spectra = new TNebulaSpectra(m_NumberOfBars);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TNebulaPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfBars; ++i) {
    Cal->AddParameter("NEBULA_ID"+ NPL::itoa(i+1)+"_T");
    Cal->AddParameter("NEBULA_ID"+ NPL::itoa(i+1)+"_Y");
  }
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Nebula",  true );
  inputChain->SetBranchAddress("Nebula", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Nebula", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree("Nebula");
  outputTree->Branch("Nebula", "TNebulaPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TNebulaPhysics::Construct() {
  return (NPL::VDetector*) new TNebulaPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Nebula{
  public:
    proxy_Nebula(){
      NPL::DetectorFactory::getInstance()->AddToken("NEBULA","Nebula");
      NPL::DetectorFactory::getInstance()->AddDetector("NEBULA",TNebulaPhysics::Construct);
    }
};

proxy_Nebula p_Nebula;
}

