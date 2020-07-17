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
 *  This class hold Strasse Treated  data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TStrassePhysics.h"

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

ClassImp(TStrassePhysics)
  ///////////////////////////////////////////////////////////////////////////
  TStrassePhysics::TStrassePhysics(){
    EventMultiplicity = 0;
    m_EventData = new TStrasseData;
    m_PreTreatedData = new TStrasseData;
    m_EventPhysics = this;
    m_Spectra = NULL;
    m_E_RAW_Threshold = 0; // adc channels
    m_E_Threshold = 0;     // MeV
    m_NumberOfDetectors = 0;
    m_MaximumStripMultiplicityAllowed = 10;
    m_StripEnergyMatching = 0.050;
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TStrassePhysics::AddDetector(TVector3){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::AddDetector(double R, double Theta, double Phi){
  m_NumberOfDetectors++;

  //double Height = 118; // mm
  double Height = 61.8; // mm
  //double Base = 95; // mm
  double Base = 78.1; // mm
  double NumberOfStrips = 128;
  double StripPitchHeight = Height / NumberOfStrips; // mm
  double StripPitchBase = Base / NumberOfStrips; // mm


  // Vector U on detector face (parallel to Y strips) Y strips are along X axis
  TVector3 U;
  // Vector V on detector face (parallel to X strips)
  TVector3 V;
  // Vector W normal to detector face (pointing to the back)
  TVector3 W;
  // Vector C position of detector face center
  TVector3 C;

  C = TVector3(R*sin(Theta)*cos(Phi),
        R*sin(Theta)*sin(Phi),
        Height*0.5+R*cos(Theta));

  TVector3 P = TVector3(cos(Theta)*cos(Phi),
      cos(Theta)*sin(Phi),
      -sin(Theta));

  W = C.Unit();
  U = W.Cross(P);
  V = W.Cross(U);

  U = U.Unit();
  V = V.Unit();

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  double X, Y, Z;

  // Moving C to the 1.1 Corner;
  TVector3 Strip_1_1;
  Strip_1_1 = C - (0.5*Base*U + 0.5*Height*V) + U*(StripPitchBase / 2.) + V*(StripPitchHeight / 2.);

  TVector3 StripPos;
  for(int i=0; i<NumberOfStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<NumberOfStrips; j++){
      StripPos = Strip_1_1 + i*U*StripPitchBase + j*V*StripPitchHeight;
      lineX.push_back(StripPos.X());
      lineY.push_back(StripPos.Y());
      lineZ.push_back(StripPos.Z());
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_StripPositionX.push_back(OneDetectorStripPositionX);
  m_StripPositionY.push_back(OneDetectorStripPositionY);
  m_StripPositionZ.push_back(OneDetectorStripPositionZ);
} 

///////////////////////////////////////////////////////////////////////////
TVector3 TStrassePhysics::GetPositionOfInteraction(const int i){
  TVector3 Position = TVector3(GetStripPositionX(DetectorNumber[i], StripX[i], StripY[i]),
      GetStripPositionY(DetectorNumber[i], StripX[i], StripY[i]),
      GetStripPositionZ(DetectorNumber[i], StripX[i], StripY[i]));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TStrassePhysics::GetDetectorNormal(const int i){
  TVector3 U = TVector3(GetStripPositionX(DetectorNumber[i],128,1),
      GetStripPositionY(DetectorNumber[i],128,1),
      GetStripPositionZ(DetectorNumber[i],128,1))

    -TVector3(GetStripPositionX(DetectorNumber[i],1,1),
      GetStripPositionY(DetectorNumber[i],1,1),
      GetStripPositionZ(DetectorNumber[i],1,1));

  TVector3 V = TVector3(GetStripPositionX(DetectorNumber[i],128,128),
      GetStripPositionY(DetectorNumber[i],128,128),
      GetStripPositionZ(DetectorNumber[i],128,128))

    -TVector3(GetStripPositionX(DetectorNumber[i],128,1),
      GetStripPositionY(DetectorNumber[i],128,1),
      GetStripPositionZ(DetectorNumber[i],128,1));

  TVector3 Normal = U.Cross(V);

  return (Normal.Unit());
}
///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  if(1 /*CheckEvent() == 1*/){
    vector<TVector2> couple = Match_X_Y();

    EventMultiplicity = couple.size();
    for(unsigned int i=0; i<couple.size(); i++){
      int N = m_PreTreatedData->GetInner_TE_DetectorNbr(couple[i].X());
      int X = m_PreTreatedData->GetInner_TE_StripNbr(couple[i].X());
      int Y = m_PreTreatedData->GetInner_LE_StripNbr(couple[i].Y());

      double TE = m_PreTreatedData->GetInner_TE_Energy(couple[i].X());
      double LE = m_PreTreatedData->GetInner_LE_Energy(couple[i].Y());
      DetectorNumber.push_back(N);
      StripX.push_back(X);
      StripY.push_back(Y);
      DE.push_back(TE);
      
      PosX.push_back(GetPositionOfInteraction(i).x());
      PosY.push_back(GetPositionOfInteraction(i).y());
      PosZ.push_back(GetPositionOfInteraction(i).z());

      int OuterMult = m_PreTreatedData->GetOuterMultTEnergy();
      for(unsigned int j=0; j<OuterMult; j++){
        if(m_PreTreatedData->GetOuter_TE_DetectorNbr(j)==N){
          double XDE = m_PreTreatedData->GetOuter_TE_Energy(j);
          double YDE = m_PreTreatedData->GetOuter_LE_Energy(j);

          E.push_back(XDE);
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TStrassePhysics::Match_X_Y(){
  vector<TVector2> ArrayOfGoodCouple;

  static unsigned int m_TEMult, m_LEMult;
  m_TEMult = m_PreTreatedData->GetInnerMultTEnergy();
  m_LEMult = m_PreTreatedData->GetInnerMultLEnergy();

  if(m_TEMult>m_MaximumStripMultiplicityAllowed || m_LEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_TEMult; i++){
    for(unsigned int j=0; j<m_LEMult; j++){

      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetInner_TE_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetInner_LE_DetectorNbr(j);

      // if same detector check energy
      if(XDetNbr == YDetNbr){
        // Declaration of variable for clarity
        double TE = m_PreTreatedData->GetInner_TE_Energy(i);
        double LE = m_PreTreatedData->GetInner_LE_Energy(i);
        double XStripNbr = m_PreTreatedData->GetInner_TE_StripNbr(i);
        double YStripNbr = m_PreTreatedData->GetInner_LE_StripNbr(i);

        // look if energy matches
        if(abs(TE-LE)/2.<m_StripEnergyMatching){
          ArrayOfGoodCouple.push_back(TVector2(i,j));
        }
      }
    }
  }

  return ArrayOfGoodCouple;
}

///////////////////////////////////////////////////////////////////////////
int TStrassePhysics::CheckEvent(){
  // Check the size of the different elements
  if(m_PreTreatedData->GetInnerMultTEnergy() == m_PreTreatedData->GetInnerMultLEnergy() )
    return 1;

  else
    return -1;
}


///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  //////
  // First Stage Energy
  unsigned int sizeFront = m_EventData->GetInnerMultTEnergy();
  for (UShort_t i = 0; i < sizeFront ; ++i) {
    if (m_EventData->GetInner_TE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetInner_TE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetInner_TE_DetectorNbr(i)),m_EventData->GetInner_TE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetInnerTE(m_EventData->GetInner_TE_DetectorNbr(i), m_EventData->GetInner_TE_StripNbr(i), Energy);
      }
    }
  }
  unsigned int sizeBack = m_EventData->GetInnerMultTEnergy();
  for (UShort_t i = 0; i < sizeBack ; ++i) {
    if (m_EventData->GetInner_LE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetInner_LE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetInner_LE_DetectorNbr(i)),m_EventData->GetInner_LE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetInnerLE(m_EventData->GetInner_LE_DetectorNbr(i), m_EventData->GetInner_LE_StripNbr(i), Energy);
      }
    }
  }

  //////
  // Second Stage Energy
  sizeFront = m_EventData->GetOuterMultTEnergy();
  for (UShort_t i = 0; i < sizeFront ; ++i) {
    if (m_EventData->GetOuter_TE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetOuter_TE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetOuter_TE_DetectorNbr(i)),m_EventData->GetOuter_TE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetOuterTE(m_EventData->GetOuter_TE_DetectorNbr(i), m_EventData->GetOuter_TE_StripNbr(i), Energy);
      }
    }
  }
  sizeBack = m_EventData->GetOuterMultTEnergy();
  for (UShort_t i = 0; i < sizeBack ; ++i) {
    if (m_EventData->GetOuter_LE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetOuter_LE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetOuter_LE_DetectorNbr(i)),m_EventData->GetOuter_LE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetOuterLE(m_EventData->GetOuter_LE_DetectorNbr(i), m_EventData->GetOuter_LE_StripNbr(i), Energy);
      }
    }
  }

}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigStrasse.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigStrasse.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigStrasse.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigStrasse.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigStrasse";
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
void TStrassePhysics::Clear() {
  EventMultiplicity = 0;

  // Position Information
  PosX.clear();
  PosY.clear();
  PosZ.clear();

  // DSSD
  DetectorNumber.clear();
  E.clear();
  StripX.clear();
  StripY.clear();
  DE.clear();
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Strasse");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");

      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");

      AddDetector(R, Theta, Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitSpectra() {
  m_Spectra = new TStrasseSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TStrassePhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Strasse", "D"+ NPL::itoa(i+1)+"_ENERGY","Strasse_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("Strasse", "D"+ NPL::itoa(i+1)+"_TIME","Strasse_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Strasse",  true );
  inputChain->SetBranchAddress("Strasse", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Strasse", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Strasse", "TStrassePhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TStrassePhysics::Construct() {
  return (NPL::VDetector*) new TStrassePhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_Strasse{
    public:
      proxy_Strasse(){
        NPL::DetectorFactory::getInstance()->AddToken("Strasse","Strasse");
        NPL::DetectorFactory::getInstance()->AddDetector("Strasse",TStrassePhysics::Construct);
      }
  };

  proxy_Strasse p_Strasse;
}

