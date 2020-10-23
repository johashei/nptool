// nptool headers
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "NPDetectorManager.h"
#include "TComptonTelescopeData.h"
#include "TComptonTelescopePhysics.h"

// root headers

// custom headers
#include "DecodeR.h"

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


int main()
{
  ///////////////////////////////////////////////////////////////////////////
  // configure option manager
//   NPOptionManager::getInstance()->Destroy();

//  string arg = "-D ./ComptonCAM.detector -C Calibration.txt -GH -E Example2.reaction -P %i --circular",port);
  string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular";
  NPOptionManager::getInstance(arg);  

  // open ROOT output file
  RootOutput::getInstance("OnlineTree.root", "OnlineTree");
  // get tree pointer
  auto m_OutputTree = RootOutput::getInstance()->GetTree();

  // configure detector manager
  string detectorfileName = NPOptionManager::getInstance()->GetDetectorFile();
  NPL::DetectorManager* m_NPDetectorManager = new NPL::DetectorManager();
  m_NPDetectorManager->ReadConfigurationFile(detectorfileName);
  m_NPDetectorManager->InitializeRootOutput();

  // instantiate raw ComptonCAM data pointer
  auto ccamData = new TComptonTelescopeData();
  ccamData->Dump();
  // connect raw CCAM data pointer to physics class
  auto ccamPhys = (TComptonTelescopePhysics*) m_NPDetectorManager->GetDetector("ComptonTelescope");
  ccamPhys->SetRawDataPointer(ccamData);

  // read data file/flux and fill ccamData object
  std::cout << "Reading data\n";
  // instantiate DecodeR object reading calorimeter data flux
  DecodeR* D = new DecodeR(false);
  D->Dump();

  // Load a file
  std::ifstream is;
  is.open("./mfm.bin", std::ios::binary);
  is.seekg (0, std::ios::end);
  int length = is.tellg();
  is.seekg (0, ios::beg);
  char* buffer = new char [length];
  is.read(buffer, length);
  is.close();

  // Read from file
  D -> setRaw(buffer);
  D -> decodeRawMFM(); // get rid of the first two (empty) events
  D -> decodeRawMFM();
  D -> decodeRawMFM();
  D -> Dump();

  int c = 0;
  const int pixelNumber = 64;
  while (D -> getCursor() < length) 
  {
     // Clear raw data and physics objects
     m_NPDetectorManager->ClearEventPhysics();
     m_NPDetectorManager->ClearEventData();

     // Read the actual data
     D -> decodeRawMFM();
     //D -> Dump();//Optionnal print

     // Set ccamData (a better way is envisionned)
     for (int i = 0; i < pixelNumber; ++i) {
       ccamData -> SetCTCalorimeterTTowerNbr( 1 );
       ccamData -> SetCTCalorimeterTDetectorNbr( 1 );//Triggered ASIC number
       ccamData -> SetCTCalorimeterTChannelNbr( D -> getPixelNumber() );//ASIC's channel number
       ccamData -> SetCTCalorimeterTTime( D -> getTime() );
       ccamData -> SetCTCalorimeterETowerNbr(1);
       ccamData -> SetCTCalorimeterEDetectorNbr( 1 );
       ccamData -> SetCTCalorimeterEChannelNbr( i );//PMT pixel number
       ccamData -> SetCTCalorimeterEEnergy( D -> getData()[i] );
     }
//     ccamData -> Dump();

     // Fill object in output ROOT file
     m_OutputTree->Fill();

     c++;
  }
  delete D;
  delete [] buffer;

  std::cout << "test compil\n";

  // Fill spectra
  m_NPDetectorManager->WriteSpectra();

  // Essential
  #if __cplusplus > 199711L && NPMULTITHREADING
   m_NPDetectorManager->StopThread();
  #endif
  RootOutput::Destroy();

  return 0;
}
