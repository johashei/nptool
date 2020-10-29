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
  //string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction";
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
  
  int i = 0;
  int c = 0;
  const int pixelNumber = 64;
  while (true)
//  while (i<1)
  {
    // Load a file(s)
    std::ifstream is;
    i = 0;
    switch (i % 6) {
      case 0: is.open("./mfm.bin", std::ios::binary); break;
      case 3: is.open("./133Ba.bin", std::ios::binary); break;
      case 2: is.open("./241Am.bin", std::ios::binary); break;
      case 1: is.open("./207Bi.bin", std::ios::binary); break;
      case 4: is.open("./241Am-1.bin", std::ios::binary); break;
      case 5: is.open("./241Am-2.bin", std::ios::binary); break;
    }
    is.seekg (0, std::ios::end);
    int length = is.tellg();
    is.seekg (0, ios::beg);
    char* buffer = new char [length];
    is.read(buffer, length);
    is.close();
    i++;
  
    // Read from file(s)
    D -> setRaw(buffer);

    D -> decodeRawMFM(); // get rid of the first two (empty) events
    D -> decodeRawMFM();
  
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
         ccamData -> SetCTCalorimeterETowerNbr(1);
         ccamData -> SetCTCalorimeterEDetectorNbr( 1 );
         ccamData -> SetCTCalorimeterEChannelNbr( i );//PMT pixel number
         ccamData -> SetCTCalorimeterEEnergy( D -> getData()[i] );
       }
       ccamData -> SetCTCalorimeterTTowerNbr( 1 );
       ccamData -> SetCTCalorimeterTDetectorNbr( 1 );//Triggered ASIC number
       ccamData -> SetCTCalorimeterTChannelNbr( D -> getPixelNumber() );//Pixel that triggered
       ccamData -> SetCTCalorimeterTTime( D -> getTime() );
  //     ccamData -> Dump();
  
       // Build physical event
       m_NPDetectorManager->BuildPhysicalEvent();
  
       // Fill object in output ROOT file
       m_OutputTree->Fill();
  
       c++;
       //usleep(10000);//Simulated 100Hz count rate
    }
  
    //std::cout << "test compil\n";
  
    delete [] buffer;

  }// End of main loop

  // Fill spectra
  m_NPDetectorManager->WriteSpectra();

  // Essential
  delete D;

  #if __cplusplus > 199711L && NPMULTITHREADING
   m_NPDetectorManager->StopThread();
  #endif
  RootOutput::Destroy();

  return 0;
}


