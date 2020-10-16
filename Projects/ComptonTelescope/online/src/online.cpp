// nptool headers
//#include "NPOptionManager.h"
//#include "RootOutput.h"
//#include "NPDetectorManager.h"
//#include "TComptonTelescopeData.h"
//#include "TComptonTelescopePhysics.h"

// root headers

// custom headers
#include "DecodeR.h"

// C++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

int main()
{
   ///////////////////////////////////////////////////////////////////////////
   // configure option manager
//   NPOptionManager::getInstance()->Destroy();

//   char arg[1000];
//   sprintf(arg,"-D ./ComptonCAM.detector -C Calibration.txt -GH -E Example2.reaction -P %i --circular",port);
//   sprintf(arg,"-D ./ComptonCAM.detector -C calibrations.txt -GH -E Example2.reaction --circular");
//   sprintf(arg,"-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular");
/*   string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular";
   NPOptionManager::getInstance(arg);  

   // ROOT output file name
   RootOutput::getInstance("OnlineTree.root", "OnlineTree");

   // configure detector manager
   string detectorfileName = NPOptionManager::getInstance()->GetDetectorFile();
   cout << "detector file name from NPOptionManager: " << detectorfileName << "\n";
   NPL::DetectorManager* m_NPDetectorManager = new NPL::DetectorManager();
   m_NPDetectorManager->ReadConfigurationFile(detectorfileName);*/
/*   m_NPDetectorManager->InitializeRootOutput();*/

   // instantiate raw ComptonCAM data pointer
//   auto ccamData = new TComptonTelescopeData();
   // connect raw CCAM data pointer to physics class
/*   auto ccamPhys = (TComptonTelescopePhysics*) m_NPDetectorManager->GetDetector("ComptonTelescope");
   ccamPhys->SetRawDataPointer(ccamData);*/
   
   // read data file/flux and fill ccamData object

   cout << "Reading data" << endl;
   DecodeR* D = new DecodeR();
   D->Dump();
//   bool verbose = false;
/*   
   // Load a file
   ifstream is;
   is.open("decodeR/mfm.bin", ios::binary);
   is.seekg (0, ios::end);
   int length = is.tellg();
   is.seekg (0, ios::beg);
   char* buffer = new char [length];
   is.read(buffer, length);
   is.close();

   // Read from file
   D -> setRaw(buffer);
   D -> decodeRawMFM(verbose);//Get rid of the first two (empty) events
   D -> decodeRawMFM(verbose);

   int c = 0;
   int i = 0;

   while (D -> getCursor() < length)
   {
      // Read the actual data
      D -> decodeRawMFM(verbose);
      //D -> Dump();//Optionnal print

      // Set ccamData (a better way is envisionned)
      for (int i=0; i<64; i++) {
        ccamData -> SetCTCalorimeterTTowerNbr( 1 );
        ccamData -> SetCTCalorimeterTDetectorNbr( D -> getSource() );//Triggered ASIC number
        ccamData -> SetCTCalorimeterTChannelNbr( D -> getChannel() );//ASIC's channel number
        ccamData -> SetCTCalorimeterTTime( D -> getTime() );
        ccamData -> SetCTCalorimeterETowerNbr(1);
        ccamData -> SetCTCalorimeterEDetectorNbr( 1 );
        ccamData -> SetCTCalorimeterEChannelNbr( i );//PMT pixel number
        ccamData -> SetCTCalorimeterEEnergy( D -> getData()[i] );
      }
      ccamData -> Dump();
      ccamData -> Clear();
      c++;
   }
   delete D;
   delete [] buffer;
*/
   // test zone...
/*
   ccamData->SetCTTrackerFrontETowerNbr(1);
   ccamData->SetCTTrackerFrontEDetectorNbr(1);
   ccamData->SetCTTrackerFrontEStripNbr(12);
   ccamData->SetCTTrackerFrontEEnergy(480);
   ccamData->Dump();
   ccamData->Clear();
   ccamData->Dump();
*/

   std::cout << "test compil\n";

   return 1;
}
