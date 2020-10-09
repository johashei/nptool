// nptool headers
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "NPDetectorManager.h"
#include "TComptonTelescopeData.h"
#include "TComptonTelescopePhysics.h"

// root headers

// C++ headers
#include <iostream>
#include <string>
using namespace std;

void online()
{
   ///////////////////////////////////////////////////////////////////////////
   // configure option manager
   NPOptionManager::getInstance()->Destroy();

//   char arg[1000];
//   sprintf(arg,"-D ./ComptonCAM.detector -C Calibration.txt -GH -E Example2.reaction -P %i --circular",port);
//   sprintf(arg,"-D ./ComptonCAM.detector -C calibrations.txt -GH -E Example2.reaction --circular");
//   sprintf(arg,"-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular");
   string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular";
   NPOptionManager::getInstance(arg);  

   // ROOT output file name
   RootOutput::getInstance("OnlineTree.root", "OnlineTree");

   // configure detector manager
   string detectorfileName = NPOptionManager::getInstance()->GetDetectorFile();
   cout << "detector file name from NPOptionManager: " << detectorfileName << "\n";
/*   NPL::DetectorManager* m_NPDetectorManager = new NPL::DetectorManager();
   m_NPDetectorManager->ReadConfigurationFile(detectorfileName);
   m_NPDetectorManager->InitializeRootOutput();

   // instantiate raw ComptonCAM data pointer
   auto ccamData = new TComptonTelescopeData();
   // connect raw CCAM data pointer to physics class
   auto ccamPhys = (TComptonTelescopePhysics*) m_NPDetectorManager->GetDetector("ComptonTelescope");
   ccamPhys->SetRawDataPointer(ccamData);
   
   // read data file/flux and fill ccamData object
   
   // test zone...
   ccamData->SetCTTrackerFrontETowerNbr(1);
   ccamData->SetCTTrackerFrontEDetectorNbr(1);
   ccamData->SetCTTrackerFrontEStripNbr(12);
   ccamData->SetCTTrackerFrontEEnergy(480);
   ccamData->Dump();
*/
}
