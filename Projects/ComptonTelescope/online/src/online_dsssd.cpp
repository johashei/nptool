// nptool headers
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "NPDetectorManager.h"
#include "TComptonTelescopeData.h"
#include "TComptonTelescopePhysics.h"

// root headers
#include "TFile.h"

// custom headers
#include "DecodeD.h"

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

//--//--//--//--//--//--//--//--//--//--//--//--//--//--//--//
int main()
{

  ///////////////////////////////////////////////////////////////////////////
  // configure option manager
//   NPOptionManager::getInstance()->Destroy();

  //string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular";
  string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction";
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

  // configure spectra server
  m_NPDetectorManager->SetSpectraServer();

  // instantiate raw ComptonCAM data pointer
  auto ccamData = new TComptonTelescopeData();
  // connect raw CCAM data pointer to physics class
  auto ccamPhys = (TComptonTelescopePhysics*) m_NPDetectorManager->GetDetector("ComptonTelescope");
  ccamPhys->SetRawDataPointer(ccamData);

  // read data file/flux and fill ccamData object
  std::cout << "Reading data\n";
  DecodeD* DD = new DecodeD(true); // Instantiates DecodeD object reading DSSSD(s) data flux
  newframe_t* event;
  DD -> setTree("../data/20200128_10h44_bi207_conv.root");
  //DD -> setTree("../data/bb7_3309-7_bi207_20210209_12h50_run10_conv.root");
  int dlen = DD -> getLength();

  while (DD -> getCursor() < 49970)
    //while (DD -> getCursor() < dlen)
  {
    // cout number of entries treated
    if (DD -> getCursor() % 10000 == 0) {
      cout << "\rEntry " << DD -> getCursor();
      cout << flush;
    }

    // Decode event
    DD->Clear();
    DD->decodeEventFinal();

    // Clear raw and physics data
    m_NPDetectorManager->ClearEventPhysics();
    m_NPDetectorManager->ClearEventData();

    // Fill data
    //cout << "event size " << DD->getEventSize() << endl;      
    for (int i = 0; i < DD->getEventSize(); i++) {
      //cout << i << endl;
      if (DD -> getFaceType(i) == 0) { // front
        ccamData->SetFrontE(1, DD->getDetNbr(i)+1, DD->getStripNbr(i), DD->getEnergy(i));
        ccamData->SetFrontT(1, DD->getDetNbr(i)+1, 33, DD->getTime());
        //cout << "face " << DD -> getFaceType(i) << " det " << DD->getDetNbr(i) << " strip " << DD->getStripNbr(i) << " E " << DD->getEnergy(i) << " T " << DD->getTime() << endl;      
      }
      else if (DD -> getFaceType(i) == 1) { // back
        ccamData->SetBackE(1, DD->getDetNbr(i)+1, DD->getStripNbr(i), DD->getEnergy(i));
        ccamData->SetBackT(1, DD->getDetNbr(i)+1, 33, DD->getTime());
        //cout << "face " << DD -> getFaceType(i) << " det " << DD->getDetNbr(i) << " strip " << DD->getStripNbr(i) << " E " << DD->getEnergy(i) << " T " << DD->getTime() << endl;      
      }
    }

    // Build physical event
    m_NPDetectorManager->BuildPhysicalEvent();

    // Fill object in output ROOT file
    m_OutputTree->Fill();

    // check spectra
    m_NPDetectorManager->CheckSpectraServer();
  }

  // fill spectra
  m_NPDetectorManager -> WriteSpectra();


  // Essential
#if __cplusplus > 199711L && NPMULTITHREADING
  m_NPDetectorManager->StopThread();
#endif
  RootOutput::Destroy();

  return 0;
}

