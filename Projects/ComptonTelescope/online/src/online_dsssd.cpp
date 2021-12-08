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
#include <cstdlib>
using namespace std;

//--//--//--//--//--//--//--//--//--//--//--//--//--//--//--//
int main(int argc, char *argv[])
{
  ///////////////////////////////////////////////////////////////////////////
  // configure option manager
//   NPOptionManager::getInstance()->Destroy();

  string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular";
  //string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction";
  NPOptionManager::getInstance(arg);  

  // open ROOT output file
  RootOutput::getInstance("OnlineTree.root", "OnlineTree");
  //RootOutput::getInstance("OnlineTree_DSSSD.root", "OnlineTree");
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

  // deal with file name and number of events to treat from command line
  std::cout << "\n";
  int nevents = -1;
  string fileName = "bb7_3309-7_cs137_20210304_14h35_conv.root";
  if (argc == 1) {
     std::cout << "Name of file to analyse should be provided as a command line argument\n";
     std::cout << "Default file " << fileName << " is considered\n";
  }
  else {
     fileName = argv[1];
     if (argc == 3) nevents = std::atoi(argv[2]);
  }

  // read data file/flux and fill ccamData object
  std::cout << "Reading data from " << fileName << "\n";
  DecodeD* DD = new DecodeD(true); // Instantiates DecodeD object reading DSSSD(s) data flux
  newframe_t* event;
  fileName = "../data/" + fileName;
  DD -> setTree(fileName.c_str());
//  DD -> setTree("../data/20200128_10h44_bi207_conv.root");
//  DD -> setTree("../data/bb7_3309-7_cs137_20210304_14h35_conv.root");
  //DD -> setTree("../data/bb7_3309-7_bi207_20210209_12h50_run10_conv.root");
  int dlen = DD -> getLength();
  // read a limited  number of events
  if (nevents > 0) dlen = nevents;
  std::cout << "Reading the first " << dlen << " entries\n";

  while (DD -> getCursor() < dlen)
  {
    // cout number of entries treated
    if (DD -> getCursor() % 10000 == 0) {
      cout << "\rEntry " << DD -> getCursor();
      cout << flush;
    }

    // Decode event
    DD->Clear();
    DD->decodeEvent();

    // Clear raw and physics data
    m_NPDetectorManager->ClearEventPhysics();
    m_NPDetectorManager->ClearEventData();

    // Fill data
    //cout << "event size " << DD->getEventSize() << endl;
    //ccamData -> Dump();
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

