// nptool headers
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "NPDetectorManager.h"
#include "TComptonTelescopeData.h"
#include "TComptonTelescopePhysics.h"

// root headers
#include "TCutG.h"
//#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"

// custom headers
#include "DecodeR.h"
#include "DecodeD.h"
#include "DecodeT.h"

#define __RUN__ 3

#define __1DET__

#define __TEST_ZONE__
#undef __TEST_ZONE__

#define __CIRCULAR_TREE__
#undef __CIRCULAR_TREE__

#define __USE_CUTG__
#undef __USE_CUTG__

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


//--//--// One-line setter for DSSSD(s) //--//--//
void setCTTracker(TComptonTelescopeData* ccamData, DecodeD* DD)
{
  for (int i = 0; i < DD->getEventSize(); i++) {
    if (DD -> getFaceType(i) == 0) { // front
      ccamData->SetFrontE(1, DD->getDetNbr(i)+1, DD->getStripNbr(i), DD->getEnergy(i));
      ccamData->SetFrontT(1, DD->getDetNbr(i)+1, 33, DD->getTime());
    }
    else if (DD -> getFaceType(i) == 1) { // back
      ccamData->SetBackE(1, DD->getDetNbr(i)+1, DD->getStripNbr(i), DD->getEnergy(i));
      ccamData->SetBackT(1, DD->getDetNbr(i)+1, 33, DD->getTime());
    }
  }
}

//--//--//--//--//--//--//--//--//--//--//--//--//--//--//--//
int main(int argc, char** argv)
{
  int resetCountSearch = 3;
  int timestampDiffSearch = 1000;
  int timestampNBins = 100;
  Option_t* tfoption;
  if (argc > 1 and std::string(argv[1]) == "-r") {tfoption = "recreate";}
  else {tfoption = "update";}
  auto fout = new TFile("timesearch.root", tfoption);
  auto bidim = new TH2F("bidim", "bidim", 
#ifdef __TEST_ZONE__
      150, 0, 150000,
#else
      2*resetCountSearch+1, -resetCountSearch-.5, resetCountSearch+.5, 
#endif
      timestampNBins, -timestampDiffSearch, timestampDiffSearch);
  auto deltaT = new TH1F("deltaT", "deltaT", timestampNBins, -timestampDiffSearch, timestampDiffSearch);

#ifdef __USE_CUTG__
  //TFile* fcut = new TFile("/disk/proto-data/data/CUT_Compton.root");
  TFile* fcut = new TFile("../data/CUT_Compton.root");
  TCutG* mcut = (TCutG*) fcut -> Get("CUT_Compton");
  fcut -> Close();
  cout << fcut << endl;
  cout << mcut << endl;
#endif

  ///////////////////////////////////////////////////////////////////////////
  // configure option manager
//   NPOptionManager::getInstance()->Destroy();
#ifdef __CIRCULAR_TREE__
  string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction --circular";
#else
  string arg = "-D ./ComptonCAM.detector -C calibrations.txt -GH -E ./10He.reaction";
#endif
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

#ifdef __1DET__
  ifstream is;
  is.open("/disk/proto-data/data/20210510_Bi207/mfm_rdd_rosmap_04_mfm_rosmap_04_2021-05-10_07_41_50.raw");
  is.seekg(0, ios::end);
  int len = is.tellg();
  is.seekg(0, ios::beg);
  char* buff = new char[len];
  is.read(buff, len);
  is.close();
  DecodeR* DR = new DecodeR(false, buff);
  while(DR -> getCursor() < len) {
    DR -> decodeBlobMFM();
    m_NPDetectorManager->ClearEventPhysics();
    m_NPDetectorManager->ClearEventData();
    ccamData -> SetCTCalorimeter(1, 4, DR->getPixelNumber(), DR->getTime(), DR->getData(), 64);
    m_NPDetectorManager->BuildPhysicalEvent();
    m_OutputTree->Fill();
    m_NPDetectorManager->CheckSpectraServer();
  }
  delete DR;
  delete [] buff;
  m_NPDetectorManager -> WriteSpectra();
#else

  // read data file/flux and fill ccamData object
  std::cout << "Reading data\n";
  DecodeR* DR = new DecodeR(false); // Instantiates DecodeR object reading calorimeter data flux
  DecodeT* DT = new DecodeT(false); // Instantiates DecodeT object reading trigger data flux
  DecodeD* DD = new DecodeD(false); // Instantiates DecodeD object reading DSSSD(s) data flux
//  newframe_t* event;
  #if __RUN__ == 0
  DD -> setTree("../data/20210210_run1/bb7_3309-7_cs137-20210210_11h05_coinc_run1_conv.root");
  #elif __RUN__ == 1
  DD -> setTree("/disk/proto-data/data/20210210_run1/bb7_3309-7_cs137-20210210_11h05_coinc_run1_conv.root");
  #elif __RUN__ == 2
  DD -> setTree("/disk/proto-data/data/20210304_run2/bb7_3309-7_cs137_20210304_14h35_conv.root");
  #elif __RUN__ == 3
  DD -> setTree("/disk/proto-data/data/20210305_run3/bb7_3309-7_cs137_20210305_14h53_conv.root");
  #elif __RUN__ == 4
  DD -> setTree("/disk/proto-data/data/20210407_run4/bb7_3309-7_cs137_20210407_14h53_conv.root");
  #endif
  int dlen = DD -> getLength();

  int i = 0;// ROSMAP files loop counter
  int c = 0;// Event counter
  int cc = 0;
  // Set some constants
  const int pixelNumber = 64;

  // Open data files
  ifstream iros, itrig;
  cout << "Loading data files " << std::flush;

  #if __RUN__ == 0
  itrig.open("../data/20210210_run1/mfm_trigger_202102101104.raw", ios::binary);
  #elif __RUN__ == 1
  itrig.open("/disk/proto-data/data/20210210_run1/mfm_trigger_202102101104.raw", ios::binary);
  #elif __RUN__ == 2
  #elif __RUN__ == 3
  itrig.open("/disk/proto-data/data/20210305_run3/mfm_trigger_20210305_run3.raw", ios::binary);
  #elif __RUN__ == 4
  itrig.open("/disk/proto-data/data/20210407_run4/mfm_trigger_20210407_run4.raw", ios::binary);
  #endif
  itrig.seekg(0, ios::end);
  int tlen = itrig.tellg();
  itrig.seekg(0, ios::beg);
  char* tbuff = new char[tlen];
  itrig.read(tbuff, tlen);
  itrig.close();
  cout << "... " << std::flush;

  #if __RUN__ == 0
  iros.open("../data/20210210_run1/mfm_rdd_rosmap_04_mfm_rosmap_04_2021-02-10_10_04_59.raw.0001", ios::binary);
  #elif __RUN__ == 1
  iros.open("/disk/proto-data/data/20210210_run1/mfm_rdd_rosmap_04_mfm_rosmap_04_2021-02-10_10_04_59.raw.0001", ios::binary);
  #elif __RUN__ == 2
  #elif __RUN__ == 3
  iros.open("/disk/proto-data/data/20210305_run3/mfm_rosmap_20210305_run3.raw", ios::binary);
  #elif __RUN__ == 4
  iros.open("/disk/proto-data/data/20210407_run4/mfm_rosmap_20210407_run4.raw", ios::binary);
  #endif
  iros.seekg(0, ios::end);
  int rlen = iros.tellg();
  iros.seekg(0, ios::beg);
  char* rbuff = new char[rlen];
  iros.read(rbuff, rlen);
  iros.close();
  cout << "Done" << endl;

  // Search for reset count in trigger data
  DT -> setRaw(tbuff);
  DT -> decodeBlobMFM();
  int resetCount = DT -> getResetCount();
  //while (not(DT->hasTrigged(2))) {
  while (not(DT->hasTrigged(0))) {
    DT -> decodeBlobMFM();
  }
  resetCount = DT->getResetCount() - resetCount;
  resetCount = -resetCount; // T B C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  cout << "Found reset count: " << resetCount << endl;

  if (resetCount == 0) {
    string ans;
    cout << "reset count is 0. Continue ? y/[n]";
    cin >> ans;
    if (ans != "y" and ans != "Y") {
      DT -> setRaw(tbuff);
      DT -> decodeBlobMFM();
      resetCount = DT -> getResetCount();
      while (not(DT->hasTrigged(2))) {
        DT -> decodeBlobMFM();
      }
      resetCount = DT->getResetCount() - resetCount;
      cout << "Found reset count: " << resetCount << endl;
    }
  }
  int cr, cd, tr, td;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (argc > 1 and std::string(argv[1]) == "-r") {
    cout << "Now looking for a potential error in reset number by filling the time bidim" << endl;
  
    // Fill control bidim
    for (int reset=-resetCountSearch; reset<resetCountSearch+1; reset++)
    {
      cr = resetCount+reset;
      cout << "Filling coincidence bidim with reset #" << cr << endl;
      DD -> rewind();
  
      DR -> setRaw(rbuff);
      DR -> decodeBlobMFM();
      DD -> decodeEvent();
  
      cd = 0;
      c = 0;
      i = 0;
      tr = DR -> getTime();
      td = DD -> getTime();
      while(c < 1000)
      {
  //    cout << DR -> getTime() << " " << DD -> getTime() << endl; 
        if (cr == cd) {
          if (abs(td-tr) < timestampDiffSearch) {
            bidim->Fill(reset, td-tr);
            c++;
          }
          if (td < tr) {
            DD -> decodeEvent();
          } else {
            DR -> decodeBlobMFM();
          }
        } else if (cr < cd) {
          DR -> decodeBlobMFM();
        } else {
          DD -> decodeEvent();
        }
        if (DR -> getTime() < tr) {
          cr++;
        }
        tr = DR -> getTime();
        if (DD -> getTime() < td) {
          cd ++;
        }
        td = DD -> getTime();
        i++; 
      } // End of event loop
    } // End of for loop
  
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    c = 0;
    cout << "Now filling deltaT histogram with all available events" << endl;

    DR -> setRaw(rbuff);
    DR -> decodeBlobMFM();
    DD -> rewind();
    DD -> decodeEvent();

    cr = resetCount; cd = 0;
    tr = DR -> getTime();
    td = DD -> getTime();

    while(DR -> getCursor() < rlen and DD -> getCursor() < dlen)
    {
      if (cr == cd) {
        if (abs(td-tr) < timestampDiffSearch) {
          deltaT -> Fill(td-tr);
          c++;
        }
        if (td < tr) {
          DD -> decodeEvent();
        } else {
          DR -> decodeBlobMFM();
        }
      } else if (cr < cd) {
        DR -> decodeBlobMFM();
      } else {
        DD -> decodeEvent();
      }
      if (DR -> getTime() < tr) {
        cr++;
      }
      tr = DR -> getTime();
      if (DD -> getTime() < td) {
        cd ++;
      }
      td = DD -> getTime();
      i++;
    } // End of event loop
    cout << "Out of " << i << " events treated, " << c << " coincidences were found and added to the histogram" << endl;


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  } else {// non reset search mode

    c = 0;
    
    DR -> setRaw(rbuff);
    DR -> decodeBlobMFM();

//#ifdef __TEST_ZONE__
//  cout << "Entering test zone." << endl;

//#else
    DD -> rewind();
    DD -> decodeEvent();

    cr = resetCount; cd = 0;
    tr = DR -> getTime();
    td = DD -> getTime();

    while(DR -> getCursor() < rlen and DD -> getCursor() < dlen and c < 10000)
    {
      if (cr == cd) {
  #ifndef __TEST_ZONE__
        if (td-tr > 20 and td-tr < 120) { // That one is the real one
        //if (td-tr > -1000 and td-tr < 1000) { // That one is the real one
  #else
        if (td-tr > -1000  and td-tr < 0) {
  #endif
          //DR -> Dump();
          //DD -> Dump();
          c++;
          // Clear raw and physics data
          m_NPDetectorManager->ClearEventPhysics();
          m_NPDetectorManager->ClearEventData();
  
          // Fill data
          ccamData -> SetCTCalorimeter(1, 4, DR->getPixelNumber(), DR->getTime(), DR->getData(), pixelNumber);
          setCTTracker(ccamData, DD);// -> getEvent(), &nb_asic, &chain, stripNumber);
          ccamData -> SetResetCount(cr);
  
          // Build physical event
          m_NPDetectorManager->BuildPhysicalEvent();
  
          // Fill object in output ROOT file
          if (ccamPhys->EventMultiplicity > 0) {
  #ifdef __USE_CUTG__
            if (mcut->IsInside(ccamPhys->Half_Energy[0], ccamPhys->Calor_E[0])) {
              cc++;
              m_OutputTree->Fill();
            }
  #else
            m_OutputTree->Fill();
            cc++;
  #endif
            cout << cc << " " << c /*<< "(" << cr << ", " << cd << ") : " << tr << " " << td*/ << " |\t";
          }
  
          // check spectra
          m_NPDetectorManager->CheckSpectraServer();
          //cout << "d" << endl;
        }
        if (td < tr) {
          DD -> decodeEvent();
        } else {
          DR -> decodeBlobMFM();
        }
      } else if (cr < cd) {
        DR -> decodeBlobMFM();
      } else {
        DD -> decodeEvent();
      }
      if (DR -> getTime() < tr) {
        cr++;
      }
      tr = DR -> getTime();
      if (DD -> getTime() < td) {
        cd ++;
      }
      td = DD -> getTime();
    
      i++;
    
    } // End of main loop
//#endif
  } // End of mode if

  delete DR;
  delete DT;
  delete DD;
  delete [] rbuff;
  delete [] tbuff;

  // fill spectra
  m_NPDetectorManager -> WriteSpectra();

  if (argc > 1 and std::string(argv[1]) == "-r") {
    fout->cd();
    bidim->Write();
    deltaT->Write();
  } 
  fout->Close();
#endif
  // Essential
#if __cplusplus > 199711L && NPMULTITHREADING
  m_NPDetectorManager->StopThread();
#endif
  RootOutput::Destroy();

  return 0;
}


