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

#define __TEST_ZONE__
#undef __TEST_ZONE__

#define __CIRCULAR_TREE__
#undef __CIRCULAR_TREE__

#define __USE_CUTG__
#undef __USE_CUTG__

#define __RESET_SEARCH__
//#undef __RESET_SEARCH__

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
int main()
{
#if defined __RESET_SEARCH__
  int resetCountSearch = 3;
  int timestampDiffSearch = 500;
  int timestampNBins = 100;
  auto fout = new TFile("pipo.root", "recreate");
  auto bidim = new TH2F("bidim", "bidim", 
#ifdef __TEST_ZONE__
      150, 0, 150000,
#else
      2*resetCountSearch+1, -resetCountSearch-.5, resetCountSearch+.5, 
#endif
      timestampNBins, -timestampDiffSearch, timestampDiffSearch);
#endif
#ifdef __USE_CUTG__
  TFile* fcut = new TFile("/disk/proto-data/data/CUT_Compton.root");
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

  // read data file/flux and fill ccamData object
  std::cout << "Reading data\n";
  DecodeR* DR = new DecodeR(false); // Instantiates DecodeR object reading calorimeter data flux
  DecodeT* DT = new DecodeT(false); // Instantiates DecodeT object reading trigger data flux
  DecodeD* DD = new DecodeD(false); // Instantiates DecodeD object reading DSSSD(s) data flux
  newframe_t* event;
  DD -> setTree("/disk/proto-data/data/20210210_run1/bb7_3309-7_cs137-20210210_11h05_coinc_run1_conv.root");
  int dlen = DD -> getLength();

  int i = 0;// ROSMAP files loop counter
  int c = 0;// Event counter
  int cc = 0;
  // Set some constants
  const int pixelNumber = 64;
  //const int stripNumber = 32;

#ifdef __CIRCULAR_TREE__

  while (DD -> getCursor() < dlen)
  {
    DD -> Clear();
    DD -> decodeEvent();

    // Clear raw and physics data
    m_NPDetectorManager->ClearEventPhysics();
    m_NPDetectorManager->ClearEventData();

    // Fill data
    setCTTracker(ccamData, DD);// -> getEvent(), &nb_asic, &chain, stripNumber);

    // Build physical event
    m_NPDetectorManager->BuildPhysicalEvent();

    // Fill object in output ROOT file
    m_OutputTree->Fill();

    // check spectra
    m_NPDetectorManager->CheckSpectraServer();
  }
#else
  // Open data files
  ifstream iros, itrig;
  cout << "Loading data files ";

  itrig.open("/disk/proto-data/data/20210210_run1/mfm_trigger_202102101104.raw", ios::binary);
  itrig.seekg(0, ios::end);
  int tlen = itrig.tellg();
  itrig.seekg(0, ios::beg);
  char* tbuff = new char[tlen];
  itrig.read(tbuff, tlen);
  itrig.close();
  cout << "... ";

  iros.open("/disk/proto-data/data/20210210_run1/mfm_rdd_rosmap_04_mfm_rosmap_04_2021-02-10_10_04_59.raw.0001", ios::binary);
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
  while (not(DT->hasTrigged(2))) {
    DT -> decodeBlobMFM();
  }
  resetCount = DT->getResetCount() - resetCount;
  cout << "Found reset count: " << resetCount << endl;

#if defined __RESET_SEARCH__ && !defined __TEST_ZONE__
  // Fill control bidim
  for (int reset=-resetCountSearch; reset<resetCountSearch+1; reset++)
  {
  int cr = resetCount+reset;
  cout << "Filling coincidence bidim with reset #" << cr << endl;
  DD -> rewind();
#else
  int cr = resetCount;
#endif

  DR -> setRaw(rbuff);
  DR -> decodeBlobMFM();
  DD -> decodeEvent();

  int cd = 0;
  c = 0;
  i = 0;
  int tr = DR -> getTime();
  int td = DD -> getTime();
//  int dt = 100;
#if defined __RESET_SEARCH__
//  while(DR -> getCursor() < rlen and DD -> getCursor() < dlen)
  while(c < 1000)
  {
//    cout << DR -> getTime() << " " << DD -> getTime() << endl; 
    if (cr == cd) {
      if (abs(td-tr) < timestampDiffSearch) {
#ifdef __TEST_ZONE__
        bidim->Fill(c, td-tr);
#else
        bidim->Fill(reset, td-tr);
#endif
        c++;
      }
#else
  while(DR -> getCursor() < rlen and DD -> getCursor() < dlen)
  {
    if (cr == cd) {
#ifndef __TEST_ZONE__
      if (td-tr > 0 and td-tr < 110) { // That one is the real one
#else
      if (td-tr > -1000  and td-tr < 1000) {
#endif
        //DR -> Dump();
        //DD -> Dump();
        c++;
        cout << cc << " " << c << "(" << cr << ", " << cd << ") : " << tr << " " << td << endl;
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
            //cout << "c" << endl;
            cc++;
            m_OutputTree->Fill();
          }
#else
          m_OutputTree->Fill();
          cc++;
#endif
        }
        //cout << "c" << endl;

        // check spectra
        m_NPDetectorManager->CheckSpectraServer();
        //cout << "d" << endl;
      }
#endif
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
#if defined __RESET_SEARCH__ && !defined __TEST_ZONE__
  }
#endif

  delete DR;
  delete DT;
  delete DD;
  delete [] rbuff;
  delete [] tbuff;
#endif

  // fill spectra
  m_NPDetectorManager -> WriteSpectra();

#if defined __RESET_SEARCH__ || defined __
  fout->cd();
  bidim->Write();
  fout->Close();
#endif

  // Essential
#if __cplusplus > 199711L && NPMULTITHREADING
  m_NPDetectorManager->StopThread();
#endif
  RootOutput::Destroy();

  return 0;
}


/*
  const bool loopForever = false;
  //while (loopForever or i<12) // for Am data
  while (loopForever or i<3) // for Bi data quick analysis
  //while (loopForever or i<489) // for Bi data all events
  { 
    // Load a root file and setup DecodeD
    //DD -> setTree("../data/20200128_11h58_am241_conv.root");
    DD -> setTree("../data/bb7_3309-7_bi207_20210203_16h25_run8_conv.root");
    int dlength = DD -> getLength();

    //while (DD -> getCursor() < dlength and (loopForever or i<12))
    while (DD -> getCursor() < dlength and (loopForever or i<3))
    {
      // Load a ROSMAP file
      std::ifstream is;
//      i = 1;
      switch (i % 3) {
        case 3: is.open("./mfm.bin", std::ios::binary); break;
        case 4: is.open("./133Ba.bin", std::ios::binary); break;
        case 5: is.open("./241Am.bin", std::ios::binary); break;
        case 0: is.open("./207Bi.bin", std::ios::binary); break;
        case 1: is.open("./241Am-1.bin", std::ios::binary); break;
        case 2: is.open("./241Am-2.bin", std::ios::binary); break;
      }
      is.seekg (0, std::ios::end);
      int length = is.tellg();
      is.seekg (0, ios::beg);
      char* buffer = new char [length];
      is.read(buffer, length);
      is.close();
      i++;

      // Setup DecodeR
      DR -> setRaw(buffer);
      DR -> decodeRawMFM(); // get rid of the first two (empty) events
      DR -> decodeRawMFM();

      // Loop on some events
      while (DD -> getCursor() < dlength and DR -> getCursor() < length)
      {
        // Clear raw and physics data
        m_NPDetectorManager->ClearEventPhysics();
        m_NPDetectorManager->ClearEventData();

        // Fill calorimeter data
        DR -> decodeRawMFM();
        setCTCalorimeter(ccamData, DR, pixelNumber);

        // Fill DSSD data
        DD -> decodeEvent();
        event = DD -> getEvent();
        setCTTracker(ccamData, event, &nb_asic, &chain, stripNumber);

        // Build physical event
        m_NPDetectorManager->BuildPhysicalEvent();

        // Fill object in output ROOT file
        m_OutputTree->Fill();

        // check spectra
        m_NPDetectorManager->CheckSpectraServer();

        c++;
        //usleep(100);//Simulated 10kHz count rate

      }// End of loop on ROSMAP events

    delete [] buffer;
  
    }// End of loop on DSSSD data

  }// End of main loop

  // fill spectra
  m_NPDetectorManager->WriteSpectra();

  // delete Decoders
  delete DR;
  delete DD;

  // Essential
  #if __cplusplus > 199711L && NPMULTITHREADING
   m_NPDetectorManager->StopThread();
  #endif
  RootOutput::Destroy();

  return 0;
}

*/




/*
    while (DD -> getCursor() < dlength)
    {
      // Clear raw and physics data
      m_NPDetectorManager->ClearEventPhysics();
      m_NPDetectorManager->ClearEventData();

      //Read some data
      DD -> decodeEvent();
      event = DD -> getEvent();
      //Fill TComptonTelescopeData here (if possible)
      for (vector<int>::iterator itchain = chain.begin(); itchain != chain.end(); ++itchain) {//Iterates on 2 faces
        for (vector<int>::iterator itasic = nb_asic.begin(); itasic != nb_asic.end(); ++itasic) {//Iterates on 1 DSSSD
          if (event->chip_data[*itchain][*itasic]) { // Test if data is present
            switch (*itchain) {
              case 0://Assuming 0 is front - to be checked
                ccamData -> SetCTTrackerFrontTTowerNbr(1);
                ccamData -> SetCTTrackerFrontTDetectorNbr(*itasic+1);
                ccamData -> SetCTTrackerFrontTStripNbr(33);
                ccamData -> SetCTTrackerFrontTTime(event->timestamp);
                for (int k = 0; k < stripNumber; k++) {//Loop on strips
                  ccamData -> SetCTTrackerFrontETowerNbr(1);
                  ccamData -> SetCTTrackerFrontEDetectorNbr(*itasic+1);
                  ccamData -> SetCTTrackerFrontEStripNbr(k+1);
                  ccamData -> SetCTTrackerFrontEEnergy(event->sample[*itchain][*itasic][k]);
                }//End of loop on strips
                break;
              case 1://Assuming 1 is back - to be checked
                ccamData -> SetCTTrackerBackTTowerNbr(1);
                ccamData -> SetCTTrackerBackTDetectorNbr(*itasic+1);
                ccamData -> SetCTTrackerBackTStripNbr(33);
                ccamData -> SetCTTrackerBackTTime(event->timestamp);
                for (int k = 0; k < stripNumber; k++) {//Loop on strips
                  ccamData -> SetCTTrackerBackETowerNbr(1);
                  ccamData -> SetCTTrackerBackEDetectorNbr(*itasic+1);
                  ccamData -> SetCTTrackerBackEStripNbr(k+1);
                  ccamData -> SetCTTrackerBackEEnergy(event->sample[*itchain][*itasic][k]);
                }//End of loop on strips
            }//End switch
          }//End if
        }//End for
      }//End for

      // Build physical event
      m_NPDetectorManager->BuildPhysicalEvent();
      //Fill output tree
      m_OutputTree->Fill();
      //check spectra
      m_NPDetectorManager->CheckSpectraServer();

      c++;

    }//End while
    cout << "Read " << c << " DSSSD events from file" << endl;

    // Read from file
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
  
       // check spectra
       m_NPDetectorManager->CheckSpectraServer();

       c++;
       usleep(100);//Simulated 10kHz count rate
    }*/
