// nptool headers
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "NPDetectorManager.h"
#include "TComptonTelescopeData.h"
#include "TComptonTelescopePhysics.h"

// root headers

// custom headers
#include "DecodeR.h"
#include "DecodeD.h"

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

//--//--// One-line setter for calorimeter //--//--//

void setCTCalorimeter(TComptonTelescopeData* ccamData, DecodeR* D, const int pixelNumber)
{
  ccamData -> SetCTCalorimeterTTowerNbr( 1 );
  ccamData -> SetCTCalorimeterTDetectorNbr( 1 );//Triggered ASIC number
  ccamData -> SetCTCalorimeterTChannelNbr( D -> getPixelNumber() );//Pixel that triggered
  ccamData -> SetCTCalorimeterTTime( D -> getTime() );
  for (int i = 0; i < pixelNumber; ++i) {//Loop on pixels
    ccamData -> SetCTCalorimeterETowerNbr(1);
    ccamData -> SetCTCalorimeterEDetectorNbr( 1 );
    ccamData -> SetCTCalorimeterEChannelNbr( i );//PMT pixel number
    ccamData -> SetCTCalorimeterEEnergy( D -> getData()[i] );
  }//End of loop on pixels
}

//--//--// One-line setter for DSSSD(s) //--//--//

// One-line setter for the Front of one DSSD
void setCTTrackerFront(TComptonTelescopeData* ccamData, newframe_t* event, int detNbr, int faceNbr, const int stripNumber)
{
  //detNbr and faceNbr are > 0
  ccamData -> SetCTTrackerFrontTTowerNbr(1);
  ccamData -> SetCTTrackerFrontTDetectorNbr(detNbr);
  ccamData -> SetCTTrackerFrontTStripNbr(33);
  ccamData -> SetCTTrackerFrontTTime(event->timestamp);
  for (int k = 0; k < stripNumber; k++) {//Loop on strips
    ccamData -> SetCTTrackerFrontETowerNbr(1);
    ccamData -> SetCTTrackerFrontEDetectorNbr(detNbr);
    ccamData -> SetCTTrackerFrontEStripNbr(k+1);
    ccamData -> SetCTTrackerFrontEEnergy(event->sample[detNbr-1][faceNbr-1][k]);
  }//End of loop on strips
}

// One-line setter for the Back of one DSSD
void setCTTrackerBack(TComptonTelescopeData* ccamData, newframe_t* event, int detNbr, int faceNbr, const int stripNumber)
{
  ccamData -> SetCTTrackerBackTTowerNbr(1);
  ccamData -> SetCTTrackerBackTDetectorNbr(detNbr);
  ccamData -> SetCTTrackerBackTStripNbr(33);
  ccamData -> SetCTTrackerBackTTime(event->timestamp);
  for (int k = 0; k < stripNumber; k++) {//Loop on strips
    ccamData -> SetCTTrackerBackETowerNbr(1);
    ccamData -> SetCTTrackerBackEDetectorNbr(detNbr);
    ccamData -> SetCTTrackerBackEStripNbr(k+1);
    ccamData -> SetCTTrackerBackEEnergy(event->sample[detNbr-1][faceNbr-1][k]);
  }//End of loop on strips
}

// One-line setter for DSSSD(s)
void setCTTracker(TComptonTelescopeData* ccamData, newframe_t* event, vector<int>* nb_asic, vector<int>* chain, const int stripNumber)
{
  for (vector<int>::iterator itchain = chain->begin(); itchain != chain->end(); ++itchain) {//Iterates on 2 faces
    for (vector<int>::iterator itasic = nb_asic->begin(); itasic != nb_asic->end(); ++itasic) {//Iterates on 1 DSSSD
      if (event->chip_data[*itchain][*itasic]) { // Test if data is present
        switch (*itchain) {
          case 0://Assuming 0 is front - to be checked
            setCTTrackerFront(ccamData, event, *itasic+1, *itchain+1, stripNumber);
            break;
          case 1://Assuming 1 is back - to be checked
            setCTTrackerBack(ccamData, event, *itasic+1, *itchain+1, stripNumber);
        }//End switch
      }//End if
    }//End for
  }//End for
}

//--//--//--//--//--//--//--//--//--//--//--//--//--//--//--//
int main()
{
  ///////////////////////////////////////////////////////////////////////////
  // configure option manager
//   NPOptionManager::getInstance()->Destroy();

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
  DecodeD* DD = new DecodeD(true); // Instantiates DecodeD object reading DSSSD(s) data flux
  newframe_t* event;

  //Sets where to look for data in DSSSD root frames
  vector <int> chain ({0, 1});//Two faces
  vector <int> nb_asic ({0});//First DSSSD
  
  int i = 0;// ROSMAP files loop counter
  int c = 0;// Event counter
  // Set some constants
  const int pixelNumber = 64;
  const int stripNumber = 32;
  const bool loopForever = true;
  while (loopForever or i<1)
  { 
    // Load a root file and setup DecodeD
    DD -> setTree("20200128_10h44_bi207_conv.root");
    int dlength = DD -> getLength();

    while (DD -> getCursor() < dlength and (loopForever or i<1))
    {
      // Load a ROSMAP file
      std::ifstream is;
      i = 1;
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
