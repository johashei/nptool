// nptool headers
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "NPDetectorManager.h"
#include "TComptonTelescopeData.h"
#include "TComptonTelescopePhysics.h"

// root headers

// Custom headers
//#include "DecodeR.h"
#include <iostream>
using namespace std;

#define _NPIXELS 64

class DecodeR
{
  private:
    // raw is a pointer to the data flux, scanned by cursor
    char* raw;
    long int cursor;
    // The most usefull data
    long int time;
    int* mat;
    // Some other data
    char sourceID;
    char channelID;
    // About the order of the pixels
    int F[_NPIXELS];
    void setF();//Order of pixels hard-coded here

  public:
    // Instanciation with all fields set to 0 (especially the cursor) except optionally the buffer dataBlock
    DecodeR();
    DecodeR(char* dataBlock);
    // Deletion deletes the data array mat
    ~DecodeR();
    // Setters and getters
    //void setTime(long int timestamp);
    //void setData(int* data);
    void setRaw(char* dataBlock);//also sets the cursor to 0
    long int getTime();
    int* getData();
    char* getRaw();
    long int getCursor();
    char getSource();
    char getChannel();
    // Utility method
    long int combineBytes(int length);
    void orderPixels();
    // Decode methods with and without MFM header (the second calls the first that fills the fields)
    void decodeRaw(bool verbose);
    void decodeRawMFM(bool verbose);
    // A dumping method
    void Dump();
};

DecodeR::DecodeR()
{
  raw = NULL;
    cursor = 0;
      time = 0;
        mat = new int [_NPIXELS];
        for (int i=0; i<_NPIXELS; i++) {
          mat[i] = 0;
        }
        setF();
}

DecodeR::DecodeR(char* dataBlock)
{
  raw = dataBlock;
  cursor = 0;
  time = 0;
  mat = new int [_NPIXELS];
  for (int i=0; i<_NPIXELS; i++) {
    mat[i] = 0;
  }
  setF();
}

DecodeR::~DecodeR()
{
  delete [] mat;
}

long int DecodeR::combineBytes(int length)
{ // Endianness fixed here. Also, probably not optimised for real time operations
  long int n = 0;
  for (int i = 0; i<length; i++) {
    n = n + (((unsigned char) raw[cursor+i]) << 8*(length-i-1));
  }
  cursor = cursor + length;
  return n;
}

void DecodeR::setF()
{
  int temp[_NPIXELS] = {31, 30, 29, 28, 23, 22, 21, 20, 15, 14, 13, 12,  7,  6,  5,  4,  0, 1,  2,  3,  8,  9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 63, 62, 61, 60, 55, 54, 53, 52, 47, 46, 45, 44, 39, 38, 37, 36, 32, 33, 34, 35, 40, 41, 42, 43, 48, 49, 50, 51, 56, 57, 58, 59};
  for (int i=0; i<_NPIXELS; i++) {
    F[i] = temp[i];
  }
}

void DecodeR::orderPixels()
{
  int data[_NPIXELS];
  for (int i=0; i<_NPIXELS; i++) {
    data[i] = mat[i];
  }
  for (int i=0; i<_NPIXELS; i++) {
    mat[F[i]] = data[i];
  }
}

void DecodeR::setRaw(char* dataBlock)
{
  raw = dataBlock;
  cursor = 0;
}

long int DecodeR::getTime()
{
  return time;
}

int* DecodeR::getData()
{
  return mat;
}

char* DecodeR::getRaw()
{
  return raw;
}

long int DecodeR::getCursor()
{
  return cursor;
}

char DecodeR::getSource()
{
  return sourceID;
}

char DecodeR::getChannel()
{
  return channelID;
}

void DecodeR::decodeRaw(bool verbose)
{
  const bool packetType0xd4Autorized = true;
  // Version and system number
  char version = raw[cursor]/32;
  char sysNumber = raw[cursor]%32;
  cursor++;
  // Packet type
  unsigned char packetType = raw[cursor];
  cursor++;
  if (verbose)
  {cout << "packetType: " << int(packetType) << endl;}
  // Packet sequence
  int seqFlag = raw[cursor]/64;
  int packetCount = combineBytes(2)%16384;
  if (verbose)
  {cout << "seqFlag: " << seqFlag << endl << "packetN: " << packetCount << endl;}
  // Timestamp
  time = combineBytes(4);
  // Data length
  int dataLength = combineBytes(2);
  if (verbose)
  {cout << int(version) << endl << int(sysNumber) << endl << int(packetType) << endl << seqFlag << endl << packetCount << endl << time << endl << dataLength << endl;}
  if (packetType == 0xd5) {
    // Source ID
    sourceID = raw[cursor];
    cursor++;
    // TriggerType
    char triggerType = raw[cursor];
    cursor++;
    // ChannelID
    channelID = raw[cursor];
    cursor++;
    // Hold Delay
    int holdDelay = combineBytes(2);
    // Number of samples
    int nofSamples = combineBytes(2);
    long int cursorAtDataBeginning = cursor;
    if (nofSamples == _NPIXELS) {
      for (int i=0; i<_NPIXELS; i++) {
        mat[i] = combineBytes(2);//In channel/asic order - not fitted for NN use !!
      }
      orderPixels();// Now it should be better
    } else
    { cout << "Number of samples " << nofSamples << " do not match the number of pixels (" << _NPIXELS << ") -> data not filled." << endl; }
    if (verbose) {
      cursor = cursorAtDataBeginning;
      cout << nofSamples << endl;
      for (int i = 0; i<nofSamples; i++) {
        cout << combineBytes(2) << " ";
      }
      cout << endl;
    }
  } else if (packetType == 0xd4 && packetType0xd4Autorized) {//Ex: reading pedestals
    // Number of events in frame
    int nofEvents = raw[cursor];
    cursor++;
    // Number of samples per event
    int nofSamples = combineBytes(2);
    if (verbose)
    { cout << "Number of events: " << nofEvents << endl << "Number of samples " << nofSamples << endl; }
    // Going through events
    for (int i=0; i<nofEvents; i++) {
      //Event timestamp
      int timestamp = combineBytes(4);
      // Trigger Type
      char triggerType = raw[cursor];
      cursor++;
      // Source ID
      sourceID = raw[cursor];
      cursor++;
      // Channel ID
      channelID = raw[cursor];
      cursor++;
      if (verbose)
      { cout << timestamp << endl << int(triggerType) << endl << int(sourceID) << " " << int(channelID) << endl; }
      // Going through samples
      for (int j=0; j<nofSamples; j++) {
        mat[j] = combineBytes(2);
      }
    }
  } else
  { cout << "Packet Type is not 0xd5 nor 0xd4 -> not analysed." << endl; }
}

void DecodeR::decodeRawMFM(bool verbose)
{
  if (raw) {
    if (verbose)
    { cout << "Decoding in progress ..." << endl; }
    long int cursorAtBeginning = cursor;
    // metaType byte
    char metaType = raw[cursor];//Some info can be gathered here
    //int unitBlockSize = pow(2, metaType%16);
    int unitBlockSize = 4;
    cursor++;
    // Size of frame
    long int frameSize = combineBytes(3);
    // Data source
    char dataSource = raw[cursor];
    cursor++;
    // Frame Type
    int frameType = combineBytes(2);
    // Revision
    char revision = raw[cursor];
    cursor++;
    // Header Size
    int headerSize = combineBytes(2);
    // Item Size
    int itemSize = combineBytes(2);
    // Number of items
    long int nItems = combineBytes(4);
    // Event number
    long int eventN = combineBytes(4);
    if (verbose)
    {cout << int(metaType) << endl <<frameSize << endl << int(dataSource) << endl << frameType << endl << int(revision) << endl << headerSize << endl << itemSize << endl << nItems << endl << eventN << endl;}
    // Giving ROSMAP data frame to decodeRaw if event number > 10
    if (eventN > 10) {
      // Going through other asics overhead (could be generalised)
      if (verbose)
      { cout << endl << "Curseur: " << cursor << endl; }
      cursor = cursorAtBeginning + (headerSize+7)*unitBlockSize;
      if (verbose)
      { cout << endl << "Curseur: " << cursor << endl; }
      // Data block length
      long int dataBlockLength = combineBytes(4);
      // Reading the actual data
      decodeRaw(verbose);
      if (verbose)
      { cout << endl << "Curseur: " << cursor << endl; }
      // Rosmap adress in ASCII
      int adressLength = raw[cursor];
      cursor++;
      char* rosmapAdress = new char [adressLength];
      for (int i=0; i<adressLength; i++)
      { rosmapAdress[i]=raw[cursor+i]; }
      if (verbose) {
        cout << "ROSMAP adress: ";
        for (int i=0; i<adressLength; i++)
        { cout << rosmapAdress[i]; }
        cout << endl;
      }
      cursor = cursor + adressLength;
      // PC Timestamp
      int timestampSize = raw[cursor];
      cursor++;
      long int pctimestamp = combineBytes(timestampSize);
      if (verbose)
      {cout << "PC timestamp: " <<  pctimestamp << endl;}
      cout << "Event number : " << eventN << " done." << endl;
    } else {
      // Event not good
      cout << "Event number : " << eventN << " -> not done." << endl;
    }
    //Skipping data reserve and ending block and going to the next event (or the end of the file)
    cursor = cursorAtBeginning + frameSize*unitBlockSize;
  } else
  { cout << "Error: Raw not set" << endl; }
}

void DecodeR::Dump()
{
  if (raw) {
    cout << "Raw adress: " << raw << endl;
    cout << "Raw value: " << raw[cursor] << endl;
  } else
  { cout << "Raw not set" << endl; }
  cout << "Time: " << time << endl;
  cout << "Data: " << endl;
  for (int i=0; i<_NPIXELS; i++) {
    if ( (i+1) % 8 == 0 ) {
      cout << mat[i] << endl;
    } else {
      cout << mat[i] << "\t";
    }
  }
}


// C++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

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
   auto ccamData = new TComptonTelescopeData();
   // connect raw CCAM data pointer to physics class
/*   auto ccamPhys = (TComptonTelescopePhysics*) m_NPDetectorManager->GetDetector("ComptonTelescope");
   ccamPhys->SetRawDataPointer(ccamData);*/
   
   // read data file/flux and fill ccamData object

   cout << "Reading data" << endl;
   DecodeR* D = new DecodeR();
   bool verbose = false;
   
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
}
