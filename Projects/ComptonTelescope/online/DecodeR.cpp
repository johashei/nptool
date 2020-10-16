#include "DecodeR.h"

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

/*void DecodeR::setTime(long int timestamp)
{
  time = timestamp;
}

void DecodeR::setData(int* data)
{
  mat = data;
}*/

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
    { cout << "Error: Raw not set" << endl;
  }
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
