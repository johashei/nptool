#include "DecodeT.h"

DecodeT::DecodeT(bool v)
{
  raw = NULL;
  cursor = 0;
  verbose = v;
  clear();
  coincidence = false;
  for (int i = 0; i < _NDETECTORS; i++) {
    triggCount[i] = 0;
  }
  resetCount = 0;
  coincWindow = 0;
}

DecodeT::DecodeT(bool v, char* data)
{
  raw = data;
  cursor = 0;
  verbose = v;
  clear();
  coincidence = false;
  for (int i = 0; i < _NDETECTORS; i++) {
    triggCount[i] = 0;
  }
  resetCount = 0;
  coincWindow = 0;
}

DecodeT::~DecodeT()
{
//  delete [] triggeredDetectors;
//  delete [] triggCount;
}

long int DecodeT::combineBytes(int length)
{
  long int n = 0;
  for (int i = 0; i<length; i++) {
    n = n + (((unsigned char) raw[cursor+i]) << 8*(length-i-1));
  }
  cursor = cursor + length;
  return n;
}

void DecodeT::setRaw(char* data)
{
  raw = data;
  cursor = 0;
}

void DecodeT::setVerbosity(bool v)
{
  verbose = v;
}

long int DecodeT::getTime()
{
  return time;
}

char* DecodeT::getRaw()
{
  return raw;
}

long int DecodeT::getCursor()
{
  return cursor;
}

void DecodeT::clear()
{
  time = 0;
  for (int i = 0; i < _NDETECTORS; i++) {
    triggeredDetectors[i] = false;
  }
}

bool DecodeT::getCoinc()
{
  return coincidence;
}

int DecodeT::getCoincW()
{
  return coincWindow;
}

bool DecodeT::hasTrigged(int i)
{
  if (i >= 0 && i < _NDETECTORS) {
    return triggeredDetectors[i];
  } else {
    return false;
  }
}

int DecodeT::getTriggCount(int i)
{
  if (i >= 0 && i < _NDETECTORS) {
    return triggCount[i];
  } else {
    return 0;
  }
}

int DecodeT::getResetCount()
{
  return resetCount;
}

void DecodeT::decodeRaw()
{
//  clear();
  // Coincidence mode
  coincidence = raw[cursor] == 1;
  cursor ++;
  // Triggered detectors
  for (int i = 0; i < _NDETECTORS; i++) {
    triggeredDetectors[i] = raw[cursor+i] == 1;
  }
  cursor += _NDETECTORS;
  // Timestamp
  time = combineBytes(6);
  // Coincidence window
  coincWindow = combineBytes(2);
  // Trigger counters
  for (int i = 0; i < _NDETECTORS; i++) {
    triggCount[i] = combineBytes(3);
  }
  // Reset Counter
  resetCount = combineBytes(3);
  // Unused fields
  cursor += 20;
  if (verbose) {
    cout << coincidence << " " << time << " " << coincWindow << " " << resetCount << endl;
    for (int i = 0; i < _NDETECTORS; i++) 
      { cout << triggeredDetectors[i] << "\t"; }
    cout << endl;
    for (int i = 0; i < _NDETECTORS; i++)
      { cout << triggCount[i] << "\t"; }
    cout << endl;
  }
}

void DecodeT::decodeBlobMFM()
{
//  clear();
  if (raw) {
    if (verbose)
      {cout << "Decoding trigger data (blob) in progress ..." << endl;}
    long int cursorAtBeginning = cursor;
    // metaType byte
    char metaType = raw[cursor];
    cursor++;
    // Size of frame
    long int frameSize = combineBytes(3);
    // Data source
    char detNbr = raw[cursor];
    cursor++;
    // Frame Type
    int frameType = combineBytes(2);
    // Revision
    char revision = raw[cursor];
    cursor++;
    if (verbose)
      {cout << int(metaType) << endl << frameSize << endl << int(detNbr) << endl << frameType << endl << int(revision) << endl; }
    // Reading the actual data
    if (verbose)
      { cout << endl << "Curseur: " << cursor << endl; }
    decodeRaw();
    if (verbose)
      { cout << endl << "Curseur: " << cursor << endl; }
    cursor = cursorAtBeginning + frameSize;
  }
}

void DecodeT::Dump()
{
  cout << "Cursor:           " << cursor << endl;
  cout << "Timestamp:        " << time << " is " << double(time)/1e8 << " s." << endl;
  cout << "Coincidence mode: " << coincidence << endl;
  cout << "Coincidence wind: " << coincWindow << endl;
  cout << "Reset count:      " << resetCount << endl;
  cout << "Triggered detectors: " << endl;
  for (int i = 0; i < _NDETECTORS; i++) {
    cout << i << ": " << triggeredDetectors[i] << " " << triggCount[i] << endl;
  }
}




























