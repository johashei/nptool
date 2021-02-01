#ifndef DECODER_T
#define DECODER_T

//#include <algorithm>
//#include <vector>
//#include <fstream>
#include <iostream>
using namespace std;

#define _NDETECTORS 5

class DecodeT
{
  private:
    // raw is a pointer to the data flux, scanned by cursor
    char* raw;
    long int cursor;
    bool verbose;
    // The most usefull data
    long int time;
    bool coincidence;
    bool triggeredDetectors[_NDETECTORS];
    // Some other data
    int triggCount[_NDETECTORS];
    int resetCount;
    int coincWindow;
    // Utility methods
    long int combineBytes(int length);

  public:
    // Instanciation with all fields set to 0 (especially the cursor) except optionally the buffer dataBlock, and sets verbosity to v
    DecodeT(bool v);
    DecodeT(bool v, char* data);
    // Deletion deletes the data array mat
    ~DecodeT();
    // Setters and getters
    //void setTime(long int timestamp);
    //void setData(int* data);
    void setRaw(char* data);//also sets the cursor to 0
    void setVerbosity(bool v);
    long int getTime();
    char* getRaw();
    long int getCursor();
    void clear();// reset time and hasTrigged
    // Coincidence
    bool getCoinc();
    int getCoincW();
    // Triggers
    bool hasTrigged(int i);
    int getTriggCount(int i);
    int getResetCount();
    // Decode methods with and without MFM header (the second calls the first that fills the fields)
    void decodeRaw();
    void decodeBlobMFM();
    // A dumping method
    void Dump();
};

#endif

