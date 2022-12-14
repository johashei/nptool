#ifndef DECODER_H
#define DECODER_H

#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

#define _NPIXELS 64

class DecodeR
{
  private:
    // raw is a pointer to the data flux, scanned by cursor
    char* raw;
    long unsigned int cursor;
    bool verbose;
    // The most usefull data
    long unsigned int time;
    int* mat;
    // Some other data
    unsigned char sourceID;
    unsigned char channelID;
    unsigned char detNbr;
    // About the order of the pixels
    int F[_NPIXELS];
    void setF();//Order of pixels hard-coded here

    // Utility methods
    long unsigned int combineBytes(int length);
    void orderPixels();

  public:
    // Instanciation with all fields set to 0 (especially the cursor) except optionally the buffer dataBlock, and sets verbosity to v
    DecodeR(bool v);
    DecodeR(bool v, char* dataBlock);
    // Deletion deletes the data array mat
    ~DecodeR();
    // Setters and getters
    //void setTime(long int timestamp);
    //void setData(int* data);
    void setRaw(char* dataBlock);//also sets the cursor to 0
    void setVerbosity(bool v);
    long int getTime();
    int* getData();
    char* getRaw();
    long int getCursor();
    char getSource();//Deprecated
    char getChannel();//Deprecated
    char getPixelNumber();
    char getDetectorNumber();
    // Decode methods with and without MFM header (the second calls the first that fills the fields)
    void decodeRaw();
    void decodeBlobMFM();
    void decodeRawMFM();
    // A dumping method
    void Dump();
};

#endif

