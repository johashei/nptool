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

#endif

