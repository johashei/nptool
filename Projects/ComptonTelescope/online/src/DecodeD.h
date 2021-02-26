#ifndef DECODED_H
#define DECODED_H

// General C++ librairies
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

// Root headers
#include "TFile.h"
#include "TTree.h"

#define NBDETECTORS 1
#define NBSTRIPS  32

enum Datatype {D_NONE = 0, D_ROOT, D_MFM};

typedef struct {
  uint8_t chain;            // Chain Number (0 or 1)
  uint8_t nb_asic;          // ASIC Number
  uint8_t chip_data;        // 0=Empty Frame;
  uint8_t analog_trigger;
  uint8_t seu;
  uint32_t ch_status;
  uint16_t ref_channel;
  uint16_t sample[32];
  uint16_t cm_data;
  uint32_t timestamp;
} frame_t;

typedef struct
{
  uint8_t chip_data[3][8];
  uint8_t analog_trigger[3][8];
  uint8_t seu[3][8];
  uint32_t ch_status[3][8];
  uint16_t ref_channel[3][8];
  uint16_t sample[3][8][32];
  uint16_t cm_data[3][8];
  uint32_t timestamp;
} newframe_t;

class DecodeD
{
  private:
    bool verbose;
    Datatype datatype;
    long int cursor;
    newframe_t event;

    // For root data
    TTree* t1;
    long int length;
    uint8_t** chip_data;
    uint8_t** analog_trigger;
    uint8_t** seu;
    uint32_t** ch_status;
    uint16_t** ref_channel;
    uint16_t*** sample;
    uint16_t** cm_data;
    uint32_t* timestamp;

    // For online data
    vector<int> FaceType;
    vector<int> DetNbr;
    vector<int> StripNbr;
    vector<double> Energy;

  public:
    DecodeD(bool v);
    ~DecodeD();

    void setTree(const char* filename);
    void setRaw();
    void rewind();

    long int getCursor();
    long int getLength();
    int getEventSize();
    int getFaceType(const int i);
    int getDetNbr(const int i);
    int getStripNbr(const int i);
    double getEnergy(const int i);
    long int getTime();
    newframe_t* getEvent();
    // One may add a few getters here and deprecate getEvent to avoid requiring the class user to know the newframe_t struct
    
    void decodeEvent();
    void decodeEventFinal();

    void Clear();
    void Dump();
};

#endif

