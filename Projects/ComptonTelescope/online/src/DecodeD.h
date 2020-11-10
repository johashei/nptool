#ifndef DECODED_H
#define DECODED_H

// General C++ librairies
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

// Root headers
#include "TFile.h"
#include "TTree.h"


class DecodeD
{
  private:
    bool verbose;

  public:
    DecodeD(bool v);
    ~DecodeD();
};

#endif
