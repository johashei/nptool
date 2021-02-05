#include "DecodeD.h"

DecodeD::DecodeD(bool v)
{
  verbose = v;
  cursor = 0;
  datatype = D_NONE;
  t1 = 0;
  length = 0;
}

DecodeD::~DecodeD()
{
  switch (datatype) {
    case D_ROOT:
      delete t1;
      break;
    case D_MFM:
      break;//TBR
  }
}

void DecodeD::setTree(const char* filename)
{
  datatype = D_ROOT;
  cursor = 0;
  TFile* f = new TFile(filename);
  t1 = (TTree*)f->Get("Events");
  length = t1 -> GetEntries();
  if (verbose) { cout << "Length of loaded tree: " << length << endl; }
  t1->SetBranchAddress("chip_data", &event.chip_data);
  t1->SetBranchAddress("analog_trigger", &event.analog_trigger);
  t1->SetBranchAddress("seu", &event.seu);
  t1->SetBranchAddress("ch_status", &event.ch_status);
  t1->SetBranchAddress("ref_channel", &event.ref_channel);
  t1->SetBranchAddress("sample", &event.sample);
  t1->SetBranchAddress("cm_data", &event.cm_data);
  t1->SetBranchAddress("timestamp", &event.timestamp);
}

void DecodeD::setRaw()
{
  datatype = D_MFM;
}

void DecodeD::rewind()
{
  cursor = 0;
}

long int DecodeD::getCursor()
{
  return cursor;
}

long int DecodeD::getLength()
{
  return length;
}

long int DecodeD::getTime()
{
  return event.timestamp;
}

newframe_t* DecodeD::getEvent()
{
  return &event;
}

void DecodeD::decodeEvent()
{
  switch (datatype) {
    case D_ROOT:
      if (cursor < length) {
        t1->GetEntry(cursor);
//        for (int i = 0; i < 3; i++) {
  //        for (int j = 0; i < 8; j++) {
/*  int i = 0; int j = 0;
            event.chip_data[i][j] = chip_data[i][j];
            event.analog_trigger[i][j] = analog_trigger[i][j];
            event.seu[i][j] = seu[i][j];
            event.ch_status[i][j] = ch_status[i][j];
            event.ref_channel[i][j] = ref_channel[i][j];
            event.cm_data[i][j] = cm_data[i][j];
            for (int k = 0; k<32; k++) {
              event.sample[i][j][k] = sample[i][j][k];
            }
//          }
  //      }
        event.timestamp = *timestamp;*/
        cursor++;
      }
      break;
    case D_MFM:
      break;
    case D_NONE:
      cout << "No data has been set to decode" << endl;
  }
}

void DecodeD::Dump()
{
  cout << "Datatype: " << datatype << " means ";
  switch (datatype) {
    case D_ROOT:
      cout << "ROOT" << endl;
      break;
    case D_MFM:
      cout << "MFM" <<  endl;
      break;
    case D_NONE:
      cout << "no data has been set" << endl;
      return;
  }
  cout << "Timestamp: " << event.timestamp << endl;
  cout << "Chip data\tanalog trigger\tseu\tchannel status\tref channel\tcm data" << endl;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 8; j++) {
      cout << (int)event.chip_data[i][j] << "\t" << (int)event.analog_trigger[i][j] << "\t" << (int)event.seu[i][j] << "\t" << (int)event.ch_status[i][j] << "\t" << (int)event.ref_channel[i][j] << "\t" << (int)event.cm_data[i][j] << endl;
      cout << "Samples: ";
      for (int k = 0; k<32; k++) {
        cout << event.sample[i][j][k] << " ";
      }
      cout << endl;
    }
  }
}



