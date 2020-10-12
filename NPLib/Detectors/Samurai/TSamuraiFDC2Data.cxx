#include "TSamuraiFDC2Data.h"
#include <iostream>

TSamuraiFDC2Data::TSamuraiFDC2Data(){};
TSamuraiFDC2Data::~TSamuraiFDC2Data(){};

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Data::SetData(const int& Det, const int& Layer, const int& Wire, const double& Time, const int& Edge){
  fDC_DetectorNbr.push_back(Det);
  fDC_LayerNbr.push_back(Layer); 
  fDC_WireNbr.push_back(Wire); 
  fDC_Time.push_back(Time); 
  fDC_Edge.push_back(Edge); 
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Data::Clear(){
  fDC_DetectorNbr.clear();
  fDC_LayerNbr.clear();
  fDC_WireNbr.clear();
  fDC_Time.clear();
  fDC_Edge.clear();
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Data::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  cout << "   - Multiplicity: " << Mult() << endl;

}
////////////////////////////////////////////////////////////////////////////////
unsigned int TSamuraiFDC2Data::MultLayer(unsigned int det , unsigned int layer, int edge){
  unsigned int mult=0;
  unsigned int size = fDC_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fDC_DetectorNbr[i]==det)
      if(fDC_LayerNbr[i]==layer)
        if(fDC_Edge[i]==edge && edge!=-1) // edge type is specified (0 or 1)
          mult++;
        else if(edge==-1)// edge type is not specified
          mult++;
  }
  return mult;

}
////////////////////////////////////////////////////////////////////////////////
std::vector<int> TSamuraiFDC2Data::GetWire(unsigned int det , unsigned int layer){
  std::vector<int> wires;
  unsigned int size = fDC_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fDC_DetectorNbr[i]==det)
      if(fDC_LayerNbr[i]==layer)
        wires.push_back(fDC_WireNbr[i]);
  }
  return wires;
}
ClassImp(TSamuraiFDC2Data); 
