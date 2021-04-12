#include "TBigRIPSPPACData.h"
#include <iostream>

TBigRIPSPPACData::TBigRIPSPPACData(){};
TBigRIPSPPACData::~TBigRIPSPPACData(){};

/*
void TBigRIPSPPACData::SetData(const int& FP, const int& ID, const int& Edge, const int& value,const int& VariableType){
  fPPAC_FP.push_back(FP);
  fPPAC_ID.push_back(ID);
  fPPAC_Edge.push_back(Edge);
  switch(VariableType){
      case 0: fPPAC_TX1.push_back(value); break; 
      case 1: fPPAC_TX2.push_back(value); break; 
      case 2: fPPAC_TY1.push_back(value); break; 
      case 3: fPPAC_TY2.push_back(value); break; 
      case 4: fPPAC_TA.push_back(value); break; 
      case 5: fPPAC_QX1.push_back(value); break; 
      case 6: fPPAC_QX2.push_back(value); break; 
      case 7: fPPAC_QY1.push_back(value); break; 
      case 8: fPPAC_QY2.push_back(value); break; 
      case 9: fPPAC_QA.push_back(value); break; 
  }
}
*/
////////////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACData::Clear(){
 // fPPAC_FP.clear();
 // fPPAC_ID.clear();
 // fPPAC_Edge.clear();
  fPPAC_TX1.clear(); 
  fPPAC_TX2.clear(); 
  fPPAC_TY1.clear(); 
  fPPAC_TY2.clear();
  fPPAC_TA.clear(); 
  fPPAC_QX1.clear(); 
  fPPAC_QX2.clear(); 
  fPPAC_QY1.clear(); 
  fPPAC_QY2.clear();
  fPPAC_QA.clear(); 

  fPPAC_TX1_ID.clear(); 
  fPPAC_TX2_ID.clear(); 
  fPPAC_TY1_ID.clear(); 
  fPPAC_TY2_ID.clear();
  fPPAC_TA_ID.clear(); 
  fPPAC_QX1_ID.clear(); 
  fPPAC_QX2_ID.clear(); 
  fPPAC_QY1_ID.clear(); 
  fPPAC_QY2_ID.clear();
  fPPAC_QA_ID.clear(); 
}

////////////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACData::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  //cout << "   - Multiplicity: " << Mult() << endl;

}
////////////////////////////////////////////////////////////////////////////////
/*
unsigned int TBigRIPSPPACData::MultLayer(unsigned int det , unsigned int layer, int edge){
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
std::vector<int> TBigRIPSPPACData::GetWire(unsigned int det , unsigned int layer){
  std::vector<int> wires;
  unsigned int size = fDC_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fDC_DetectorNbr[i]==det)
      if(fDC_LayerNbr[i]==layer)
        wires.push_back(fDC_WireNbr[i]);
  }
  return wires;
}
*/
ClassImp(TBigRIPSPPACData); 
