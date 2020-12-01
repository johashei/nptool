/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Cyril Lenain   contact address: lenain@lpcaen.in2p3.fr *
 *                                                                           *
 * Creation Date   : November 2020                                           *
 * Last update     : November 2020                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class find multi-neutrons in an neutron detection array by          *
 *    rejecting "crosstalk"                                                  *
 *****************************************************************************/

#include"NPCrossTalk.h"
#include "Math/Factory.h"
#include "TError.h"
#include "TGraph.h"

#include <iostream>
using namespace std;
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
CrossTalk::CrossTalk(){
  // Radius of a cluster = Coef * dXdYdZ
  coef = 3;
  // this avoid error 
  gErrorIgnoreLevel = kError;
  /* FirstHit = -1; */

}
////////////////////////////////////////////////////////////////////////////////
CrossTalk::~CrossTalk(){
}

////////////////////////////////////////////////////////////////////////////////

void CrossTalk::AddHitVector(const vector<double>& X, const vector<double>& Y, const vector<double>& Z,const vector<double>& dX, const vector<double>& dY, const vector<double>& dZ, const vector<double> &T){

  HitX = &X;
  HitY = &Y;
  HitZ = &Z;

  HitdX = &dX;
  HitdY = &dY;
  HitdZ = &dZ;

  HitT = &T;
  sizeHit = X.size();
}

vector<double> CrossTalk::ComputeCrossTalk(){

  FirstHit = -1;
  // take back the first hit, automaticaly identified as the first Neutron -> probably useless with this new method
  for(int i = 0; i < sizeHit; i++){
    static double minT;
    minT = 1000000;
    if((*HitT)[i] < minT){
      minT = (*HitT)[i];
      FirstHit = i;
    }
  }

  static double x1,y1,z1,dx1,dy1,dz1,t1;
  static double x2,y2,z2,dx2,dy2,dz2,t2;
  static double Dist, dR1, dR2, v1;

  // A different Cluster number (starting at 1) is assigned to each hit
  static vector<int> ID_ClustHit;
  ID_ClustHit.clear();
  for(int i =0 ; i < sizeHit; i++){
    ID_ClustHit.push_back(i+1);
  }
  
  /*
  // Identify Hits around the 1st neutron, 1st cluster
  x1 = (*HitX)[FirstHit], y1 = (*HitY)[FirstHit], z1 = (*HitZ)[FirstHit], dx1 = (*HitdX)[FirstHit], dy1 = (*HitdY)[FirstHit], dz1 = (*HitdZ)[FirstHit], t1 = (*HitT)[FirstHit];   
  dR1 = sqrt(x1*dx1 + dy1*dy1 + dz1*dz1);
  v1 = sqrt(x1*x1 + y1*y1 + z1*z1)/minT; 

  for(int i = 0; i < sizeHit; i++){
    if(i == FirstHit){
      ID_ClustHit.push_back(1);
      continue;
    }
    x2 = (*HitX)[i], y2 = (*HitY)[i], z2 = (*HitZ)[i], dx2 = (*HitdX)[i], dy2 = (*HitdX)[i], dz2 = (*HitdZ)[i], t2 = (*HitT)[i];   
    Dist = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
    if(Dist < coef*dR1){ //coef*dR1 : Radius of the cluster
      ID_ClustHit.push_back(1); // this hit is part of the cluster 1
    }
    else(){
      ID_ClustHit.push_back(-1);
    }
  }
  */
  
  //Test each Dist(n-n) to find clusters
  //When 2 neutrons are part of a the same cluster, change their ID_ClustHit for the lowest
  for(int j = 0; j < sizeHit-1; j++){
    x1 = (*HitX)[j], y1 = (*HitY)[j], z1 = (*HitZ)[j], dx1 = (*HitdX)[j], dy1 = (*HitdY)[j], dz1 = (*HitdZ)[j], t1 = (*HitT)[j];   
    dR1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1);  
    for(int jj = j+1; jj < sizeHit ; jj++){
      x2 = (*HitX)[jj], y2 = (*HitY)[jj], z2 = (*HitZ)[jj];   
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      if(Dist < coef*dR1){
        if(ID_ClustHit[jj] < ID_ClustHit[j]){
          ID_ClustHit[j] == ID_ClustHit[jj]; 
        }
        else if(ID_ClustHit[jj] > ID_ClustHit[j]){
          ID_ClustHit[jj] == ID_ClustHit[j];
        }
      }
    }
  }
  //Create a map to hold the Hit_ID for each cluster 
  static map<unsigned int, vector<unsigned int>> mapOfClust;
  mapOfClust.clear();
  for(unsigned int i = 0; i < sizeHit; i++){
    mapOfClust[ID_ClustHit[i]].push_back(i);
  }
  
  //Put first hit in time of each cluster in a new map mapOfFirstN
  //And also find the first hit of all : FirstN
  static map< unsigned int, unsigned int > mapOfFirstN;
  static unsigned int NbrOfClust, FirstN;
  NbrOfClust = 0;
  NbrOfClust = mapOfClust.size();
  static double GeneralTmin; 
  GeneralTmin = 1000000;
  for(int i = 1; i <= NbrOfClust; i++){
    static unsigned int clustSize;
    static double TminClust;
    TminClust = 10000000;
    clustSize = mapOfClust[i].size();
    for(int j = 0; j < clustSize; j++){
      if((*HitT)[mapOfClust[i][j]] < TminClust){
        TminClust = (*HitT)[mapOfClust[i][j]];
        mapOfFirstN[i] = mapOfClust[i][j];
      }
      if((*HitT)[mapOfClust[i][j]] < GeneralTmin){
        GeneralTmin = (*HitT)[mapOfClust[i][j]];
        FirstN = mapOfClust[i][j];
      }
    }
  }
   
  FirstHit = FirstN;

  static vector<double> m_Neutrons;
  m_Neutrons.clear();
  return m_Neutrons;
}

int CrossTalk::GetFirstN(){
  return FirstHit;
}

