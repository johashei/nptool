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

bool cmp(pair<int,double>& a, pair<int,double>&b){
    return a.second < b.second;
}

vector<int> CrossTalk::ComputeCrossTalk(){

  FirstHit = -1;

  static double x1,y1,z1,dx1,dy1,dz1,t1;
  static double x2,y2,z2,dx2,dy2,dz2,t2;
  static double Dist, dR1, dR2;

  //Attribute a new index based on the time arrive from the first to the last hit
  //Using vector of pairs (ID,Time)
  static vector<pair<int, double>> pairSortedID;
  pairSortedID.clear();
  for(int i = 0; i < sizeHit; i++){
    pairSortedID.emplace_back(i, (*HitT)[i]);
  }
  //Sort pair vector (ID,Time) in Time
  sort(pairSortedID.begin(), pairSortedID.end(), cmp);
  static vector<unsigned int> SortedID;
  SortedID.clear();
  //Put new ID sorted in a vector
  for(int i = 0; i < sizeHit; i++){
    SortedID.push_back(pairSortedID[i].first);
  }

  // A different Cluster number (starting at 1) is assigned to each hit
  static vector<int> ID_ClustHit;
  ID_ClustHit.clear();
  for(int i =0 ; i < sizeHit; i++){
    ID_ClustHit.push_back(i+1);
  }
  
  FirstHit = SortedID[0];
  
  //Test each Dist(n-n) to find clusters
  //When 2 neutrons are part of a the same cluster, change their ID_ClustHit for the lowest
  //Create a map to hold the Hit_ID for each cluster 
  static map<unsigned int, vector<unsigned int>> mapOfClust;
  mapOfClust.clear();
  for(int j = 0; j < sizeHit-1; j++){
    x1 = (*HitX)[SortedID[j]], y1 = (*HitY)[SortedID[j]], z1 = (*HitZ)[SortedID[j]], dx1 = (*HitdX)[SortedID[j]], dy1 = (*HitdY)[SortedID[j]], dz1 = (*HitdZ)[SortedID[j]], t1 = (*HitT)[SortedID[j]];   
    dR1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1);  
    for(int jj = j+1; jj < sizeHit ; jj++){
      x2 = (*HitX)[SortedID[jj]], y2 = (*HitY)[SortedID[jj]], z2 = (*HitZ)[SortedID[jj]];   
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      if(Dist < coef*dR1){
        if(ID_ClustHit[jj] < ID_ClustHit[j]){
          ID_ClustHit[j] = ID_ClustHit[jj]; 
        }
        else if(ID_ClustHit[jj] > ID_ClustHit[j]){
          ID_ClustHit[jj] = ID_ClustHit[j];
        }
      }
    }
  }

  for(unsigned int i = 0; i < sizeHit; i++){
    mapOfClust[ID_ClustHit[i]].push_back(SortedID[i]);
  }

  //Put first hit in time of each cluster in a new map mapOfHead and remake the numbering of clusters (1,2,3...)
  static map< unsigned int, unsigned int > mapOfHead;
  static unsigned int NbrOfClust;
  mapOfHead.clear();
  NbrOfClust = mapOfClust.size();
 
  static unsigned int count;
  count=1;
  for(auto itr = mapOfClust.begin(); itr != mapOfClust.end(); itr++){
    mapOfHead[count] = itr->second[0];
    count++;
  } 

  static double v_n, Dmax;
  static vector<double> CrossTalk;
  CrossTalk.clear();

  //Test now the "Friend" clusters
  for(int i = 1; i < NbrOfClust-1; i++){ 
    x1 = (*HitX)[mapOfHead[i]], y1 = (*HitY)[mapOfHead[i]], z1 = (*HitZ)[mapOfHead[i]], t1 = (*HitT)[mapOfHead[i]]; 
    v_n = sqrt(x1*x1+y1*y1+z1*z1)/t1;
    for(int j = i+1; i < NbrOfClust; i++){
      x2 = (*HitX)[mapOfHead[i]], y2 = (*HitY)[mapOfHead[i]], z2 = (*HitZ)[mapOfHead[i]], t2 = (*HitT)[mapOfHead[i]]; 
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      Dmax = (t2-t1)*v_n;
      if(Dist < Dmax){
        CrossTalk.push_back(j);
      }
    }
  }

  //elimate potential CrossTalk in mapOfHead
  static unsigned int sizeCT;
  sizeCT = CrossTalk.size();
  for(unsigned int i = 0; i < sizeCT; i++  ){
    mapOfHead.erase(CrossTalk[i]);
  }
  //return vector of real neutrons IDs 
  m_Neutrons.clear();
  for(auto itr = mapOfHead.begin(); itr != mapOfHead.end(); itr++){
    m_Neutrons.push_back(itr->second);
  }

  return m_Neutrons;

}

int CrossTalk::GetFirstN(){
  return FirstHit;
}

