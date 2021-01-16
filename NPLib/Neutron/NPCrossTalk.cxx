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
#include "TPad.h"
#include "TGraph.h"

#include <iostream>
using namespace std;
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
CrossTalk::CrossTalk(){
  // Radius of a cluster = Coef * dXdYdZ
  coef = 5;
  // this avoid error 
  gErrorIgnoreLevel = kError;

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

  const double C_light = 299.792458;
  static double x1,y1,z1,dx1,dy1,dz1,t1;
  static double x2,y2,z2,dx2,dy2,dz2,t2;
  static double Dist, dR1, dR2;

  //Attribute a new index based on the time arrive from the first to the last hit
  //Using vector of pairs (ID,Time)
  static vector<pair<int, double>> pair_SortedID;
  pair_SortedID.clear();
  for(int i = 0; i < sizeHit; i++){
    pair_SortedID.emplace_back(i, (*HitT)[i]);
  }
  //Sort pair vector (ID,Time) in Time
  sort(pair_SortedID.begin(), pair_SortedID.end(), cmp);
  
  m_SortedID.clear();
  //Put new ID sorted in a vector
  for(int i = 0; i < sizeHit; i++){
    m_SortedID.push_back(pair_SortedID[i].first);
  }

  // A different Cluster number (starting at 1) is assigned to each hit
  static vector<int> ID_ClustHit;
  ID_ClustHit.clear();
  for(int i =0 ; i < sizeHit; i++){
    ID_ClustHit.push_back(i+1);
  }
  
  //Test each Dist(n-n) to find clusters
  //When 2 neutrons are part of a the same cluster, change their ID_ClustHit for the lowest
  //Create a map to hold the Hit_ID for each cluster 
  mapOfClust.clear();
  for(int j = 0; j < sizeHit-1; j++){
    x1 = (*HitX)[m_SortedID[j]], y1 = (*HitY)[m_SortedID[j]], z1 = (*HitZ)[m_SortedID[j]], dx1 = (*HitdX)[m_SortedID[j]], dy1 = (*HitdY)[m_SortedID[j]], dz1 = (*HitdZ)[m_SortedID[j]], t1 = (*HitT)[m_SortedID[j]];   
    dR1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1);  
    for(int jj = j+1; jj < sizeHit ; jj++){
      x2 = (*HitX)[m_SortedID[jj]], y2 = (*HitY)[m_SortedID[jj]], z2 = (*HitZ)[m_SortedID[jj]];   
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
    mapOfClust[ID_ClustHit[i]].push_back(m_SortedID[i]);
  }

  //Put first hit in time of each cluster in a new map mapOfHead and remake 
  //the numbering of clusters (1,2,3...)
  static map< unsigned int, unsigned int > mapOfHead;
  static unsigned int NbrOfClust;
  mapOfHead.clear();
  NbrOfClust = mapOfClust.size();

  static unsigned int count, count2;
  count=1,count2=9;  m_HeadClust.clear();
  Clusters_3D->Clear();
  /* Clusters_3D->SetTitle("Hit Clusters in NEBULA & NeuLAND"); */
  /* Clusters_3D->SetFillColor(29); */
  /* Clusters_3D->SetMarkerSize(10); */
  /* Clusters_3D->SetMarkerStyle(20); */
  /* Clusters_3D->SetLineWidth(2); */
  
  Clusters_3D->SetPoint(1,10000,-2500,0);
  Clusters_3D->SetPoint(2,10000,-2500,2000);
  Clusters_3D->SetPoint(3,10000,2500,0);
  Clusters_3D->SetPoint(4,10000,2500,2000);
  Clusters_3D->SetPoint(5,12500,-2500,0);
  Clusters_3D->SetPoint(6,12500,-2500,2000);
  Clusters_3D->SetPoint(7,12500,2500,0);
  Clusters_3D->SetPoint(8,12500,2500,2000);
  Clusters_3D->SetPointError(1,0,0,0);
  Clusters_3D->SetPointError(2,0,0,0);
  Clusters_3D->SetPointError(3,0,0,0);
  Clusters_3D->SetPointError(4,0,0,0);
  Clusters_3D->SetPointError(5,0,0,0);
  Clusters_3D->SetPointError(6,0,0,0);
  Clusters_3D->SetPointError(7,0,0,0);
  Clusters_3D->SetPointError(8,0,0,0);
  for(auto itr = mapOfClust.begin(); itr != mapOfClust.end(); itr++){
    mapOfHead[count] = itr->second[0];
    m_HeadClust.push_back(itr->second[0]);
    Clusters_3D->SetMarkerColor(count);
    Clusters_3D->SetLineColor(count);
    count++;   
    for(int j =0; j < itr->second.size();j++){
      /* XY_Clusters->Set(count2, (*HitX)[itr->second[j]], (*HitY)[itr->second[j]]); */
      /* XZ_Clusters->Set(count2, (*HitX)[itr->second[j]], (*HitZ)[itr->second[j]]); */
      /* YZ_Clusters->Set(count2, (*HitZ)[itr->second[j]], (*HitY)[itr->second[j]]); */
      /* cout << (*HitX)[itr->second[j]] << " " <<  (*HitdX)[itr->second[j]] << " " << (*HitdY)[itr->second[j]]<< ""<< (*HitdZ)[itr->second[j]]<<  endl; */
      Clusters_3D->SetPoint(count2, (*HitZ)[itr->second[j]], (*HitX)[itr->second[j]], (*HitY)[itr->second[j]]);
      Clusters_3D->SetPointError(count2, (*HitdZ)[itr->second[j]], (*HitdX)[itr->second[j]], (*HitdY)[itr->second[j]]);
      count2++;
    }
  } 
  
  /* Clusters_3D->Draw("err"); */
/*   Clusters_2D.push_back(XY_Clusters); */
/*   Clusters_2D.push_back(XZ_Clusters); */
/*   Clusters_2D.push_back(YZ_Clusters); */

  static double v_n, Dmax, beta1, beta12;
  static vector<double> CrossTalk;
  CrossTalk.clear();
  //Test now the "Friend" clusters
  for(int i = 1; i < NbrOfClust; i++){ 
    x1 = (*HitX)[mapOfHead[i]], y1 = (*HitY)[mapOfHead[i]], z1 = (*HitZ)[mapOfHead[i]], t1 = (*HitT)[mapOfHead[i]]; 
    v_n = sqrt(x1*x1+y1*y1+z1*z1)/t1;
    beta1 = v_n/C_light;
    for(int j = i+1; j <= NbrOfClust; j++){
      x2 = (*HitX)[mapOfHead[j]], y2 = (*HitY)[mapOfHead[j]], z2 = (*HitZ)[mapOfHead[j]], t2 = (*HitT)[mapOfHead[j]]; 
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      Dmax = (t2-t1)*v_n;
      beta12 = Dist/(t2-t1)/C_light;
      if(beta1/beta12 > 1){
      /* if(Dist < Dmax){ */
        CrossTalk.push_back(j);
        /* cout << "HIHIHIHI" << endl; */
      }
    }
  }

  //elimate potential CrossTalk in mapOfHead
  static unsigned int sizeCT;
  sizeCT = CrossTalk.size();
  /* cout << "Event " <<  sizeCT<< endl; */
  for(unsigned int i = 0; i < sizeCT; i++  ){
    mapOfHead.erase(CrossTalk[i]);
    /* cout << "CT " << CrossTalk[i] << endl; */
  }
  
  for (auto itr = mapOfHead.begin(); itr != mapOfHead.end(); itr ++){
    /* cout << "neutrons " << itr->first << endl; */
  }

  //return vector of real neutrons IDs 
  m_Neutrons.clear();
  for(auto itr = mapOfHead.begin(); itr != mapOfHead.end(); itr++){
    m_Neutrons.push_back(itr->second);
  }
  return m_Neutrons;
}

vector<int> CrossTalk::ComputeCrossTalk_HyperSphere(){

  const double C_light = 299.792458;
  static double x1,y1,z1,dx1,dy1,dz1,t1;
  static double x2,y2,z2,dx2,dy2,dz2,t2;
  static double Dist, dR1, dR2, HyperSphere, Beta1;

  //Attribute a new index based on the time arrive from the first to the last hit
  //Using vector of pairs (ID,Time)
  static vector<pair<int, double>> pair_SortedID;
  pair_SortedID.clear();
  for(int i = 0; i < sizeHit; i++){
    pair_SortedID.emplace_back(i, (*HitT)[i]);
  }
  //Sort pair vector (ID,Time) in Time
  sort(pair_SortedID.begin(), pair_SortedID.end(), cmp);
  
  m_SortedID.clear();
  //Put new ID sorted in a vector
  for(int i = 0; i < sizeHit; i++){
    m_SortedID.push_back(pair_SortedID[i].first);
  }

  // A different Cluster number (starting at 1) is assigned to each hit
  static vector<int> ID_ClustHit;
  ID_ClustHit.clear();
  for(int i =0 ; i < sizeHit; i++){
    ID_ClustHit.push_back(i+1);
  }
  
  //Test each Dist(n-n) to find clusters
  //When 2 neutrons are part of a the same cluster, change their ID_ClustHit for the lowest
  //Create a map to hold the Hit_ID for each cluster 
  mapOfClust.clear();
  for(int j = 0; j < sizeHit-1; j++){
    x1 = (*HitX)[m_SortedID[j]], y1 = (*HitY)[m_SortedID[j]], z1 = (*HitZ)[m_SortedID[j]], dx1 = (*HitdX)[m_SortedID[j]], dy1 = (*HitdY)[m_SortedID[j]], dz1 = (*HitdZ)[m_SortedID[j]], t1 = (*HitT)[m_SortedID[j]];   
    dR1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1);  
    for(int jj = j+1; jj < sizeHit ; jj++){
      x2 = (*HitX)[m_SortedID[jj]], y2 = (*HitY)[m_SortedID[jj]], z2 = (*HitZ)[m_SortedID[jj]];   
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      Beta1 = sqrt(x1*x1 + y1*y1 + z1*z1)/t1/C_light;
      HyperSphere = ((dx1/Dist)*(dx1/Dist)+(dy1/Dist)*(dy1/Dist)+(dz1/Dist)*(dz1/Dist)+((t2-t1)/(Dist/(Beta1*C_light)))*((t2-t1)/(Dist/(Beta1*C_light))))/4;
      if(HyperSphere < 2){
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
    mapOfClust[ID_ClustHit[i]].push_back(m_SortedID[i]);
  }

  //Put first hit in time of each cluster in a new map mapOfHead and remake 
  //the numbering of clusters (1,2,3...)
  static map< unsigned int, unsigned int > mapOfHead;
  static unsigned int NbrOfClust;
  mapOfHead.clear();
  NbrOfClust = mapOfClust.size();

  static unsigned int count, count2;
  count=1,count2=9;  m_HeadClust.clear();
  for(auto itr = mapOfClust.begin(); itr != mapOfClust.end(); itr++){
    mapOfHead[count] = itr->second[0];
    m_HeadClust.push_back(itr->second[0]);
    count++;   
  } 
  
  static double v_n, Dmax, beta1, beta12;
  static vector<double> CrossTalk;
  CrossTalk.clear();
  //Test now the "Friend" clusters
  for(int i = 1; i < NbrOfClust; i++){ 
    x1 = (*HitX)[mapOfHead[i]], y1 = (*HitY)[mapOfHead[i]], z1 = (*HitZ)[mapOfHead[i]], t1 = (*HitT)[mapOfHead[i]]; 
    v_n = sqrt(x1*x1+y1*y1+z1*z1)/t1;
    beta1 = v_n/C_light;
    for(int j = i+1; j <= NbrOfClust; j++){
      x2 = (*HitX)[mapOfHead[j]], y2 = (*HitY)[mapOfHead[j]], z2 = (*HitZ)[mapOfHead[j]], t2 = (*HitT)[mapOfHead[j]]; 
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      Dmax = (t2-t1)*v_n;
      beta12 = Dist/(t2-t1)/C_light;
      if(beta1/beta12 > 1){
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
  
  for (auto itr = mapOfHead.begin(); itr != mapOfHead.end(); itr ++){
  }

  //return vector of real neutrons IDs 
  m_Neutrons.clear();
  for(auto itr = mapOfHead.begin(); itr != mapOfHead.end(); itr++){
    m_Neutrons.push_back(itr->second);
  }
  return m_Neutrons;
}

vector<int> CrossTalk::GetSortedHits(){
  return m_SortedID;
}
vector<TGraph> CrossTalk::GetClusters(){
  return Clusters_2D;
}
TGraph2DErrors *CrossTalk::Get3DClusters(){
  return Clusters_3D;
}
vector<int> CrossTalk::GetHeadClust(){
  return m_HeadClust;
}


