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
  // take back the first hit, automaticaly identified as the first Neutron
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
  static double Dist, dR1, dR2;

  // Identify Hits around the 1st neutron, 1st cluster
  x1 = (*HitX)[FirstHit], y1 = (*HitY)[FirstHit], z1 = (*HitZ)[FirstHit], dx1 = (*HitdX)[FirstHit], dy1 = (*HitdY)[FirstHit], dz1 = (*HitdZ)[FirstHit], t1 = (*HitT)[FirstHit];   
  dR1 = sqrt(x1*dx1 + dy1*dy1 + dz1*dz1);
  for(int i = 0; i < sizeHit; i++){
    if(i == FirstHit){
      ClustHit.push_back(1);
      continue;
    }
    x2 = (*HitX)[i], y2 = (*HitY)[i], z2 = (*HitZ)[i], dx2 = (*HitdX)[i], dy2 = (*HitdX)[i], dz2 = (*HitdZ)[i], t2 = (*HitT)[i];   
    Dist = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))-dR1-dR2;
    if(Dist < coef*dR1){ 
      ClustHit.push_back(1); // this hit is part of the cluster 1
    }
  }

  static vector<double> m_Neutrons;
  m_Neutrons.clear();
  return m_Neutrons;
}

int CrossTalk::GetFirstN(){
  return FirstHit;
}

