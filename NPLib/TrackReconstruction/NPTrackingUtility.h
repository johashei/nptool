#ifndef NPUTILITY_H
#define NPUTILITY_H
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta  contact address: matta@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : July 2020                                               *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with finding all the track event by event                *
 *****************************************************************************/

namespace NPL{
  // return the minimum distance between v and w defined respectively by points 
  // v1,v2 and w1 w1
  // Also compute the best crossing position BestPosition, i.e. average position
  // at the minimum distance.
  double MinimumDistance(const TVector3& v1,const TVector3& v2, const TVector3& w1, const TVector3& w2, TVector3& BestPosition, TVector3& delta){
  TVector3 v = v2-v1;
  TVector3 w = w2-w1;
  // Finding best position
  TVector3 e = v1-w1;
  double A = -(v.Mag2()*w.Mag2()-(v.Dot(w)*v.Dot(w)));
  double s = (-v.Mag2()*(w.Dot(e))+(v.Dot(e))*(w.Dot(v)))/A;
  double t = (w.Mag2()*(v.Dot(e))-(w.Dot(e)*w.Dot(v)))/A;
  double d = sqrt((e+v*t-w*s).Mag2());
 
  BestPosition = 0.5*(v1+t*v+w1+s*w);
  delta = (v1+t*v-w1-s*w);
  return d;
  }
}


#endif
