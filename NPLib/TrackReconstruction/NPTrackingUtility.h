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
  double MinimunDistance(TVector3 v1,TVector3 v2, TVector3 w1, TVector3 w2){
  TVector3 v=v1-v2;
  TVector3 w = w2-w1;
  // Minimum distance
  // let be n perpendicular to both line
  TVector3 n = v.Cross(w);

  double d = n.Dot(v1-w1)/n.Mag();
  return d;
  }
}


#endif
