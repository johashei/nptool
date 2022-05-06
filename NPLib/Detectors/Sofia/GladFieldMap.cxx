/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@ceafr *
 * Creation Date  : May 2022                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold basic field map and brho reconstruction                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "GladFieldMap.h"
#include "NPPhysicalConstants.h"
#include "Math/Factory.h"
using namespace NPUNITS;

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(GladFieldMap)


//////////////////////////////////////////////////////////////////////
GladFieldMap::GladFieldMap() {
  m_BrhoScan=NULL;
  m_Bfield = TVector3(0,1.5/1000,0);
  m_Z_Glad = 4.434*m;
  m_Leff = 2.*m;
  m_Tilt = 14.*deg;

  m_Zmax = 9.*m;
  m_Limit = 1000;
  m_dt = 0.1*nanosecond;

  m_CentralTheta = 20.*deg;
  m_X_MWPC3 = -1436.;
  m_Z_MWPC3 = 8380.;

  m_InitPos = TVector3(0,0,0);
  m_InitDir = TVector3(0,0,1);
}



//////////////////////////////////////////////////////////////////////
GladFieldMap::~GladFieldMap() {
}

//////////////////////////////////////////////////////////////////////
double GladFieldMap::FindBrho(TVector3 Pos_init, TVector3 Dir_init, TVector3 Pos_final){
  double Brho;
  m_InitPos  = Pos_init;
  m_InitDir  = Dir_init;
  m_FinalPos = Pos_final;

  m_InitDir = m_InitDir.Unit();

  if(!m_BrhoScan)
    BrhoScan(1,11,0.2);

  Brho = m_BrhoScan->Eval(m_FinalPos.X());

  return Brho;
}
//////////////////////////////////////////////////////////////////////
TGraph* GladFieldMap::BrhoScan(double Brho_min, double Brho_max, double Brho_step){
  if(m_BrhoScan)
    delete m_BrhoScan;

  m_BrhoScan = new TGraph();

  int i=0;
  for(double Brho=Brho_min; Brho<Brho_max; Brho+=Brho_step){
    vector<TVector3> vPos = Propagate(Brho, m_InitPos, m_InitDir);

    TVector3 M_inter = CalculateIntersectionPoint(vPos);

    m_BrhoScan->SetPoint(i,M_inter.X(),Brho);
    i++;
  }

  return m_BrhoScan;
}

//////////////////////////////////////////////////////////////////////
vector<TVector3> GladFieldMap::Propagate(double Brho, TVector3 Pos, TVector3 Dir){
  vector<TVector3> vPos;

  TVector3 xk1, xk2, xk3, xk4;
  TVector3 pk1, pk2, pk3, pk4;

  static NPL::Particle N("90Zr");
  N.SetBrho(Brho);
  double T = N.GetEnergy();
  double M = N.Mass();
  double E = T + M;
  double P = sqrt(T*T + 2*M*T)/c_light;

  Dir = Dir.Unit();

  double px = P*Dir.X();
  double py = P*Dir.Y();
  double pz = P*Dir.Z();
  TVector3 Imp = TVector3(px,py,pz);

  int count=0;
  while(Pos.Z()<m_Zmax && count<m_Limit){
    func(N, Pos, Imp, xk1, pk1);
    func(N, Pos + (m_dt/2)*xk1, Imp + (m_dt/2)*pk1, xk2, pk2);
    func(N, Pos + (m_dt/2)*xk2, Imp + (m_dt/2)*pk2, xk3, pk3);
    func(N, Pos + m_dt*xk3, Imp + m_dt*pk3, xk4, pk4);

    Pos += (m_dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4);
    Imp += (m_dt/6)*(pk1 + 2*pk2 + 2*pk3 + pk4);

    vPos.push_back(Pos);
    count++;
  }

  return vPos;

}

//////////////////////////////////////////////////////////////////////
void GladFieldMap::func(NPL::Particle& N, TVector3 Pos, TVector3 Imp, TVector3& xk, TVector3& pk){
  double px, py, pz;
  px = Imp.X();
  py = Imp.Y();
  pz = Imp.Z();

  double M = N.Mass();
  double T = N.GetEnergy();
  double E = T + M;

  double vx, vy, vz;
  vx = px*pow(c_light,2)/E;
  vy = py*pow(c_light,2)/E;
  vz = pz*pow(c_light,2)/E;
  xk.SetX(vx);
  xk.SetY(vy);
  xk.SetZ(vz);

  TVector3 B = LoadMap(Pos);
  double Bx, By, Bz;
  Bx = B.X();
  By = B.Y();
  Bz = B.Z();
  double q = N.GetZ()*eplus;
  pk.SetX(q*(vy*Bz - vz*By));
  pk.SetY(q*(vz*Bx - vx*Bz));
  pk.SetZ(q*(vx*By - vy*Bx));

  return;
}

//////////////////////////////////////////////////////////////////////
TVector3 GladFieldMap::CalculateIntersectionPoint(vector<TVector3> vPos){
 
  unsigned int size = vPos.size();

  // Track equation Z = a0*X + b0
  double a0, b0;
  a0 = (vPos[size-1].Z() - vPos[size-100].Z()) / (vPos[size-1].X() - vPos[size-100].X());
  b0 = vPos[size-1].Z() - a0*vPos[size-1].X();

  // MWPC3 equation Z_MWPC = a1*X_MWPC + b1
  double a1, b1;
  a1 = tan(m_CentralTheta);
  b1 = m_Z_MWPC3 - m_X_MWPC3*tan(m_CentralTheta);

  double Mx, My, Mz;
  Mx = (b0 - b1) / (a1 -a0);
  My = 0;
  Mz = a0*Mx + b0;

  TVector3 M = TVector3(Mx,My,Mz);
  return M;
}

//////////////////////////////////////////////////////////////////////
TVector3 GladFieldMap::LoadMap(TVector3 Pos) {

  double Bx = 0;
  double By = 0;
  double Bz = 0;

  double Zinf = m_Z_Glad + Pos.X()*tan(m_Tilt) - m_Leff/2;
  double Zsup = m_Z_Glad + Pos.X()*tan(m_Tilt) + m_Leff/2;
  
  if(Pos.Z()>Zinf && Pos.Z()<Zsup){
    Bx = m_Bfield.X();
    By = m_Bfield.Y();
    Bz = m_Bfield.Z();
  }
  
  TVector3 B = TVector3(Bx,By,Bz);
  return B;
}

