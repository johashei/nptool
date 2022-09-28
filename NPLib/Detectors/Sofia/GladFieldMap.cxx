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
    m_min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    m_func = ROOT::Math::Functor(this,&GladFieldMap::Delta,1);
    m_min->SetFunction(m_func);
    m_min->SetPrintLevel(-1);
    m_bin = 50;
    m_Scale = -2135./3583.81;
    m_Z_Glad = 2724.;
    m_Leff = 2.*m;
    m_Tilt = 14.*deg;

    m_Zmax = 8.5*m;
    m_Limit = 1000;
    m_dt = 0.1*nanosecond;
    m_x_min=1e6;
    m_x_max=-1e6;
    m_y_min=1e6;
    m_y_max=-1e6;
    m_z_min=1e6;
    m_z_max=-1e6;

    m_CentralTheta = 20.*deg;
    m_X_MWPC3 = -1436.;
    m_Z_MWPC3 = 8380.;
    m_Angle_MWPC3 = 20.*deg;
    m_R_MWPC3 = 4199.; 

    m_InitPos = TVector3(0,0,0);
    m_InitDir = TVector3(0,0,1);
  }



//////////////////////////////////////////////////////////////////////
GladFieldMap::~GladFieldMap() {
}

//////////////////////////////////////////////////////////////////////
double GladFieldMap::Delta(const double* parameter){
  static TVector3 diff;
  static vector<TVector3> pos;

  pos = Propagate(parameter[0], m_InitPos, m_InitDir, false);
  //pos.back().RotateY(-m_Angle_MWPC3+m_Tilt);
  
  diff = pos.back() - m_FinalPos;
  //cout << pos.back().X() << " " << pos.back().Z() << endl;
  return diff.Mag2();
}

//////////////////////////////////////////////////////////////////////
double GladFieldMap::FindBrho(TVector3 Pos_init, TVector3 Dir_init, TVector3 Pos_final){
    
  if(!m_BrhoScan)
    BrhoScan(6.,11,1);

  m_InitPos  = Pos_init;
  m_InitDir  = Dir_init;
  m_FinalPos = Pos_final;

  if(m_FinalPos.X()>-2000){
    m_InitDir = m_InitDir.Unit();

    double param[1];
    param[0] = m_BrhoScan->Eval(m_FinalPos.X());
    
    return param[0];
    
    m_min->Clear();
    m_min->SetPrecision(1e-6);
    m_min->SetMaxFunctionCalls(1000);
    m_min->SetLimitedVariable(0,"B",param[0],0.1,6,11);
    m_min->Minimize();
    
    return m_min->X()[0];
  }

  else
    return -1;
}

//////////////////////////////////////////////////////////////////////
TGraph* GladFieldMap::BrhoScan(double Brho_min, double Brho_max, double Brho_step){
  if(m_BrhoScan)
    delete m_BrhoScan;

  m_BrhoScan = new TGraph();
  unsigned int size = (Brho_max-Brho_min)/Brho_step;
  m_BrhoScan->Set(size);

  unsigned int i=0;
  TVector3 pos = TVector3(0,0,0);
  TVector3 dir = TVector3(0,0,1);
  //pos.RotateY(m_Tilt);
  //dir.RotateY(m_Tilt);
  for(double Brho=Brho_min; Brho<Brho_max; Brho+=Brho_step){
    vector<TVector3> vPos = Propagate(Brho, pos, dir, false);
    //vPos.back().RotateY(-m_Angle_MWPC3);

    m_BrhoScan->SetPoint(i++,vPos.back().X(),Brho);
    //cout << vPos.back().X() << " " << Brho << endl;
  }

  m_BrhoScan->Sort();

  return m_BrhoScan;
}

//////////////////////////////////////////////////////////////////////
TVector3 GladFieldMap::PropagateToMWPC(TVector3 pos, TVector3 dir){
  // go to MWPC3 frame reference
  //pos.RotateY(-m_Angle_MWPC3);
  //dir.RotateY(-m_Angle_MWPC3);

  double deltaZ = m_Z_MWPC3 - pos.Z();
  dir*=deltaZ/dir.Z();
  pos+=dir;
  pos.SetX(pos.X());
  //pos.RotateY(m_Angle_MWPC3);

  return pos;

}
//////////////////////////////////////////////////////////////////////
vector<TVector3> GladFieldMap::Propagate(double Brho, TVector3 Pos, TVector3 Dir, bool store){
  //Pos.RotateY(m_Tilt);
  //Dir.RotateY(m_Tilt);
  static NPL::Particle N("90Zr");
  N.SetBrho(Brho);
  
  // track result
  static std::vector<TVector3> track;
  track.clear();

  // starting point of the track
  if(store){
    //Pos.RotateY(-m_Tilt);
    track.push_back(Pos);
    //Pos.RotateY(-m_Tilt);
  }

  static double r;
  r = sqrt(Pos.X()*Pos.X() + Pos.Z()*Pos.Z());

  // Propagate to m_R_max with one line
  /*while(r>m_R_max && count<m_Limit){
    Pos+=(r-m_R_max)/cos(Dir.Theta())*Dir.Unit();
    r = 1.01*sqrt(Pos.X()*Pos.X() + Pos.Z()*Pos.Z());
  }

  if(r<=m_R_max){ // success
    if(store){
      //Pos.RotateY(-m_Tilt);
      track.push_back(Pos);
      //Pos.RotateY(m_Tilt);
    }
  }
  else{
    cout << "Fail" << endl;
    return track;
  }*/

  static TVector3 xk1, xk2, xk3, xk4;
  static TVector3 pk1, pk2, pk3, pk4;
  static double T, M, E, P;
  static double px, py, pz;
  static TVector3 Imp;
  
  T = N.GetEnergy();
  M = N.Mass();
  E = T + M;
  P = sqrt(T*T + 2*M*T)/c_light;
 
  Dir = Dir.Unit();

  px = P*Dir.X();
  py = P*Dir.Y();
  pz = P*Dir.Z();
  Imp = TVector3(px,py,pz);

  unsigned int count=0;
  while(Pos.Z()<=m_Zmax && count<m_Limit){
    func(N, Pos, Imp, xk1, pk1);
    func(N, Pos + (m_dt/2)*xk1, Imp + (m_dt/2)*pk1, xk2, pk2);
    func(N, Pos + (m_dt/2)*xk2, Imp + (m_dt/2)*pk2, xk3, pk3);
    func(N, Pos + m_dt*xk3, Imp + m_dt*pk3, xk4, pk4);

    Pos += (m_dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4);
    Imp += (m_dt/6)*(pk1 + 2*pk2 + 2*pk3 + pk4);

    if(store){
      //Pos.RotateY(-m_Tilt);
      track.push_back(Pos);
      //Pos.RotateY(m_Tilt);
    }
    r = sqrt(Pos.X()*Pos.X() + Pos.Z()*Pos.Z());
    count++;
  }

  Imp=Imp.Unit();
  Pos = PropagateToMWPC(Pos,Imp);
  //Pos.RotateY(-m_Tilt);
  track.push_back(Pos);

  return track;

}

//////////////////////////////////////////////////////////////////////
void GladFieldMap::func(NPL::Particle& N, TVector3 Pos, TVector3 Imp, TVector3& xk, TVector3& pk){
  static vector<double> B;
  static double px, py, pz;
  static double Bx, By, Bz;
  static double vx, vy, vz;
  static double M, T, E;
  static double q;

  px = Imp.X();
  py = Imp.Y();
  pz = Imp.Z();

  M = N.Mass();
  T = N.GetEnergy();
  E = T + M;

  vx = px*pow(c_light,2)/E;
  vy = py*pow(c_light,2)/E;
  vz = pz*pow(c_light,2)/E;
  xk.SetX(vx);
  xk.SetY(vy);
  xk.SetZ(vz);

  //vector<double> position = {Pos.X(),Pos.Y(),Pos.Z()};
  B = InterpolateB(Pos);
  Bx = B[0];
  By = B[1];
  Bz = B[2];

  q = N.GetZ()*eplus;
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
void GladFieldMap::LoadMap(string filename) {
  ifstream ifile(filename.c_str());
  if(!ifile.is_open()){
    cout << "Error: failed to load Glad field map: " << filename << endl;
    exit(1);
  }

  cout << "///////// Loading ASCII Glad field map: " << filename << endl;
  double x,y,z,Bx,By,Bz;

  unsigned int count=0;
  //while(!ifile.eof()){
  for(unsigned int i=0; i<401841; i++){
    ifile >> x >> y >> z >> Bx >> By >> Bz;
    if(++count%10000==0)
      cout << "\r - Loading " << count << " values " << flush;

    x = x*10;
    y = y*10;
    z = z*10;

    z = z + x*sin(m_Tilt);
    z += m_Z_Glad;
    
    Bx *= m_Scale;
    By *= m_Scale;
    Bz *= m_Scale;

    vector<double> p = {x,y,z};
    Bx*=tesla;
    By*=tesla;
    Bz*=tesla;
    vector<double> B = {Bx,By,Bz};
    m_field[p] = B;

    if(x<m_x_min)
      m_x_min=x;
    if(y<m_y_min)
      m_y_min=y;
    if(z<m_z_min)
      m_z_min=z;
    if(x>m_x_max)
      m_x_max=x;
    if(y>m_y_max)
      m_y_max=y;
    if(z>m_z_max)
      m_z_max=z;

  }

  m_R_max = m_z_max;
  cout << endl;
  cout << "///////// ASCII file loaded"<< endl;
  cout << "m_field size= " << m_field.size() << endl;
  cout << "m_x_min= " << m_x_min << endl;
  cout << "m_x_max= " << m_x_max << endl;
  cout << "m_y_min= " << m_y_min << endl;
  cout << "m_y_max= " << m_y_max << endl;
  cout << "m_z_min= " << m_z_min << endl;
  cout << "m_z_max= " << m_z_max << endl;
  cout << "/////////"<< endl;
  
  ifile.close();
}

//////////////////////////////////
vector<double> GladFieldMap::InterpolateB(const vector<double>& pos)
{
  static vector<double> nullv = {0,0,0};

  double x,y,z;
  x=pos[0];
  y=pos[1];
  z=pos[2];

  // out of bound
  if(x<m_x_min || x>m_x_max)
    return nullv;
  if(y<m_y_min || y>m_y_max)
    return nullv;
  if(z<m_z_min || z>m_z_max)
    return nullv;

  double x0 = (double)((int)(x)/m_bin)*m_bin;
  if(x<=x0)
    x0=(double)((int)(x-m_bin)/m_bin)*m_bin; 

  double x1 = (double)((int)(x)/m_bin)*m_bin;
  if(x>=x1)
    x1=(double)((int)(x+m_bin)/m_bin)*m_bin; 

  double y0 = (double)((int)(y)/m_bin)*m_bin;
  if(y<=y0)
    y0=(double)((int)(y-m_bin)/m_bin)*m_bin; 

  double y1 = (double)((int)(y)/m_bin)*m_bin;
  if(y>=y1)
    y1=(double)((int)(y+m_bin)/m_bin)*m_bin; 

  double z0 = (double)((int)(z)/m_bin)*m_bin;
  if(z<=z0)
    z0=(double)((int)(z-m_bin)/m_bin)*m_bin; 

  double z1 = (double)((int)(z)/m_bin)*m_bin;
  if(z>=z1)
    z1=(double)((int)(z+m_bin)/m_bin)*m_bin; 

  //vector<double> X={xm,ym,zm};
  vector<double> X000={x0,y0,z0};
  vector<double> X111={x1,y1,z1};
  vector<double> X100={x1,y0,z0};
  vector<double> X010={x0,y1,z0};
  vector<double> X001={x0,y0,z1};
  vector<double> X101={x1,y0,z1};
  vector<double> X011={x0,y1,z1};
  vector<double> X110={x1,y1,z0};

  vector<double> C000;
  vector<double> C111;
  vector<double> C100;
  vector<double> C010;
  vector<double> C001;
  vector<double> C101;
  vector<double> C011;
  vector<double> C110;

  if(m_field.lower_bound(X000)!=m_field.end())
    C000=m_field.lower_bound(X000)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X111)!=m_field.end())
    C111=m_field.lower_bound(X111)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X100)!=m_field.end())
    C100=m_field.lower_bound(X100)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X010)!=m_field.end())
    C010=m_field.lower_bound(X010)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X001)!=m_field.end())
    C001=m_field.lower_bound(X001)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X101)!=m_field.end())
    C101=m_field.lower_bound(X101)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X011)!=m_field.end())
    C011=m_field.lower_bound(X011)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X110)!=m_field.end())
    C110=m_field.lower_bound(X110)->second;
  else 
    return nullv;

  double xd = (x-x0)/(x1-x0);    
  double yd = (y-y0)/(y1-y0);    
  double zd = (z-z0)/(z1-z0);    
  double alphaX = 1-xd;
  double alphaY = 1-yd;
  double alphaZ = 1-zd;
  // X 
  vector<double> C00 = {C000[0]*alphaX+C100[0]*xd,C000[1]*alphaX+C100[1]*xd,C000[2]*alphaX+C100[2]*xd};
  vector<double> C01 = {C001[0]*alphaX+C101[0]*xd,C001[1]*alphaX+C101[1]*xd,C001[2]*alphaX+C101[2]*xd};
  vector<double> C10 = {C010[0]*alphaX+C110[0]*xd,C010[1]*alphaX+C110[1]*xd,C010[2]*alphaX+C110[2]*xd};
  vector<double> C11 = {C011[0]*alphaX+C111[0]*xd,C011[1]*alphaX+C111[1]*xd,C011[2]*alphaX+C111[2]*xd};
  // Y
  vector<double> C0  = {C00[0] *alphaY+C10[0] *yd,C00[1] *alphaY+C10[1] *yd,C00[2] *alphaY+C10[2] *yd};
  vector<double> C1  = {C01[0] *alphaY+C11[0] *yd,C01[1] *alphaY+C11[1] *yd,C01[2] *alphaY+C11[2] *yd};
  // Z
  vector<double> res = {C0[0]  *alphaZ+C1[0]  *zd,C0[1]  *alphaZ+C1[1]  *zd,C0[2]  *alphaZ+C1[2]  *zd};

  return res;
}

//////////////////////////////////
vector<double> GladFieldMap::GetB(const vector<double>& pos)
{
  static vector<double> nullv = {0,0,0};

  auto it = m_field.find(pos);
  if(it!=m_field.end()){
    return it->second;
  }
  else
    return nullv;

}

