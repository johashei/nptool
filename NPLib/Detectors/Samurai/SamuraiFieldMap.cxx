/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : May  2021                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Samurai field map data                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "SamuraiFieldMap.h"
#include "NPPhysicalConstants.h"
using namespace NPUNITS;

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(SamuraiFieldMap);

////////////////////////////////////////////////////////////////////////////////
double SamuraiFieldMap::FindBrho(TVector3& p_fdc0,TVector3& d_fdc0,TVector3& p_fdc2,TVector3& d_fdc2){
  if(!m_BrhoScan)
    BrhoScan(2.5,10,0.1);

  // do a first guess based on fdc2 pos
  double b = m_BrhoScan->Eval(p_fdc2.X()); 
  double b0 = b;
  vector<TVector3> pos=Propagate(3000,b,p_fdc0,d_fdc0);
  double step = 1;
  double d = (pos.back()-p_fdc2).Mag();
  double dd=d;
  short sign =1;
  unsigned int limit =1000;
  unsigned count=0;
  while(step>1e-6 && count<limit){
   dd=d;
   b+=sign*step;
   pos=Propagate(3000,b,p_fdc0,d_fdc0); 
   d = (pos.back()-p_fdc2).Mag();
   if(d>=dd){
    step/=10;
    sign=-sign;
   }
   count++;
  }
  return b-sign*0.5*step;
}

////////////////////////////////////////////////////////////////////////////////
TGraph* SamuraiFieldMap::BrhoScan(double min, double max,double step){
  if(m_BrhoScan)
    delete m_BrhoScan;
  m_BrhoScan=new TGraph;
  unsigned int size = (max-min)/step;
  m_BrhoScan->Set(size);
  unsigned int i=0;
  for(double b = min ; b < max ; b+=step){
   vector<TVector3> pos= Propagate(3000,b,TVector3(0,0,0),TVector3(0,0,1));
   pos.back().RotateY(-m_fdc2angle);
   m_BrhoScan->SetPoint(i++,pos.back()[0],b); 
  }
  return m_BrhoScan;
}

////////////////////////////////////////////////////////////////////////////////
TVector3 SamuraiFieldMap::PropagateToFDC2(TVector3 pos, TVector3 dir){
  // go to FDC2 frame reference
  pos.RotateY(-m_fdc2angle);
  dir.RotateY(-m_fdc2angle);

  double deltaZ=m_fdc2R-pos.Z();
  dir*=deltaZ/dir.Z();
  pos+=dir;
  pos.SetX(pos.X());
  pos.RotateY(m_fdc2angle);
  return pos;
}

////////////////////////////////////////////////////////////////////////////////
std::vector< TVector3 > SamuraiFieldMap::Propagate(double rmax, double Brho, TVector3 pos, TVector3 dir){
  pos.RotateY(m_angle);
  dir.RotateY(m_angle);
  dir=dir.Unit();
  // Property of a particle with the correct Brho:
  // We assume a 4He to compute v
  // The choice of the particle is of no importance
  static NPL::Particle N("4He");
  N.SetBrho(Brho);

  // track result
  std::vector< TVector3 > track;
  
  // starting point of the track
  pos.RotateY(-m_angle);
  track.push_back(pos);
  pos.RotateY(m_angle);

  dir=dir.Unit();
  double r = sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
  // number of step taken
  unsigned int count = 0;
  // maximum number of state before giving up
  unsigned int limit = 1000;

  // First propagate to r_max with one line
  while(r>rmax && count<limit){
    pos+=(r-rmax)/cos(dir.Theta())*dir.Unit();
    r= 1.01*sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
  }

  if(r<=rmax){ // success
    pos.RotateY(-m_angle);
    track.push_back(pos);
    pos.RotateY(m_angle);
  }
  else {// failure
    //cout << "Fail" << endl;
    return track;
  }

  TVector3 xk1,xk2,xk3,xk4; // position
  TVector3 pk1,pk2,pk3,pk4; // impulsion

  double K = N.GetEnergy(); // kinetic energy
  double m = N.Mass(); // mc2
  double P = sqrt(K*K+2*K*m)/c_light; // P
  double px = P*dir.X();//px
  double py = P*dir.Y();//py
  double pz = P*dir.Z();//pz
  TVector3 imp = P*dir;
  double h = 0.1*nanosecond; // typical step length  ~1mm at beta 0.6
  while(r<=rmax && count < limit){
    func(N, pos           , imp            , xk1, pk1);
    func(N, pos+xk1*(h/2.), imp+pk1*(h/2.) , xk2, pk2);
    func(N, pos+xk2*(h/2.), imp+pk2*(h/2.) , xk3, pk3);
    func(N, pos+xk3*h     , imp+pk3*h      , xk4, pk4);
    pos +=(xk1+2*xk2+2*xk3+xk4)*(h/6.); 
    imp +=(pk1+2*pk2+2*pk3+pk4)*(h/6.); 
    pos.RotateY(-m_angle);
    track.push_back(pos);
    pos.RotateY(m_angle);
    r = sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
    count++;
  }
  imp=imp.Unit();
  pos = PropagateToFDC2(pos, imp);
  pos.RotateY(-m_angle);
  track.push_back(pos);
  
  return track;

}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::func(NPL::Particle& N, TVector3 pos, TVector3 imp, TVector3& new_pos, TVector3& new_imp){
  double px,py,pz;
  px=imp.X(); 
  py=imp.Y();
  pz=imp.Z();

  double P2,D,m2c4;
  P2=imp.Mag2(); // P2
  m2c4 = N.Mass()*N.Mass();
  D=sqrt(m2c4+P2*c_squared); // sqrt(m2c4+P2c2)
  double vx=px*c_squared/D;// pxc * c / D = pxc2/D
  double vy=py*c_squared/D;
  double vz=pz*c_squared/D;
  new_pos.SetX(vx);
  new_pos.SetY(vy);
  new_pos.SetZ(vz);
  vector<float> B = InterpolateB(pos);
  double Bx= B[0]; 
  double By= B[1];
  double Bz= B[2];
  double q = N.GetZ()*eplus; // issue with the tesla/coulomb definition
  new_imp.SetX(q*(vy*Bz-vz*By));// q*pyc2*Bz/D -q*pzc2*By/D
  new_imp.SetY(q*(vz*Bx-vx*Bz));
  new_imp.SetZ(q*(vx*By-vy*Bx));
}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadMap(double angle,std::string file,unsigned int bin){
  m_bin=bin;
  m_angle=angle;
  if(file.find(".bin")!=std::string::npos)
    LoadBinary(file);
  else
    LoadAscii(file);
}
////////////////////////////////////////////////////////////////////////////////
std::vector<float> SamuraiFieldMap::GetB(std::vector<float>& pos){
  static vector<float> nullv ={0,0,0};
  // the map is only 1/4 of the detector so we apply symetrie:
  double x,y,z ;
  
  if(pos[0]<0)
    pos[0] = -pos[0];

  if(pos[2]<0)
    pos[2] = -pos[2];

  auto it=m_field.find(pos);
  if(it!=m_field.end()){
    return it->second;
  }
  else 
    return nullv;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<float> SamuraiFieldMap::InterpolateB(const std::vector<float>& pos){
  static vector<float> nullv ={0,0,0};
  // the map is only 1/4 of the detector so we apply symetrie:
  double x,y,z ;
  
  if(pos[0]>0)
    x = pos[0];
  else
    x = -pos[0];

  y = pos[1];

  if(pos[2]>0)
    z = pos[2];
  else
    z = -pos[2];

 // out of bound 
  if(x<m_x_min || x>m_x_max)
    return nullv;
  if(y<m_y_min || y>m_y_max)
    return nullv;
  if(z<m_z_min || z>m_z_max)
    return nullv;
 
   

  float xm = (float)((int)x/m_bin*m_bin);
  float ym = (float)((int)y/m_bin*m_bin);
  float zm = (float)((int)z/m_bin*m_bin);

  vector<float> p0={xm,ym,zm};
  vector<float> p1={xm+m_bin,ym,zm};
  vector<float> p2={xm,ym+m_bin,zm};
  vector<float> p3={xm,ym,zm+m_bin};
  vector<float> p4={xm-m_bin,ym,zm};
  vector<float> p5={xm,ym-m_bin,zm};
  vector<float> p6={xm,ym,zm-m_bin};

  vector<map<vector<float>,vector<float>>::iterator> it=
  { m_field.lower_bound(p0),
    m_field.lower_bound(p1),m_field.lower_bound(p2),m_field.lower_bound(p3),
    m_field.lower_bound(p4),m_field.lower_bound(p5),m_field.lower_bound(p6)};

  float Bx=0;
  float By=0;
  float Bz=0;
  float totalW=0;
  auto end=m_field.end();
  unsigned int size = it.size();
  for(unsigned int i = 0 ; i < size; i++){
    if(it[i]!=end){
      double d = 1e-6+sqrt( (x-it[i]->first[0])*(x-it[i]->first[0])+
          (y-it[i]->first[1])*(y-it[i]->first[1])+
          (z-it[i]->first[2])*(z-it[i]->first[2]));

      Bx+=it[i]->second[0]/(d*d);
      By+=it[i]->second[1]/(d*d);
      Bz+=it[i]->second[2]/(d*d);
      totalW+=1./(d*d);
    }
  }
  vector<float> res = {Bx/totalW,By/totalW,Bz/totalW};
  return res;
}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadAscii(std::string file){
  ifstream in(file.c_str());
  if(!in.is_open()){
    cout << "Error: failed to load samurai field map " << file << endl;
    exit(1);
  }

  cout << "//////// Loading Ascii Samurai field map " << file << endl; 
  float x,y,z,Bx,By,Bz;

  m_x_max=m_y_max=m_z_max=-1e32;
  m_x_min=m_y_min=m_z_min=1e32;
  unsigned int  count =0 ;

  // ignore 8 first line 
  string buffer;
  for(unsigned int i = 0 ; i < 8 ; i++){
    getline(in,buffer);
  }

  while(in >> x >> y >> z >> Bx >> By >> Bz){
    if(++count%50000==0)
      cout << "\r  - Loading " << count << " values " << flush; 
    vector<float> p = {x,y,z};
    Bx*=tesla;
    By*=tesla;
    Bz*=tesla;
    vector<float> B = {Bx,By,Bz};
    m_field[p]=B;
    if(x<m_x_min)
      m_x_min=x;
    if(x>m_x_max)
      m_x_max=x;  
    if(y<m_y_min)
      m_y_min=y;
    if(y>m_y_max)
      m_y_max=y;  
    if(z<m_z_min)
      m_z_min=z;
    if(z>m_z_max)
      m_z_max=z;  
  }

  cout << "\r  - " << count << " values loaded" << endl; 
  cout << "  - min(" << m_x_min <<";"<< m_y_min <<";" << m_z_min<< ") max(" << m_x_max <<";"<< m_y_max <<";" << m_z_max<< ")" << endl; 
  in.close();
}
////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadBinary(std::string file){
  ifstream in(file.c_str(),std::ifstream::binary);
  if(!in.is_open()){
    cout << "Error: failed to load samurai field map " << file << endl;
    exit(1);
  }

  cout << "//////// Loading Binary Samurai field map " << file << endl; 
  float x,y,z,Bx,By,Bz;

  m_x_max=m_y_max=m_z_max=-1e32;
  m_x_min=m_y_min=m_z_min=1e32;
  unsigned int  count =0 ;
  while(!in.eof()){

    if(++count%50000==0)
      cout << "\r  - Loading " << count << " values " << flush; 

    in.read((char*)&x,sizeof(x));
    in.read((char*)&y,sizeof(y));
    in.read((char*)&z,sizeof(z));
    in.read((char*)&Bx,sizeof(Bx));
    in.read((char*)&By,sizeof(By));
    in.read((char*)&Bz,sizeof(Bz));

    vector<float> p = {x,y,z};
    Bx*=tesla;
    By*=tesla;
    Bz*=tesla;
    vector<float> B = {Bx,By,Bz};
    m_field[p]=B;
    if(x<m_x_min)
      m_x_min=x;
    if(x>m_x_max)
      m_x_max=x;  
    if(y<m_y_min)
      m_y_min=y;
    if(y>m_y_max)
      m_y_max=y;  
    if(z<m_z_min)
      m_z_min=z;
    if(z>m_z_max)
      m_z_max=z;  
  }
  cout << "\r  - " << count << " values loaded" << endl; 
  cout << "  - min(" << m_x_min <<";"<< m_y_min <<";" << m_z_min<< ") max(" << m_x_max <<";"<< m_y_max <<";" << m_z_max<< ")" << endl; 
  in.close();
}
