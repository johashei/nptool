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

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(SamuraiFieldMap);

////////////////////////////////////////////////////////////////////////////////
SamuraiFieldMap::SamuraiFieldMap(std::string file){
  LoadMap(file);
}
////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadMap(std::string file){
  if(file.find(".bin")!=std::string::npos)
    LoadBinary(file);
  else
    LoadAscii(file);
}
////////////////////////////////////////////////////////////////////////////////
std::vector<float> SamuraiFieldMap::GetB(std::vector<float>& pos){
  static vector<float> nullv ={0,0,0};
  auto it=m_field.find(pos);
  if(it!=m_field.end()){
    return it->second;
  }
  else 
    return nullv;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<float> SamuraiFieldMap::InterpolateB(std::vector<float>& pos){
  static vector<float> nullv ={0,0,0};
  double x = pos[0];
  double y = pos[1];
  double z = pos[2];

  // out of bound 
  if(x<m_x_min || x>m_x_max)
    return nullv;
  if(y<m_y_min || y>m_y_max)
    return nullv;
  if(z<m_z_min || z>m_z_max)
    return nullv;

  m_bin = 10;
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
      double d = 1e-6+ sqrt( (pos[0]-it[i]->first[0])*(pos[0]-it[i]->first[0])+
          (pos[1]-it[i]->first[1])*(pos[1]-it[i]->first[1])+
          (pos[2]-it[i]->first[2])*(pos[2]-it[i]->first[2]));

      Bx+=it[i]->second[0]/d;
      By+=it[i]->second[1]/d;
      Bz+=it[i]->second[2]/d;
      totalW+=1./d;
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
