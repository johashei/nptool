#ifndef SamuraiFieldMap_h
#define SamuraiFieldMap_h
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

#include<string>
#include<vector>
#include<map>
#include"TObject.h"
#include"TGraph.h"
#include"TVector3.h"
#include "NPParticle.h"
class SamuraiFieldMap{

  public:
    SamuraiFieldMap(){m_BrhoScan=NULL;};
    SamuraiFieldMap(std::string file);
    ~SamuraiFieldMap(){};
  
  public: // Map reading
    void LoadMap(double angle, std::string file, unsigned int bin);

  private:
    void LoadAscii(std::string file);
    void LoadBinary(std::string file);

  private:
    // map[Pos]=B;
    std::map<std::vector<float>,std::vector<float>> m_field;
    float m_x_max,m_y_max,m_z_max,m_x_min,m_y_min,m_z_min;
    int m_bin;
    double m_angle;

  public: // getting the field at a point in space
    // return B at an existing point
    std::vector<float> GetB(std::vector<float>& pos); 
    inline std::vector<float> GetB(float x,float y ,float z){
      std::vector<float> pos = {x,y,z};
      return GetB(pos);
    };
    
    // interpolate B witin volume (0 outside volume)
    std::vector<float> InterpolateB(const std::vector<float>& pos);
    // interpolate B witin volume (0 outside volume)
    inline std::vector<float> InterpolateB(const TVector3& pos){
      std::vector<float> p={(float)pos.X(),(float)pos.Y(),(float)pos.Z()};
      return InterpolateB(p);
    };
 
  public: // Propagation of a particule in the field
    // return a 3D track of the particle in the field
    std::vector< TVector3 > Propagate(double rmax, double Brho, TVector3 pos, TVector3 dir);
    void func(NPL::Particle& N, TVector3 pos, TVector3 dir, TVector3& new_pos, TVector3& new_dir);
  private:
    double m_fdc2angle;
    double m_fdc2R;
  public:
    void SetFDC2Angle(double angle){m_fdc2angle=angle;};
    void SetFDC2R(double R){m_fdc2R=R;};
    TVector3 PropagateToFDC2(TVector3 pos, TVector3 dir);

  public:
    TGraph* BrhoScan(double min,double max,double step);
    double  FindBrho(TVector3& p_fdc0,TVector3& d_fdc0,TVector3& p_fdc2,TVector3& d_fdc2);

  private:
    TGraph* m_BrhoScan;
    //
    ClassDef(SamuraiFieldMap,1);
};

#endif
