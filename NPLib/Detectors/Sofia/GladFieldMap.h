#ifndef GladFieldMap_h
#define GladFieldMap_h
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2022                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold a basic field map and is used for Brho reconstruction    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
#include <map>
using namespace std;

// ROOT
#include "TObject.h"
#include "TGraph.h"
#include "TVector3.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"

#include "NPParticle.h"

#define NX 80;
#define NY 40;
#define NZ 120;

class GladFieldMap{
  public: 
    GladFieldMap();
    ~GladFieldMap();
    
  private:
    // GLAD parameters
    map<vector<double>, vector<double>> m_field;
    vector<double> m_Bx;
    vector<double> m_By;
    vector<double> m_Bz;
    double m_Z_Glad;
    double m_Leff;
    double m_Tilt;
    int m_Nx;
    int m_Ny;
    int m_Nz;
    double m_x_min;
    double m_y_min;
    double m_z_min;
    double m_x_max;
    double m_y_max;
    double m_z_max;
    double m_R_max;
    double m_Scale;
    int m_bin;
  private:
    // MWPC3 paramters
    double m_CentralTheta;
    double m_X_MWPC3;
    double m_Z_MWPC3;
    double m_R_MWPC3;
    double m_Angle_MWPC3;
  private:
    // Runge-Kunta 4 paramaters
    double m_dt;
    int m_Limit;
    double m_Zmax;

  private:
    TVector3 m_InitPos;
    TVector3 m_InitDir;
    TVector3 m_FinalPos;

  public:
    void SetZGlad(double val) {m_Z_Glad = val;}
    void SetLeff(double val) {m_Leff = val;}
    void SetGladTiltAngle(double val) {m_Tilt = val;}
    void SetScale(double val) {m_Scale = val;}

    void SetCentralTheta(double val) {m_CentralTheta = val;}
    void SetX_MWPC3(double val) {m_X_MWPC3 = val;}
    void SetZ_MWPC3(double val) {m_Z_MWPC3 = val;}
  
    void SetPropagationTimeInterval(double val) {m_dt = val;}
    void SetLimit(int val) {m_Limit = val;}
    void SetPropagationMaxZ(double val) {m_Zmax = val;}

    void SetInitPos(TVector3 Pos) {m_InitPos = Pos;}
    void SetInitDir(TVector3 Dir) {m_InitPos = Dir;}

  public:
    void LoadMap(string filename);
    vector<double> InterpolateB(const vector<double>& pos);
    inline vector<double> InterpolateB(const TVector3& pos){
      vector<double> position={(double)pos.X(),(double)pos.Y(),(double)pos.Z()};
      return InterpolateB(position);
    };
    double GetB(TVector3 localpoint, string field_component);
    TGraph* BrhoScan(double Brho_min, double Brho_max, double Brho_step, TVector3 pos, TVector3 dir);
    TVector3 CalculateIntersectionPoint(vector<TVector3> vPos);
    vector<TVector3> Propagate(double Brho, TVector3 Pos, TVector3 Dir, bool store);
    TVector3 PropagateToMWPC(TVector3 pos, TVector3 dir);
    void func(NPL::Particle& N, TVector3 Pos, TVector3 Imp, TVector3& xk, TVector3& pk);
    double FindBrho(TVector3 Pos_init, TVector3 Dir_init, TVector3 Pos_final);

  private:
    TGraph* m_BrhoScan;
    ROOT::Math::Minimizer* m_min;
    ROOT::Math::Functor m_func;
    double Delta(const double* parameter);
    bool IsInside(TVector3 localpoint, int& ix, int& iy, int& iz, double& dx, double& dy, double& dz);
    double Interpolate(double dx, double dy, double dz);

    double m_C0[2][2][2]; //!
    double m_C1[2][2]; //!
    double m_C2[2]; //!

  private:

    ClassDef(GladFieldMap,1)  
};

#endif
