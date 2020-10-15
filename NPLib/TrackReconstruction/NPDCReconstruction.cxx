#include"NPDCReconstruction.h"
#include "Math/Factory.h"
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta   contact address: matta@lpcaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : October 2020                                            *
 * Last update     : October 2020                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class have all the method needed to analyse Drift Chambers          *
 *****************************************************************************/

using namespace std;
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
DCReconstruction::DCReconstruction(){
  m_min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  m_func=ROOT::Math::Functor(this,&NPL::DCReconstruction::SumD,2); 
  }
////////////////////////////////////////////////////////////////////////////////
DCReconstruction::~DCReconstruction(){
  delete m_min;
  }

////////////////////////////////////////////////////////////////////////////////
void DCReconstruction::BuildTrack2D(const vector<double>& X,const vector<double>& Z,const vector<double>& R,double& X0,double& X100,double& a, double& b ){
  fitX=&X;
  fitZ=&Z;
  fitR=&R;
  // assume all X,Z,R of same size
  unsigned int size = X.size();
  // Define the starting point of the fit: a straight line passing through the 
  // the first and last wire
  // z = ax+b -> x=(z-b)/a
  double ai = (Z[size-1]-Z[0])/(X[size-1]-R[size-1]-X[0]-R[0]);
  double bi = Z[0]-ai*(X[0]+R[0]);
  double parameter[2]={ai,bi};

  m_min->SetFunction(m_func);
  m_min->SetVariable(0,"a",parameter[0],1000);
  m_min->SetVariable(1,"b",parameter[1],1000);
  m_min->SetTolerance(0.1);
  // Perform minimisation
  m_min->Minimize(); 
  // access set of parameter that produce the minimum
  const double *xs = m_min->X();
  a=xs[0];
  b=xs[1];
  X0=-b/a;
  X100=(100-b)/a;
}

////////////////////////////////////////////////////////////////////////////////
void DCReconstruction::ResolvePlane(const TVector3& L,const double& ThetaU ,const TVector3& H, const double& ThetaV, TVector3& PosXY){
  // direction of U and V wire
  TVector3 u = TVector3(0,1,0);
  u.RotateZ(ThetaU);

  TVector3 v = TVector3(0,1,0);
  v.RotateZ(ThetaV);


  // Compute the coeff of the two line of vecotr u (v) going through H (L)
  // dv : y = av*x+bv
  double av = v.Y()/v.X();
  double bv = H.Y() - av*H.X();

  // du : y = au*x+bv
  double au = u.Y()/u.X();
  double bu = L.Y() - au*L.X();

  // We look for M(xM, yM) that intersect du and dv:
  double xM,yM;
  if(!isinf(au) && !isinf(av)){ // au and av are not inf, i.e. not vertical line
    xM = (bv-bu)/(au-av);
    yM = au*xM+bu;
  }
  else if(isinf(av)){// av is inf, so v is along Y axis, H is direct measure of X
    xM = H.X();
    yM = au*xM+bu;
  }
  else if (isinf(au)){//au is inf, so u is along Y axis, L is direct measure of X
    xM = L.X();
    yM = av*xM+bv;
  }
  else{ // all is lost
    xM=-10000;
    yM=-10000;
  }

  PosXY=TVector3(xM,yM,0);

}


////////////////////////////////////////////////////////////////////////////////
double DCReconstruction::SumD(const double* parameter ){
  //cout << " "<<parameter[0] << " h " << parameter[1] ;
  unsigned int size = fitX->size();
  // Compute the sum P of the distance between the circle and the track
  double P = 0;
  double a = parameter[0];
  double b = parameter[1];
  double ab= a*b;
  double a2=a*a;
  double c,d,x,z,r;

  for(unsigned int i = 0 ; i < size ; i++){
    c = (*fitX)[i];
    d = (*fitZ)[i];
    r = (*fitR)[i];
    x = (a*d-ab+c)/(1+a2);
    z = a*x+b;
    P+= abs( (x-c)*(x-c)+(z-d)*(z-d)-r*r)/r;
  }
  return P;

}


