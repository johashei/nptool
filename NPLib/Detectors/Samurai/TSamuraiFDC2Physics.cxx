/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiFDC2 treated data                                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TSamuraiFDC2Physics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "NPOptionManager.h"
#include "NPDetectorFactory.h"
#include "NPSystemOfUnits.h"
//   ROOT
#include "TChain.h"
#include "TVector2.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
using namespace NPUNITS;
///////////////////////////////////////////////////////////////////////////

ClassImp(TSamuraiFDC2Physics)
  ///////////////////////////////////////////////////////////////////////////
  TSamuraiFDC2Physics::TSamuraiFDC2Physics(){
    m_EventData         = new TSamuraiFDC2Data ;
    m_PreTreatedData    = new TSamuraiFDC2Data ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
    ToTThreshold = 180;
  }

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::BuildPhysicalEvent(){
  PreTreat();
  RemoveNoise();

  // Map[detector & plane angle, vector of spatial information]
  static map<std::pair<unsigned int,double>, vector<double> > X ; 
  static map<std::pair<unsigned int,double>, vector<double> > Z ; 
  static map<std::pair<unsigned int,double>, vector<double> > R ; 
  X.clear();Z.clear();R.clear();
  unsigned int size = Detector.size();
  for(unsigned int i = 0 ; i < size ; i++){
    if(DriftLength[i]>0.1){
      int det = Detector[i];
      int layer = Layer[i];
      int wire = Wire[i]; 
      SamuraiDCIndex idx(det,layer,wire);
      std::pair<unsigned int, double> p(det,Wire_Angle[idx]);
      X[p].push_back(Wire_X[idx]); 
      Z[p].push_back(Wire_Z[idx]); 
      R[p].push_back(DriftLength[i]); 
    }
  }

  // Reconstruct the vector for each of the plane of each of the detector
  static double X0,X100;
  static map<std::pair<unsigned int,double>, TVector3 > VX0 ;  
  static map<std::pair<unsigned int,double>, TVector3 > VX100 ;  
  static map<std::pair<unsigned int,double>, int > MultPlane ;  
  VX0.clear();VX100.clear();
  for(auto it = X.begin();it!=X.end();++it){
    double a = Track2D(X[it->first],Z[it->first],R[it->first],X0,X100); 
    // very small a means track perpendicular to the chamber, what happen when there is pile up

    if(abs(a)>1000)
      PileUp++;

    MultPlane[it->first] = X[it->first].size() ;
    Mult+=X[it->first].size();
    // Position at z=0
    TVector3 P(X0,0,0);
    P.RotateZ(it->first.second);
    VX0[it->first]=P;
    // Direction of the vector in the plane
    TVector3 D = TVector3(X100,0,0);
    D.RotateZ(it->first.second);
    VX100[it->first]=D;

  }
  // Reconstruct the central position (z=0) for each detector
  static map<unsigned int,vector<TVector3> > C ;  
  C.clear();
  TVector3 P;

  for(auto it1 = VX0.begin();it1!=VX0.end();++it1){
    for(auto it2 = it1;it2!=VX0.end();++it2){
      if(it1!=it2 && it1->first.first==it2->first.first){// different plane, same detector
        ResolvePlane(it1->second,it1->first.second,it2->second,it2->first.second,P);
        if(P.X()!=-10000)
          C[it1->first.first].push_back(P);
      }
    }
  }
  // Reconstruct the position at z=100 for each detector
  static map<unsigned int,vector<TVector3> > C100 ;  
  C100.clear();
  for(auto it1 = VX100.begin();it1!=VX100.end();++it1){
    for(auto it2 = it1;it2!=VX100.end();++it2){
      if(it1!=it2 && it1->first.first==it2->first.first){// different plane, same detector
        ResolvePlane(it1->second,it1->first.second,it2->second,it2->first.second,P);

        if(P.X()!=-10000)
          C100[it1->first.first].push_back(P);
      }
    }
  }
  // Build the Reference position by averaging all possible pair 
  size = C[2].size();
  double PosX100,PosY100;
  if(size){
    PosX=0;
    PosY=0;
    PosX100=0;
    PosY100=0;
    for(unsigned int i = 0 ; i < size ; i++){
      PosX+= C[2][i].X(); 
      PosY+= C[2][i].Y(); 
      PosX100+= C100[2][i].X(); 
      PosY100+= C100[2][i].Y(); 
    //       cout << C[2][i].X() << " (" << C[2][i].Y() << ") ";
    } 
   // cout << endl;
    MultMean=size;
    // Mean position at Z=0
    PosX=PosX/size; 
    PosY=PosY/size; 
    // Mean position at Z=100
    PosX100=PosX100/size; 
    PosY100=PosY100/size; 
    // Compute ThetaX, angle between the Direction vector projection in XZ with
    // the Z axis
    ThetaX=atan((PosX100-PosX)/100.);
    // Compute PhiY, angle between the Direction vector projection in YZ with
    // the Z axis
    PhiY=atan((PosY100-PosY)/100.);
    Dir=TVector3(PosX100-PosX,PosY100-PosY,100).Unit();
  }

  return;
}
////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::ResolvePlane(const TVector3& L,const double& ThetaU ,const TVector3& H, const double& ThetaV, TVector3& PosXY){
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
};
///////////////////////////////////////////////////////////////////////////
double TSamuraiFDC2Physics::Track2D(const vector<double>& X,const vector<double>& Z,const vector<double>& R,double& X0,double& X100 ){
  fitX=X;
  fitZ=Z;
  fitR=R;
  // assume all X,Z,R of same size
  unsigned int size = X.size();
/*
    ofstream out("data.txt");
      for(unsigned int i = 0 ; i < size ; i++){
      out << fitX[i] << " " << fitZ[i] << " "<< fitR[i] << endl;
      }
      out.close();
  */      
  // Define the starting point of the fit: a straight line passing through the 
  // the first and last wire
  // z = ax+b -> x=(z-b)/a
  double a = (fitZ[size-1]-fitZ[0])/(fitX[size-1]-fitR[size-1]-fitX[0]-fitR[0]);
  double b = fitZ[0]-a*(fitX[0]+fitR[0]);
  double parameter[2]={a,b};
/*    out.open("line.txt");
      out << a << " " << b <<endl;
      out.close();
      */
  //cout << a << " - " << b << " " << SumD(parameter) <<endl;

  static ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  static ROOT::Math::Functor f(this,&TSamuraiFDC2Physics::SumD,2); 
  min->SetFunction(f);
  min->SetVariable(0,"a",parameter[0],1000);
  min->SetVariable(1,"b",parameter[1],1000);
  min->SetTolerance(0.1);
  // cout <<"start " << endl;
  min->Minimize(); 
  //cout << "start" << endl;
  const double *xs = min->X();
  X0=-xs[1]/xs[0];
  X100=(100-xs[1])/xs[0];
/*  
     out.open("fit.txt");
     out << xs[0] << " " << xs[1] <<endl;
     out.close();
     if(X0==0){
     exit(0);
     }
  */   
  return xs[0];
  //  cout << " done " << SumD(xs) <<endl;  
  //cout << xs[0] << " / " << xs[1] << " " << SumD(xs) <<endl;
}
////////////////////////////////////////////////////////////////////////////////
double TSamuraiFDC2Physics::SumD(const double* parameter){
  //cout << " "<<parameter[0] << " h " << parameter[1] ;
  unsigned int size = fitX.size();
  // Compute the sum P of the distance between the circle and the track
  double P = 0;
  double a = parameter[0];
  double b = parameter[1];
  double ab= a*b;
  double a2=a*a;
  double c,d,x,z;

  for(unsigned int i = 0 ; i < size ; i++){
    c = fitX[i];
    d = fitZ[i];
    x = (a*d-ab+c)/(1+a2);
    z = a*x+b;
    //P+= sqrt(abs( (x-c)*(x-c)+(z-d)*(z-d)-fitR[i]*fitR[i]))/fitR[i];
    P+= abs( (x-c)*(x-c)+(z-d)*(z-d)-fitR[i]*fitR[i])/fitR[i];
    //  cout << sqrt(abs((x-c)*(x-c)+(z-d)*(z-d)-fitR[i]*fitR[i] ))/fitR[i]<< " ";
  }
  //cout << endl;
  return P;
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::PreTreat(){
  ClearPreTreatedData();
  static CalibrationManager* Cal = CalibrationManager::getInstance();
  static string channel;
  // one map per detector
  map<unsigned int, vector<double> > X ; 
  map<unsigned int, vector<double> > Z ; 
  map<unsigned int, vector<double> > R ; 

  unsigned int size = m_EventData->Mult();
  for(unsigned int i = 0 ; i < size ; i++){
    // EDGE=1 is the leading edge, IE, the real time.
    // EDGE=0 is the trailing edge, so it helps build Tot
    if(m_EventData->GetEdge(i)==1){
      int det   = m_EventData->GetDetectorNbr(i); 
      int layer = m_EventData->GetLayerNbr(i); 
      int wire  = m_EventData->GetWireNbr(i); 
      double time = m_EventData->GetTime(i);
      double etime = 0;
      // look for matching trailing edge   
      for(unsigned int j = 0 ; j < size ; j++){
        if(m_EventData->GetEdge(j)==0){
          int edet   = m_EventData->GetDetectorNbr(j); 
          int elayer = m_EventData->GetLayerNbr(j); 
          int ewire  = m_EventData->GetWireNbr(j); 
          // same wire
          if(wire==ewire && layer==elayer && det==edet){
            etime = m_EventData->GetTime(j); 
          }    
        }
        if(etime && etime>time)
          break;
        else
          etime=0;
      }
      // a valid wire must have an edge
      if(etime && time && etime-time>ToTThreshold){
        Detector.push_back(det);
        Layer.push_back(layer);       
        Wire.push_back(wire);
        Time.push_back(time);
        ToT.push_back(etime-time);
        channel="SamuraiFDC2/L" + NPL::itoa(layer);
        DriftLength.push_back(Cal->ApplySigmoid(channel,etime));
      }
    }

  }
  return;
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::RemoveNoise(){
  // Remove the noise by looking if a matching wire exist in the adjacent layer
  // this done by looking at the closest plane with the same orientation
  unsigned int size = Detector.size(); 
  for(unsigned int i = 0 ; i < size ; i++){
    bool match=false;
    int det = Detector[i];
    int layer = Layer[i];
    int wire = Wire[i];
    // look for matching adjacent wire   

    for(unsigned int j = 0 ; j < size ; j++){
      int adet = Detector[j];
      int alayer = Layer[j];
      int awire = Wire[j];
      bool blayer = false;
      if(layer%2==0){
        if(layer+1==alayer)
          blayer=true;
      }

      else{
        if(layer-1==alayer)
          blayer=true;
      }

      if(det==adet && blayer && abs(wire-awire)<=1){
        match=true;
        break;
      }
    }

    if(match)
      Matched.push_back(true);
    else
      Matched.push_back(false);
  }
  return;
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::Clear(){
  MultMean=0;
  PileUp=0;
  Mult=0;
  PosX=PosY=-10000;
  DriftLength.clear();
  Detector.clear();
  Layer.clear();
  Wire.clear();
  Time.clear();
  ToT.clear();
  ParticleDirection.clear();
  MiddlePosition.clear();
  Matched.clear();
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SAMURAIFDC2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    cout << endl << "////  Samurai FDC2 (" << i+1 << ")" << endl;
    string xmlpath = blocks[i]->GetString("XML");
    NPL::XmlParser xml;
    xml.LoadFile(xmlpath);
    AddDC("SAMURAIFDC2",xml);
  }
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::AddDC(string name, NPL::XmlParser& xml){
  std::vector<NPL::block*> b = xml.GetAllBlocksWithName(name);  
  // FDC2 case
  if(name=="SAMURAIFDC2"){
    unsigned int det=2;
    unsigned int size = b.size();
    for(unsigned int i = 0 ; i < size ; i++){
      unsigned int layer = b[i]->AsInt("layer"); 
      unsigned int wire  = b[i]->AsInt("wireid"); 
      double X = b[i]->AsDouble("wirepos");  
      double Z = b[i]->AsDouble("wirez");  
      string sDir = b[i]->AsString("anodedir");
      double T=0;
      if(sDir=="X")
        T=0*deg;
      else if(sDir=="Y")
        T=90*deg;
      else if(sDir=="U")
        T=30*deg;
      else if(sDir=="V")
        T=-30*deg;
      else{
        cout << "ERROR: Unknown layer orientation for Samurai FDC2"<< endl;
        exit(1);
      }
      SamuraiDCIndex idx(det,layer,wire);
      Wire_X[idx]=X;
      Wire_Z[idx]=Z;
      Wire_Angle[idx]=T;
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitSpectra(){  
  //m_Spectra = new TSamuraiFDC2Spectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::FillSpectra(){  
  //  m_Spectra -> FillRawSpectra(m_EventData);
  //  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::CheckSpectra(){  
  //  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TSamuraiFDC2Physics::GetSpectra() {
  /*  if(m_Spectra)
      return m_Spectra->GetMapHisto();
      else{
      map< string , TH1*> empty;
      return empty;
      }*/
  map< string , TH1*> empty;
  return empty;

} 

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::WriteSpectra(){
  // m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::AddParameterToCalibrationManager(){
  CalibrationManager* Cal = CalibrationManager::getInstance();

  // each layer
  for( int l = 0 ; l < 14 ; ++l){
    Cal->AddParameter("SamuraiFDC2", "L"+ NPL::itoa(l),"FDC2_L"+ NPL::itoa(l));
  }

}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "SamuraiFDC2" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fDC_*"))
    inputChain->SetBranchStatus( "fDC_*",true);
  inputChain->SetBranchAddress( "SamuraiFDC2" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitializeRootInputPhysics(){
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch( "SamuraiFDC2" , "TSamuraiFDC2Physics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSamuraiFDC2Physics::Construct(){
  return (NPL::VDetector*) new TSamuraiFDC2Physics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_samuraiFDC2{
    public:
      proxy_samuraiFDC2(){
        NPL::DetectorFactory::getInstance()->AddToken("SAMURAIFDC2","Samurai");
        NPL::DetectorFactory::getInstance()->AddDetector("SAMURAIFDC2",TSamuraiFDC2Physics::Construct);
      }
  };

  proxy_samuraiFDC2 p_samuraiFDC2;
}

