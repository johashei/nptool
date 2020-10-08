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
using namespace NPUNITS;
///////////////////////////////////////////////////////////////////////////

ClassImp(TSamuraiFDC2Physics)
  ///////////////////////////////////////////////////////////////////////////
  TSamuraiFDC2Physics::TSamuraiFDC2Physics(){
    m_EventData         = new TSamuraiFDC2Data ;
    m_PreTreatedData    = new TSamuraiFDC2Data ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
    ToTThreshold = 0;
  }

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::BuildPhysicalEvent(){
  PreTreat();
  RemoveNoise();
 
  // one map per detector
  map<unsigned int, vector<double> > X ; 
  map<unsigned int, vector<double> > Z ; 
  map<unsigned int, vector<double> > R ; 

  unsigned int size = Detector.size();
  for(unsigned int i = 0 ; i < size ; i++){
      
    
    }

  return;
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::Track2D(const vector<double>& X,const vector<double>& Z,const vector<double>& R,double& dirX, double& dirZ,double& refX ){
  
  }

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::PreTreat(){
  ClearPreTreatedData();
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
        SamuraiDCIndex idx(det,layer,wire);
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
     Wire_T[idx]=T;
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

