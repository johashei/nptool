/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F, Flavigny   contact address: flavigny@lpccaen.in2p3.fr *
 *                                                                           *
 * Creation Date  : April 2021                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold RIBF PPAC treated data                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TBigRIPSPPACPhysics.h"

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
using namespace NPUNITS;
///////////////////////////////////////////////////////////////////////////

ClassImp(TBigRIPSPPACPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TBigRIPSPPACPhysics::TBigRIPSPPACPhysics(){
    m_EventData         = new TBigRIPSPPACData ;
    m_PreTreatedData    = new TBigRIPSPPACData ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
  }

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::BuildPhysicalEvent(){
  PreTreat();
  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::PreTreat(){

  ClearPreTreatedData();

  //static map<std::pair<int,int>, std::vector<int> > data ; 
  static map<int, BigRIPSPPACVariables > Fdata ; 
  //static map<std::pair<int,int>, bool > multiHit ; 

  //pair of detector ID and variable type (TX1,TX2,TY1,TY2,TA)=(0,1,2,3,4)
  std::pair<unsigned int, double> pair_id_type; 
  int id, TimeRaw;

  //pair_id_type.second = 0; //TX1
  unsigned int sizeTX1 = m_EventData->GetTX1Mult();
  for(unsigned int i = 0 ; i < sizeTX1 ; i++){
        id = m_EventData->GetTX1ID(i);
        TimeRaw = m_EventData->GetTX1(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
/*
            pair_id_type.first = id;
            if(data[pair_id_type].size()==0){
                data[pair_id_type].push_back(TimeRaw); 
            }else {
                multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
*/
            if(Fdata[id].FTX1.size()==0){
                Fdata[id].FTX1.push_back(TimeRaw*ch2ns_TX1[id]); 
            }else {
               // multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
        }
  }

//  pair_id_type.second = 1; //TX2
  unsigned int sizeTX2 = m_EventData->GetTX2Mult();
  for(unsigned int i = 0 ; i < sizeTX2 ; i++){
        id = m_EventData->GetTX2ID(i);
        TimeRaw = m_EventData->GetTX2(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
/*
            pair_id_type.first = id;
            if(data[pair_id_type].size()==0){
                data[pair_id_type].push_back(TimeRaw); 
            }else {
                multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
*/
            if(Fdata[id].FTX2.size()==0){
                Fdata[id].FTX2.push_back(TimeRaw*ch2ns_TX2[id]); 
            }else {
               // multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
        }
  }

//  pair_id_type.second = 2; //TY1
  unsigned int sizeTY1 = m_EventData->GetTY1Mult();
  for(unsigned int i = 0 ; i < sizeTY1 ; i++){
        id = m_EventData->GetTY1ID(i);
        TimeRaw = m_EventData->GetTY1(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
/*
            pair_id_type.first = id;
            if(data[pair_id_type].size()==0){
                data[pair_id_type].push_back(TimeRaw); 
            }else {
                multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
*/
            if(Fdata[id].FTY1.size()==0){
                Fdata[id].FTY1.push_back(TimeRaw*ch2ns_TY1[id]); 
            }else {
               // multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
        }
  }

// pair_id_type.second = 3; //TY2
  unsigned int sizeTY2 = m_EventData->GetTY2Mult();
  for(unsigned int i = 0 ; i < sizeTY2 ; i++){
        id = m_EventData->GetTY2ID(i);
        TimeRaw = m_EventData->GetTY2(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
/*
            pair_id_type.first = id;
            if(data[pair_id_type].size()==0){
                data[pair_id_type].push_back(TimeRaw); 
            }else {
                multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
*/
            if(Fdata[id].FTY2.size()==0){
                Fdata[id].FTY2.push_back(TimeRaw*ch2ns_TY2[id]); 
            }else {
               // multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
        }
  }

// pair_id_type.second = 4; //TA
  unsigned int sizeTA = m_EventData->GetTAMult();
  for(unsigned int i = 0 ; i < sizeTA ; i++){
        id = m_EventData->GetTAID(i);
        TimeRaw = m_EventData->GetTA(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
/*
            pair_id_type.first = id;
            if(data[pair_id_type].size()==0){
                data[pair_id_type].push_back(TimeRaw); 
            }else {
                multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }
*/
            if(Fdata[id].FTA.size()==0){
                Fdata[id].FTA.push_back(TimeRaw*ch2ns_TA[id]); 
            }else {
               // multiHit[pair_id_type]=true; //multiple value for same TDC ch.
            }

        }
  }

/*  
  for(auto it = data.begin();it!=data.end();++it){
    pair_id_type = it->first;
    cout<< "ID:"<<pair_id_type.first<<"\tType:"<<pair_id_type.second<<" \tT:"<<data[pair_id_type][0]<<endl; 
  }
*/
  int j=0;
    BigRIPSPPACVariables toto;
  for(auto it = Fdata.begin();it!=Fdata.end();++it){
    id = it->first;
    ID.push_back(id);
    FP.push_back(FPL[id]);
    toto.Clear();
    toto = it->second;
    //cout<< "ID:"<<id<<endl;
    //cout<< "ch2nsX1:"<<ch2ns_TX1[id]<<endl;
    //cout<< "ch2nsX2:"<<ch2ns_TX2[id]<<endl;
    //cout<< "ch2nsY1:"<<ch2ns_TY1[id]<<endl;
    //cout<< "ch2nsY2:"<<ch2ns_TY2[id]<<endl;
    //cout<< "ch2nsA:"<<ch2ns_TA[id]<<endl;
    //toto.Print(); 
    for (int i=0; i<toto.FTX1.size(); i++) TX1.push_back(toto.FTX1[i]);
    for (int i=0; i<toto.FTX2.size(); i++) TX2.push_back(toto.FTX2[i]);
    for (int i=0; i<toto.FTY1.size(); i++) TY1.push_back(toto.FTY1[i]);
    for (int i=0; i<toto.FTY2.size(); i++) TY2.push_back(toto.FTY2[i]);
    for (int i=0; i<toto.FTA.size(); i++) TA.push_back(toto.FTA[i]);
    //TX1=toto.FTX1;
    //TX2=toto.FTX2;
    //TY1=toto.FTY1;
    //TY2=toto.FTY2;
    //TA=toto.FTA;

    //if(toto.HasTXs()) {TDiffX.push_back(TX1[j]-TX2[j]);}
    if(toto.HasTXs()) {TDiffX.push_back(toto.FTX1[0]-toto.FTX2[0]);}
    else {TDiffX.push_back(-99999);}
    //if(toto.HasTXs() && toto.HasTA()) {TSumX.push_back(TX1[j]+TX2[j]-2*TA[j]);}
    if(toto.HasTXs() && toto.HasTA()) {TSumX.push_back(toto.FTX1[0]+toto.FTX2[0]-2*toto.FTA[0]);}
    else {TSumX.push_back(-99999);}
    //if(toto.HasTYs()) {TDiffY.push_back(TY1[j]-TY2[j]);}
    if(toto.HasTYs()) {TDiffY.push_back(toto.FTY1[0]-toto.FTY2[0]);}
    else {TDiffY.push_back(-99999);}
    //if(toto.HasTYs() && toto.HasTA()) {TSumY.push_back(TY1[j]+TY2[j]-2*TA[j]);}
    if(toto.HasTYs() && toto.HasTA()) {TSumY.push_back(toto.FTY1[0]+toto.FTY2[0]-2*toto.FTA[0]);}
    else {TSumY.push_back(-99999);}
    j++;
  }
  //Print(); 
  Fdata.clear();

  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::Clear(){
    TX1.clear();
    TX2.clear();
    TY1.clear();
    TY2.clear();
    TA.clear();
    TSumX.clear();
    TDiffX.clear();
    TSumY.clear();
    TDiffY.clear();
    ID.clear();
    FP.clear();
    //Data.clear();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::Print(){
    cout << "XXXXXXXXXXXXXXXXXXXXXXXX PPAC Physics Event XXXXXXXXXXXXXXXXX" << endl;
    cout << "TX1_Mult = " << TX1.size();
    for (UShort_t i = 0; i < TX1.size(); i++){cout << "\tTX1: " << TX1[i] << endl;}
    cout << "TX2_Mult = " << TX2.size();
    for (UShort_t i = 0; i < TX2.size(); i++){cout << "\tTX2: " << TX2[i] << endl;}
    cout << "TY1_Mult = " << TY1.size();
    for (UShort_t i = 0; i < TY1.size(); i++){cout << "\tTY1: " << TY1[i] << endl;}
    cout << "TY2_Mult = " << TY2.size();
    for (UShort_t i = 0; i < TY2.size(); i++){cout << "\tTY2: " << TY2[i] << endl;}
    cout << "TA_Mult = " << TA.size();
    for (UShort_t i = 0; i < TA.size(); i++){cout << "\tTA: " << TA[i] << endl;}
    cout << "TDiffX_Mult = " << TDiffX.size();
    for (UShort_t i = 0; i < TDiffX.size(); i++){cout << "\tTDiffX: " << TDiffX[i] << endl;}
    cout << "TDiffY_Mult = " << TDiffY.size();
    for (UShort_t i = 0; i < TDiffY.size(); i++){cout << "\tTDiffY: " << TDiffY[i] << endl;}
    cout << "TSumX_Mult = " << TSumX.size();
    for (UShort_t i = 0; i < TSumX.size(); i++){cout << "\tTSumX: " << TSumX[i] << endl;}
    cout << "TSumY_Mult = " << TSumY.size();
    for (UShort_t i = 0; i < TSumY.size(); i++){cout << "\tTSumY: " << TSumY[i] << endl;}
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("BigRIPSPPAC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " PPAC detector file found " << endl; 

  vector<string> token= {"XML"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    cout << endl << "////  BigRIPSPPAC file (" << i+1 << ")" << endl;
    string xmlpath = blocks[i]->GetString("XML");
    NPL::XmlParser xml;
    xml.LoadFile(xmlpath);
    AddPPACs("BigRIPSPPAC",xml);
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::AddPPACs(string name, NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName(name);  

  if(name=="BigRIPSPPAC"){
    unsigned int size = b.size();
    for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = b[i]->AsInt("ID"); 
      //string sDir = b[i]->AsString("anodedir");
      RawLowerLimit[ID] = b[i]->AsInt("tdc_underflow"); 
      RawUpperLimit[ID] = b[i]->AsInt("tdc_overflow"); 
      FPL[ID] = b[i]->AsInt("FPL"); 
      ch2ns_TX1[ID] = b[i]->AsDouble("x1_ch2ns"); 
      ch2ns_TX2[ID] = b[i]->AsDouble("x2_ch2ns"); 
      ch2ns_TY1[ID] = b[i]->AsDouble("y1_ch2ns"); 
      ch2ns_TY2[ID] = b[i]->AsDouble("y2_ch2ns"); 
      ch2ns_TA[ID] = b[i]->AsDouble("a_ch2ns"); 
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitSpectra(){  
  //m_Spectra = new TBigRIPSPPACSpectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::FillSpectra(){  
  //  m_Spectra -> FillRawSpectra(m_EventData);
  //  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::CheckSpectra(){  
  //  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TBigRIPSPPACPhysics::GetSpectra() {
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
void TBigRIPSPPACPhysics::WriteSpectra(){
  // m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::AddParameterToCalibrationManager(){
  //CalibrationManager* Cal = CalibrationManager::getInstance();

  // each layer
  //for( int l = 0 ; l < 14 ; ++l){
  //  Cal->AddParameter("SamuraiFDC2", "L"+ NPL::itoa(l),"FDC2_L"+ NPL::itoa(l));
  //}

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "BigRIPSPPAC" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fPPAC_*"))
    inputChain->SetBranchStatus( "fPPAC_*",true);
  inputChain->SetBranchAddress( "BigRIPSPPAC" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitializeRootInputPhysics(){
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch( "PPAC" , "TBigRIPSPPACPhysics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TBigRIPSPPACPhysics::Construct(){
  return (NPL::VDetector*) new TBigRIPSPPACPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_bigripsPPAC{
    public:
      proxy_bigripsPPAC(){
        NPL::DetectorFactory::getInstance()->AddToken("BigRIPSPPAC","BigRIPS");
        NPL::DetectorFactory::getInstance()->AddDetector("BigRIPSPPAC",TBigRIPSPPACPhysics::Construct);
      }
  };

  proxy_bigripsPPAC p_bigripsPPAC;
}

