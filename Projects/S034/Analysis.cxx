/*****************************************************************************
 * Copyright (C) 2009-2021    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta contact address: matta@lpccaen.in2p3.fr         *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  S034 analysis project                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/



#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"RootInput.h"
#include"RootOutput.h"
#include"NPParticle.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  Minos= (TMinosPhysics*) m_DetectorManager->GetDetector("Minos");
  Nebula = (TNebulaPhysics*) m_DetectorManager->GetDetector("NEBULA");
  BDC = (TSamuraiBDCPhysics*) m_DetectorManager->GetDetector("SAMURAIBDC");
  FDC0 = (TSamuraiFDC0Physics*) m_DetectorManager->GetDetector("SAMURAIFDC0");
  FDC2 = (TSamuraiFDC2Physics*) m_DetectorManager->GetDetector("SAMURAIFDC2");
  Hodo = (TSamuraiHodoscopePhysics*) m_DetectorManager->GetDetector("SAMURAIHOD");

  FragmentTarget = NPL::EnergyLoss("He6_LH2.G4table","G4Table",1000 );
  m_field.LoadMap(30*deg,"field_map/180702-2,40T-3000.table.bin",10);
  m_field.SetFDC2Angle((59.930-90.0)*deg);
  m_field.SetFDC2R(FDC2->GetOffset().Z());
  InitOutputBranch();
  InitInputBranch();
  // for fdc/bdc alignement
  //file.open("Calibration/Pos/bdc.txt");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  Clear();
  
  static NPL::Particle He6("6He");
  static NPL::Particle He8("8He");
  static NPL::Particle n("neutron");
  static double mhe6 = He6.Mass();
  static double mn   = n.Mass();
  static double sumM= mhe6+mn;
  static TLorentzVector LVHe6;
  static TLorentzVector LVn;
  He8.SetBeta(0.5152);
  Trigger=Trigger&0x00ff;
  if( FDC2->PosX>-1500 && FDC2->PosX<1000 
      && FDC2->PosY>-500 && FDC2->PosY<500 
      && FDC0->PosX>-80 && FDC0->PosX<80 
      && FDC0->PosY>-80 && FDC0->PosY<80 // both FDC ok
      && Minos->Tracks_P0.size()>0     ) {  // p,pn or p,2p
    // Compute ThetaX and PhiY using Minos vertex and FDC0 X
    // Check if both BDC are reconstructed
    TVector3 BDC1=BDC->GetPos(1);
    TVector3 BDC2=BDC->GetPos(2);

    if( BDC1.Z()!=-10000 && BDC2.Z()!=-10000){
      TVector3 Vertex,delta;
      TVector3 P1 = Minos->Tracks_P0[0]+Minos->Tracks_Dir[0];

      // Vertex in Samurai Reference frame 
      MinimumDistanceTwoLines(BDC1,BDC2, 
          Minos->Tracks_P0[0], P1,
          Vertex, delta) ;
     
      TVector3 FDC0_Dir= FDC0->GetPos()-Vertex;
      FDC0_Dir=FDC0_Dir.Unit();
      TVector3 BDCDir=BDC2-BDC1;
      BDCDir=BDCDir.Unit();
      BDCDir*=(Vertex.Z()-BDC2.Z())/BDCDir.Z();
      BDCX=(BDC2+BDCDir).X();
      BDCY=(BDC2+BDCDir).Y();

      // XYZ relative to Minos entrance
      X=Vertex.X();
      Y=Vertex.Y();
      Z=Vertex.Z()+4650;

      // Fragment analysis 
      if(FDC0_Dir.Z()>0.6 && Z>0 && Z<150 && sqrt(X*X+Y*Y)<15){
        double FDC0_ThetaX = atan((FDC0->PosX-Vertex.X())/(1254.39-Z));
        double FDC0_PhiY   = atan((FDC0->PosY-Vertex.Y())/(1254.39-Z));
        Brho=m_field.FindBrho(FDC0->GetPos(),FDC0_Dir,FDC2->GetPos(),TVector3(0,0,1));
        He6.SetBrho(Brho);
        double Energy = He6.GetEnergy();
        if(Energy){
          Energy=FragmentTarget.EvaluateInitialEnergy(Energy,150-Z,FDC0_Dir.Angle(TVector3(0,0,1)));
          He6.SetKineticEnergy(Energy);
          Beta_f=He6.GetBeta();
          LVHe6.SetVectM(TVector3(0,0,0),mhe6);
          LVHe6.Boost(Beta_f*FDC0_Dir.Unit());
        }
        
      }
      
      // Neutron analysis
      if(Nebula->Charge.size()>0 && !Nebula->HasVeto()){
        // Index of first neutron hit
        unsigned int first = Nebula->GetFirstHit();
        TVector3 Pfirst = (Nebula->GetPos(first)-Vertex);
        double L = Pfirst.Mag();
        double TSBT= (Vertex.Z()+7377.56)/He8.GetVelocity();
        double TOF = Nebula->TOF[first]-TSBT;
        Beta_n = (L/TOF)/NPUNITS::c_light;
        LVn.SetVectM(TVector3(0,0,0),mn);
        LVn.Boost(Beta_n*Pfirst.Unit());
      }
      
      // a full event
      if(Beta_n&&Beta_f){
        He6.SetEnergyImpulsion(LVHe6);
        n.SetEnergyImpulsion(LVn);
        Erel = (LVHe6+LVn).Mag()-sumM;
      }
      // Calib//////////////////////////////////////////////////////////////////
 /*     static int count=0;
      if(Minos->Delta_Vertex < 5 && FDC2->PosX-252.55>0&&FDC0->GetPos().X()>-10000 && FDC0->Dir.Z()>0.9 && Minos->Z_Vertex>0&& sqrt(Minos->X_Vertex*Minos->X_Vertex+Minos->Y_Vertex*Minos->Y_Vertex)<15){
        file << FDC0->GetPos().X()   <<" " << FDC0->GetPos().Y() << " " << FDC0->GetPos().Z() <<" " ;
        file << Minos->X_Vertex      <<" " << Minos->Y_Vertex    << " " << Minos->Z_Vertex    << " " ;
        file << FDC2->GetPos().X()   <<" " << FDC2->GetPos().Y() << " " << FDC2->GetPos().Z() <<" " << FDC2->Dir.X() <<" " << FDC2->Dir.Y() << " " << FDC2->Dir.Z()<< endl;
        count ++;
      }
      if(count>1000)
        exit(1);
        */
    /*  static int count=0;
      if(Minos->Delta_Vertex < 5 && sqrt(Minos->X_Vertex*Minos->X_Vertex+Minos->Y_Vertex*Minos->Y_Vertex)<15 && Minos->Z_Vertex>-4650){
        file << BDC1.X()  <<" " << BDC1.Y()<< " " << BDC1.Z() <<" " ;
        file << BDC2.X()  <<" " << BDC2.Y()<< " " << BDC2.Z() <<" " ;
        file << Minos->X_Vertex      <<" " << Minos->Y_Vertex    << " " << Minos->Z_Vertex    << endl ;
        count ++;
      }
      if(count>10000)
        exit(1);
      */  

    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Clear(){
  Brho=-1000;
  BDCX=-1000;
  BDCY=-1000;
  Beta_f=-1000;
  Beta_n=-1000;
  X=Y=Z=-1000;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Brho",&Brho,"Brho/D");
  RootOutput::getInstance()->GetTree()->Branch("BDCX",&BDCX,"BDCX/D");
  RootOutput::getInstance()->GetTree()->Branch("BDCY",&BDCY,"BDCY/D");
  RootOutput::getInstance()->GetTree()->Branch("X",&X,"X/D");
  RootOutput::getInstance()->GetTree()->Branch("Y",&Y,"Y/D");
  RootOutput::getInstance()->GetTree()->Branch("Z",&Z,"Z/D");
  RootOutput::getInstance()->GetTree()->Branch("Erel",&Erel,"Erel/D");
  RootOutput::getInstance()->GetTree()->Branch("Beta_f",&Beta_f,"Beta_f/D");
  RootOutput::getInstance()->GetTree()->Branch("Beta_n",&Beta_n,"Beta_n/D");
  RootOutput::getInstance()->GetTree()->Branch("Trigger",&Trigger,"Trigger/I");
} 
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchAddress("Trigger",&Trigger);
}
  ////////////////////////////////////////////////////////////////////////////////
  //            Construct Method to be pass to the AnalysisFactory              //
  ////////////////////////////////////////////////////////////////////////////////
  NPL::VAnalysis* Analysis::Construct(){
    return (NPL::VAnalysis*) new Analysis();
  }

  ////////////////////////////////////////////////////////////////////////////////
  //            Registering the construct method to the factory                 //
  ////////////////////////////////////////////////////////////////////////////////
  extern "C"{
    class proxy_analysis{
      public:
        proxy_analysis(){
          NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
        }
    };

  proxy_analysis p_analysis;
}
