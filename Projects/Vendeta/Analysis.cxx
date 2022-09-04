/*****************************************************************************
	* Copyright (C) 2009-2016    this file is part of the NPTool Project        *
	*                                                                           *
	* For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
	* For the list of contributors see $NPTOOL/Licence/Contributors             *
	*****************************************************************************/

/*****************************************************************************
	* Original Author: XAUTHORX  contact address: XMAILX                        *
	*                                                                           *
	* Creation Date  : XMONTHX XYEARX                                           *
	* Last update    :                                                          *
	*---------------------------------------------------------------------------*
	* Decription:                                                               *
	*  This class describe  Vendeta analysis project                       *
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
#include"NPOptionManager.h"
#include"RootInput.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
		InitOutputBranch();

		Vendeta= (TVendetaPhysics*) m_DetectorManager->GetDetector("Vendeta");
		FC= (TFissionChamberPhysics*) m_DetectorManager->GetDetector("FissionChamber");

		neutron = new NPL::Particle("1n");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

		ReInitValue();

		unsigned int FC_mult = FC->AnodeNumber.size();
		unsigned int HF_mult = FC->Time_HF.size();

		double GammaOffset[11] = {971.37, 970.67, 972.73, 972.59, 989.33, 982.03, 999.59, 985.13, 980.03, 976.83, 958.46};

		double incomingDT=0;
		double incomingE=0;
		double flight_path = 21500.;
		for(unsigned int j=0; j<HF_mult; j++){
				for(unsigned int i=0; i<FC_mult; i++){
						incomingDT = FC->Time[i] - FC->Time_HF[j] - GammaOffset[FC->AnodeNumber[i]-1];
						if(incomingDT<0){
							incomingDT += 1790;
						}
						double length = flight_path;// + 6*FC->AnodeNumber[i];	
						neutron->SetBeta((length/incomingDT) / c_light);
						incomingE = neutron->GetEnergy();

						inToF.push_back(incomingDT);
						inEnergy.push_back(incomingE);
				}
		}

		if(FC_mult==1 && incomingDT<1790 && incomingDT>0){

				int anode = FC->AnodeNumber[0];
				double Time_FC = FC->Time[0];
				bool isFake = FC->isFakeFission[0];
				double FC_DT = FC->DT_FC[0];

				Vendeta->SetAnodeNumber(anode);
				Vendeta->BuildPhysicalEvent();
				//FC->BuildPhysicalEvent();

				// VENDETA LG 
				unsigned int Vendeta_LG_mult = Vendeta->LG_DetectorNumber.size();
				for(unsigned int i=0; i<Vendeta_LG_mult; i++){

						int DetNbr          = Vendeta->LG_DetectorNumber[i];
						double Time_Vendeta = Vendeta->LG_Time[i];
						double Rdet         = Vendeta->GetDistanceFromTarget(DetNbr);
						TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);

						double DT = Time_Vendeta - Time_FC;// + ToF_Shift_Vendlg[DetNbr-1];

						if(DT>-500){

								double DeltaTheta = atan(63.5/Rdet);
								double Theta_Vendeta = DetPos.Theta();
								double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
								//cout << DetNbr << " " << Rdet << endl;
								//neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
								//neutron->SetTimeOfFlight(DT*1e-9/(0.55));
								neutron->SetBeta(  (Rdet/DT) / c_light); 

								double En = neutron->GetEnergy();

								// Filling output tree
								LG_DT.push_back(FC_DT);
								LG_Tof.push_back(DT);
								LG_ID.push_back(DetNbr);
								LG_Anode_ID.push_back(anode);
								LG_ELab.push_back(En);
								LG_ThetaLab.push_back(Theta_random);
								LG_Q1.push_back(Vendeta->LG_Q1[i]);
								LG_Q2.push_back(Vendeta->LG_Q2[i]);
								LG_Qmax.push_back(Vendeta->LG_Qmax[i]);
								LG_FakeFission.push_back(isFake);
								FC_Q1.push_back(FC->Q1[0]);
								FC_Q2.push_back(FC->Q2[0]);
						}
				}

				// VENDETA HG 
				unsigned int Vendeta_HG_mult = Vendeta->HG_DetectorNumber.size();
				for(unsigned int i=0; i<Vendeta_HG_mult; i++){
						int DetNbr          = Vendeta->HG_DetectorNumber[i];
						double Time_Vendeta = Vendeta->HG_Time[i];
						double Rdet         = Vendeta->GetDistanceFromTarget(DetNbr);
						TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);

						double DT = Time_Vendeta - Time_FC;// + ToF_Shift_Vendhg[DetNbr-1];

						if(DT>-500){
								double DeltaTheta = atan(63.5/Rdet);
								double Theta_Vendeta = DetPos.Theta();
								double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
								//cout << DetNbr << " " << Rdet << endl;
								//neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
								//neutron->SetTimeOfFlight(DT*1e-9/(0.55));
								neutron->SetBeta( (Rdet/DT) / c_light); 
								double En = neutron->GetEnergy();

								// Filling output tree
								HG_DT.push_back(FC_DT),
								HG_ID.push_back(DetNbr);
								HG_Anode_ID.push_back(anode);
								HG_Tof.push_back(DT);
								HG_ELab.push_back(En);
								HG_ThetaLab.push_back(Theta_random);
								HG_Q1.push_back(Vendeta->HG_Q1[i]);
								HG_Q2.push_back(Vendeta->HG_Q2[i]);
								HG_Qmax.push_back(Vendeta->HG_Qmax[i]);
								HG_FakeFission.push_back(isFake);
								FC_Q1.push_back(FC->Q1[0]);
								FC_Q2.push_back(FC->Q2[0]);
						}
				}

				//Process coincidences signals in VENDETA LG / HG
				if(HG_Tof.size() > 0 && LG_Tof.size() > 0 ){
						for(int j = 0; j < LG_Tof.size();j++){
								for(int k = 0; k < HG_Tof.size(); k++){
										if(abs(HG_Tof[k]-LG_Tof[j]) < 2 && HG_ID[k] == LG_ID[j]){
												if( HG_Q2[k]>120000){
														HG_ID[k] = -1;
														HG_Tof[k] = - 100000;
														HG_ELab[k] = - 100000;
														HG_ThetaLab[k] = - 100000;
														HG_Q1[k] = - 100000;
														HG_Q2[k] = - 100000;
														HG_Qmax[k] = - 100000;
												}  
												else if( HG_Q2[k]<120000){
														LG_ID[j] = -1;
														LG_Tof[j] = - 100000;
														LG_ELab[j] = - 100000;
														LG_ThetaLab[j] = - 100000;
														LG_Q1[j] = - 100000;
														LG_Q2[j] = - 100000;
														LG_Qmax[j] = - 100000;

												}
										}
								}
						}
				} // if LG && HG*/

		}// if FC = 1
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
		// Incoming neutron
		RootOutput::getInstance()->GetTree()->Branch("inToF",&inToF);
		RootOutput::getInstance()->GetTree()->Branch("inEnergy",&inEnergy);

		RootOutput::getInstance()->GetTree()->Branch("LG_DT",&LG_DT);
		RootOutput::getInstance()->GetTree()->Branch("LG_ID",&LG_ID);
		RootOutput::getInstance()->GetTree()->Branch("LG_Anode_ID",&LG_Anode_ID);
		RootOutput::getInstance()->GetTree()->Branch("LG_ThetaLab",&LG_ThetaLab);
		RootOutput::getInstance()->GetTree()->Branch("LG_ELab",&LG_ELab);
		RootOutput::getInstance()->GetTree()->Branch("LG_Tof",&LG_Tof);
		RootOutput::getInstance()->GetTree()->Branch("LG_Q1",&LG_Q1);
		RootOutput::getInstance()->GetTree()->Branch("LG_Q2",&LG_Q2);
		RootOutput::getInstance()->GetTree()->Branch("LG_Qmax",&LG_Qmax);
		RootOutput::getInstance()->GetTree()->Branch("LG_FakeFission",&LG_FakeFission);

		RootOutput::getInstance()->GetTree()->Branch("HG_DT",&HG_DT);
		RootOutput::getInstance()->GetTree()->Branch("HG_ID",&HG_ID);
		RootOutput::getInstance()->GetTree()->Branch("HG_Anode_ID",&HG_Anode_ID);
		RootOutput::getInstance()->GetTree()->Branch("HG_ThetaLab",&HG_ThetaLab);
		RootOutput::getInstance()->GetTree()->Branch("HG_ELab",&HG_ELab);
		RootOutput::getInstance()->GetTree()->Branch("HG_Tof",&HG_Tof);
		RootOutput::getInstance()->GetTree()->Branch("HG_Q1",&HG_Q1);
		RootOutput::getInstance()->GetTree()->Branch("HG_Q2",&HG_Q2);
		RootOutput::getInstance()->GetTree()->Branch("HG_Qmax",&HG_Qmax);
		RootOutput::getInstance()->GetTree()->Branch("FC_Q1",&FC_Q1);
		RootOutput::getInstance()->GetTree()->Branch("FC_Q2",&FC_Q2);
		RootOutput::getInstance()->GetTree()->Branch("HG_FakeFission",&HG_FakeFission);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
		inToF.clear();
		inEnergy.clear();

		LG_ThetaLab.clear();
		LG_ELab.clear();
		LG_Tof.clear();
		LG_ID.clear();
		LG_Anode_ID.clear();
		LG_Q1.clear();
		LG_Q2.clear();
		LG_Qmax.clear();
		LG_FakeFission.clear();
		LG_DT.clear();

		HG_ThetaLab.clear();
		HG_ELab.clear();
		HG_Tof.clear();
		HG_ID.clear();
		HG_Anode_ID.clear();
		HG_Q1.clear();
		HG_Q2.clear();
		HG_Qmax.clear();
		HG_FakeFission.clear();
		HG_DT.clear();

		FC_Q1.clear();
		FC_Q2.clear();

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
		return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
		class proxy{
				public:
						proxy(){
								NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
						}
		};

		proxy p;
}

