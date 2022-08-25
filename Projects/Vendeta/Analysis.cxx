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
	if(FC_mult==1 ){
		int anode = FC->AnodeNumber[0];
		double Time_FC = FC->Time[0];

		Vendeta->SetAnodeNumber(anode);
		Vendeta->BuildPhysicalEvent();
		FC->BuildPhysicalEvent();

		// VENDETA LG 
		unsigned int Vendeta_LG_mult = Vendeta->LG_DetectorNumber.size();
		for(unsigned int i=0; i<Vendeta_LG_mult; i++){

			int DetNbr          = Vendeta->LG_DetectorNumber[i];
			double Time_Vendeta = Vendeta->LG_Time[i];
			double Rdet         = Vendeta->GetDistanceFromTarget(DetNbr);
			TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);

			double DT = Time_Vendeta - Time_FC;// + ToF_Shift_Vendlg[DetNbr-1];

			if(DT>0){

				double DeltaTheta = atan(63.5/Rdet);
				double Theta_Vendeta = DetPos.Theta();
				double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
				//cout << DT << " " << Rdet << endl;
				//neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
				neutron->SetTimeOfFlight(DT*1e-9/(0.55));
				//neutron->SetBeta( (LoF_Vend[DetNbr-1] / DT )/ c_light); 

				double En = neutron->GetEnergy();

				// Filling output tree
				LG_Tof.push_back(DT);
				LG_ID.push_back(DetNbr);
				LG_ELab.push_back(En);
				LG_ThetaLab.push_back(Theta_random);
				LG_Q1.push_back(Vendeta->LG_Q1[i]);
				LG_Q2.push_back(Vendeta->LG_Q2[i]);
				LG_Qmax.push_back(Vendeta->LG_Qmax[i]);
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

			if(DT>0){
				double DeltaTheta = atan(63.5/Rdet);
				double Theta_Vendeta = DetPos.Theta();
				double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
				//cout << DT << " " << Rdet << endl;
				//neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
				neutron->SetTimeOfFlight(DT*1e-9/(0.55));
				//neutron->SetBeta( (LoF_Vend[DetNbr-1] / DT )/ c_light); 
				double En = neutron->GetEnergy();

				// Filling output tree
				HG_ID.push_back(DetNbr);
				HG_Tof.push_back(DT);
				HG_ELab.push_back(En);
				HG_ThetaLab.push_back(Theta_random);
				HG_Q1.push_back(Vendeta->HG_Q1[i]);
				HG_Q2.push_back(Vendeta->HG_Q2[i]);
				HG_Qmax.push_back(Vendeta->HG_Qmax[i]);
			}
		}

		//Process coincidences signals in VENDETA LG / HG

		/*if(HG_Tof.size() > 0 && LG_Tof.size() > 0 ){
			for(int j = 0; j < LG_Tof.size();j++){
				for(int k = 0; k < HG_Tof.size(); k++){
					if(abs(HG_Tof[k]-LG_Tof[j]) < 2 && HG_ID[k] == LG_ID[j]){
						if( HG_Q2[k]>120000){
							//  HG_ID[k] = 
							HG_Tof[k] = - 100000;
							HG_ELab[k] = - 100000;
							HG_ThetaLab[k] = - 100000;
							HG_Q1[k] = - 100000;
							HG_Q2[k] = - 100000;
							HG_Qmax[k] = - 100000;
						}  
						else if( HG_Q2[k]<120000){
							// HG_ID[k] = 
							LG_Tof[k] = - 100000;
							LG_ELab[k] = - 100000;
							LG_ThetaLab[k] = - 100000;
							LG_Q1[k] = - 100000;
							LG_Q2[k] = - 100000;

						}
					}
				}
			}
		} // if LG && HG*/
	}// if FC = 1


	}

	////////////////////////////////////////////////////////////////////////////////
	void Analysis::InitOutputBranch(){
		RootOutput::getInstance()->GetTree()->Branch("LG_ID",&LG_ID);
		RootOutput::getInstance()->GetTree()->Branch("LG_ThetaLab",&LG_ThetaLab);
		RootOutput::getInstance()->GetTree()->Branch("LG_ELab",&LG_ELab);
		RootOutput::getInstance()->GetTree()->Branch("LG_Tof",&LG_Tof);
		RootOutput::getInstance()->GetTree()->Branch("LG_Q1",&LG_Q1);
		RootOutput::getInstance()->GetTree()->Branch("LG_Q2",&LG_Q2);
		RootOutput::getInstance()->GetTree()->Branch("LG_Qmax",&LG_Qmax);

		RootOutput::getInstance()->GetTree()->Branch("HG_ID",&HG_ID);
		RootOutput::getInstance()->GetTree()->Branch("HG_ThetaLab",&HG_ThetaLab);
		RootOutput::getInstance()->GetTree()->Branch("HG_ELab",&HG_ELab);
		RootOutput::getInstance()->GetTree()->Branch("HG_Tof",&HG_Tof);
		RootOutput::getInstance()->GetTree()->Branch("HG_Q1",&HG_Q1);
		RootOutput::getInstance()->GetTree()->Branch("HG_Q2",&HG_Q2);
		RootOutput::getInstance()->GetTree()->Branch("HG_Qmax",&HG_Qmax);

	}

	////////////////////////////////////////////////////////////////////////////////
	void Analysis::ReInitValue(){
		LG_ThetaLab.clear();
		LG_ELab.clear();
		LG_Tof.clear();
		LG_ID.clear();
		LG_Q1.clear();
		LG_Q2.clear();
		LG_Qmax.clear();

		HG_ThetaLab.clear();
		HG_ELab.clear();
		HG_Tof.clear();
		HG_ID.clear();
		HG_Q1.clear();
		HG_Q2.clear();
		HG_Qmax.clear();

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

