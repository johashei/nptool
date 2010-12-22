#include "ObjectManager.hh"

using namespace std;

int main(int argc,char** argv)
{	
	NPOptionManager* myOptionManager = NPOptionManager::getInstance(argc,argv)  ;
	string detectorfileName 		= myOptionManager->GetDetectorFilePath()	      ;
	string reactionfileName 	  = myOptionManager->GetCalibrationFilePath()	    ;
	string calibrationfileName 	= myOptionManager->GetCalibrationFilePath()	    ;
	string runToTreatFileName 	= myOptionManager->GetRunToReadFilePath()       ;
  
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//	First of All instantiate RootInput and Output
	//	Detector will be attached later
	RootInput:: getInstance(runToTreatFileName)	;
	RootOutput::getInstance("Analysis/60Fe_AnalyzedData", "AnalyzedTree")	;
	
	//	Instantiate some Reaction
	NPL::Reaction*  e530Reaction = new Reaction								;
	e530Reaction	->	ReadConfigurationFile(reactionfileName)	;

	//	Instantiate the Calibration Manger using a file (WARNING:prior to the detector instantiation)
	CalibrationManager* myCalibration = CalibrationManager::getInstance(calibrationfileName) ;

	//	Instantiate the detector using a file 
	NPA::DetectorManager* myDetector = new DetectorManager 		;
	myDetector	->	ReadConfigurationFile(detectorfileName)		;
	
	//	Ask the detector manager to load the parameter added by the detector in the calibrationfileName
	myCalibration->LoadParameterFromFile() ;
	/////////////////////////////////////////////////////////////////////////////////////////////////////

  //	Attach more branch to the output
  double ELab=0,ThetaLab=0,ExcitationEnergy=0;
	RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab,"ELab/D") 							;
	RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D") 	;
	RootOutput::getInstance()->GetTree()->Branch("ExcitationEnergy",&ExcitationEnergy,"ExcitationEnergy/D") ;
  //	Get the formed Chained Tree and Treat it
  TChain* Chain = RootInput:: getInstance() -> GetChain()	;

  TMust2Physics* M2 		= (TMust2Physics*) 	myDetector -> m_Detector["MUST2"] 	;
  cout <<  " ///////// Starting Analysis ///////// "<< endl << endl ;
	
  int i ,N=Chain -> GetEntries();
	
  //N = 1000;
  cout << " Number of Event to be treated : " << N << endl ;
	
  clock_t begin=clock();
  clock_t end=begin;
  for ( i = 0 ; i < N ; i ++ )
    {
      // Minimum code
    	if( i%10000 == 0 && i!=0) 	{	
				cout.precision(5);
				end=clock();										
				double TimeElapsed = (end-begin)/CLOCKS_PER_SEC;
				double percent = (double)i/N ;
				double TimeToWait = (TimeElapsed/percent) - TimeElapsed	;					
				cout	<< "\r Progression:" << percent*100 
	     				<< " % \t | \t Remaining time : ~" 
	     				<<  TimeToWait <<"s"<< flush;
      }	
										
      else if (i==N-1) 	cout << "\r Progression:" 
			     << " 100% " <<endl;
					
      Chain -> GetEntry(i);
      // Clear Previous Event
      myDetector -> ClearEventPhysics()			;
      // Build the new one
      myDetector -> BuildPhysicalEvent()		;
      ////
		
			
			// Must 2
			for(int hit = 0; hit < M2 -> Si_E.size() ; hit ++)
				{
					ELab = -1 ; ThetaLab = -1;
					//	Get Hit Direction
					TVector3 HitDirection  = M2 -> GetPositionOfInteraction(hit) - TVector3(0,0,-40);
					// Angle between beam and particle
					ThetaLab  = ThetaCalculation ( HitDirection , TVector3(0,0,1)   ) ;	
					ELab = M2 -> Si_E[hit] + M2 -> SiLi_E[hit]	;
			  }
			
      RootOutput::getInstance()->GetTree()->Fill()	;
    }

  cout << " A total of " << i << " event has been annalysed " << endl ;
  cout << endl << " ///////////////////////////////////// "<< endl<< endl ;
  RootOutput::getInstance()->Destroy();
  return 0	;
}

double ThetaCalculation (TVector3 A , TVector3 B)
{
  double Theta = acos( (A.Dot(B)) / (A.Mag()*B.Mag()) ) ;
  return Theta*rad ;
}
