#include "MinimizeBeamSpot.h"
////////////////////////////////////////////////////////////////////////////////
double devE(const double * parameter) {
  //Beam energy: 7.7 [MeV/A] * 47 [A] = 361.9 [MeV]
  static NPL::Reaction reaction("47K(d,p)48K@355");

  //Beam spot offset
  TVector3 offset(parameter[0], parameter[1], parameter[2]);

  //Other variable initilizations
  unsigned int size = pos.size();
  TVector3 MugastNormal;
  double dE, Theta;
  TVector3 dir;

  //Initilize histogram
  h -> Reset();
  h1 -> Reset();
  h2 -> Reset();
  h3 -> Reset();
  h4 -> Reset();
  h5 -> Reset();
  h7 -> Reset();
  hT -> Reset();
  hEL -> Reset();

  //Now that initial range is wide, crop to single peak
  if(refE==0.143){
    h->SetAxisRange( -1.0, +2.6, "X");
    h1->SetAxisRange(-1.0, +2.6, "X");
    h2->SetAxisRange(-1.0, +2.6, "X");
    h3->SetAxisRange(-1.0, +2.6, "X");
    h4->SetAxisRange(-1.0, +2.6, "X");
    h5->SetAxisRange(-1.0, +2.6, "X");
    h7->SetAxisRange(-1.0, +2.6, "X");
  }

  //Initilize results array
  //    7 => Sum in 0 and them MG's in 1-6
  //    7 => Mean, MeanErr, StdDev, StdDevErr, Chi2/NDF, 
  //    Mean2, StdDev2, Mean3, StdDev3
  double FitResultMatrix[7][9];

  //Loop over events
  for (unsigned int i = 0; i < size; i++) {
    //Particle path vector
    dir = * (pos[i]) - offset;

    //Set final beam energy (added 05July)
    double FinalBeamEnergy = BeamTarget.Slow(
//		    reaction.GetBeamEnergy(), 
                    354.75,                        //post-CATS beam energy
		    0.5*parameter[3]*micrometer,   //thickness
		    0);                            //angle, probably
    reaction.SetBeamEnergy(FinalBeamEnergy);

    //Define normal vector for the MG# of detection
    DetectorSwitch(detnum[i], MugastNormal);

    //Angle leaving target, angle entering MUGAST & energy deposited in MUGAST
    double ThetaTarget = dir.Angle(TVector3(0.0, 0.0, 1.0));
    double ThetaMugast = dir.Angle(MugastNormal);
    double Energy = energy[i];
//cout << "@ angle " << ThetaTarget/deg << "   " << Energy << " -> ";
    //Energy loss in Al
    Energy = Al.EvaluateInitialEnergy(
      Energy,           //energy  Al 
      0.4 * micrometer, //0.4 * micrometer, //thickness of Al
      ThetaMugast       //angle impinging on MUGAST
    );

//cout << Energy << " -> ";
    //Energy loss in target
    Energy = CD2.EvaluateInitialEnergy(
      Energy,                          //energy after leaving target
      0.5 * parameter[3] * micrometer, //pass through half target
      ThetaTarget                      //angle leaving target
    );
    
//cout << Energy << endl;
    //Final value of Ex
    double Ex = reaction.ReconstructRelativistic(Energy, ThetaTarget);

    //Fill histograms with
//if(ThetaTarget/deg<130.){ // TESTING
//if(ThetaTarget/deg>130.){ // TESTING
    if(allButMG3){
      if(detnum[i]!=3){
        h -> Fill(Ex);
        cout << "removed MG3 from sum!!" << endl;
      }
    } else {
      h -> Fill(Ex);
    }
//}
    DetectorSwitch(detnum[i], Ex);
    hT -> Fill(ThetaTarget/deg,Ex);
    hEL -> Fill(ThetaTarget/deg,Energy-energy[i]);//ELoss in target and Al
  }

  //Initilise, Draw & Fit histograms
  InitiliseCanvas(FitResultMatrix);

  //Write vals to screen
  if (flagDraw) {cout << "==================================================" << endl;}

  /*** Minimize by one peak ***/
/**/
  double multiplier = 0.05; //0.08;
  double metric = abs(FitResultMatrix[mgSelect][0]-refE) + abs(multiplier*FitResultMatrix[mgSelect][2]);
  //double metric = abs(FitResultMatrix[mgSelect][0]-0.143) + abs(multiplier*FitResultMatrix[mgSelect][2]) + abs(FitResultMatrix[mgSelect][5]-1.410) + abs(multiplier*FitResultMatrix[mgSelect][6]);
  //double metric = abs(FitResultMatrix[mgSelect][0]-0.143) + abs(multiplier*FitResultMatrix[mgSelect][2]) + abs(FitResultMatrix[mgSelect][5]-1.410) + abs(multiplier*FitResultMatrix[mgSelect][6])+ abs(FitResultMatrix[mgSelect][7]-1.980) + abs(multiplier*FitResultMatrix[mgSelect][8]) ;
  //double metric = abs(FitResultMatrix[mgSelect][7]-1.981) + abs(multiplier*FitResultMatrix[mgSelect][8]); ;

/**/

  /*** Minimize by all peaks ***/
/**
  //double multiplier = 0.125;
  double multiplier = 0.05;
  double metric = 
	            (1.0/6.0)*abs(FitResultMatrix[1][0]-refE) 
  	          + (1.0/6.0)*abs(FitResultMatrix[2][0]-refE) 
  	          + (1.0/6.0)*abs(FitResultMatrix[3][0]-refE)
 	          + (1.0/6.0)*abs(FitResultMatrix[4][0]-refE)
 	          + (1.0/6.0)*abs(FitResultMatrix[5][0]-refE)
  	          + (1.0/6.0)*abs(FitResultMatrix[6][0]-refE)
  	          //+ 1.0*abs(FitResultMatrix[0][0]-refE) 
		  + multiplier*FitResultMatrix[1][2]
		  + multiplier*FitResultMatrix[2][2]
		  + multiplier*FitResultMatrix[3][2]
		  + multiplier*FitResultMatrix[4][2]
		  + multiplier*FitResultMatrix[5][2]
		  + multiplier*FitResultMatrix[6][2]
                  ;  
**/

  /*** Minimize by variation in peaks ***/
/**
  vector<double> r = {FitResultMatrix[1][0],
                      FitResultMatrix[2][0], 
  	              FitResultMatrix[3][0],
 	              FitResultMatrix[4][0],
 	              FitResultMatrix[5][0],
  	              FitResultMatrix[6][0]};
  double metric = (*max_element(r.begin(),r.end()) - *min_element(r.begin(),r.end()));
**/

  WriteToCout(FitResultMatrix[mgSelect], parameter[3], metric);
  if (flagDraw) { WriteToFile(FitResultMatrix[mgSelect], parameter, metric); }
  
  return metric;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MinimizeBeamSpot() {

  //filename = "XYZE_GateOn1838_11Apr22.txt"; refE = 1.981; // the energy of the selected states
  filename = "XYZE_GateOn0143_11Apr22.txt"; refE = 0.143; // the energy of the selected states


  //Read data in
  LoadFile();

  //Output formatting
  cout << fixed << showpoint << setprecision(6) << showpos;

  //Read in
  cout << "==================================================" << endl;
  cout << "=--------- SELECT TELESCOPE TO MINIMIZE ---------=" << endl;
  cout << "= Type MG# of telescope metric to use, or:       =" << endl;
  cout << "=   0 = sum of all MG's                          =" << endl;
  cout << "=   8 = sum of all MG's without MG3              =" << endl;
  cout << "==================================================" << endl;
  cin >> mgSelect;
  if (mgSelect == 7) {
    mgSelect = 6;
  } // Correct the input for MG7
  if (mgSelect == 8) {
    allButMG3=true;
    mgSelect = 0;
  }

  cout << "==================================================" << endl;

  if (mgSelect >= 7 ) {
   cout << "ERROR!! INVALID SELECTION" << endl;
    return;
  }

  //Open output file
  file.open("gridMinResults.txt", ios::out);
  file << "minimX  \tminimY  \tminimZ  \tminimT  \tMetric  \tMean    \tMeanErr \tStdDev  \tStdDvErr\tChi2/NDF" << endl;
  file << setprecision(6) << fixed;

  //Start timer
  auto start = high_resolution_clock::now();


  //Start with beam (0,0,0) and 4.76um 0.5mg/c2 target
  double parameter[4] = {
    //0.0, 0.0, 0.0, 4.76 
    //-3.9164, +0.0550, 1.3958, 1.3008

    //-3.9164, +0.0550, +1.0, +2.4
    //-3.9164, +0.0550, +0.681611, +2.256769
    //-3.9164, +0.0550, +0.786721, +2.113878

    //-3.9164, 0.0550, 0.5, 2.6
    //-3.9164, 0.0550, +0.146952, +2.517784

    //-4.675715, +0.143686, +0.115426, +2.586484
    //-4.636827, -0.206620, +0.115426, +2.586484 


    //-4., +1.0, +0.115426, +2.586484 
    //-5.088554, +1.476001, +0.342552, +2.586484 
    //-5.088554, +1.476001, +0.236014, +2.523850 
    //-3.989139, 0.99116, 0.115426, 2.586484 
    //-3.989139, +0.991160, +0.587660, +2.586484

    //-3.989139, +0.991160, +0.587660, +2.586484
    //-3.9164, +0.0550, 0.0, 2.6
    //-3.9164, +0.0550, 0.0, +2.838327
    //-3.912110, +0.083580, +0.000000, +2.838327
    //-3.912110, +0.083580, +0.000000, +2.6
    //-3.9164, +0.0550, +0.0, 2.3
    //-3.9164, +0.0550, +0.0, +2.591512
    -5.334864, -0.698900, +0.0, +2.591512

  };

  //Don't draw iterations of minimizer
  flagDraw = 0;
  gROOT -> SetBatch(kTRUE);

  //Initial pass through
  devE(parameter);

  //Function with 4 parameter XYZ and Target thickness
  auto func = ROOT::Math::Functor( & devE, 4);

  //Initilise minimizer
  auto minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minim -> SetMaxFunctionCalls(100000000); // used by Minuit and Minuit2 
  minim -> SetMaxIterations(100000000); // used by GSL
  minim -> SetPrintLevel(3);
  minim -> SetPrecision(1e-06);

  //Set minimizer function
  minim -> SetFunction(func);

  //Assign variable limits
  //minim -> SetLimitedVariable(0, "X", parameter[0], 0.10, -8.0, -0.0);
  minim -> SetFixedVariable(0, "X", parameter[0]); 
  //minim -> SetLimitedVariable(1, "Y", parameter[1], 0.10, -6.0, +6.0);
  minim -> SetFixedVariable(1, "Y", parameter[1]);
  //minim -> SetLimitedVariable(2, "Z", parameter[2], 0.05, 0.90, +1.10);
  //minim -> SetLimitedVariable(2, "Z", parameter[2], 0.10, -2.00, +2.00);
  minim -> SetFixedVariable(2, "Z", parameter[2]);
  //minim -> SetLimitedVariable(3, "T", parameter[3], 0.05, +0.9, +1.5); // ELASTICS, 1.2(3)
  //minim -> SetLimitedVariable(3, "T", parameter[3], 0.05, +2.05, +3.25); // ELASTICS, 2.65(59)
  minim -> SetFixedVariable(3, "T", parameter[3]);


  //Don't draw iterations of minimizer
  flagDraw = 0;
  gROOT -> SetBatch(kTRUE);

  //Shrink it, babeyyy
  minim -> Minimize();

  //Draw minimal value
  flagDraw = 1;
  gROOT -> SetBatch(kFALSE);

  //Pull values from minimizer
  const double * x = minim -> X();
  cout << "==================================================" << endl;
  cout << "=---------------- FINAL PEAK FITS ---------------=" << endl;
  cout << "==================================================" << endl;
    devE(x);
  cout << "==================================================" << endl;
  cout << "=------------ RESULTS OF MINIMIZATION -----------=" << endl;
    if (mgSelect == 6) { mgSelect = 7; } // Correct the input for MG7
  cout << "=------------------- USING MG " << mgSelect << " -----------------=" << endl;
  cout << "==================================================" << endl;
  cout << "\t\tX = " << x[0] << " mm" << endl;
  cout << "\t\tY = " << x[1] << " mm" << endl;
  cout << "\t\tZ = " << x[2] << " mm" << endl;
  cout << "\t\tT = " << x[3] << " um" << endl;
  cout << "==================================================" << endl;

  //Stop timer
  auto stop = high_resolution_clock::now();

  //Close output file
  file.close();

  //Program runtime
  auto duration = duration_cast < seconds > (stop - start);
  cout << " Runtime = " << duration.count() << " s" << endl;
}
////////////////////////////////////////////////////////////////////////////////
