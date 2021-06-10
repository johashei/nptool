#include "MinimizeBeamSpot.h"
////////////////////////////////////////////////////////////////////////////////
double devE(const double * parameter) {
  //Beam energy: 7.7 [MeV/A] * 47 [A] = 361.9 [MeV]
  static NPL::Reaction reaction("47K(d,p)48K@362");

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

  //Initilize results array
  //    7 => Sum in 0 and them MG's in 1-6
  //    5 => Mean, MeanErr, StdDev, StdDevErr, Chi2/NDF
  double FitResultMatrix[7][5];

  //Loop over events
  for (unsigned int i = 0; i < size; i++) {
    //Particle path vector
    dir = * (pos[i]) - offset;

    //Define normal vector for the MG# of detection
    DetectorSwitch(detnum[i], MugastNormal);

    //Angle leaving target, angle entering MUGAST & energy deposited in MUGAST
    double ThetaTarget = dir.Angle(TVector3(0, 0, 1));
    double ThetaMugast = dir.Angle(MugastNormal);
    double Energy = energy[i];

    //Energy loss in Al
    Energy = Al.EvaluateInitialEnergy(
      Energy,           //energy after Al 
      0.4 * micrometer, //thickness of Al
      ThetaMugast       //angle impinging on MUGAST
    );

    //Energy loss in target
    Energy = CD2.EvaluateInitialEnergy(
      Energy,                          //energy after leaving target
      0.5 * parameter[3] * micrometer, //pass through half target
      ThetaTarget                      //angle leaving target
    );

    //Final value of Ex
    double Ex = reaction.ReconstructRelativistic(Energy, ThetaTarget);

    //Fill histograms with Ex
    h -> Fill(Ex);
    DetectorSwitch(detnum[i], Ex);
  }

  //Initilise, Draw & Fit histograms
  InitiliseCanvas(FitResultMatrix);

  //Write vals to screen
  if (flagDraw) {
    cout << "==================================================" << endl;
  }

  //Adapt the metric as needed
  double multiplier = 0.07; //0.1;
  double metric = sqrt(pow(FitResultMatrix[mgSelect][0] - refE, 2) 
		     + pow(multiplier * FitResultMatrix[mgSelect][2], 2));

  WriteToCout(FitResultMatrix[mgSelect], parameter[3], metric);

  //  double multiplier = 0.2;
  //  double metric = abs(FitResultMatrix[1][0]-refE) + multiplier*FitResultMatrix[1][2]
  //	        + abs(FitResultMatrix[2][0]-refE) + multiplier*FitResultMatrix[2][2]
  //	        + abs(FitResultMatrix[3][0]-refE) + multiplier*FitResultMatrix[3][2]
  //	        + abs(FitResultMatrix[4][0]-refE) + multiplier*FitResultMatrix[4][2]
  //	        + abs(FitResultMatrix[5][0]-refE) + multiplier*FitResultMatrix[5][2]
  //	        + abs(FitResultMatrix[6][0]-refE) + multiplier*FitResultMatrix[6][2];

  if (flagDraw) { WriteToFile(FitResultMatrix[mgSelect], parameter, metric); }

  return metric;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void MinimizeBeamSpot() {

  //Read data in
  LoadFile();

  //Output formatting
  cout << fixed << showpoint << setprecision(6) << showpos;

  //Read in
  cout << "==================================================" << endl;
  cout << "=--------- SELECT TELESCOPE TO MINIMIZE ---------=" << endl;
  cout << "= Type MG# of telescope metric to use, or type 0 =" << endl;
  cout << "= to use the sum of all MG's                     =" << endl;
  cout << "==================================================" << endl;
  cin >> mgSelect;
  if (mgSelect == 7) {
    mgSelect = 6;
  } // Correct the input for MG7
  cout << "==================================================" << endl;

  if (mgSelect >= 7) {
   cout << "ERROR!! INVALID SELECTION" << endl;
    return;
  }

  //Open output file
  file.open("gridMinResults.txt", ios::out);
  file << "minimX  \tminimY  \tminimZ  \tminimT  \tMetric  \tMean    \tMeanErr \tStdDev  \tStdDvErr\tChi2/NDF" << endl;
  file << setprecision(6) << fixed;

  //Start timer
  auto start = high_resolution_clock::now();

  //TESTING: Grid method of avoiding local minima
  for (int x = 0; x < 2; x++) {
    for (int y = 0; y < 2; y++) {
      for (int z = 0; z < 2; z++) {
        for (int t = 0; t < 5; t++) {

          //Start with beam (0,0,0) and 4.76um 0.5mg/c2 target
          //double parameter[4] = {0.0, 0.0, 0.0, 4.76};   
          double parameter[4] = {
            9.0 - (9.0 * x),
            9.0 - (9.0 * y),
            9.0 - (9.0 * z),
            (1.5 - (0.25 * t)) * 4.75
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
          minim -> SetMaxFunctionCalls(1000000); // used by Minuit and Minuit2 
          minim -> SetMaxIterations(1000000); // used by GSL
          minim -> SetPrintLevel(3);
          minim -> SetPrecision(1e-10);

          //Set minimizer function
          minim -> SetFunction(func);

          //Assign variable limits
          minim -> SetLimitedVariable(0, "X", parameter[0], 0.01, -10, 10);
          minim -> SetLimitedVariable(1, "Y", parameter[1], 0.01, -10, 10);
          minim -> SetLimitedVariable(2, "Z", parameter[2], 0.01, -5, 5);
          minim -> SetLimitedVariable(3, "T", parameter[3], 0.01, 4.76 * 0.5, 4.76 * 1.5);

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
        }
      }
    }
  }

  //Stop timer
  auto stop = high_resolution_clock::now();

  //Close output file
  file.close();

  //Program runtime
  auto duration = duration_cast < seconds > (stop - start);
  cout << " Runtime = " << duration.count() << " s" << endl;
}
////////////////////////////////////////////////////////////////////////////////
