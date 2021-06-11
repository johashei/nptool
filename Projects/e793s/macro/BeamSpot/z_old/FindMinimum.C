/***************************************************************************
 *   Charlie Paxman (cp00474@surrey.ac.uk)                                 *
 ***************************************************************************/

//-------------------------------
//C++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
//Root
//#include <TVector3.h>
//NPTool
//#include "NPEnergyLoss.h"
//#include "NPReaction.h"
//#include "NPSystemOfUnits.h"
using namespace std;
//-------------------------------
double Minimize(vector<double>, vector<vector<double>>, int, int);
//-------------------------------

void FindMinimum(){
  ifstream metricsFile;
  vector<double> T, X, Y, Z;
  vector<vector<double>> mean, sigma, metric;
  /* loop temps */
  double tT, tX, tY, tZ;
  double tmn, tsg, tmt;
  vector<double> tmnV, tsgV, tmtV;
  int lines = 0;

  const char* metricsFilePath = "output_Full_02June_coarse_metrics.txt";

  metricsFile.open(metricsFilePath, ios::in);
  if(!metricsFile){
    cout << "ERROR: File not opened." << endl;
    return 1;
  }

  // Read in the output from BeamSpot.C
  while(true){
    // Reset temps
    tT = tX = tY = tZ = tmn = tsg = tmt = 0.0;
    tmnV.clear(); tsgV.clear(); tmtV.clear();

    // Read from file
    metricsFile >> tT >> tX >> tY >> tZ;
    for(int i=0; i<7; i++){
      tmn = tsg = tmt = 0.0;
      metricsFile >> tmn >> tsg >> tmt;

      tmnV.push_back(tmn);
      tsgV.push_back(tsg);
      tmtV.push_back(tmt);
    }

    // Push temps to vectors
    T.push_back(tT);
    X.push_back(tX);
    Y.push_back(tY);
    Z.push_back(tZ);
    mean.push_back(tmnV);
    sigma.push_back(tsgV);
    metric.push_back(tmtV);

    // Count lines
    lines++;
    // Break out
    if(metricsFile.eof()){break;}
  }  

  cout << "==================================================" << endl;
  cout << "= Reading " << metricsFilePath << endl;
  cout << "==================================================" << endl;
  cout << "=--------- SELECT TELESCOPE TO MINIMIZE ---------=" << endl;
  cout << "= Type MG# of telescope metric to use, or type 0 =" << endl;
  cout << "= to use the sum of all MG's                     =" << endl;
  cout << "==================================================" << endl;

  unsigned int selectVar;
  cin >> selectVar;

  if(selectVar==7){selectVar=6;} // Correct the input for MG7

  if(selectVar>=7){
    cout << " FAIL! Invalid selection" << endl;
    cout << " Exiting..." << endl;
    return 1;
  }

  cout << fixed << showpoint;
  cout << ">   Thickness: " << setprecision(6) 
	                    << Minimize(T, metric, lines, selectVar) 
			    << " mm" 
			    << endl;
  cout << ">   X:         " << setprecision(6) 
	                    << Minimize(X, metric, lines, selectVar) 
			    << " mm" 
			    << endl;
  cout << ">   Y:         " << setprecision(6) 
	                    << Minimize(Y, metric, lines, selectVar) 
			    << " mm" 
			    << endl;
  cout << ">   Z:         " << setprecision(6) 
	                    << Minimize(Z, metric, lines, selectVar) 
			    << " mm" 
			    << endl;


  //TO DO:
  // - Automatically open the histogram that corresponds to the values 
  //   above, to save time and check for unreasonable results



}


/* * * * * * * * * * * * */
/* * * * * * * * * * * * */
/* * * * * * * * * * * * */

double Minimize(vector<double> V, vector<vector<double>> metric, int lines, int metSelect){
  vector<double> minVector;
  double minVal;
  unsigned int minIndex;

  //this is pretty inefficient, but there are no speed issues so oh well
  for(int i=0; i<lines-1; i++){
    minVector.push_back(metric[i][metSelect]);
  }

  minVal = *min_element(minVector.begin(), minVector.end());
  minIndex = min_element(minVector.begin(), minVector.end()) - minVector.begin();

  return V[minIndex];
}
/* * * * * * * * * * * * */

