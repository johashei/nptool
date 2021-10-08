#include "NPTrackingUtility.h"
#include "NPSystemOfUnits.h"

using namespace NPL;
using namespace NPUNITS;

void Z(){
  double res_z= 1;
  double res_theta= 1*deg;
  double angle=87.7*deg;
  TRandom3 rand;
  auto h = new TH2D("dd","dd",250,-50,250,250,70,105);
  for(unsigned int i = 0 ; i < 10000; i++){
    // Create to tracks with a fixe angle
    TVector3 T1(0,0,1);
    TVector3 T2(0,0,1);
    // Blurred angle 
    T1.SetTheta(rand.Gaus(angle*0.5,res_theta)); 
    T2.SetTheta(rand.Gaus(-angle*0.5,res_theta)); 
    // Build direction
    TVector3 P1A,P2A,P1B,P2B; 

    // Blurred Z
    double z = rand.Uniform(0,150);
    TVector3 dz1= TVector3(0,0,rand.Gaus(z,res_z));
    P1A=15*T1+dz1;  
    P1B=30*T1+dz1;  

    TVector3 dz2= TVector3(0,0,rand.Gaus(z,res_z));
    P2A=15*T2+dz2;  
    P2B=30*T2+dz2;  

    // Compute vertex
    TVector3 vertex,delta;
    MinimumDistanceTwoLines(P1A,P1B,P2A,P2B,vertex,delta);
  //  cout <<vertex.Z() << " " << ((P1B-P1A).Theta()+(P2B-P2A).Theta())/deg << endl;
    h->Fill(vertex.Z(),((P1B-P1A).Theta()+(P2B-P2A).Theta())/deg);
  }
  h->Draw("colz");
}
