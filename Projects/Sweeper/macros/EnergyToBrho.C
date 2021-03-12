#include "mass.h"

Double_t EnergyToBrho(Double_t Energy, Int_t A, Int_t Z){

  Double_t c = 3e8; //m/s
  Double_t u = 1.66e-27; //Kg
  Double_t e = 1.60217662e-19; //C 


  Double_t mass =  Nuke_Mass_Tab[A][Z];

  Double_t gamma = Energy/mass + 1;
  Double_t beta = TMath::Sqrt(1-1/(gamma*gamma));

  cout<<"mass: "<<mass<<" gamma: "<<gamma<<" beta: "<<beta<<endl;

  Double_t Brho = A*u*c*beta*gamma/(Z*e); 

  cout<< "Brho : "<< Brho << " Tm " << endl;
  
  return Brho;
}
