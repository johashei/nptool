#include "Garfield/GeometrySimple.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewMedium.hh"
		 
void GARFDRIFT(double energy, double t_start, G4ThreeVector start, G4ThreeVector delta_pos, double R, int Pad_start, double ID,
	  double ScoLength, double SegLength) {


  const double rWire = 1.2;
  const double rTube = 5.;
  const double lTube = 25.19;
  //For R6 = 2, R7 = 1
  //const double vWire = -375.; //For V_total = 750V
  //const double vWire = -645.3312; //For V_total = 1200V //SHOULD BE 583.741 !?
  //const double vWire = -535.0959; //For V_total = 1100V
  //const double vWire = -875.6116; //For V_total = 1800V 
  const double vWire =  -681.0312; //For V_total = 1400V
  const double vTube = 0.;
  double x_start, y_start, z_start;
  ofstream file;
  double production_bias = 1.e-01;
  double t_end, x_end, y_end, z_end;
  double Pad_end;
  //double w_value = 35.; //for HeCO2 --guess
  double w_value = 26.31; //for argon
  int electrons  = (int)(energy/w_value*production_bias);//*production_bias; //Need to carry excess energy over later.
  
  Garfield::MediumMagboltz* gas = new Garfield::MediumMagboltz();
  Garfield::ViewMedium* mediumView = new Garfield::ViewMedium();
  Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();
  Garfield::SolidTube* tube = new Garfield::SolidTube(0, 0, 0, 0, rTube, lTube/2);  
  Garfield::ComponentAnalyticField* cmp = new Garfield::ComponentAnalyticField();
  Garfield::Sensor* sensor = new Garfield::Sensor();
  Garfield::AvalancheMC* drift = new Garfield::AvalancheMC();
  
  cout << "\n" << electrons << "\n" << endl;
  
  if(R>rWire && electrons>0) {

    gas->LoadGasFile("ar_90_ch4_10.gas"); //atomic fraction in this case (P10).
    //gas->LoadGasFile("he_90_co2_10_500mbar.gas"); //mass fraction in this case
   
    gas->Initialise(true);
    
    mediumView->SetMedium(gas);
    geo->AddSolid(tube, gas);
    cmp->SetGeometry(geo);
    cmp->AddWire(0, 0, 2*rWire, vWire, "c");
    cmp->AddTube(rTube, vTube, 0, "a");
    sensor->AddComponent(cmp);
    drift->SetSensor(sensor);
      
    for(int e=0;e<electrons;e++) {
      
      double randomize = (double)std::rand() / (double)RAND_MAX;

      x_start = start.x() + delta_pos.x()*randomize;
      y_start = start.y() + delta_pos.y()*randomize;
      z_start = start.z() + delta_pos.z()*randomize;

      drift->DriftElectron(x_start, y_start, z_start, t_start);
      int np = drift->GetNumberOfDriftLinePoints();
      drift->GetDriftLinePoint(np-1, x_end, y_end, z_end, t_end); //np-1 is last DriftLineEndPoint
      
      Pad_end = (int)((z_end + ScoLength / 2.) / SegLength ) + 1; //new Pad number 
      
      file.open("signal.dat", std::ios::app);
      file <<  Pad_end << "\t" << t_end << endl;
      file.close();
      
    }
  }

  delete gas; delete mediumView; delete geo; delete tube; delete cmp; delete sensor; delete drift; //otherwise these are held in memory and cause a killed: 9 crash.
  
}
