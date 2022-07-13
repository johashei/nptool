void CalcEn(double PosY, double PosZ, double TOF){
  double c_light = 299.795;
  double mass_neutron = 939.57;
  double distance = sqrt(pow(PosY,2)+pow(PosZ,2));
  double speed = distance/TOF;
  double gamma = 1/sqrt(1-pow(speed,2)/pow(c_light,2));
  double energy = (gamma-1)*mass_neutron;
  std::cout << "Energy is : " << energy << " MeV" << std::endl;
 }
