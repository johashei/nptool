#include"NPPhysicalConstants.h"
double CS(double Z2, double Ei, double Ef, double Mt, double theta){
  static double a2 = NPUNITS::fine_structure_const*NPUNITS::fine_structure_const;  
  double sin4 = pow(sin(theta*0.5),4);
  double sin2 = sin(theta*0.5)*sin(theta*0.5);
  double cos2 = cos(theta*0.5)*cos(theta*0.5);
  return (M_PI*a2)/(Ei*Ei*sin4)*cos2/(1+(2*Ef/Mt)*sin2);
  }

void Mott(){
 NPL::Reaction r;
 r.ReadConfigurationFile("Sn132.reac");

 double Ei=r.GetBeamEnergy(); 
 double Ef,Thetaf;// electron
 double HE,ThetaE;// heavy ion
 unsigned int size = 180;
 double step = 180./size;
 double Z2 = r.GetParticle2()->GetZ()*r.GetParticle2()->GetZ();
 double Mt = r.GetParticle2()->Mass();
 vector<double> x,y;
 ofstream out("mott.txt");
 for(unsigned int i = 0 ; i < size ; i++){
   r.SetThetaCM(i*step);  
   r.KineRelativistic(Thetaf,Ef,ThetaE,HE);
   if(Thetaf){
     double val = CS(Z2,Ei,Ef,Mt,Thetaf); 
    y.push_back(val);
    x.push_back(Thetaf/NPUNITS::deg);
    out << i*step << " " << val << endl; 
     }
   } 
 
 auto g = new TGraph(x.size(),&x[0],&y[0]);
 g->Draw("ap");
 out.close();
  }
