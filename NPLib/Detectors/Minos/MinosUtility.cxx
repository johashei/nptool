#include "MinosUtility.h"
#include "Math/Factory.h"
#include <algorithm>
#include <fstream>
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
MinosUtility::MinosUtility(){
  // Setting up for signal fitting 
  m_signal_offset=0;
  m_signal_min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  m_signal_func=ROOT::Math::Functor(this,&NPL::MinosUtility::signal_chi2,4); 
  m_signal_min->SetFunction(m_signal_func);

  m_signal_parameter[0]=1000;
  m_signal_parameter[1]=250;
  m_signal_parameter[2]=15;
  m_signal_parameter[3]=0;
  m_Dp1=8; // mean distance between t0 and T[0]
  m_p2=15;// mean tau
  m_Ap0=20;// mean ratio between q max and A

  m_signal_min->SetLimitedVariable(0,"a",m_signal_parameter[0],1,0,100000);
  m_signal_min->SetLimitedVariable(1,"b",m_signal_parameter[1],1,0,512);
  m_signal_min->SetLimitedVariable(2,"c",m_signal_parameter[2],1,1,50);
  m_signal_min->SetLimitedVariable(3,"d",m_signal_parameter[3],1,-50,50);

  m_counter=0;

}

////////////////////////////////////////////////////////////////////////////////
double MinosUtility::signal_shape(const double& x, const double* p){
  static double a;
  a=(x-p[1])/p[2];
  if(x > p[1] && x < 512.) 
    return (p[0]*exp(-3.*a)*sin(a)*(a*a*a) + p[3]);
  else 
    return p[3];
}

////////////////////////////////////////////////////////////////////////////////
double MinosUtility::signal_chi2(const double* p){
  static double chi2, diff,expected;
  static unsigned int counter;
  counter=0;
  chi2=0;
  for(unsigned int i = 0 ; i < m_signal_size ; i++ ){
    expected = signal_shape((*m_fitSignalT)[i],p);
    diff  = (*m_fitSignalQ)[i] - expected ; 
    chi2 += diff*diff; 
    counter++;
    //if the track is too long, only the first part of the signal is taken
    if((*m_fitSignalT)[i]>t0+50)
      break;
  }
  return chi2/counter;
}

////////////////////////////////////////////////////////////////////////////////
double MinosUtility::Calibrate(const std::vector<double>& T, const std::vector<double>& Q, double& time, double& charge){
  // setting up for the minisation
  m_fitSignalT=&T;
  m_fitSignalQ=&Q;
  m_signal_size=T.size();   
  static double delta,qmax,Qmax,Dp1,Ap0;

  if(!m_signal_size){
    time=-10000;    
    charge=-10000;
    return -10000;
  }

  t0=T[0];
  qmax=*(std::max_element(Q.begin(),Q.end()));
  Qmax = qmax*m_Ap0;
  m_signal_min->SetLimitedVariable(0,"a",Qmax,1,Qmax*0.8,Qmax*1.2);
  m_signal_min->SetLimitedVariable(1,"b",t0+m_Dp1,1,t0-10,t0+10);
  m_signal_min->SetLimitedVariable(2,"c",m_p2,1,1,50);

  // Perform minimisation
  m_signal_min->Minimize(); 

  // access set of parameter that produce the minimum
  const double *xs = m_signal_min->X();
  charge= xs[0];
  time=xs[1];


  // Self optimisation of the starting parameter
  m_counter++;
  Dp1=time-t0;
  m_Dp1 +=(Dp1-m_Dp1)/(m_counter) ;
  Ap0=charge/qmax;
  m_Ap0 += (Ap0-m_Ap0)/(m_counter);
  m_p2 +=(xs[2]-m_p2)/(m_counter);


 /* 
     using namespace std;
     ofstream out("signal.txt");
     for(unsigned int i = 0 ; i < m_signal_size ; i++){
     out << T[i] << " " << Q[i] << endl;
     }
     out.close();
     out.open("fit.txt");
     out << xs[0] << " " << xs[1] << " " << xs[2] << " " << xs[3] << endl; 

     cout << xs[1]-T[0] << " " << xs[1] << " " << xs[2] << " " << xs[3] << " " << delta << " " << endl;
     static int first=0;

     if(first==40)
     exit(0);

     first++;
*/


  return m_signal_min->MinValue() ;

}
