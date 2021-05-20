#ifndef TBigRIPSReco_H
#define TBigRIPSReco_H
#include "TObject.h"
#include <vector>
#include <iostream>
class TBigRIPSReco: public TObject{

  public:
    TBigRIPSReco();
    ~TBigRIPSReco();

  private:
    double aoq;
    double beta;
    double delta;
    double brho;
    double angle;
    double z;
  public:
    void Clear();
    void Init();
    void Print();
    void RecBrho(std::vector<double>,std::vector<double>,std::vector<std::vector<double>>,double);
    void RecAoqOne(double, double);
    void RecZet(double,double, std::vector<double>);
    const double c_mm_ns = 299.7792458; // in mm/ns
    const double mnucleon = 931.49432;  // MeV
  
  public:
    void SetAoq(double value){aoq=value;}
    void SetBeta(double value){beta=value;}
    void SetDelta(double value){delta=value;}
    void SetBrho(double value){brho=value;}
    void SetZ(double value){z=value;}

    ClassDef(TBigRIPSReco,1); 
};

#endif
