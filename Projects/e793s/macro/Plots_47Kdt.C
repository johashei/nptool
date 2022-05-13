#include "GausFit.h"
//#include "KnownPeakFitter.h"
#include "DrawPlots.h"

//#include "CS2.h"
//#include "ThreeBodyBreakup.h"
//#include "ThreeBodyBreakup_FitPhaseSpace.h"

/* MAIN FUNCTION */

void Plots_47Kdt(){

  LoadChain47Kdt();
  gStyle->SetOptStat("nemMrRi");

  tCentre = 2750;  tRange = 350;
  timegate = "abs(T_MUGAST_VAMOS-" + to_string(tCentre) + ")<" + to_string(tRange);
  det_gate = "MUST2.TelescopeNumber>0 && MUST2.TelescopeNumber<5";
  
  cout << "==============================================" << endl;
  cout << "=============== CROSS SECTIONS ===============" << endl;
  cout << "==============================================" << endl;

}
