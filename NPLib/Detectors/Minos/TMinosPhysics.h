#ifndef TMinosPHYSICS_H
#define TMinosPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Cyril Lenain  contact address: lenain@lpccaen.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : 2019                                                     *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// C++ headers 
#include <vector>
#include <map>
#include <string>
using namespace std;

// ROOT headers
#include "TObject.h"
#include "TH1.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMinuit.h"


// NPTool headers
#include "TMinosData.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "Tracking.h"
#include "TMinosResult.h"

#include "TClonesArray.h"

using namespace NPL;
// forward declaration

class TMinosPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TMinosPhysics();
    ~TMinosPhysics() {};
  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   

    void Clear(const Option_t*) {};
    
    /* static  void SumDistance(int &, double *, double & sum, double * par,  int); */
  
    //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    
    vector<double>  X_Pad;  
    vector<double>  Y_Pad;  
    vector<double>  Z_Pad;
    vector<double>  Q_Pad;
    vector<double>  T_Pad;

  /// A usefull method to bundle all operation to add a detector
    void AddDetector(TVector3 POS, string shape); 
    void AddDetector(double R, double Theta,TVector3 POS, double Phi, double TargetLenght); 
  
  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);
    // add parameters to the CalibrationManger
    void AddParameterToCalibrationManager();
    // method called event by event, aiming at extracting the 
    // physical information from detector
    void BuildPhysicalEvent();
    // same as BuildPhysicalEvent() method but with a simpler
    // treatment
    void BuildSimplePhysicalEvent();
    // same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};
    // activate raw data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputRaw();
    // activate physics data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputPhysics();
    // create branches of output ROOT file
    void InitializeRootOutput();
    // clear the raw and physical data objects event by event
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   
    // declare list of histograms

  //////////////////////////////////////////////////////////////
  // specific methods to Minos array
  public:
    
    NPL::Tracking*      Tracking_functions; //!
    
    /* int NclusterFit; //! */ 
    // fit function for the Q(t) signals at each pad
    static  double conv_fit(double *x, double *p);
    static double distance2(double x,double y,double z, double *p);
    static void SumDistance(int &, double *, double & sum, double * par,  int);

    
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {
    
    m_PreTreatedData->Clear();
    Q_Pad.clear();
    X_Pad.clear();
    Y_Pad.clear();
    Z_Pad.clear();
    T_Pad.clear();
    }

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TMinosData object to TMinosPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TMinosData* rawDataPointer) {m_EventData = rawDataPointer;}
     
  // objects are not written in the TTree
  private:
    TMinosData*         m_EventData;        //!
    TMinosData*         m_PreTreatedData;   //!
    TMinosPhysics*      m_EventPhysics;     //!
    TMinuit*            min;//!
  
TClonesArray data_result;//!
TMinosResult *minosdata_result;//!

  TClonesArray fitdata;//!

  // getters for raw and pre-treated data object
  public:
    TMinosData* GetRawData()        const {return m_EventData;}
    TMinosData* GetPreTreatedData() const {return m_PreTreatedData;}

    vector<double> GetPad_X()  {return X_Pad;} 
    vector<double> GetPad_Y()  {return Y_Pad;} 
    vector<double> GetPad_Z()  {return Z_Pad;} 
    vector<double> GetParFit1()  {return parFit1;} 
    vector<double> GetParFit2()  {return parFit2;} 
    vector<double> GetParFit3()  {return parFit3;} 
    vector<double> GetParFit4()  {return parFit4;} 
    double GetVertexZ()  {return Zvertex[0];} 

    /* TClonesArray fitdata; */
  
  /* fitdata.SetClass("TMinosClust"); */
  
  // parameters used in the analysis
  private:
   

    bool SimulationBool=false;//!
    
    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!

    // used in PreTreatedData()
    vector<int> clusterringboolTemp;//!
    vector<int> clusterringbool;//!
    vector<int> clusternbr;//!
    vector<int> clusterpads;//!
    
    vector<double> Xpad;//!
    vector<double> Ypad;//!
    vector<double> Qpad;//!
    vector<double> XpadNew;
    vector<double> YpadNew;
    vector<double> QpadNew;//!
    vector<double> ZpadNew;//!
    
    vector<double> XpadTemp;//!
    vector<double> YpadTemp;//!
    vector<double> QpadTemp;//!
    vector<double> ZpadTemp;//!
    
    double x_mm;//!
    double y_mm;//!
    double z_mm;//!
    double q_pad;//!
    double t_pad;//!
    double Chi2;//!
    double maxCharge;//!
    double ChargeBin;//!
    
    double hfit_max; //!
    double hfit_max_T; //!
    double T_min; //!
    double T_max; //!
    
    
    int Iteration;///
    int filter_result;//!
    int filled;//!
    int fit2DStatus;//!
    int indexfill;//!
    int trackNbr;//!
   
    TF1 *fit_function;//!
    /* TH1F *hfit = new TH1F("hfit","hfit",512,0,512);//! */
    TH1F *hfit ;//!
    
    int npoint_temp=0;//!
    int cluster_temp;//!
    int ringsum;//!
    int ringtouch[19];//!

    TGraph *gryztmp;//!
    TGraph* grxztmp;//!
    
    double zmax;//!
    vector<double> xin;//!
    vector<double> yin;//!
    vector<double> zin;//!
    vector<double> qin;//!
    vector<double> xout;//!
    vector<double> yout;//!
    vector<double> zout;//!
    vector<double> qout;//!
    vector<double> xoutprime;//!
    vector<double> youtprime;//!
    vector<double> zoutprime;//!
    vector<double> qoutprime;//!

    vector<double> TOTxoutprime;
    vector<double> TOTyoutprime;
    vector<double> TOTzoutprime;
    vector<double> TOTqout;
 
    vector<double> Xvertex;
    vector<double> Yvertex;
    vector<double> Zvertex;
    
    vector<double> sumTheta;
    vector<double> Theta1;
    vector<double> Theta2;

    vector<int> trackclusternbr;//!
    vector<int> tracknbr;//!
    
    double Tot_tracks=0; //!
    int trackNbr_FINAL;//!
    vector<TGraph> gryz;//!
    vector<TGraph> grxz;//!
    
    int npoint;//!
    int array_final;//!
    TVector3 point;//!

    // Variables only filled when trackNbr_FINAL>=1
    int allevt_2pfiltered;//!
    Double_t pStart[4];//!
    double parFit_temp[4];//!
    double err_temp[4];//!
    Double_t chi[2];//!
    Int_t fitStatus[2];//!
    Double_t arglist[10];//!
    Int_t iflag;//!
    int nvpar;//!
    int nparx;//!
    double amin;//!
    double edm;//!
    double errdef;//!
      
    vector<double> lenght;
    vector<double> chargeTot;
    vector<double> parFit1;
    vector<double> parFit2;
    vector<double> parFit3;
    vector<double> parFit4;
    vector<double> errFit1_local;
    vector<double> errFit2_local;
    vector<double> errFit3_local;
    vector<double> errFit4_local;
  

    double xv;//!
    double yv;//!
    double zv;//!
    
    double Dist_min;//!
    double Theta_tr1;//!
    double Theta_tr2;//!
  // number of detectors
  private:
    int m_NumberOfDetectors;  //!

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TMinosPhysics,1)  // MinosPhysics structure
};

#endif
