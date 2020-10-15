#ifndef TSAMURAIFDC2PHYSICS_H
#define TSAMURAIFDC2PHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiFDC2 treated data                                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <vector>
#include <map>

// NPL
#include "TSamuraiFDC2Data.h"
//#include "TSamuraiFDC2Spectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"
#include "NPDCReconstruction.h"
// ROOT 
#include "TVector3.h" 
// Forward declaration
//class TSamuraiFDC2Spectra;


using namespace std ;
// little class to index each of the DC wire
class SamuraiDCIndex{
  public:
   SamuraiDCIndex(){};  
   ~SamuraiDCIndex(){};  
   SamuraiDCIndex(unsigned int det, unsigned int layer, unsigned int wire){
     m_det=det;
     m_layer=layer;
     m_wire=wire;
     m_norme=Norme();
   };  
  
  unsigned int m_det;
  unsigned int m_layer;
  unsigned int m_wire;
  unsigned int m_norme;
    
  int Norme() const {return (m_det*1000000000+m_layer*1000000+m_wire);} ;
  bool operator<(const SamuraiDCIndex i2){
    return this->Norme()<i2.Norme();
    }
  
  friend bool operator<(const SamuraiDCIndex i1,const SamuraiDCIndex i2){
    return i1.Norme()<i2.Norme();
    }
  
  friend bool operator==(const SamuraiDCIndex i1,const SamuraiDCIndex i2){
   return i1.Norme()==i2.Norme();
   }
  };


class TSamuraiFDC2Physics : public TObject, public NPL::VDetector{
  public:
    TSamuraiFDC2Physics();
    ~TSamuraiFDC2Physics() {};

  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    //   Provide Physical Multiplicity
    vector<double> DriftLength;
    vector<int> Detector;
    vector<int> Layer;
    vector<int> Wire;
    vector<double> Time;
    vector<double> ToT;
    vector<bool>   Matched;
    // Computed variable
    vector<TVector3> ParticleDirection;
    vector<TVector3> MiddlePosition;

    double PosX;
    double PosY;
    double ThetaX;
    double PhiY;
    TVector3 Dir;
    int Mult;
    int MultMean;
    int PileUp;

  public:
    // Projected position at given Z plan
    TVector3 ProjectedPosition(double Z);

  private: // Charateristic of the DC 
    void AddDC(string name, NPL::XmlParser&);//! take the XML file and fill in Wire_X and Layer_Angle
    map<SamuraiDCIndex,double> Wire_X;//! X position of the wires
    map<SamuraiDCIndex,double> Wire_Z;//! Z position of the wires
    map<SamuraiDCIndex,double> Wire_Angle;//! Wire Angle (0 for X, 90 for Y, U and V are typically at +/-30)
  
  private: // Analysis
    double ToTThreshold;//! a ToT threshold to remove noise
    void RemoveNoise();
    // Construct the 2D track and ref position at Z=0 and Z=100 based on X,Z and Radius provided

    // Object use to perform the DC reconstruction
    NPL::DCReconstruction m_reconstruction;//!

  public: //   Innherited from VDetector Class

    // Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
    void ReadConfiguration(NPL::InputParser) ;


    // Add Parameter to the CalibrationManger
    void AddParameterToCalibrationManager() ;      

    // Activated associated Branches and link it to the private member DetectorData address
    // In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
    void InitializeRootInputRaw() ;

    // Activated associated Branches and link it to the private member DetectorPhysics address
    // In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ;

    // Create associated branches and associated private member DetectorPhysics address
    void InitializeRootOutput() ;

    // This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
    void BuildPhysicalEvent() ;

    // Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
    // This method aimed to be used for analysis performed during experiment, when speed is requiered.
    // NB: This method can eventually be the same as BuildPhysicalEvent.
    void BuildSimplePhysicalEvent() ;

    // Same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    // Those two method all to clear the Event Physics or Data
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra();
    // Fill the spectra hold by the spectra class
    void FillSpectra();
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra();
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra();
    // Write Spectra to file
    void WriteSpectra();

  public:      //   Specific to SamuraiFDC2 Array

    //   Clear The PreTeated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    // Retrieve raw and pre-treated data
    TSamuraiFDC2Data* GetRawData()        const {return m_EventData;}
    TSamuraiFDC2Data* GetPreTreatedData() const {return m_PreTreatedData;}

  private:   //   Root Input and Output tree classes
    TSamuraiFDC2Data*         m_EventData;//!
    TSamuraiFDC2Data*         m_PreTreatedData;//!
    TSamuraiFDC2Physics*      m_EventPhysics;//!


  private: // Spectra Class
   // TSamuraiFDC2Spectra* m_Spectra; // !

  public: // Spectra Getter
    map< string , TH1*> GetSpectra(); 

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    ClassDef(TSamuraiFDC2Physics,1)  // SamuraiFDC2Physics structure
};

#endif
