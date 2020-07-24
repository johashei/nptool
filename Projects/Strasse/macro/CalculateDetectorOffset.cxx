#include <iostream>
#include <ctime>
#include <cstdlib>
using namespace std;

// ROOT headers
#include "TString.h"

// nptool headers
#include "NPInputParser.h"

using namespace NPL;


void CalculateDetectorOffset(const char * fname = "strasse_optimized"){

  // Open output ROOT file from NPTool simulation run
  string path = "";
  string inFileName = fname;
  inFileName += ".detector";


  InputParser *inParser = new InputParser(inFileName,true);
  vector<NPL::InputBlock*> blocks_info = inParser->GetAllBlocksWithTokenAndValue("Strasse","Info");

  if(blocks_info.size()>1){
    cout << "ERROR: can only accepte one info block, " << blocks_info.size() << " info block founds." << endl; 
    exit(1); 
  }


  vector<string> info = {
    "Inner_Wafer_Length",         
    "Inner_Wafer_Width",          
    "Inner_Wafer_Thickness",     
    "Inner_Wafer_AlThickness",    
    "Inner_Wafer_PADExternal",    
    "Inner_Wafer_PADInternal",  
    "Inner_Wafer_GuardRing",    
    "Inner_PCB_PortWidth",      
    "Inner_PCB_StarboardWidth", 
    "Inner_PCB_BevelAngle",     
    "Inner_PCB_UpstreamWidth",  
    "Inner_PCB_DownstreamWidth",
    "Inner_PCB_MidWidth",       
    "Inner_PCB_Thickness",      
    "Inner_Wafer_TransverseStrips",
    "Inner_Wafer_LongitudinalStrips",
    "Outer_Wafer_Length",       
    "Outer_Wafer_Width",        
    "Outer_Wafer_Thickness",    
    "Outer_Wafer_AlThickness",  
    "Outer_Wafer_PADExternal",  
    "Outer_Wafer_PADInternal",  
    "Outer_Wafer_GuardRing",    
    "Outer_PCB_PortWidth",      
    "Outer_PCB_StarboardWidth", 
    "Outer_PCB_BevelAngle",     
    "Outer_PCB_UpstreamWidth",  
    "Outer_PCB_DownstreamWidth",
    "Outer_PCB_MidWidth",       
    "Outer_PCB_Thickness",      
    "Outer_Wafer_TransverseStrips",
    "Outer_Wafer_LongitudinalStrips",
    "Chamber_Thickness",
    "Chamber_Cylinder_Length",
    "Chamber_Radius",
    "Chamber_ExitTube_Radius",
    "Chamber_ExitTube_Length",
    "Chamber_Flange_Inner_Radius",
    "Chamber_Sphere_Radius",
    "Chamber_Sphere_Shift"
  };

  ////////////////////
  // Inner Detector //
  ////////////////////
  // Wafer parameter
  double Inner_Wafer_Length=-999;
  double Inner_Wafer_Width=-999;
  double Inner_Wafer_Thickness=-999;
  double Inner_Wafer_AlThickness=-999;
  double Inner_Wafer_PADExternal=-999;
  double Inner_Wafer_PADInternal=-999;
  double Inner_Wafer_GuardRing=-999;

  // PCB parameter
  double Inner_PCB_PortWidth=-999;
  double Inner_PCB_StarboardWidth=-999;
  double Inner_PCB_BevelAngle=-999;
  double Inner_PCB_UpstreamWidth=-999;
  double Inner_PCB_DownstreamWidth=-999;
  double Inner_PCB_MidWidth=-999;
  double Inner_PCB_Thickness=-999;
  double Inner_Wafer_TransverseStrips=-999;
  double Inner_Wafer_LongitudinalStrips=-999;

  ////////////////////
  // Outer Detector //
  ////////////////////
  // Wafer parameter
  double Outer_Wafer_Length=-999;
  double Outer_Wafer_Width=-999;
  double Outer_Wafer_Thickness=-999;
  double Outer_Wafer_AlThickness=-999;
  double Outer_Wafer_PADExternal=-999;
  double Outer_Wafer_PADInternal=-999;
  double Outer_Wafer_GuardRing=-999;

  // PCB parameter
  double Outer_PCB_PortWidth=-999;
  double Outer_PCB_StarboardWidth=-999;
  double Outer_PCB_BevelAngle=-999;
  double Outer_PCB_UpstreamWidth=-999;
  double Outer_PCB_DownstreamWidth=-999;
  double Outer_PCB_MidWidth=-999;
  double Outer_PCB_Thickness=-999;
  double Outer_Wafer_TransverseStrips=-999;
  double Outer_Wafer_LongitudinalStrips=-999;

  // Vacuum Chamber //
  double Chamber_Thickness=-999;
  double Chamber_Cylinder_Length=-999;
  double Chamber_Radius=-999;
  double Chamber_ExitTube_Radius=-999;
  double Chamber_ExitTube_Length=-999;
  double Chamber_Flange_Inner_Radius=-999;
  double Chamber_Sphere_Radius=-999;
  double Chamber_Sphere_Shift=-999;

  if(blocks_info[0]->HasTokenList(info)){
    cout << endl << "////  Strasse info block" <<  endl;
    Inner_Wafer_Length = blocks_info[0]->GetDouble("Inner_Wafer_Length","mm");
    Inner_Wafer_Width = blocks_info[0]->GetDouble("Inner_Wafer_Width","mm");          
    Inner_Wafer_Thickness = blocks_info[0]->GetDouble("Inner_Wafer_Thickness","micrometer");      
    Inner_Wafer_AlThickness = blocks_info[0]->GetDouble("Inner_Wafer_AlThickness","micrometer");     
    Inner_Wafer_PADExternal = blocks_info[0]->GetDouble("Inner_Wafer_PADExternal","mm");     
    Inner_Wafer_PADInternal = blocks_info[0]->GetDouble("Inner_Wafer_PADInternal","mm");   
    Inner_Wafer_GuardRing = blocks_info[0]->GetDouble("Inner_Wafer_GuardRing","mm");     
    Inner_Wafer_TransverseStrips = blocks_info[0]->GetInt("Inner_Wafer_TransverseStrips");        
    Inner_Wafer_LongitudinalStrips = blocks_info[0]->GetInt("Inner_Wafer_LongitudinalStrips");       
    Inner_PCB_PortWidth = blocks_info[0]->GetDouble("Inner_PCB_PortWidth","mm");       
    Inner_PCB_StarboardWidth = blocks_info[0]->GetDouble("Inner_PCB_StarboardWidth","mm");  
    Inner_PCB_BevelAngle = blocks_info[0]->GetDouble("Inner_PCB_BevelAngle","mm");      
    Inner_PCB_UpstreamWidth = blocks_info[0]->GetDouble("Inner_PCB_UpstreamWidth","mm");   
    Inner_PCB_DownstreamWidth = blocks_info[0]->GetDouble("Inner_PCB_DownstreamWidth","mm"); 
    Inner_PCB_MidWidth = blocks_info[0]->GetDouble("Inner_PCB_MidWidth","mm");        
    Inner_PCB_Thickness = blocks_info[0]->GetDouble("Inner_PCB_Thickness","mm");       
    Outer_Wafer_Length = blocks_info[0]->GetDouble("Outer_Wafer_Length","mm");        
    Outer_Wafer_Width = blocks_info[0]->GetDouble("Outer_Wafer_Width","mm");         
    Outer_Wafer_Thickness = blocks_info[0]->GetDouble("Outer_Wafer_Thickness","mm");     
    Outer_Wafer_AlThickness = blocks_info[0]->GetDouble("Outer_Wafer_AlThickness","micrometer");   
    Outer_Wafer_PADExternal = blocks_info[0]->GetDouble("Outer_Wafer_PADExternal","mm");   
    Outer_Wafer_PADInternal = blocks_info[0]->GetDouble("Outer_Wafer_PADInternal","mm");   
    Outer_Wafer_GuardRing = blocks_info[0]->GetDouble("Outer_Wafer_GuardRing","mm");     
    Outer_Wafer_TransverseStrips = blocks_info[0]->GetInt("Outer_Wafer_TransverseStrips");        
    Outer_Wafer_LongitudinalStrips = blocks_info[0]->GetInt("Outer_Wafer_LongitudinalStrips");       
    Outer_PCB_PortWidth = blocks_info[0]->GetDouble("Outer_PCB_PortWidth","mm");       
    Outer_PCB_StarboardWidth = blocks_info[0]->GetDouble("Outer_PCB_StarboardWidth","mm");  
    Outer_PCB_BevelAngle = blocks_info[0]->GetDouble("Outer_PCB_BevelAngle","deg");      
    Outer_PCB_UpstreamWidth = blocks_info[0]->GetDouble("Outer_PCB_UpstreamWidth","mm");   
    Outer_PCB_DownstreamWidth = blocks_info[0]->GetDouble("Outer_PCB_DownstreamWidth","mm"); 
    Outer_PCB_MidWidth = blocks_info[0]->GetDouble("Outer_PCB_MidWidth","mm");        
    Outer_PCB_Thickness = blocks_info[0]->GetDouble("Outer_PCB_Thickness","mm");       
    Chamber_Thickness= blocks_info[0]->GetDouble("Chamber_Thickness","mm"); 
    Chamber_Cylinder_Length= blocks_info[0]->GetDouble("Chamber_Cylinder_Length","mm");        
    Chamber_Radius= blocks_info[0]->GetDouble("Chamber_Radius","mm");       
    Chamber_ExitTube_Radius=blocks_info[0]->GetDouble("Chamber_ExitTube_Radius","mm");
    Chamber_ExitTube_Length=blocks_info[0]->GetDouble("Chamber_ExitTube_Length","mm");
    Chamber_Flange_Inner_Radius=blocks_info[0]->GetDouble("Chamber_Flange_Inner_Radius","mm");
    Chamber_Sphere_Radius=blocks_info[0]->GetDouble("Chamber_Sphere_Radius","mm");
    Chamber_Sphere_Shift=blocks_info[0]->GetDouble("Chamber_Sphere_Shift","mm");
  }




  vector<NPL::InputBlock*> starget = inParser->GetAllBlocksWithToken("Target");

    double TargetThickness = -999; 
    double TargetX = -999; 
    double TargetY = -999; 
    double TargetZ = -999; 

  if(starget.size()==1){
    cout << "////       TARGET      ////" << endl;
    cout << "//// Solid Target found " << endl;
    vector<string> token = {"Thickness","Radius","Material","Angle","X","Y","Z"};
    if(starget[0]->HasTokenList(token)){
      TargetThickness= starget[0]->GetDouble("Thickness","mm");
      TargetX=starget[0]->GetDouble("X","mm");
      TargetY=starget[0]->GetDouble("Y","mm");
      TargetZ=starget[0]->GetDouble("Z","mm");
    }
    else{
      cout << "ERROR: Target token list incomplete, check your input file" << endl;
      exit(1);
    }
  }
 

//////////////////////////////////////////////////


double d_TargetFront_InnerActive = 15; //mm
double d_TargetFront_OuterActive = 43; //mm

double InnerTotalLength = 
                      Inner_PCB_UpstreamWidth
                    + Inner_Wafer_Length * 2
                    + Inner_PCB_MidWidth
                    + Inner_PCB_DownstreamWidth;

double d_TargetCenter_InnerCenter = 
                    - TargetThickness/2. 
                    + d_TargetFront_InnerActive
                    + InnerTotalLength/2.
                    - Inner_Wafer_GuardRing
                    - Inner_PCB_UpstreamWidth;

double OuterTotalLength = 
                      Outer_PCB_UpstreamWidth
                    + Outer_Wafer_Length * 2
                    + Outer_PCB_MidWidth
                    + Outer_PCB_DownstreamWidth;

double d_TargetCenter_OuterCenter = 
                    - TargetThickness/2. 
                    + d_TargetFront_OuterActive
                    + OuterTotalLength/2.
                    - Outer_Wafer_GuardRing
                    - Outer_PCB_UpstreamWidth;


cout<<endl;
cout<< "----------  INPUT DISTANCES (mm)  -------------"<<endl;
cout<<endl;
cout<<"Target Thickness : "<<TargetThickness<<endl;
cout<<"Beginning Target - Beginning Inner Active : "<<d_TargetFront_InnerActive<<endl;
cout<<"Beginning Target - Beginning Outer Active : "<<d_TargetFront_OuterActive<<endl;
cout<<endl;
cout<< "--------- CALCULATED DISTANCES (mm) -----------"<<endl;

cout << "InnerTotalLength = "<<InnerTotalLength<<endl; 
cout << "d_TargetCenter_InnerCenter = "<<d_TargetCenter_InnerCenter<<endl; 
cout<<endl;

cout << "OuterTotalLength = "<<OuterTotalLength<<endl; 
cout << "d_TargetCenter_OuterCenter = "<<d_TargetCenter_OuterCenter<<endl; 
cout<<endl;
cout<< "---------------------------------- -----------"<<endl;
cout<<"Remark: this calculation assumes that the center of target is at (0,0,0)"<<endl;
}


