/*****************************************************************************
 * Copyright (C) 2009-2021    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta contact address: matta@lpccaen.in2p3.fr         *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  S034 analysis project                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/



#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"RootInput.h"
#include"RootOutput.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
   Minos= (TMinosPhysics*) m_DetectorManager->GetDetector("Minos");
   //Nebula = (TNebulaPhysics*) m_DetectorManager->GetDetector("NEBULA");
   BDC = (TSamuraiBDCPhysics*) m_DetectorManager->GetDetector("SAMURAIBDC");
   FDC0 = (TSamuraiFDC0Physics*) m_DetectorManager->GetDetector("SAMURAIFDC0");
   FDC2 = (TSamuraiFDC2Physics*) m_DetectorManager->GetDetector("SAMURAIFDC2");
   Hodo = (TSamuraiHodoscopePhysics*) m_DetectorManager->GetDetector("SAMURAIHOD");
   m_field.LoadMap("field_map/180702-2,40T-3000.table.bin");

   InitOutputBranch();
   InitInputBranch();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  Clear();
  //cout << Trigger << " " ; 
  Trigger=Trigger&0x00ff;
//cout << Trigger << endl;
  // Compute Brho 
  if(FDC2->PosX>-10000 && FDC0->PosX>-10000 ){ // if both are correctly build
   // Compute ThetaX and PhiY using Minos vertex and FDC0 XY
   double FDC0_ThetaX = FDC0->ThetaX;
   double FDC0_PhiY   = FDC0->PhiY;

   if(Minos->Z_Vertex>0){
    FDC0_ThetaX = atan((FDC0->PosX-Minos->X_Vertex)/(1283.7-Minos->Z_Vertex));
    FDC0_PhiY   = atan((FDC0->PosY-Minos->Y_Vertex)/(1283.7-Minos->Z_Vertex));
   } 
    double brho_param[6]={FDC0->PosX/*+1.77*/, FDC0->PosY, tan(FDC0_ThetaX), tan(FDC0_PhiY), FDC2->PosX/*-252.55*/, FDC2->ThetaX};
    Brho=r_fit(brho_param);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Clear(){
  Brho=-1000;
  Beta_f=-1000;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Brho",&Brho,"Brho/D");
  RootOutput::getInstance()->GetTree()->Branch("Beta_f",&Beta_f,"Beta_f/D");
  RootOutput::getInstance()->GetTree()->Branch("Trigger",&Trigger,"Trigger/I");
} 
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
    RootInput::getInstance()->GetChain()->SetBranchAddress("Trigger",&Trigger);
}

////////////////////////////////////////////////////////////////////////////////
//
static int    e_gNVariables    = 6;
static int    e_gNCoefficients = 56;
static double e_gDMean         = 140.96;
// Assignment to mean vector.
static double e_gXMean[] = {
  0.0025732, -0.284362, -0.566756, -0.946845, 0.361982, -0.52295 };

// Assignment to minimum vector.
static double e_gXMin[] = {
  -331.189, -196.066, -257.522, -304.613, -333.804, -292.129 };

// Assignment to maximum vector.
static double e_gXMax[] = {
  329.322, 351.904, 931.885, 11176.9, 1147.99, 933.394 };

// Assignment to coefficients vector.
static double e_gCoefficient[] = {
  -1915.5,
  -3376.77,
  9001.24,
  1517.31,
  6.5778,
  -296.117,
  -3507.57,
  102.821,
  1180.97,
  366.41,
  4244.46,
  1292.38,
  -6569.85,
  48.611,
  -1145.28,
  6475.24,
  -461.55,
  618.848,
  8140.28,
  1052.09,
  149.209,
  761.819,
  -184.68,
  43.8596,
  2527.57,
  -12.1075,
  -4.13165,
  -6.7228,
  -2.09925,
  43.2175,
  13.6207,
  -1162.11,
  2455.15,
  309.178,
  -4.67617,
  0.881522,
  -185.94,
  -21.6162,
  -55.7686,
  -1034.58,
  -23.6721,
  2179.29,
  740.309,
  8385.12,
  2721.25,
  -2174.5,
  159.758,
  1176.06,
  -4487.32,
  1365.38,
  5.61114,
  177.229,
  268.174,
  -45.8185,
  -85.2077,
  -45.3608
 };

// Assignment to error coefficients vector.
static double e_gCoefficientRMS[] = {
  164.828,
  294.513,
  506.942,
  116.759,
  111.515,
  20.2561,
  184.244,
  25.4792,
  228.82,
  34.9984,
  275.86,
  98.3259,
  349.254,
  189.385,
  85.2052,
  437.683,
  54.1842,
  46.0126,
  447.321,
  62.0962,
  36.481,
  177.681,
  19.1419,
  49.0136,
  155.81,
  65.5158,
  2.72524,
  4.55078,
  1.88398,
  48.9402,
  14.3597,
  100.324,
  135.04,
  59.2191,
  0.820262,
  0.701362,
  21.6537,
  7.42057,
  20.6061,
  86.9029,
  18.9372,
  164.455,
  54.5583,
  549.607,
  184.047,
  158.979,
  48.9787,
  334.83,
  254.413,
  85.1526,
  6.5432,
  22.5547,
  83.4628,
  25.3156,
  84.8766,
  48.3489
 };

// Assignment to powers vector.
// The powers are stored row-wise, that is
//  p_ij = e_gPower[i * NVariables + j];
static int    e_gPower[] = {
  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  2,  1,
  1,  2,  1,  1,  1,  1,
  1,  1,  1,  2,  1,  1,
  1,  1,  2,  1,  1,  1,
  1,  1,  1,  1,  1,  3,
  1,  1,  1,  1,  2,  2,
  1,  3,  1,  1,  1,  1,
  2,  1,  2,  1,  1,  1,
  1,  2,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  1,  2,  1,  1,  1,  2,
  1,  1,  3,  1,  1,  1,
  2,  1,  1,  1,  1,  2,
  3,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  1,  2,  1,  2,  1,
  2,  1,  1,  1,  2,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  2,  1,  1,  3,  1,  1,
  1,  2,  1,  3,  1,  1,
  3,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  2,  1,
  1,  3,  1,  1,  2,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  1,  2,  1,  2,  2,  1,
  1,  1,  1,  3,  2,  1,
  2,  1,  1,  1,  3,  1,
  1,  2,  1,  1,  3,  1,
  1,  1,  2,  1,  3,  1,
  3,  2,  1,  1,  1,  1,
  1,  2,  3,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  2,  1,  3,  1,  1,  1,
  2,  2,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  2,
  3,  1,  2,  1,  1,  1,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  2,  2,
  1,  1,  1,  1,  3,  2,
  2,  3,  1,  1,  1,  1,
  1,  2,  1,  1,  1,  3,
  2,  1,  2,  1,  2,  1,
  1,  3,  2,  1,  1,  1,
  2,  2,  2,  1,  1,  1,
  2,  1,  1,  1,  1,  3
};

// 
// The function   double e_fit(double *x)
// 
double e_fit(double *x) {
  double returnValue = e_gDMean;
  int    i = 0, j = 0, k = 0;
  for (i = 0; i < e_gNCoefficients ; i++) {
    // Evaluate the ith term in the expansion
    double term = e_gCoefficient[i];
    for (j = 0; j < e_gNVariables; j++) {
      // Evaluate the polynomial in the jth variable.
      int power = e_gPower[e_gNVariables * i + j]; 
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / (e_gXMax[j] - e_gXMin[j]) * (x[j] - e_gXMax[j]);
      // what is the power to use!
      switch(power) {
      case 1: r = 1; break; 
      case 2: r = v; break; 
      default: 
        p2 = v; 
        for (k = 3; k <= power; k++) { 
          p3 = p2 * v;
          p3 = 2 * v * p2 - p1; 
          p1 = p2; p2 = p3; 
        }
        r = p3;
      }
      // multiply this term by the poly in the jth var
      term *= r; 
    }
    // Add this term to the final result
    returnValue += term;
  }
  return returnValue;
}

// EOF for ../Geant4/samurai/simtree_energy.C
// -*- mode: c++ -*-
// 
// File ../Geant4/samurai/simtree_momentum.C generated by TMultiDimFit::MakeRealCode
// on Tue Jun  9 23:41:17 2020
// ROOT version 5.32/00
//
// This file contains the function 
//
//    double  p_fit(double *x); 
//
// For evaluating the parameterization obtained
// from TMultiDimFit and the point x
// 
// See TMultiDimFit class documentation for more information 
// 
//
// Static data variables
//
static int    p_gNVariables    = 6;
static int    p_gNCoefficients = 57;
static double p_gDMean         = 524.342;
// Assignment to mean vector.
static double p_gXMean[] = {
  0.0025732, -0.284362, -0.566756, -0.946845, 0.361982, -0.52295 };

// Assignment to minimum vector.
static double p_gXMin[] = {
  -331.189, -196.066, -257.522, -304.613, -333.804, -292.129 };

// Assignment to maximum vector.
static double p_gXMax[] = {
  329.322, 351.904, 931.885, 11176.9, 1147.99, 933.394 };

// Assignment to coefficients vector.
static double p_gCoefficient[] = {
  -4150.11,
  -7349.66,
  -1376.07,
  726.65,
  -604.755,
  157.175,
  -7402.61,
  137.906,
  1266,
  607.703,
  9432.11,
  2336.63,
  -12820,
  -981.145,
  -2219,
  -1146.1,
  1069.02,
  -9607.95,
  -135.98,
  236.78,
  -5477.97,
  2213.85,
  6282.03,
  89.2823,
  -10.974,
  -104.594,
  -62.7333,
  -9.89892,
  -21.8523,
  -18.8223,
  16.9863,
  3546.43,
  -461.655,
  -17.9853,
  -2026.01,
  -134.173,
  5565.29,
  18164.2,
  5239.56,
  -4192.18,
  22.9582,
  -15479.1,
  -114.231,
  2974.51,
  -10402.5,
  2516.63,
  25.4881,
  303.465,
  17.2719,
  1202.72,
  50.8832,
  -1.17495,
  1284.76,
  17.0545,
  -0.894684,
  -136.872,
  -134.371
 };

// Assignment to error coefficients vector.
static double p_gCoefficientRMS[] = {
  5.0702e-12,
  9.69107e-12,
  3.89769e-11,
  1.13183e-11,
  1.53099e-11,
  5.35431e-12,
  8.94257e-12,
  1.12034e-11,
  2.16374e-11,
  6.32208e-12,
  6.87455e-11,
  2.70028e-11,
  1.70926e-11,
  2.9263e-11,
  1.4199e-11,
  5.24502e-12,
  6.53991e-12,
  4.1161e-11,
  1.19525e-11,
  1.99626e-11,
  3.47403e-11,
  9.32059e-11,
  1.19932e-10,
  6.39059e-12,
  2.49985e-11,
  4.91273e-11,
  1.92969e-11,
  2.95574e-10,
  1.39773e-11,
  6.12734e-11,
  3.16965e-11,
  1.26652e-10,
  1.15348e-11,
  1.5818e-11,
  1.00252e-11,
  1.09154e-10,
  2.29228e-10,
  1.31396e-10,
  5.16128e-11,
  2.71397e-11,
  2.1153e-10,
  7.86728e-11,
  9.25091e-12,
  1.78175e-10,
  6.64139e-11,
  1.25002e-11,
  8.61292e-11,
  3.38293e-11,
  1.16568e-11,
  1.64392e-10,
  1.42658e-11,
  4.97717e-11,
  1.20839e-11,
  4.82003e-11,
  1.9685e-11,
  4.28747e-11,
  1.11506e-11
 };

// Assignment to powers vector.
// The powers are stored row-wise, that is
//  p_ij = p_gPower[i * NVariables + j];
static int    p_gPower[] = {
  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  2,  1,
  1,  2,  1,  1,  1,  1,
  1,  1,  1,  2,  1,  1,
  1,  1,  2,  1,  1,  1,
  1,  1,  1,  1,  1,  3,
  1,  1,  1,  1,  2,  2,
  1,  3,  1,  1,  1,  1,
  2,  1,  2,  1,  1,  1,
  1,  2,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  1,  2,  1,  1,  1,  2,
  1,  1,  3,  1,  1,  1,
  3,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  1,  2,  1,  2,  1,
  1,  2,  1,  1,  2,  1,
  2,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  1,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  2,  1,  1,  3,  1,  1,
  1,  2,  1,  3,  1,  1,
  2,  2,  1,  1,  2,  1,
  1,  3,  1,  1,  2,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  3,  1,
  3,  2,  1,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  2,  1,  3,  1,  1,  1,
  2,  2,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  2,
  2,  2,  2,  1,  1,  1,
  2,  1,  1,  2,  1,  2,
  3,  1,  2,  1,  1,  1,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  2,  2,
  1,  1,  1,  1,  3,  2,
  2,  1,  1,  1,  1,  3,
  1,  2,  1,  1,  1,  3,
  3,  1,  1,  1,  2,  1,
  2,  1,  2,  1,  2,  1,
  1,  1,  1,  3,  2,  1,
  2,  1,  1,  1,  3,  1,
  1,  3,  1,  1,  1,  2,
  2,  3,  1,  1,  1,  1,
  1,  2,  1,  1,  3,  1,
  1,  2,  3,  1,  1,  1,
  1,  3,  2,  1,  1,  1
};

// 
// The function   double p_fit(double *x)
// 
double p_fit(double *x) {
  double returnValue = p_gDMean;
  int    i = 0, j = 0, k = 0;
  for (i = 0; i < p_gNCoefficients ; i++) {
    // Evaluate the ith term in the expansion
    double term = p_gCoefficient[i];
    for (j = 0; j < p_gNVariables; j++) {
      // Evaluate the polynomial in the jth variable.
      int power = p_gPower[p_gNVariables * i + j]; 
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / (p_gXMax[j] - p_gXMin[j]) * (x[j] - p_gXMax[j]);
      // what is the power to use!
      switch(power) {
      case 1: r = 1; break; 
      case 2: r = v; break; 
      default: 
        p2 = v; 
        for (k = 3; k <= power; k++) { 
          p3 = p2 * v;
          p3 = 2 * v * p2 - p1; 
          p1 = p2; p2 = p3; 
        }
        r = p3;
      }
      // multiply this term by the poly in the jth var
      term *= r; 
    }
    // Add this term to the final result
    returnValue += term;
  }
  return returnValue;
}

// EOF for ../Geant4/samurai/simtree_momentum.C
// -*- mode: c++ -*-
// 
// File ../Geant4/samurai/simtree_rigidity.C generated by TMultiDimFit::MakeRealCode
// on Tue Jun  9 23:41:17 2020
// ROOT version 5.32/00
//
// This file contains the function 
//
//    double  r_fit(double *x); 
//
// For evaluating the parameterization obtained
// from TMultiDimFit and the point x
// 
// See TMultiDimFit class documentation for more information 
// 
//
// Static data variables
//
static int    r_gNVariables    = 6;
static int    r_gNCoefficients = 57;
static double r_gDMean         = 5.24705;
// Assignment to mean vector.
static double r_gXMean[] = {
  0.0025732, -0.284362, -0.566756, -0.946845, 0.361982, -0.52295 };

// Assignment to minimum vector.
static double r_gXMin[] = {
  -331.189, -196.066, -257.522, -304.613, -333.804, -292.129 };

// Assignment to maximum vector.
static double r_gXMax[] = {
  329.322, 351.904, 931.885, 11176.9, 1147.99, 933.394 };

// Assignment to coefficients vector.
static double r_gCoefficient[] = {
  -41.5298,
  -73.5475,
  -13.7703,
  7.27153,
  -6.05174,
  1.57283,
  -74.0774,
  1.38001,
  12.6688,
  6.08124,
  94.3864,
  23.3825,
  -128.289,
  -9.81824,
  -22.2054,
  -11.4689,
  10.6976,
  -96.146,
  -1.36074,
  2.36944,
  -54.8176,
  22.1538,
  62.8638,
  0.893442,
  -0.109816,
  -1.04666,
  -0.627767,
  -0.0990577,
  -0.218674,
  -0.188353,
  0.16998,
  35.4888,
  -4.61974,
  -0.179978,
  -20.2741,
  -1.34266,
  55.6915,
  181.768,
  52.4319,
  -41.9508,
  0.229741,
  -154.898,
  -1.1431,
  29.7657,
  -104.097,
  25.1837,
  0.255057,
  3.03675,
  0.172838,
  12.0355,
  0.509184,
  -0.0117577,
  12.8565,
  0.170664,
  -0.00895304,
  -1.36967,
  -1.34464
 };

// Assignment to error coefficients vector.
static double r_gCoefficientRMS[] = {
  5.0702e-12,
  9.69107e-12,
  3.89769e-11,
  1.13183e-11,
  1.53099e-11,
  5.35431e-12,
  8.94257e-12,
  1.12034e-11,
  2.16374e-11,
  6.32208e-12,
  6.87455e-11,
  2.70028e-11,
  1.70926e-11,
  2.9263e-11,
  1.4199e-11,
  5.24502e-12,
  6.53991e-12,
  4.1161e-11,
  1.19525e-11,
  1.99626e-11,
  3.47403e-11,
  9.32059e-11,
  1.19932e-10,
  6.39059e-12,
  2.49985e-11,
  4.91273e-11,
  1.92969e-11,
  2.95574e-10,
  1.39773e-11,
  6.12734e-11,
  3.16965e-11,
  1.26652e-10,
  1.15348e-11,
  1.5818e-11,
  1.00252e-11,
  1.09154e-10,
  2.29228e-10,
  1.31396e-10,
  5.16128e-11,
  2.71397e-11,
  2.1153e-10,
  7.86728e-11,
  9.25091e-12,
  1.78175e-10,
  6.64139e-11,
  1.25002e-11,
  8.61292e-11,
  3.38293e-11,
  1.16568e-11,
  1.64392e-10,
  1.42658e-11,
  4.97717e-11,
  1.20839e-11,
  4.82003e-11,
  1.9685e-11,
  4.28747e-11,
  1.11506e-11
 };

// Assignment to powers vector.
// The powers are stored row-wise, that is
//  p_ij = r_gPower[i * NVariables + j];
static int    r_gPower[] = {
  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  2,  1,
  1,  2,  1,  1,  1,  1,
  1,  1,  1,  2,  1,  1,
  1,  1,  2,  1,  1,  1,
  1,  1,  1,  1,  1,  3,
  1,  1,  1,  1,  2,  2,
  1,  3,  1,  1,  1,  1,
  2,  1,  2,  1,  1,  1,
  1,  2,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  1,  2,  1,  1,  1,  2,
  1,  1,  3,  1,  1,  1,
  3,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  1,  2,  1,  2,  1,
  1,  2,  1,  1,  2,  1,
  2,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  1,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  2,  1,  1,  3,  1,  1,
  1,  2,  1,  3,  1,  1,
  2,  2,  1,  1,  2,  1,
  1,  3,  1,  1,  2,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  3,  1,
  3,  2,  1,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  2,  1,  3,  1,  1,  1,
  2,  2,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  2,
  2,  2,  2,  1,  1,  1,
  2,  1,  1,  2,  1,  2,
  3,  1,  2,  1,  1,  1,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  2,  2,
  1,  1,  1,  1,  3,  2,
  2,  1,  1,  1,  1,  3,
  1,  2,  1,  1,  1,  3,
  3,  1,  1,  1,  2,  1,
  2,  1,  2,  1,  2,  1,
  1,  1,  1,  3,  2,  1,
  2,  1,  1,  1,  3,  1,
  1,  3,  1,  1,  1,  2,
  2,  3,  1,  1,  1,  1,
  1,  2,  1,  1,  3,  1,
  1,  2,  3,  1,  1,  1,
  1,  3,  2,  1,  1,  1
};

// 
// The function   double r_fit(double *x)
// 
double r_fit(double *x) {
  double returnValue = r_gDMean;
  int    i = 0, j = 0, k = 0;
  for (i = 0; i < r_gNCoefficients ; i++) {
    // Evaluate the ith term in the expansion
    double term = r_gCoefficient[i];
    for (j = 0; j < r_gNVariables; j++) {
      // Evaluate the polynomial in the jth variable.
      int power = r_gPower[r_gNVariables * i + j]; 
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / (r_gXMax[j] - r_gXMin[j]) * (x[j] - r_gXMax[j]);
      // what is the power to use!
      switch(power) {
      case 1: r = 1; break; 
      case 2: r = v; break; 
      default: 
        p2 = v; 
        for (k = 3; k <= power; k++) { 
          p3 = p2 * v;
          p3 = 2 * v * p2 - p1; 
          p1 = p2; p2 = p3; 
        }
        r = p3;
      }
      // multiply this term by the poly in the jth var
      term *= r; 
    }
    // Add this term to the final result
    returnValue += term;
  }
  return returnValue;
}

// EOF for ../Geant4/samurai/simtree_rigidity.C
// -*- mode: c++ -*-
// 
// File ../Geant4/samurai/simtree_length.C generated by TMultiDimFit::MakeRealCode
// on Tue Jun  9 23:41:17 2020
// ROOT version 5.32/00
//
// This file contains the function 
//
//    double  l_fit(double *x); 
//
// For evaluating the parameterization obtained
// from TMultiDimFit and the point x
// 
// See TMultiDimFit class documentation for more information 
// 
//
// Static data variables
//
static int    l_gNVariables    = 6;
static int    l_gNCoefficients = 56;
static double l_gDMean         = 9302.09;
// Assignment to mean vector.
static double l_gXMean[] = {
  0.0025732, -0.284362, -0.566756, -0.946845, 0.361982, -0.52295 };

// Assignment to minimum vector.
static double l_gXMin[] = {
  -331.189, -196.066, -257.522, -304.613, -333.804, -292.129 };

// Assignment to maximum vector.
static double l_gXMax[] = {
  329.322, 351.904, 931.885, 11176.9, 1147.99, 933.394 };

// Assignment to coefficients vector.
static double l_gCoefficient[] = {
  -5063.74,
  -8210.4,
  8531.82,
  -734.412,
  -360.971,
  -8616.61,
  331.099,
  -854.689,
  13989.9,
  -14949,
  -2334.33,
  23725.3,
  -232.89,
  863.204,
  -5135.67,
  -671.923,
  193.744,
  -7892.44,
  5470.48,
  2106.04,
  220.14,
  392.627,
  -750.825,
  -13.1656,
  -15.4854,
  -15.0127,
  -692.948,
  -46.1235,
  -10.895,
  -2288.18,
  13.9869,
  -7.15341,
  0.581274,
  42.6349,
  -8.22167,
  -116.588,
  8175.36,
  1080.47,
  129.214,
  26910.7,
  7002.79,
  -4598.69,
  -14.7257,
  2944.04,
  -467.226,
  -15008,
  2105.01,
  6.06541,
  361.147,
  -1073.29,
  3745.34,
  -284.423,
  1666.15,
  -499.767,
  8786.59,
  119.522
 };

// Assignment to error coefficients vector.
static double l_gCoefficientRMS[] = {
  1092.36,
  1874.66,
  2316.91,
  992.765,
  293.062,
  1327.01,
  164.035,
  1469.78,
  1822.87,
  2291.78,
  587.934,
  2842.24,
  222.135,
  300.266,
  2772.37,
  341.844,
  284.914,
  849.513,
  1225.2,
  5295.19,
  203.834,
  313.389,
  883.734,
  23.072,
  38.4748,
  15.8304,
  233.869,
  387.214,
  119.356,
  5419.4,
  456.337,
  6.95297,
  5.94522,
  60.7779,
  168.339,
  150.975,
  1134.29,
  329.863,
  198.355,
  3547.52,
  1330.81,
  1077.78,
  656.951,
  1360.04,
  385.463,
  1653.15,
  541.211,
  56.6718,
  276.195,
  707.914,
  681.038,
  733.694,
  700.282,
  173.163,
  2230.31,
  136.351
 };

// Assignment to powers vector.
// The powers are stored row-wise, that is
//  p_ij = l_gPower[i * NVariables + j];
static int    l_gPower[] = {
  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  2,  1,
  1,  1,  1,  2,  1,  1,
  1,  1,  2,  1,  1,  1,
  1,  1,  1,  1,  1,  3,
  1,  1,  1,  1,  2,  2,
  2,  1,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  1,
  2,  1,  1,  1,  1,  2,
  3,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  1,  2,  1,  2,  1,
  1,  2,  1,  1,  2,  1,
  2,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  1,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  2,  1,  1,  3,  1,  1,
  3,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  2,  1,
  1,  3,  1,  1,  2,  1,
  1,  3,  1,  2,  1,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  1,  1,  1,  3,  2,  1,
  2,  1,  1,  1,  3,  1,
  1,  2,  1,  1,  3,  1,
  3,  2,  1,  1,  1,  1,
  1,  2,  3,  1,  1,  1,
  2,  1,  3,  1,  1,  1,
  2,  2,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  1,  3,  2,  1,  1,  1,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  2,
  2,  2,  2,  1,  1,  1,
  1,  2,  1,  2,  1,  2,
  3,  1,  2,  1,  1,  1,
  1,  2,  1,  1,  2,  2,
  1,  1,  1,  1,  3,  2,
  2,  3,  1,  1,  1,  1,
  2,  1,  1,  1,  1,  3,
  1,  2,  1,  1,  1,  1,
  1,  2,  2,  1,  1,  1,
  1,  2,  1,  3,  1,  1,
  2,  1,  2,  1,  2,  1,
  1,  1,  2,  1,  3,  1,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  1,  3
};

// 
// The function   double l_fit(double *x)
// 
double l_fit(double *x) {
  double returnValue = l_gDMean;
  int    i = 0, j = 0, k = 0;
  for (i = 0; i < l_gNCoefficients ; i++) {
    // Evaluate the ith term in the expansion
    double term = l_gCoefficient[i];
    for (j = 0; j < l_gNVariables; j++) {
      // Evaluate the polynomial in the jth variable.
      int power = l_gPower[l_gNVariables * i + j]; 
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / (l_gXMax[j] - l_gXMin[j]) * (x[j] - l_gXMax[j]);
      // what is the power to use!
      switch(power) {
      case 1: r = 1; break; 
      case 2: r = v; break; 
      default: 
        p2 = v; 
        for (k = 3; k <= power; k++) { 
          p3 = p2 * v;
          p3 = 2 * v * p2 - p1; 
          p1 = p2; p2 = p3; 
        }
        r = p3;
      }
      // multiply this term by the poly in the jth var
      term *= r; 
    }
    // Add this term to the final result
    returnValue += term;
  }
  return returnValue;
}

// EOF for ../Geant4/samurai/simtree_length.C
// -*- mode: c++ -*-
// 
// File ../Geant4/samurai/simtree_tof.C generated by TMultiDimFit::MakeRealCode
// on Tue Jun  9 23:41:17 2020
// ROOT version 5.32/00
//
// This file contains the function 
//
//    double  t_fit(double *x); 
//
// For evaluating the parameterization obtained
// from TMultiDimFit and the point x
// 
// See TMultiDimFit class documentation for more information 
// 
//
// Static data variables
//
static int    t_gNVariables    = 6;
static int    t_gNCoefficients = 56;
static double t_gDMean         = 63.5106;
// Assignment to mean vector.
static double t_gXMean[] = {
  0.0025732, -0.284362, -0.566756, -0.946845, 0.361982, -0.52295 };

// Assignment to minimum vector.
static double t_gXMin[] = {
  -331.189, -196.066, -257.522, -304.613, -333.804, -292.129 };

// Assignment to maximum vector.
static double t_gXMax[] = {
  329.322, 351.904, 931.885, 11176.9, 1147.99, 933.394 };

// Assignment to coefficients vector.
static double t_gCoefficient[] = {
  296.575,
  509.773,
  -1298.15,
  -177.102,
  38.3448,
  495.622,
  -7.92753,
  -88.3852,
  -44.7681,
  -204.021,
  916.644,
  158.841,
  -954.907,
  54.8928,
  -71.5493,
  -1160.79,
  -151.56,
  -6.48142,
  -20.8229,
  -95.551,
  21.8388,
  4.53363,
  -368.211,
  -11.1978,
  -0.122623,
  0.899953,
  0.664196,
  -4.26981,
  -2.76381,
  -353.964,
  -41.6085,
  -0.104552,
  28.9895,
  2.23953,
  8.78264,
  123.105,
  3.83715,
  -338.603,
  -92.1278,
  -1248.05,
  -432.176,
  299.722,
  10.2068,
  17.7278,
  -20.2405,
  -135.529,
  644.109,
  -168.691,
  -0.0126702,
  -23.4315,
  -633.876,
  -171.644,
  -48.5535,
  7.73714,
  -3.93291,
  0.0801397
 };

// Assignment to error coefficients vector.
static double t_gCoefficientRMS[] = {
  90.0633,
  153.552,
  286.58,
  57.4084,
  10.6353,
  98.8219,
  12.2663,
  110.128,
  18.1831,
  52.623,
  177.01,
  44.7459,
  220.758,
  29.2743,
  22.5181,
  280.17,
  36.8213,
  66.015,
  23.1365,
  90.5589,
  10.4606,
  23.5416,
  99.1486,
  50.5072,
  1.90998,
  3.21059,
  1.31771,
  31.4861,
  9.39573,
  67.4273,
  33.8813,
  0.497038,
  12.9849,
  4.9474,
  13.3374,
  41.0747,
  11.7911,
  86.8346,
  27.0137,
  280.477,
  98.0117,
  82.1288,
  53.8616,
  98.7791,
  31.4828,
  165.676,
  125.204,
  40.4042,
  4.27406,
  11.177,
  143.484,
  55.4336,
  52.6846,
  15.5019,
  22.798,
  0.581552
 };

// Assignment to powers vector.
// The powers are stored row-wise, that is
//  p_ij = t_gPower[i * NVariables + j];
static int    t_gPower[] = {
  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  2,  1,
  1,  1,  1,  2,  1,  1,
  1,  1,  2,  1,  1,  1,
  1,  1,  1,  1,  1,  3,
  1,  1,  1,  1,  2,  2,
  1,  3,  1,  1,  1,  1,
  1,  2,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  1,
  2,  1,  1,  1,  1,  2,
  3,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  2,  1,
  2,  1,  1,  1,  2,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  2,  1,  1,  3,  1,  1,
  1,  2,  1,  3,  1,  1,
  3,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  2,  1,
  1,  3,  1,  1,  2,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  1,  2,  1,  2,  2,  1,
  1,  1,  1,  3,  2,  1,
  1,  2,  1,  1,  3,  1,
  1,  1,  2,  1,  3,  1,
  3,  2,  1,  1,  1,  1,
  1,  2,  3,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  2,  1,  3,  1,  1,  1,
  2,  2,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  2,
  2,  2,  2,  1,  1,  1,
  1,  2,  1,  2,  1,  2,
  3,  1,  2,  1,  1,  1,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  2,  2,
  1,  1,  1,  1,  3,  2,
  2,  3,  1,  1,  1,  1,
  1,  2,  1,  1,  1,  3,
  2,  1,  2,  1,  1,  1,
  2,  2,  1,  1,  1,  1,
  2,  1,  2,  1,  2,  1,
  1,  3,  2,  1,  1,  1,
  2,  1,  1,  1,  1,  3,
  2,  1,  1,  1,  3,  1
};

// 
// The function   double t_fit(double *x)
// 
double t_fit(double *x) {
  double returnValue = t_gDMean;
  int    i = 0, j = 0, k = 0;
  for (i = 0; i < t_gNCoefficients ; i++) {
    // Evaluate the ith term in the expansion
    double term = t_gCoefficient[i];
    for (j = 0; j < t_gNVariables; j++) {
      // Evaluate the polynomial in the jth variable.
      int power = t_gPower[t_gNVariables * i + j]; 
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / (t_gXMax[j] - t_gXMin[j]) * (x[j] - t_gXMax[j]);
      // what is the power to use!
      switch(power) {
      case 1: r = 1; break; 
      case 2: r = v; break; 
      default: 
        p2 = v; 
        for (k = 3; k <= power; k++) { 
          p3 = p2 * v;
          p3 = 2 * v * p2 - p1; 
          p1 = p2; p2 = p3; 
        }
        r = p3;
      }
      // multiply this term by the poly in the jth var
      term *= r; 
    }
    // Add this term to the final result
    returnValue += term;
  }
  return returnValue;
}

// EOF for ../Geant4/samurai/simtree_tof.C


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy{
  public:
    proxy(){
      NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
    }
};

proxy p;
}

