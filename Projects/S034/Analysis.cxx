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
  // m_field.LoadMap("field_map/180702-2,40T-3000.table.bin",10);

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
    //double brho_param[6]={FDC0->PosX/*+1.77*/, FDC0->PosY, tan(FDC0_ThetaX), tan(FDC0_PhiY), FDC2->PosX/*-252.55*/, FDC2->ThetaX};

    double brho_param[6]={FDC0->PosX/*+1.77*/, FDC0->PosY, 0, 0, FDC2->PosX/*-252.55*/, FDC2->ThetaX};
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

// -*- mode: c++ -*-
// 
// File simtree_energy.C generated by TMultiDimFit::MakeRealCode
// on Fri May 28 00:25:23 2021
// ROOT version 5.32/00
//
// This file contains the function 
//
//    double  e_fit(double *x); 
//
// For evaluating the parameterization obtained
// from TMultiDimFit and the point x
// 
// See TMultiDimFit class documentation for more information 
// 
//
// Static data variables
//
static int    e_gNVariables    = 6;
static int    e_gNCoefficients = 76;
static double e_gDMean         = 69.5201;
// Assignment to mean vector.
static double e_gXMean[] = {
  0.00503955, -0.0144156, -0.101042, -0.00374735, 0.252452, -0.512868 };

// Assignment to minimum vector.
static double e_gXMin[] = {
  -64.3965, -52.0748, -0.0335722, -0.0401226, -1065.23, -299.807 };

// Assignment to maximum vector.
static double e_gXMax[] = {
  52.3459, 55.8258, 0.0432727, 0.04243, 1147.99, 931.885 };

// Assignment to coefficients vector.
static double e_gCoefficient[] = {
  -871.253,
  -1919.06,
  -788.479,
  2645.64,
  443.115,
  984.07,
  243.005,
  3552.25,
  -157.556,
  2188.2,
  618.533,
  1176.32,
  1303.87,
  -1056.05,
  -772.684,
  6.86646,
  1343.46,
  3390.39,
  871.459,
  7388.91,
  -599.788,
  -4023.35,
  -299.575,
  1670.13,
  0.620256,
  10.0888,
  -0.557305,
  -1.12306,
  9.51285,
  0.728805,
  -0.620644,
  -4.05522,
  0.633792,
  -11.1877,
  -4.56183,
  -5.71034,
  -4.45376,
  1.38591,
  1.18557,
  -4.92592,
  -5.16675,
  -4.66713,
  -1506.82,
  -306.44,
  4269.09,
  2292.62,
  19.7398,
  2612.26,
  1682.7,
  -1163.29,
  -597.13,
  -7838.7,
  14404.4,
  6607.08,
  -532.085,
  -317.801,
  681.575,
  -546.576,
  -775.207,
  -3803.59,
  -1160.61,
  0.0848953,
  4.22188,
  -1.96972,
  0.794997,
  -2260.46,
  -1.21301,
  -7409.46,
  269.548,
  -0.00925447,
  2.77788,
  1.95742,
  0.453985,
  -0.188488,
  -0.675742,
  -0.181284
 };

// Assignment to error coefficients vector.
static double e_gCoefficientRMS[] = {
  66.6599,
  130.954,
  45.1508,
  146.993,
  49.7567,
  42.2773,
  8.59505,
  193.537,
  31.5244,
  92.7045,
  67.4399,
  97.1825,
  63.1346,
  63.8809,
  37.0527,
  0.0875267,
  99.2757,
  205.578,
  127.446,
  361.809,
  107.175,
  235.193,
  35.8634,
  105.166,
  2.78933,
  2.0882,
  3.14393,
  3.06267,
  8.25249,
  0.616472,
  1.47943,
  2.11076,
  0.930677,
  7.39957,
  2.26951,
  2.50688,
  0.95993,
  2.42862,
  0.393511,
  0.458493,
  2.0148,
  2.57894,
  72.1773,
  61.4191,
  180.529,
  189.341,
  7.78812,
  193.422,
  248.274,
  208.754,
  69.8817,
  458.225,
  704.917,
  400.536,
  30.129,
  31.4199,
  27.3927,
  40.2961,
  54.0649,
  207.698,
  79.9551,
  9.4176,
  1.89081,
  0.422792,
  1.61055,
  155.773,
  1.45066,
  404.667,
  36.5963,
  0.729026,
  2.55087,
  0.355723,
  1.83536,
  0.547759,
  2.82987,
  2.05714
 };

// Assignment to powers vector.
// The powers are stored row-wise, that is
//  p_ij = e_gPower[i * NVariables + j];
static int    e_gPower[] = {
  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  2,  1,
  1,  1,  1,  2,  1,  1,
  1,  1,  2,  1,  1,  1,
  1,  1,  1,  1,  1,  3,
  1,  1,  1,  1,  2,  2,
  1,  3,  1,  1,  1,  1,
  2,  1,  2,  1,  1,  1,
  1,  1,  1,  2,  1,  2,
  1,  2,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  2,  1,  1,  1,  1,  2,
  3,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  2,  1,
  1,  1,  2,  2,  1,  1,
  2,  1,  1,  1,  2,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  2,  1,  1,  3,  1,  1,
  1,  2,  1,  3,  1,  1,
  1,  1,  2,  3,  1,  1,
  1,  1,  3,  2,  1,  1,
  1,  2,  2,  2,  1,  1,
  1,  3,  1,  1,  2,  1,
  1,  3,  1,  2,  1,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  2,  1,  1,  2,  2,  1,
  1,  2,  1,  2,  2,  1,
  1,  1,  1,  3,  2,  1,
  3,  1,  1,  2,  1,  1,
  2,  1,  1,  1,  3,  1,
  1,  1,  2,  1,  3,  1,
  3,  2,  1,  1,  1,  1,
  1,  2,  3,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  2,  2,  2,  1,  1,  1,
  2,  1,  1,  2,  1,  2,
  1,  2,  1,  2,  1,  2,
  1,  1,  2,  2,  1,  2,
  1,  1,  1,  3,  1,  2,
  2,  1,  1,  1,  2,  2,
  1,  1,  2,  1,  2,  2,
  1,  1,  1,  2,  2,  2,
  2,  1,  1,  1,  1,  3,
  1,  2,  1,  1,  1,  3,
  1,  1,  2,  1,  1,  3,
  1,  2,  1,  1,  1,  1,
  1,  2,  1,  1,  1,  2,
  1,  2,  1,  1,  2,  1,
  2,  2,  1,  1,  1,  1,
  2,  1,  2,  2,  1,  1,
  2,  2,  1,  1,  2,  1,
  1,  1,  1,  2,  3,  1,
  2,  1,  3,  1,  1,  1,
  2,  2,  1,  1,  1,  2,
  3,  1,  2,  1,  1,  1,
  1,  2,  1,  1,  2,  2,
  1,  1,  1,  2,  1,  3,
  3,  1,  1,  1,  2,  1,
  1,  1,  2,  2,  2,  1,
  1,  2,  1,  1,  3,  1,
  2,  3,  1,  1,  1,  1,
  1,  1,  3,  1,  1,  1,
  2,  1,  2,  1,  2,  1,
  1,  3,  2,  1,  1,  1
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

// EOF for simtree_energy.C
// -*- mode: c++ -*-
// 
// File simtree_momentum.C generated by TMultiDimFit::MakeRealCode
// on Fri May 28 00:25:23 2021
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
static int    p_gNCoefficients = 76;
static double p_gDMean         = 361.833;
// Assignment to mean vector.
static double p_gXMean[] = {
  0.00503955, -0.0144156, -0.101042, -0.00374735, 0.252452, -0.512868 };

// Assignment to minimum vector.
static double p_gXMin[] = {
  -64.3965, -52.0748, -0.0335722, -0.0401226, -1065.23, -299.807 };

// Assignment to maximum vector.
static double p_gXMax[] = {
  52.3459, 55.8258, 0.0432727, 0.04243, 1147.99, 931.885 };

// Assignment to coefficients vector.
static double p_gCoefficient[] = {
  -1679.5,
  -3768.96,
  -1773.35,
  5927.56,
  -1052.24,
  786.483,
  2247.73,
  553.038,
  7920.2,
  -326.613,
  1147.19,
  3010.2,
  -1520.59,
  -2398.64,
  12.0318,
  2542.51,
  7843.96,
  1874.12,
  17655.7,
  -902.702,
  -9667.32,
  -2259.82,
  -673.898,
  3744.62,
  23.4495,
  -8.41591,
  23.1301,
  23.3445,
  -0.850973,
  -21.3182,
  -5.51343,
  -13.2501,
  0.0420299,
  4.8444,
  4.76409,
  -11.8105,
  -4406.58,
  -637.169,
  8498.88,
  4695.17,
  4938.56,
  3623.63,
  -1738.8,
  -1348.34,
  34420.8,
  15289.1,
  -1177.71,
  -587.346,
  1529.61,
  450.733,
  4356.32,
  2406.82,
  0.828932,
  -1494.72,
  -7.21335,
  -4.61121,
  2.36389,
  -9.47311,
  9.7888,
  -3.50166,
  -12.2852,
  -2912.85,
  -0.277797,
  42.1753,
  -1.26653,
  -18825.4,
  -8755.82,
  3.69156,
  8.16369,
  -7.45288,
  -0.273504,
  -9.0052,
  -17056.8,
  -0.145767,
  -1.13422,
  0.865558
 };

// Assignment to error coefficients vector.
static double p_gCoefficientRMS[] = {
  1.69031e-11,
  3.29223e-11,
  6.40666e-11,
  2.89085e-11,
  6.58543e-11,
  8.69785e-11,
  7.0898e-11,
  3.57517e-11,
  5.62969e-11,
  1.89724e-11,
  1.69397e-10,
  1.38078e-10,
  1.28278e-10,
  1.24785e-10,
  2.44521e-11,
  3.47315e-10,
  1.44807e-10,
  2.73924e-10,
  1.18391e-10,
  3.51888e-10,
  1.21525e-10,
  2.55786e-10,
  1.82279e-11,
  6.11644e-11,
  7.41771e-11,
  9.78353e-11,
  1.16067e-09,
  1.31561e-09,
  3.16376e-11,
  1.03428e-09,
  6.01156e-10,
  3.13162e-11,
  9.7543e-11,
  9.29464e-11,
  8.80057e-11,
  1.16199e-10,
  4.98203e-10,
  3.6952e-11,
  4.19168e-10,
  5.29444e-10,
  6.7647e-10,
  5.33585e-10,
  6.85262e-10,
  3.55029e-11,
  2.30548e-10,
  2.81963e-10,
  1.35502e-10,
  1.39257e-10,
  1.49982e-10,
  1.83996e-10,
  2.15232e-10,
  2.71861e-10,
  1.88434e-11,
  1.92304e-11,
  7.69361e-11,
  1.02682e-10,
  3.52953e-10,
  6.53796e-10,
  5.36745e-10,
  1.27112e-10,
  7.41853e-11,
  3.7455e-11,
  7.64436e-11,
  9.64104e-10,
  8.7203e-11,
  2.36674e-10,
  1.34232e-10,
  6.85812e-11,
  5.67547e-10,
  4.46781e-10,
  3.2384e-11,
  7.34435e-11,
  2.61417e-10,
  3.21076e-11,
  8.11735e-11,
  7.22006e-11
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
  1,  1,  1,  2,  1,  2,
  1,  1,  2,  1,  1,  2,
  1,  2,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  2,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  2,  1,
  1,  1,  2,  2,  1,  1,
  2,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  1,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  1,  2,  1,  3,  1,  1,
  1,  1,  3,  2,  1,  1,
  1,  2,  2,  2,  1,  1,
  2,  1,  2,  2,  1,  1,
  1,  3,  1,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  1,  2,  1,  2,  2,  1,
  1,  1,  1,  3,  2,  1,
  3,  1,  1,  2,  1,  1,
  2,  1,  1,  1,  3,  1,
  1,  2,  1,  1,  3,  1,
  1,  1,  2,  1,  3,  1,
  2,  2,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  2,  1,  1,  2,  1,  2,
  1,  2,  1,  2,  1,  2,
  1,  1,  2,  2,  1,  2,
  1,  1,  1,  3,  1,  2,
  1,  1,  2,  1,  2,  2,
  1,  1,  1,  2,  2,  2,
  2,  1,  1,  1,  1,  3,
  1,  2,  1,  1,  1,  3,
  1,  1,  2,  1,  1,  3,
  1,  1,  1,  2,  1,  3,
  2,  1,  2,  1,  1,  1,
  1,  2,  2,  1,  1,  1,
  1,  1,  3,  1,  1,  1,
  3,  1,  1,  1,  1,  1,
  1,  1,  2,  3,  1,  1,
  1,  3,  1,  2,  1,  1,
  2,  1,  2,  1,  2,  1,
  2,  1,  1,  2,  2,  1,
  1,  1,  2,  2,  2,  1,
  1,  1,  1,  2,  3,  1,
  3,  2,  1,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  2,  1,  3,  1,  1,  1,
  2,  2,  2,  1,  1,  1,
  3,  1,  2,  1,  1,  1,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  2,  1,
  2,  1,  1,  3,  1,  1,
  2,  2,  1,  1,  2,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  1,  2,  3,  1,  1,  1,
  1,  2,  1,  1,  2,  2,
  3,  1,  1,  1,  2,  1,
  1,  3,  2,  1,  1,  1,
  2,  3,  1,  1,  1,  1
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

// EOF for simtree_momentum.C
// -*- mode: c++ -*-
// 
// File simtree_rigidity.C generated by TMultiDimFit::MakeRealCode
// on Fri May 28 00:25:23 2021
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
static int    r_gNCoefficients = 76;
static double r_gDMean         = 4.82778;
// Assignment to mean vector.
static double r_gXMean[] = {
  0.00503955, -0.0144156, -0.101042, -0.00374735, 0.252452, -0.512868 };

// Assignment to minimum vector.
static double r_gXMin[] = {
  -64.3965, -52.0748, -0.0335722, -0.0401226, -1065.23, -299.807 };

// Assignment to maximum vector.
static double r_gXMax[] = {
  52.3459, 55.8258, 0.0432727, 0.04243, 1147.99, 931.885 };

// Assignment to coefficients vector.
static double r_gCoefficient[] = {
  -22.4088,
  -50.2876,
  -23.6611,
  79.0889,
  -14.0396,
  10.4937,
  29.9905,
  7.37895,
  105.676,
  -4.35786,
  15.3065,
  40.1638,
  -20.2886,
  -32.0041,
  0.160535,
  33.9237,
  104.659,
  25.0055,
  235.572,
  -12.0444,
  -128.987,
  -30.1518,
  -8.99153,
  49.9628,
  0.312876,
  -0.11229,
  0.308614,
  0.311475,
  -0.0113542,
  -0.28444,
  -0.0735633,
  -0.17679,
  0.000560787,
  0.0646367,
  0.0635652,
  -0.157583,
  -58.795,
  -8.50147,
  113.397,
  62.6455,
  65.8931,
  48.3485,
  -23.2001,
  -17.9904,
  459.262,
  203.996,
  -15.7137,
  -7.8367,
  20.4089,
  6.01393,
  58.1245,
  32.1131,
  0.0110601,
  -19.9434,
  -0.0962446,
  -0.0615253,
  0.0315403,
  -0.126396,
  0.130608,
  -0.0467211,
  -0.163916,
  -38.8649,
  -0.00370653,
  0.562727,
  -0.0168988,
  -251.179,
  -116.825,
  0.0492549,
  0.108925,
  -0.0994405,
  -0.00364924,
  -0.120152,
  -227.582,
  -0.0019449,
  -0.0151334,
  0.0115488
 };

// Assignment to error coefficients vector.
static double r_gCoefficientRMS[] = {
  1.69031e-11,
  3.29223e-11,
  6.40666e-11,
  2.89085e-11,
  6.58543e-11,
  8.69785e-11,
  7.0898e-11,
  3.57517e-11,
  5.62969e-11,
  1.89724e-11,
  1.69397e-10,
  1.38078e-10,
  1.28278e-10,
  1.24785e-10,
  2.44521e-11,
  3.47315e-10,
  1.44807e-10,
  2.73924e-10,
  1.18391e-10,
  3.51888e-10,
  1.21525e-10,
  2.55786e-10,
  1.82279e-11,
  6.11644e-11,
  7.41771e-11,
  9.78353e-11,
  1.16067e-09,
  1.31561e-09,
  3.16376e-11,
  1.03428e-09,
  6.01156e-10,
  3.13162e-11,
  9.7543e-11,
  9.29464e-11,
  8.80057e-11,
  1.16199e-10,
  4.98203e-10,
  3.6952e-11,
  4.19168e-10,
  5.29444e-10,
  6.7647e-10,
  5.33585e-10,
  6.85262e-10,
  3.55029e-11,
  2.30548e-10,
  2.81963e-10,
  1.35502e-10,
  1.39257e-10,
  1.49982e-10,
  1.83996e-10,
  2.15232e-10,
  2.71861e-10,
  1.88434e-11,
  1.92304e-11,
  7.69361e-11,
  1.02682e-10,
  3.52953e-10,
  6.53796e-10,
  5.36745e-10,
  1.27112e-10,
  7.41853e-11,
  3.7455e-11,
  7.64436e-11,
  9.64104e-10,
  8.7203e-11,
  2.36674e-10,
  1.34232e-10,
  6.85812e-11,
  5.67547e-10,
  4.46781e-10,
  3.2384e-11,
  7.34435e-11,
  2.61417e-10,
  3.21076e-11,
  8.11735e-11,
  7.22006e-11
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
  1,  1,  1,  2,  1,  2,
  1,  1,  2,  1,  1,  2,
  1,  2,  1,  1,  1,  2,
  2,  1,  1,  1,  1,  2,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  2,  1,
  1,  1,  2,  2,  1,  1,
  2,  1,  1,  1,  2,  1,
  2,  2,  1,  1,  1,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  1,  2,  1,  3,  1,  1,
  1,  1,  3,  2,  1,  1,
  1,  2,  2,  2,  1,  1,
  2,  1,  2,  2,  1,  1,
  1,  3,  1,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  1,  2,  1,  2,  2,  1,
  1,  1,  1,  3,  2,  1,
  3,  1,  1,  2,  1,  1,
  2,  1,  1,  1,  3,  1,
  1,  2,  1,  1,  3,  1,
  1,  1,  2,  1,  3,  1,
  2,  2,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  2,  1,  1,  2,  1,  2,
  1,  2,  1,  2,  1,  2,
  1,  1,  2,  2,  1,  2,
  1,  1,  1,  3,  1,  2,
  1,  1,  2,  1,  2,  2,
  1,  1,  1,  2,  2,  2,
  2,  1,  1,  1,  1,  3,
  1,  2,  1,  1,  1,  3,
  1,  1,  2,  1,  1,  3,
  1,  1,  1,  2,  1,  3,
  2,  1,  2,  1,  1,  1,
  1,  2,  2,  1,  1,  1,
  1,  1,  3,  1,  1,  1,
  3,  1,  1,  1,  1,  1,
  1,  1,  2,  3,  1,  1,
  1,  3,  1,  2,  1,  1,
  2,  1,  2,  1,  2,  1,
  2,  1,  1,  2,  2,  1,
  1,  1,  2,  2,  2,  1,
  1,  1,  1,  2,  3,  1,
  3,  2,  1,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  2,  1,  3,  1,  1,  1,
  2,  2,  2,  1,  1,  1,
  3,  1,  2,  1,  1,  1,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  2,  1,
  2,  1,  1,  3,  1,  1,
  2,  2,  1,  1,  2,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  1,  2,  3,  1,  1,  1,
  1,  2,  1,  1,  2,  2,
  3,  1,  1,  1,  2,  1,
  1,  3,  2,  1,  1,  1,
  2,  3,  1,  1,  1,  1
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

// EOF for simtree_rigidity.C
// -*- mode: c++ -*-
// 
// File simtree_length.C generated by TMultiDimFit::MakeRealCode
// on Fri May 28 00:25:23 2021
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
static int    l_gNCoefficients = 77;
static double l_gDMean         = 9270.11;
// Assignment to mean vector.
static double l_gXMean[] = {
  0.00503955, -0.0144156, -0.101042, -0.00374735, 0.252452, -0.512868 };

// Assignment to minimum vector.
static double l_gXMin[] = {
  -64.3965, -52.0748, -0.0335722, -0.0401226, -1065.23, -299.807 };

// Assignment to maximum vector.
static double l_gXMax[] = {
  52.3459, 55.8258, 0.0432727, 0.04243, 1147.99, 931.885 };

// Assignment to coefficients vector.
static double l_gCoefficient[] = {
  -6590.16,
  -13914.6,
  -2835.19,
  11286.2,
  -2318.17,
  1506.93,
  2806.41,
  1055.27,
  15246.9,
  -582.244,
  14287.3,
  2024.12,
  5826.32,
  2883.06,
  -3211.64,
  -2565.14,
  -3508.49,
  -3945.73,
  81.3328,
  6976.65,
  15491.2,
  3340.02,
  32961.2,
  -4257.58,
  -18445.7,
  -1204.58,
  7213.73,
  1.15164,
  56.7961,
  -6.55251,
  -3.1322,
  32.3364,
  17.5881,
  2.39586,
  -17.4049,
  4.99685,
  -44.3394,
  -18.7891,
  -32.238,
  -26.3708,
  6.54775,
  7.12538,
  -7694.27,
  -1133.08,
  27850.4,
  11364,
  -5002.96,
  13573.4,
  6437.24,
  -8275.75,
  -2420.82,
  -5.63532,
  -35904.8,
  64281.3,
  30190.8,
  -1440.27,
  2620.77,
  1005.06,
  -17864.9,
  -5198.67,
  1.04599,
  -1.61611,
  -3.90405,
  -15.7528,
  -9.48741,
  -25.0997,
  -10129.3,
  2.95363,
  99.1311,
  -34800.1,
  -2335.5,
  -8.51776,
  11.7194,
  9.32626,
  -25.8219,
  2.18449,
  -0.917242
 };

// Assignment to error coefficients vector.
static double l_gCoefficientRMS[] = {
  647.441,
  1258.13,
  378.249,
  1228.53,
  338.526,
  439.055,
  444.451,
  68.299,
  1618.83,
  242.058,
  929.339,
  595.678,
  754.718,
  687.373,
  454.652,
  379.081,
  542.258,
  290.313,
  1.10705,
  744.15,
  1703.8,
  1010.19,
  2870.46,
  856.075,
  1864.98,
  296.741,
  874.288,
  33.2233,
  23.4443,
  37.7505,
  37.4312,
  100.661,
  24.1279,
  7.88562,
  26.9518,
  11.868,
  89.8789,
  28.9366,
  32.0384,
  12.1687,
  29.4468,
  4.91973,
  565.232,
  471.598,
  1808.9,
  1470.56,
  738.43,
  1449.61,
  1967.42,
  1666.8,
  577.982,
  17.9786,
  3633.48,
  5592.56,
  3319.49,
  257.515,
  239.116,
  312.322,
  1701.52,
  563.862,
  9.27568,
  114.437,
  17.2737,
  5.63337,
  5.32805,
  24.809,
  1098.33,
  25.1201,
  96.0282,
  3315.12,
  242.887,
  36.0287,
  32.4886,
  4.45102,
  31.7885,
  20.2812,
  22.29
 };

// Assignment to powers vector.
// The powers are stored row-wise, that is
//  p_ij = l_gPower[i * NVariables + j];
static int    l_gPower[] = {
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
  1,  1,  1,  2,  1,  2,
  1,  2,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  1,  2,  1,  1,  1,  2,
  1,  1,  3,  1,  1,  1,
  2,  1,  1,  1,  1,  2,
  3,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  3,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  2,  1,
  1,  1,  2,  2,  1,  1,
  2,  1,  1,  1,  2,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  2,  1,  1,  3,  1,  1,
  1,  2,  1,  3,  1,  1,
  1,  1,  2,  3,  1,  1,
  1,  1,  3,  2,  1,  1,
  1,  2,  2,  2,  1,  1,
  2,  2,  1,  1,  2,  1,
  1,  3,  1,  1,  2,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  2,  1,  1,  2,  2,  1,
  1,  2,  1,  2,  2,  1,
  1,  1,  1,  3,  2,  1,
  3,  1,  1,  2,  1,  1,
  2,  1,  1,  1,  3,  1,
  3,  1,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  2,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  2,
  2,  1,  1,  2,  1,  2,
  1,  2,  1,  2,  1,  2,
  1,  1,  2,  2,  1,  2,
  1,  1,  1,  3,  1,  2,
  3,  1,  2,  1,  1,  1,
  2,  1,  1,  1,  2,  2,
  1,  1,  2,  1,  2,  2,
  1,  1,  1,  2,  2,  2,
  1,  2,  1,  1,  1,  3,
  1,  1,  2,  1,  1,  3,
  1,  1,  1,  2,  1,  3,
  1,  2,  1,  1,  2,  1,
  2,  2,  1,  1,  1,  1,
  3,  1,  1,  1,  2,  1,
  2,  1,  2,  2,  1,  1,
  1,  3,  1,  2,  1,  1,
  1,  1,  2,  1,  3,  1,
  1,  1,  1,  2,  3,  1,
  3,  2,  1,  1,  1,  1,
  2,  2,  1,  1,  1,  2,
  1,  3,  2,  1,  1,  1,
  2,  2,  2,  1,  1,  1,
  1,  2,  1,  1,  2,  2,
  2,  1,  1,  1,  1,  3,
  2,  1,  2,  1,  2,  1,
  1,  1,  2,  2,  2,  1,
  1,  2,  1,  1,  3,  1,
  1,  2,  3,  1,  1,  1,
  2,  1,  3,  1,  1,  1,
  2,  3,  1,  1,  1,  1
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

// EOF for simtree_length.C
// -*- mode: c++ -*-
// 
// File simtree_tof.C generated by TMultiDimFit::MakeRealCode
// on Fri May 28 00:25:23 2021
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
static int    t_gNCoefficients = 77;
static double t_gDMean         = 88.0547;
// Assignment to mean vector.
static double t_gXMean[] = {
  0.00503955, -0.0144156, -0.101042, -0.00374735, 0.252452, -0.512868 };

// Assignment to minimum vector.
static double t_gXMin[] = {
  -64.3965, -52.0748, -0.0335722, -0.0401226, -1065.23, -299.807 };

// Assignment to maximum vector.
static double t_gXMax[] = {
  52.3459, 55.8258, 0.0432727, 0.04243, 1147.99, 931.885 };

// Assignment to coefficients vector.
static double t_gCoefficient[] = {
  282.419,
  600.69,
  199.214,
  -610.381,
  16.8965,
  -181.103,
  -65.4683,
  27.7219,
  -586.193,
  -221.165,
  131.422,
  257.64,
  1.87553,
  -199.501,
  -2326.17,
  1262.92,
  88.6015,
  -380.386,
  -2.50251,
  0.192506,
  -0.0682083,
  -0.0605175,
  -0.00933975,
  2.54256,
  2.13194,
  1.70986,
  -0.734912,
  0.725959,
  0.566525,
  1.41423,
  262.324,
  54.1268,
  -1143.9,
  255.544,
  -5.4915,
  -383.569,
  158.012,
  177.209,
  2458.18,
  2069.51,
  -4536.21,
  -1851.8,
  144.418,
  21.6852,
  -147.382,
  49.0783,
  -791.882,
  13.6185,
  -269.147,
  75.7361,
  134.56,
  -230.177,
  -950.155,
  1062.38,
  81.3635,
  -0.315919,
  -1.87225,
  -0.378093,
  0.88221,
  -0.437349,
  -0.563023,
  1.42746,
  -525.101,
  -447.633,
  0.297019,
  20.3115,
  195.852,
  0.424978,
  -0.802113,
  0.778162,
  381.743,
  -0.338946,
  -0.109268,
  -0.201413,
  -0.0377476,
  -0.150836,
  0.0571845
 };

// Assignment to error coefficients vector.
static double t_gCoefficientRMS[] = {
  51.6658,
  100.306,
  31.0942,
  100.473,
  36.5499,
  36.1352,
  5.48829,
  19.3199,
  74.8745,
  55.5697,
  31.1116,
  44.4306,
  0.110749,
  81.8808,
  231.073,
  150.427,
  24.4632,
  71.4537,
  2.1768,
  3.60448,
  0.792011,
  1.63035,
  3.6143,
  8.62908,
  3.21448,
  1.21207,
  0.490659,
  0.55933,
  0.532584,
  2.39502,
  45.341,
  37.6402,
  145.644,
  60.5957,
  9.27703,
  159.431,
  133.982,
  47.623,
  293.067,
  271.101,
  450.201,
  272.399,
  19.7108,
  21.145,
  19.4434,
  28.0278,
  132.355,
  49.4246,
  60.5858,
  37.5262,
  23.2998,
  58.9152,
  139.819,
  139.147,
  68.8468,
  3.6068,
  9.68357,
  10.9918,
  2.70767,
  1.19191,
  0.444182,
  3.0704,
  118.051,
  114.732,
  1.74076,
  25.8244,
  43.6475,
  3.15977,
  2.42206,
  2.89927,
  84.993,
  3.25617,
  1.97584,
  2.41757,
  0.930273,
  2.8245,
  2.14083
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
  1,  3,  1,  1,  1,  1,
  2,  1,  2,  1,  1,  1,
  1,  1,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  1,
  2,  1,  1,  1,  1,  2,
  1,  1,  1,  1,  3,  1,
  1,  2,  1,  2,  1,  1,
  1,  1,  2,  1,  2,  1,
  2,  1,  1,  1,  2,  1,
  1,  1,  1,  3,  1,  1,
  1,  1,  1,  1,  2,  3,
  1,  2,  1,  3,  1,  1,
  1,  1,  3,  2,  1,  1,
  1,  3,  1,  1,  2,  1,
  1,  3,  1,  2,  1,  1,
  2,  1,  2,  1,  2,  1,
  2,  2,  1,  2,  1,  1,
  1,  2,  1,  2,  2,  1,
  1,  1,  1,  3,  2,  1,
  2,  1,  1,  1,  3,  1,
  1,  1,  2,  1,  3,  1,
  1,  1,  1,  2,  3,  1,
  3,  2,  1,  1,  1,  1,
  3,  1,  1,  1,  1,  2,
  1,  3,  1,  1,  1,  2,
  2,  1,  2,  1,  1,  2,
  1,  1,  3,  1,  1,  2,
  2,  2,  2,  1,  1,  1,
  1,  2,  1,  2,  1,  2,
  1,  1,  2,  2,  1,  2,
  1,  1,  1,  3,  1,  2,
  2,  1,  1,  1,  2,  2,
  1,  2,  1,  1,  2,  2,
  1,  1,  2,  1,  2,  2,
  1,  1,  1,  2,  2,  2,
  2,  1,  1,  1,  1,  3,
  1,  2,  1,  1,  1,  3,
  1,  1,  2,  1,  1,  3,
  1,  2,  1,  1,  1,  1,
  1,  1,  1,  1,  2,  2,
  1,  1,  1,  2,  1,  2,
  1,  2,  2,  1,  1,  1,
  1,  2,  1,  1,  1,  2,
  3,  1,  1,  1,  1,  1,
  2,  1,  1,  2,  1,  1,
  1,  1,  1,  2,  2,  1,
  1,  2,  1,  1,  2,  1,
  1,  1,  2,  2,  1,  1,
  1,  1,  2,  3,  1,  1,
  1,  2,  2,  2,  1,  1,
  2,  1,  2,  2,  1,  1,
  1,  2,  2,  1,  2,  1,
  1,  1,  3,  1,  2,  1,
  1,  2,  1,  1,  3,  1,
  1,  2,  3,  1,  1,  1,
  1,  2,  2,  1,  1,  2,
  2,  1,  1,  2,  1,  2,
  3,  1,  2,  1,  1,  1,
  1,  1,  1,  2,  1,  3,
  2,  2,  1,  1,  1,  1,
  2,  1,  1,  3,  1,  1,
  2,  2,  1,  1,  2,  1,
  2,  1,  1,  2,  2,  1,
  2,  2,  1,  1,  1,  2,
  1,  1,  2,  2,  2,  1,
  2,  1,  3,  1,  1,  1,
  1,  3,  2,  1,  1,  1,
  3,  1,  1,  1,  2,  1,
  3,  1,  1,  2,  1,  1,
  2,  3,  1,  1,  1,  1
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

// EOF for simtree_tof.C


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

