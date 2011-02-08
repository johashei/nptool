/*****************************************************************************
 * Copyright (C) 2009-2010   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 10/06/09                                                 *
 * Last update    : 04/09/09                                                 *
 *---------------------------------------------------------------------------*
 * Decription: This class records all the information concerning the event   *
 *             generators, e.g. vertex of interaction, angles of emitted     *
 *             particles...                                                  *
 *             This class derives from TObject (ROOT) and its aim is to be   *
 *             stored in the output TTree of the G4 simulation               *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + 04/09/09: Add private members for emittance  (N. de Sereville)       *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#ifndef __INITIALCONDITIONS__
#define __INITIALCONDITIONS__

#include <vector>
#include "TObject.h"
using namespace std ;


class TInitialConditions : public TObject
{
private:
   // Incident particle properties (before interactions in the target)
   // Vertex of interaction
   vector<Double_t>   fIC_Position_X;
   vector<Double_t>   fIC_Position_Y;
   vector<Double_t>   fIC_Position_Z;
   // Theta and Phi angles for the emittance
   vector<Double_t>   fIC_Incident_Emittance_Theta;
   vector<Double_t>   fIC_Incident_Emittance_Phi;
   // Incident particle angles
   vector<Double_t>   fIC_Incident_Angle_Theta;
   vector<Double_t>   fIC_Incident_Angle_Phi;
   // Incident particle energy
   vector<Double_t>   fIC_Incident_Energy;

   // Emitted particle properties (after interactions in the target)
   vector<Double_t>   fIC_Emitted_Angle_ThetaCM;
   // Emitted particle angles in the incident frame
   vector<Double_t>   fIC_Emitted_Angle_ThetaLab_IncidentFrame;
   vector<Double_t>   fIC_Emitted_Angle_Phi_IncidentFrame;
   // Emitted particle angles in the world frame
   vector<Double_t>   fIC_Emitted_Angle_ThetaLab_WorldFrame;
   vector<Double_t>   fIC_Emitted_Angle_Phi_WorldFrame;
   // Emittedparticle energy
   vector<Double_t>   fIC_Emitted_Energy;


public:
   TInitialConditions();
   virtual ~TInitialConditions();

   void  Clear();
   void  Clear(const Option_t*) {};
   void  Dump() const;

   /////////////////////           SETTERS           ////////////////////////
   // Incident particle properties (before interactions in the target)
   // Vertex of interaction
   void SetICPositionX(Double_t PositionX)      {fIC_Position_X.push_back(PositionX);}
   void SetICPositionY(Double_t PositionY)      {fIC_Position_Y.push_back(PositionY);}
   void SetICPositionZ(Double_t PositionZ)      {fIC_Position_Z.push_back(PositionZ);}
   // Theta and Phi angles for the emittance
   void SetICIncidentEmittanceTheta(Double_t Theta)   {fIC_Incident_Emittance_Theta.push_back(Theta);}
   void SetICIncidentEmittancePhi(Double_t Phi)       {fIC_Incident_Emittance_Phi.push_back(Phi);}
   // Incident particle angles
   void SetICIncidentAngleTheta(Double_t AngleTheta) {fIC_Incident_Angle_Theta.push_back(AngleTheta);}
   void SetICIncidentAnglePhi(Double_t AnglePhi)     {fIC_Incident_Angle_Phi.push_back(AnglePhi);}
   // Incident particle energy
   void SetICIncidentEnergy(Double_t Energy)         {fIC_Incident_Energy.push_back(Energy);}
   
   // Emitted particle angles
   // Center of mass
   void SetICEmittedAngleThetaCM(Double_t AngleTheta) {fIC_Emitted_Angle_ThetaCM.push_back(AngleTheta);}
   // Angles in the incident frame
   void SetICEmittedAngleThetaLabIncidentFrame(Double_t AngleTheta)   {fIC_Emitted_Angle_ThetaLab_IncidentFrame.push_back(AngleTheta);}
   void SetICEmittedAnglePhiIncidentFrame(Double_t AnglePhi)          {fIC_Emitted_Angle_Phi_IncidentFrame.push_back(AnglePhi);}
   // Angles in the world frame
   void SetICEmittedAngleThetaLabWorldFrame(Double_t AngleTheta)   {fIC_Emitted_Angle_ThetaLab_WorldFrame.push_back(AngleTheta);}
   void SetICEmittedAnglePhiWorldFrame(Double_t AnglePhi)          {fIC_Emitted_Angle_Phi_WorldFrame.push_back(AnglePhi);}
   // Emitted particle energy
   void SetICEmittedEnergy(Double_t Energy) {fIC_Emitted_Energy.push_back(Energy);}


   /////////////////////           GETTERS           ////////////////////////
   // Incident particle properties (before interactions in the target)
   // Vertex of interaction
   Double_t GetICPositionX(Int_t i) {return fIC_Position_X.at(i);}
   Double_t GetICPositionY(Int_t i) {return fIC_Position_Y.at(i);}
   Double_t GetICPositionZ(Int_t i) {return fIC_Position_Z.at(i);}
   // Theta and Phi angles for the emittance
   Double_t GetICIncidentEmittanceTheta(Int_t i) {return fIC_Incident_Emittance_Theta.at(i);}
   Double_t GetICIncidentEmittancePhi(Int_t i)   {return fIC_Incident_Emittance_Phi.at(i);}
   // Incident particle angles
   Double_t GetICIncidentAngleTheta(Int_t i)   {return fIC_Incident_Angle_Theta.at(i);}
   Double_t GetICIncidentAnglePhi(Int_t i)     {return fIC_Incident_Angle_Phi.at(i);}
   // Incident particle energy
   Double_t GetICIncidentEnergy(Int_t i)   {return fIC_Incident_Energy.at(i);}
   
   // Emitted particle angles
   // Center of Mass
   Double_t GetICEmittedAngleThetaCM(Int_t i) {return fIC_Emitted_Angle_ThetaCM.at(i);}
   // Angles in the incident frame
   Double_t GetICEmittedAngleThetaLabIncidentFrame(Int_t i) {return fIC_Emitted_Angle_ThetaLab_IncidentFrame.at(i);}
   Double_t GetICEmittedAnglePhiIncidentFrame(Int_t i)      {return fIC_Emitted_Angle_Phi_IncidentFrame.at(i);}
   // Angles in the world frame
   Double_t GetICEmittedAngleThetaLabWorldFrame(Int_t i) {return fIC_Emitted_Angle_ThetaLab_WorldFrame.at(i);}
   Double_t GetICEmittedAnglePhiWorldFrame(Int_t i)      {return fIC_Emitted_Angle_Phi_WorldFrame.at(i);}
   // Emitted particle energy
   Double_t GetICEmittedEnergy(Int_t i) {return fIC_Emitted_Energy.at(i);}

   ClassDef(TInitialConditions, 1) // InitialConditions structure
};

#endif
