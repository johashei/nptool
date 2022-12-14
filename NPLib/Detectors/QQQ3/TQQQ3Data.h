#ifndef __QQQ3DATA__
#define __QQQ3DATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold the QQQ3 Silicon array raw data (Made for TIG64 card)   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>

// ROOT
#include "TNamed.h"

class TQQQ3Data : public TNamed {
 private:
  // QQQ3
  // Energy
  std::vector<UShort_t> fQQQ3_StripFront_DetectorNbr;
  std::vector<UShort_t> fQQQ3_StripFront_StripNbr;
  std::vector<Double_t> fQQQ3_StripFront_Energy;
  std::vector<Double_t> fQQQ3_StripFront_TimeCFD;
  std::vector<Double_t> fQQQ3_StripFront_TimeLED;
  std::vector<Double_t> fQQQ3_StripFront_Time;

  std::vector<UShort_t> fQQQ3_StripBack_DetectorNbr;
  std::vector<UShort_t> fQQQ3_StripBack_StripNbr;
  std::vector<Double_t> fQQQ3_StripBack_Energy;
  std::vector<Double_t> fQQQ3_StripBack_TimeCFD;
  std::vector<Double_t> fQQQ3_StripBack_TimeLED;
  std::vector<Double_t> fQQQ3_StripBack_Time;

  std::vector<UShort_t> fQQQ3_PAD_DetectorNbr;
  std::vector<Double_t> fQQQ3_PAD_Energy;
  std::vector<Double_t> fQQQ3_PAD_TimeCFD;
  std::vector<Double_t> fQQQ3_PAD_TimeLED;
  std::vector<Double_t> fQQQ3_PAD_Time;

 public:
  TQQQ3Data();
  virtual ~TQQQ3Data();

  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

  /////////////////////           SETTERS           ////////////////////////
  inline void SetFront_DetectorNbr(const UShort_t& DetNbr) { fQQQ3_StripFront_DetectorNbr.push_back(DetNbr); }
  inline void SetFront_StripNbr(const UShort_t& StripNbr) { fQQQ3_StripFront_StripNbr.push_back(StripNbr); }
  inline void SetFront_Energy(const Double_t& Energy) { fQQQ3_StripFront_Energy.push_back(Energy); }
  inline void SetFront_TimeCFD(const Double_t& TimeCFD) { fQQQ3_StripFront_TimeCFD.push_back(TimeCFD); }
  inline void SetFront_TimeLED(const Double_t& TimeLED) { fQQQ3_StripFront_TimeLED.push_back(TimeLED); }
  inline void SetFront_Time(const Double_t& Time) { fQQQ3_StripFront_Time.push_back(Time); }

  inline void SetBack_DetectorNbr(const UShort_t& DetNbr) { fQQQ3_StripBack_DetectorNbr.push_back(DetNbr); }
  inline void SetBack_StripNbr(const UShort_t& StripNbr) { fQQQ3_StripBack_StripNbr.push_back(StripNbr); }
  inline void SetBack_Energy(const Double_t& Energy) { fQQQ3_StripBack_Energy.push_back(Energy); }
  inline void SetBack_TimeCFD(const Double_t& TimeCFD) { fQQQ3_StripBack_TimeCFD.push_back(TimeCFD); }
  inline void SetBack_TimeLED(const Double_t& TimeLED) { fQQQ3_StripBack_TimeLED.push_back(TimeLED); }
  inline void SetBack_Time(const Double_t& Time) { fQQQ3_StripBack_Time.push_back(Time); }

  inline void SetPAD_DetectorNbr(const UShort_t& DetNbr) { fQQQ3_PAD_DetectorNbr.push_back(DetNbr); }
  inline void SetPAD_Energy(const Double_t& Energy) { fQQQ3_PAD_Energy.push_back(Energy); }
  inline void SetPAD_TimeCFD(const Double_t& TimeCFD) { fQQQ3_PAD_TimeCFD.push_back(TimeCFD); }
  inline void SetPAD_TimeLED(const Double_t& TimeLED) { fQQQ3_PAD_TimeLED.push_back(TimeLED); }
  inline void SetPAD_Time(const Double_t& Time) { fQQQ3_PAD_Time.push_back(Time); }

  inline void SetFront(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy,
                       const Double_t& TimeCFD, const Double_t& TimeLED, const Double_t& Time = 0) {
    SetFront_DetectorNbr(DetNbr);
    SetFront_StripNbr(StripNbr);
    SetFront_Energy(Energy);
    SetFront_TimeCFD(TimeCFD);
    SetFront_TimeLED(TimeLED);
    SetFront_Time(Time);
  };
  inline void SetBack(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy, const Double_t& TimeCFD,
                      const Double_t& TimeLED, const Double_t& Time = 0) {
    SetBack_DetectorNbr(DetNbr);
    SetBack_StripNbr(StripNbr);
    SetBack_Energy(Energy);
    SetBack_TimeCFD(TimeCFD);
    SetBack_TimeLED(TimeLED);
    SetBack_Time(Time);
  };
  inline void SetPAD(const UShort_t& DetNbr, const Double_t& Energy, const Double_t& TimeCFD, const Double_t& TimeLED,
                     const Double_t& Time = 0) {
    SetPAD_DetectorNbr(DetNbr);
    SetPAD_Energy(Energy);
    SetPAD_TimeCFD(TimeCFD);
    SetPAD_TimeLED(TimeLED);
    SetPAD_Time(Time);
  };

  /////////////////////           GETTERS           ////////////////////////
  inline UShort_t GetFront_DetectorNbr(const unsigned int& i) const { return fQQQ3_StripFront_DetectorNbr[i]; } //!
  inline UShort_t GetFront_StripNbr(const unsigned int& i) const { return fQQQ3_StripFront_StripNbr[i]; }       //!
  inline Double_t GetFront_Energy(const unsigned int& i) const { return fQQQ3_StripFront_Energy[i]; }           //!
  inline Double_t GetFront_TimeCFD(const unsigned int& i) const { return fQQQ3_StripFront_TimeCFD[i]; }         //!
  inline Double_t GetFront_TimeLED(const unsigned int& i) const { return fQQQ3_StripFront_TimeLED[i]; }         //!
  inline Double_t GetFront_Time(const unsigned int& i) const { return fQQQ3_StripFront_Time[i]; }               //!

  inline UShort_t GetBack_DetectorNbr(const unsigned int& i) const { return fQQQ3_StripBack_DetectorNbr[i]; } //!
  inline UShort_t GetBack_StripNbr(const unsigned int& i) const { return fQQQ3_StripBack_StripNbr[i]; }       //!
  inline Double_t GetBack_Energy(const unsigned int& i) const { return fQQQ3_StripBack_Energy[i]; }           //!
  inline Double_t GetBack_TimeCFD(const unsigned int& i) const { return fQQQ3_StripBack_TimeCFD[i]; }         //!
  inline Double_t GetBack_TimeLED(const unsigned int& i) const { return fQQQ3_StripBack_TimeLED[i]; }         //!
  inline Double_t GetBack_Time(const unsigned int& i) const { return fQQQ3_StripBack_Time[i]; }               //!

  inline UShort_t GetPAD_DetectorNbr(const unsigned int& i) const { return fQQQ3_PAD_DetectorNbr[i]; } //!
  inline Double_t GetPAD_Energy(const unsigned int& i) const { return fQQQ3_PAD_Energy[i]; }           //!
  inline Double_t GetPAD_TimeCFD(const unsigned int& i) const { return fQQQ3_PAD_TimeCFD[i]; }         //!
  inline Double_t GetPAD_TimeLED(const unsigned int& i) const { return fQQQ3_PAD_TimeLED[i]; }         //!
  inline Double_t GetPAD_Time(const unsigned int& i) const { return fQQQ3_PAD_Time[i]; }               //!

  inline unsigned int GetMultiplicityFront() const { return fQQQ3_StripFront_DetectorNbr.size(); } //!
  inline unsigned int GetMultiplicityBack() const { return fQQQ3_StripBack_DetectorNbr.size(); }   //!
  inline unsigned int GetMultiplicityPAD() const { return fQQQ3_PAD_DetectorNbr.size(); }          //!

  ClassDef(TQQQ3Data, 1) // QQQ3Data structure
};

#endif
