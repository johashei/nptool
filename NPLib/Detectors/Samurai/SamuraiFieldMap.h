#ifndef SamuraiFieldMap_h
#define SamuraiFieldMap_h
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : May  2021                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Samurai field map data                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include"TObject.h"
#include<string>
#include<vector>
#include<map>

class SamuraiFieldMap{

  public:
    SamuraiFieldMap(){};
    SamuraiFieldMap(std::string file);
    ~SamuraiFieldMap(){};
  
  public: // Map reading
    void LoadMap(std::string file);
    void LoadAscii(std::string file);
    void LoadBinary(std::string file);

  private:
    // map[Pos]=B;
    std::map<std::vector<float>,std::vector<float>> m_field;
    float m_x_max,m_y_max,m_z_max,m_x_min,m_y_min,m_z_min;

  public:
    std::vector<float>& GetB(std::vector<float>& pos);
    inline std::vector<float>& GetB(float x,float y ,float z){
      std::vector<float> pos = {x,y,z};
      return GetB(pos);
    };

    ClassDef(SamuraiFieldMap,1);
};

#endif
