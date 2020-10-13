#ifndef NPXMLPARSER_H
#define NPXMLPARSER_H
#include<string>
#include<map>
#include<set>
#include<vector>
#include "TXMLEngine.h"
namespace NPL{
  /////////////////////
  class parameter{
    public:
      parameter();
      parameter(std::string name, std::string value){
        m_name=name;
        m_value=value;
        };

     ~parameter();

    private:
      std::string m_name;
      std::string m_value;
    
    public:
      std::string GetName()  const {return m_name;};
      std::string GetValue() const {return m_value;};

    public:
    bool operator<(const parameter p2){
      return this->m_name<p2.m_name;
      }
    friend bool operator<(const parameter p1,const parameter p2){
      return p1.m_name<p2.m_name;
      }
    friend bool operator==(const parameter p1,const parameter p2){
      return p1.m_name==p2.m_name;
      };
  };

  /////////////////////
  class Channel{
    public:
      Channel();
      ~Channel();
      Channel(int device,int geo,int ch){
        m_device=device;
        m_geo=geo;
        m_ch=ch;
       };

    private:
      int m_device;
      int m_geo;
      int m_ch;
    
    public:
      int norme() const {return (m_device*1000000+m_geo*1000+m_ch);} ;

    public:
      bool operator<(const Channel p2){
        return this->norme()<p2.norme();
        }
      
      friend bool operator<(const Channel p1,const Channel p2){
        return p1.norme()<p2.norme();
        }
      
      friend bool operator==(const Channel p1,const Channel p2){
       return p1.norme()==p2.norme();
       }
  };
  /////////////////////
  class block{
    public:   
      block();
      ~block();
    
    public:
      int AsInt(std::string name);
      double AsDouble(std::string name);
      std::string AsString(std::string name); 
      void AddParameter(parameter p){ m_parameters.insert(p);};

    private:
      std::string m_name;
      std::set<parameter> m_parameters;
  };
  
  /////////////////////
  class XmlParser{
    public:
      XmlParser(){};
      ~XmlParser(){};

    public:
      void CheckFile(std::string file, std::string DTD);
      void LoadFile(std::string file);
      void DisplayNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level);
      void LoadNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level);
    
    public: // access by channel id
      block* GetBlock(Channel c){
        auto it = m_Channels.find(c);
        if(it!=m_Channels.end())
          return it->second;
        else
          return NULL;
        };
      void SetBlock(const Channel& c,block* b){
        m_Channels[c]=b; 
        } 
 
    public:
      std::vector<block*> GetAllBlocksWithName(std::string name) ;
      std::vector<std::string> GetAllBlocksName();
    private:
      std::map<std::string,std::vector<block>> m_blocks;
      std::map<Channel,block*> m_Channels;
      
  };

}

#endif
