#include"NPDetectorFactory.h"
using namespace NPA;
#include<dlfcn.h>
#include<fstream>
#include<iostream>
#include<stdlib.h>

DetectorFactory* DetectorFactory::m_Instance = 0;
#ifdef __linux__
std::string SHARED_LIB_EXTENSION = ".so";
#endif

#ifdef __FreeBSD__
std::string SHARED_LIB_EXTENSION = ".so";
#endif

#ifdef __APPLE__
std::string SHARED_LIB_EXTENSION = ".dylib";
#endif



////////////////////////////////////////////////////////////////////////////////
DetectorFactory::DetectorFactory(){

}

/////////////////////////////////////////////////////////////////////////////////
DetectorFactory::~DetectorFactory(){
  if(m_Instance!=0)
    delete m_Instance;

  m_Instance = 0 ;
}

////////////////////////////////////////////////////////////////////////////////
DetectorFactory* DetectorFactory::getInstance(){
  if(!m_Instance){
    m_Instance = new DetectorFactory();
  }

  return m_Instance;

}
////////////////////////////////////////////////////////////////////////////////
void DetectorFactory::ReadClassList(std::string FileList){
    std::ifstream InFile(FileList.c_str());

    if(!InFile.is_open()){
      cout << "ERROR: Detector Class List file " << FileList << " Not found" << endl;
      exit(1);
    }

    string Token, LibName;
    while(InFile >> Token >> LibName)
      m_TokenLib[Token] = LibName; 
}
////////////////////////////////////////////////////////////////////////////////
void DetectorFactory::CreateClassList(std::string FileList){
  ofstream outFile(FileList.c_str());
  std::map<string,string>::iterator it;

  for(it = m_TokenLib.begin();it!=m_TokenLib.end();it++){
    outFile << it->first << " " << it->second << std::endl;
  
  }


}
////////////////////////////////////////////////////////////////////////////////
NPA::VDetector* DetectorFactory::Construct(std::string Token){
  std::map<std::string,ClassDetectorFactoryFn>::iterator it;
  
  if(m_Construct.find(Token)!=m_Construct.end())
    return  m_Construct[Token]();

  else if( m_TokenLib.find(Token)!=m_TokenLib.end()){
    // Add absolute path to library name
    std::string path = getenv("NPTOOL"); 
    std::string libName = path+"NPLib/lib/"+m_TokenLib[Token];
    dlopen(libName.c_str(),RTLD_NOW);
    
    if(m_Construct.find(Token)!=m_Construct.end())
      return  m_Construct[Token]();
  
    else{
      std::cout << "Warning: Detector with Token " << Token << " has no Constructor or no Library" << std::endl;
      return NULL;
    }
  }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////
void DetectorFactory::AddDetector(std::string Token, ClassDetectorFactoryFn fn){
  m_Construct[Token] = fn;
}
////////////////////////////////////////////////////////////////////////////////
void DetectorFactory::AddToken(std::string Token, std::string LibName){
  // Make sure the lib name is correct:
  // Strip it from lib and NP and extension
  std::string remword = "lib";
  size_t pos;
  if((pos=LibName.find(remword))!=std::string::npos){
    LibName.swap(LibName.erase(pos,remword.length()));
  }

  remword = "NP";
  if((pos=LibName.find(remword))!=std::string::npos){
    LibName.swap(LibName.erase(pos,remword.length()));
  }
 
  remword = ".so";
  if((pos=LibName.find(remword))!=std::string::npos){
    LibName.swap(LibName.erase(pos,remword.length()));
  }

  remword = ".dylib";
  if((pos=LibName.find(remword))!=std::string::npos){
    LibName.swap(LibName.erase(pos,remword.length()));
  }


  // Create the full name:
  LibName = "libNP"+LibName + SHARED_LIB_EXTENSION;
  
  m_TokenLib[Token] = LibName;
}