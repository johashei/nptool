/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 21/07/09                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: This class is a singleton class which deals with the ROOT     *
 *             input file and tree both for NPSimulation and NPAnalysis.     *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <sys/stat.h>
#include <stdlib.h>
#include <dlfcn.h>

#include "RootInput.h"
#include "TAsciiFile.h"
#include "NPOptionManager.h"
using namespace std;

RootInput* RootInput::instance = 0;
////////////////////////////////////////////////////////////////////////////////
RootInput* RootInput::getInstance(std::string configFileName){
  // A new instance of RootInput is created if it does not exist:
  if (instance == 0) {
    instance = new RootInput(configFileName);
  }

  // The instance of RootInput is returned:
  return instance;
}

////////////////////////////////////////////////////////////////////////////////
void RootInput::Destroy(){
  if (instance != 0) {
    delete instance;
    instance = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void RootInput::ReadOldStyleInputFile(NPL::InputParser& parser){
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("TTreeName");
    pTreeName = blocks[0]->GetLines()[0];
    vector<NPL::InputBlock*> files = parser.GetAllBlocksWithToken("RootFileName");
    if(files.size()>0){
      vector<string> lines=files[0]->GetLines();
      unsigned int size = lines.size();
      for(unsigned int i = 0 ; i < size ; i++){
       pTreePath.push_back(lines[i]); 
      }
    }
}

////////////////////////////////////////////////////////////////////////////////
void RootInput::ReadInputFile(NPL::InputParser& parser){
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Tree");
    pTreeName = blocks[0]->GetMainValue();
    std::vector<std::string> lines=blocks[0]->GetLines();
    unsigned int size = lines.size();
    for(unsigned int i = 0 ; i < size ; i++){
        if(lines[i].find(".root")!=string::npos)
          pTreePath.push_back(lines[i]); 
        else if(lines[i].find(".tree")!=string::npos)
          ReadTreeFile(lines[i]);
    }

    vector<NPL::InputBlock*> friends = parser.GetAllBlocksWithToken("Friend");
    unsigned int sizeF = friends.size();
    for(unsigned int i = 0 ; i < sizeF ; i++){
      pFriends.insert(pair< string,vector<string> > (friends[i]->GetMainValue(),friends[i]->GetLines()));
    }
}

////////////////////////////////////////////////////////////////////////////////
void  RootInput::ReadTreeFile(std::string path){
  ifstream tree(path.c_str());
  path=path.substr(0,path.rfind("/")+1);
  path+="/";
  std::string buffer;
  bool first=true;
  while(tree>>buffer){
    if(first){
      pTreePath.push_back(path+buffer);
      first=false;
    }
    else{
      vector<string> friends={path+buffer};
      pFriends.insert(pair<string,vector<string>>(pTreeName, friends));
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
RootInput::RootInput(std::string configFileName){

  NumberOfFriend = 0;
  pRootFile  = NULL;
  std::string lastfile= NPOptionManager::getInstance()->GetLastFile();
  if(lastfile!="VOID"){
    configFileName = lastfile;
  }

  std::cout << "/////////// ROOT Input files ///////////" << std::endl;
  std::cout << "Initializing input TChain using : " << configFileName << std::endl;
  NPL::InputParser parser(configFileName);

  // Old style file
  vector<NPL::InputBlock*> old_blocks = parser.GetAllBlocksWithToken("TTreeName");
  if(old_blocks.size()==1){
    ReadOldStyleInputFile(parser);
  }
  // New style file
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Tree");
  if(blocks.size()==1){
    ReadInputFile(parser);
  }
  // If the tree come from a simulation, the InteractionCoordinates
  // and InitialConditions lib are loaded
  if(pTreeName=="SimulatedTree"){
   std::string path = getenv("NPTOOL");
   path+="/NPLib/lib/";
   std::string libName="libNPInteractionCoordinates"+NPOptionManager::getInstance()->GetSharedLibExtension();
   libName=path+libName;
   dlopen(libName.c_str(),RTLD_NOW);
   libName="libNPInitialConditions"+NPOptionManager::getInstance()->GetSharedLibExtension();
   libName=path+libName;
   dlopen(libName.c_str(),RTLD_NOW);
  }

  // Initialise the chain
  pRootChain = new TChain(pTreeName.c_str());

  // Add all the files
  unsigned int size = pTreePath.size();
  std::string firstfile;
  for(unsigned int i = 0 ; i < size ; i++){
    cout << "  - Adding file : " << pTreePath[i].c_str() << endl;
    pRootChain->Add(pTreePath[i].c_str());
    // Test if the file is a regex or a single file
    double counts;
    std::string command = "ls " + pTreePath[i] + " > .ls_return";
    counts= system(command.c_str());
    std::ifstream return_ls(".ls_return");
    std::string files;
    while(return_ls >> files){
      if(counts == 0)
        firstfile = files;
      counts++;
    }
  }
  
  // Add all the friends
  for(auto it = pFriends.begin(); it!=pFriends.end() ; it++){
    unsigned int size = it->second.size();
    for(unsigned int i = 0 ; i < size ; i++){
     cout << "  - Adding friend : " << it->second[i].c_str() << endl;
     pRootChain->AddFriend(it->first.c_str(),it->second[i].c_str());
    }
  }

  if (!pRootFile) 
    pRootFile = new TFile(firstfile.c_str());

  if( pRootChain->GetEntries() ==0){
    std::cout << "\033[1;31m**** ERROR: No entries to analyse ****\033[0m" << std::endl; 
    exit(1);
  }
 else
    std::cout << "\033[1;32mROOTInput:  " << pRootChain->GetEntries() << " entries loaded in the input chain\033[0m" << std::endl ;
}

////////////////////////////////////////////////////////////////////////////////
void RootInput::AddFriendChain(std::string RunToAdd){
  NumberOfFriend++;
  std::ostringstream suffix_buffer;
  suffix_buffer << "_" << NumberOfFriend ; 
  std::string suffix = suffix_buffer.str();
  bool CheckTreeName     = false;
  bool CheckRootFileName = false;

  // Read configuration file Buffer
  std::string lineBuffer, dataBuffer;

  // Open file
  std::ifstream inputConfigFile;
  inputConfigFile.open(RunToAdd.c_str());

  TChain* localChain = new TChain();
  
  std::cout << "/////////////////////////////////" << std::endl;
  std::cout << "Adding friend to current TChain" << std::endl;

  if (!inputConfigFile) {
    std::cout << "Run to Add file :" << RunToAdd << " not found " << std::endl; 
    return;
  }

  else {
    while (!inputConfigFile.eof()) {
      getline(inputConfigFile, lineBuffer);

      // search for token giving the TTree name
      if (lineBuffer.compare(0, 9, "TTreeName") == 0) {
        inputConfigFile >> dataBuffer;
        // adding suffix to insure uniquity of the chain name
        dataBuffer+suffix;
        // initialize localChain
        localChain->SetName(dataBuffer.c_str());
        CheckTreeName = true ;
      }

      // search for token giving the list of Root files to treat
      else if (lineBuffer.compare(0, 12, "RootFileName") == 0  &&  localChain) {
        CheckRootFileName = true ;

        while (!inputConfigFile.eof()) {
          inputConfigFile >> dataBuffer;

          // ignore comment Line 
          if (dataBuffer.compare(0, 1, "%") == 0) {
            inputConfigFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
          }

          else if (!inputConfigFile.eof()) {
            localChain->Add(dataBuffer.c_str());
            std::cout << "Adding file " << dataBuffer << " to TChain" << std::endl;
          }
        }
      }
    }
  }

  if (!CheckRootFileName || !CheckTreeName) 
    std::cout << "WARNING: Token not found for InputTree Declaration : Input Tree has not be Added to the current Chain" << std::endl;

  else
    pRootChain->AddFriend( localChain->GetName() );

  std::cout << "/////////////////////////////////" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
std::string RootInput::DumpAsciiFile(const char* type, const char* folder){
  std::string name="fail";

  std::string sfolder = folder;
  // create folder if not existing
  // first test if folder exist
  struct stat dirInfo;
  int res = stat(folder, &dirInfo);
  if (res != 0) {   // if folder does not exist, create it
    std::string cmd = "mkdir " + sfolder;
    if (system(cmd.c_str()) != 0) std::cout << "RootInput::DumpAsciiFile() problem creating directory" << std::endl;
  }

  std::string stype = type;
  if (stype == "EventGenerator") {
    TAsciiFile *aFile = (TAsciiFile*)pRootFile->Get(stype.c_str());

    if(aFile)
    {
      // build file name
      std::string title = aFile->GetTitle();
      size_t pos = title.rfind("/");
      if (pos != std::string::npos) name = sfolder + title.substr(pos);
      else name = sfolder + "/" + title;
      aFile->WriteToFile(name.c_str());
    }
  }

  else if (stype == "DetectorConfiguration") {
    TAsciiFile *aFile = (TAsciiFile*)pRootFile->Get(stype.c_str());
    if(aFile)
    {
      // build file name
      std::string title = aFile->GetTitle();
      size_t pos = title.rfind("/");
      if (pos != std::string::npos) name = sfolder + title.substr(pos);
      else name = sfolder + "/" + title;
      aFile->WriteToFile(name.c_str());
    }
  }

  else if (stype == "Calibration") {
    TAsciiFile *aFile = (TAsciiFile*)pRootFile->Get(stype.c_str());
    if(aFile)
    {
      // build file name
      std::string title = aFile->GetTitle();
      size_t pos = title.rfind("/");
      if (pos != std::string::npos) name = sfolder + title.substr(pos);
      else name = sfolder + "/" + title;
      aFile->WriteToFile(name.c_str());
    }
  }

  else if (stype == "RunToTreat") {
  }
  else {
    std::cout << "RootInput::DumpAsciiFile() unkwown keyword" << std::endl;
  }

  return name;
}

////////////////////////////////////////////////////////////////////////////////
RootInput::~RootInput(){
  // delete default directory ./.tmp
  struct stat dirInfo;
  int res = stat("./.tmp", &dirInfo);
  if (res == 0) {   // if does exist, delete it
    if (system("rm -rf ./.tmp") != 0) std::cout << "RootInput::~RootInput() problem deleting ./.tmp directory" << std::endl; 
  }
  std::cout << std::endl << "Root Input summary" << std::endl;
  std::cout << "  - Number of bites read: " << pRootFile->GetBytesRead() << std::endl;
  std::cout << "  - Number of transactions: " << pRootFile->GetReadCalls() << std::endl;
  // Close the Root file
  pRootFile->Close();
}

////////////////////////////////////////////////////////////////////////////////
TChain* MakeFriendTrees(std::string RunToRead1,std::string RunToRead2){
  RootInput:: getInstance(RunToRead1)	;
  RootInput:: getInstance()->AddFriendChain(RunToRead2);
  return RootInput:: getInstance()->GetChain();
}
