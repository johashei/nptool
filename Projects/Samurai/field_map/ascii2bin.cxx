void ascii2bin(){
  std::string file ="3T.table";
  ifstream in(file.c_str());
  if(!in.is_open()){
    cout << "Error: failed to load samurai field map " << file << endl;
    exit(1);
  }

  cout << "//////// Loading Samurai Field Map " << file << endl; 
  double x,y,z,Bx,By,Bz;
  std::string name = file+".bin";
  ofstream out(name.c_str(),std::ofstream::binary);

  unsigned int  count =0 ;

  // ignore 8 first line 
  string buffer;
  for(unsigned int i = 0 ; i < 8 ; i++){
    getline(in,buffer);
  }

  while(!in.eof()){
    if(++count%10000==0)
      cout << "\r  - Loading " << count << " values " << flush;

    if(in >> x >> y >> z >> Bx >> By >> Bz){
      out.write((char*)&x,sizeof(x));
      out.write((char*)&y,sizeof(y));
      out.write((char*)&z,sizeof(z));
      out.write((char*)&Bx,sizeof(Bx));
      out.write((char*)&By,sizeof(By));
      out.write((char*)&Bz,sizeof(Bz));
    }
    else{
    cout << x <<" " <<  y << " "<< z << " "<< Bx <<" "<< By << " " << Bz << endl;;
    if(!in.eof())
      in.clear();
    }
    //     cout << x << endl;
  }

  cout << "\r  - " << count << " values loaded" << endl; 

  out.close();
  in.close();
}
