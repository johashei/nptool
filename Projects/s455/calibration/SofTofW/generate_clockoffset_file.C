void generate_clockoffset_file()
{

  string output_filename = "ClockOffset.cal";
  ofstream ofile;
  ofile.open(output_filename.c_str());

  for(int i=0; i<28; i++){
    for(int j=0; j<2; j++){
      ofile << "SofTofW_TOFW" + to_string(i+1) + "_PMT" + to_string(j+1) + "_CLOCKOFFSET" << " " << 0 << endl;
    }
  }
  ofile.close();

}
