void transform_file(int det=1, int signal=1)
{
  string input_filename  = "VFTX_DET" + to_string(det) + "_SIGNAL" + to_string(signal) + ".r3b";
  string output_filename = "VFTX_DET" + to_string(det) + "_SIGNAL" + to_string(signal) + ".cal";

  ifstream ifile;
  ifile.open(input_filename.c_str());
  ofstream ofile;
  ofile.open(output_filename.c_str());

  string buffer;
  getline(ifile,buffer);
  ofile << buffer << " ";

  double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;
  string endofline;
  //while(!ifile.eof()){
  for(int i=0; i<100; i++){
    ifile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> endofline;
    cout << v1 << " " << v2 << " " << v3 << " " << v4 << " " << v5 << " " << v6 << " " << v7 << " " << v8 << " " << v9 << " " << v10 << " " << endofline << endl;
    ofile << v1 << " " << v2 << " " << v3 << " " << v4 << " " << v5 << " " << v6 << " " << v7 << " " << v8 << " " << v9 << " " << v10 << " ";
  }
}
