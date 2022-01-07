void transform_file()
{
  string input_filename = "VFTX_SOFTOFW.r3b";

  ifstream ifile;
  ifile.open(input_filename.c_str());

  double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;
  string endofline;
  int line=0;
  for(int i=0; i<28; i++){
    for(int j=0; j<2; j++){
      string output_filename = "VFTX_TOFW" + to_string(i+1) + "_PMT" + to_string(j+1) + ".cal";
      ofstream ofile;
      ofile.open(output_filename.c_str());

      ofile << "SofTofW_TOFW" << i+1 << "_PMT" << j+1 << "_TIME ";
      for(int p=0; p<100; p++){
        line++;
        ifile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> endofline;
        cout << v1 << " " << v2 << " " << v3 << " " << v4 << " " << v5 << " " << v6 << " " << v7 << " " << v8 << " " << v9 << " " << v10 << " " << endofline << endl;
        cout << line << endl;
        ofile << v1 << " " << v2 << " " << v3 << " " << v4 << " " << v5 << " " << v6 << " " << v7 << " " << v8 << " " << v9 << " " << v10 << " ";
      }
      ofile.close();
    }
  }


}
