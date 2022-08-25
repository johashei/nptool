void GenerateCalibFile(){
  string output_filename = "Vendeta_Time.cal";
  ofstream ofile;
  ofile.open(output_filename.c_str());

  int det_number=72;
  int anode_number=11;
  double time_offset=0;
  for(int i=1;i<det_number+1;i++){
    for(int j=1;j<anode_number+1;j++){
      TString token_lg = Form("Vendeta_DET%i_LG_ANODE%i_TIMEOFFSET",i,j);
      TString token_hg = Form("Vendeta_DET%i_HG_ANODE%i_TIMEOFFSET",i,j);

      ofile << token_lg << " " << time_offset << endl;
      ofile << token_hg << " " << time_offset << endl;
    }
  }
  ofile.close();
}
