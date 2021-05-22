void convert(){

 ifstream in("NEB_cal_T.txt");
 ofstream out("Nebula_T.txt");
 string buffer;
 int id;
 double val;
 getline(in,buffer); // ignore first line
 while(in >> id >> buffer){
    out << "NEBULA_ID" << id << "_T " << buffer  << endl;
 }
 out.close();
 in.close();
 in.open("cal_Y.txt");
 out.open("Nebula_Y.txt");

 getline(in,buffer); // ignore first line
 while(in >> id >> buffer){
    out << "NEBULA_ID" << id << "_Y " << buffer  << endl;
 }
 out.close();
 in.close();
}
