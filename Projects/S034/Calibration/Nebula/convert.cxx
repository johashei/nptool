void convert(){

 ifstream in("NEB_cal_T.txt");
 ofstream out("Nebula_T.txt");
 string buffer;
 int id;
 double val;
 getline(in,buffer); // ignore first line
 while(in >> id >> buffer){
    out << "NEBULA_T_ID" << id << " " << buffer  << endl;
 }
 out.close();
 in.close();
 in.open("cal_Y.txt");
 out.open("Nebula_Y.txt");

 getline(in,buffer); // ignore first line
 while(in >> id >> buffer){
    out << "NEBULA_Y_ID" << id << " " << buffer  << endl;
 }
 out.close();
 in.close();
}
