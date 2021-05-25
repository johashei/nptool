void dist(){

  ifstream f("distance.txt");
  double d;
  auto h = new TH1D("h","h",1000,0,100);

  while(f>>d){
    h->Fill(d);
  }

  h->Draw();

}
