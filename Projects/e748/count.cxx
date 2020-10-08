void count(){
  
  std::vector<int> run = {315,316,317,318,320,321,323,325,326,327,328,329,330,331,339,341,342,346,347,348 };
  unsigned int size = run.size();
  ofstream file("event.txt");
  int sum=0;
  for(unsigned int i = 0 ; i < size ; i++){
    
    TChain* chain = new TChain("AutoTree");
    chain->Add(Form("/data/Transfert/e748/root/run_0%i_*.root",run[i]));
    file << run[i] << " " << chain->GetEntries() << endl;
    sum+=chain->GetEntries();
    }
    file << "total " << sum << endl;
}
