#include "esercizio2.h"

//questo programma fornisce la stima migliore dei parametri mu e sigma che
//minimizzano l'energia e fa il grafico di <H> in funzione dello step del
//simulated annealing

int main()
{ 
  TApplication app("app", 0, 0);
  TGraphErrors *graph = new TGraphErrors();
  Input(); //inizializzazione
  equilibrazione();
  //inizializzo l'energia e il suo errore
  vector<double> vec = Energia();
  energy = vec[0];
  errore = vec[1];
  cout << "Energia iniziale: " << energy << endl;
  //simulated annealing: è fatto un ciclo per diverse temperature diminuendola
  //gradualmente e all'interno di questo ciclo è fatto un ulteriore ciclo
  //muovendo con la funzione Move2() l'energia locale e il suo errore
  int contatore = 0;
  for(int j = 0; j < 9; j++){
    double T = 0.05;
    T = T - (j+1)*0.0049;
    for(int i = 0; i < 20; i++){
      Move2(T);
      graph->SetPoint(contatore, (contatore+1), energy);
      graph->SetPointError(contatore, 0.0, errore);
      contatore += 1;
    }
  }
  cout << "Miglior stima energia: " << best_E << endl;
  cout << "Errore migliore stima energia: " << best_err << endl;
  cout << "Miglior stima mu: " << best_mu << endl;
  cout << "Miglior stima sigma: " << best_sigma << endl;

  TCanvas SA("SA", "SA");
  graph->SetTitle("<H>");
  graph->GetXaxis()->SetTitle("step");
  graph->GetYaxis()->SetTitle("<H>");
  graph->Draw("ALP");
  app.Run();
  return 0;
}
