#include "algoritmo_genetico.h"

//funzione che inizializza il generatore di numeri casuali
void setup(){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

int main(){
   setup();
   
   //coordinate delle città disposte su una circonferenza
   
   const vector<double> x = {
   1, 0.982973, 0.932472, 0.850217, 0.739009, 0.602635, 0.445738, 0.273663,
   0.0922684, -0.0922684, -0.273663, -0.445738, -0.602635, -0.739009,
   -0.850217, -0.932472, -0.982973, -1, -0.982973, -0.932472, -0.850217,
   -0.739009, -0.602635, -0.445738, -0.273663, -0.0922684, 0.0922684,
   0.273663, 0.445738, 0.602635, 0.739009, 0.850217, 0.932472, 0.982973
   };

   const vector<double> y = {
   0, 0.18375, 0.361242, 0.526432, 0.673696, 0.798017, 0.895163, 0.961826,
   0.995734, 0.995734, 0.961826, 0.895163, 0.798017, 0.673696, 0.526432,
   0.361242, 0.18375, 1.22465e-16, -0.18375, -0.361242, -0.526432,
   -0.673696, -0.798017, -0.895163, -0.961826, -0.995734, -0.995734,
   -0.961826, -0.895163, -0.798017, -0.673696, -0.526432, -0.361242,
   -0.18375
   }; 

   //coordinate città disposte all'interno di un quadrato 2x2 con vertici in
   //(0,0), (0,2), (2,2) e (2,0) 
   /*
   const vector<double> x = {
   0.0555556, 0.111111, 0.166667, 0.222222, 0.277778, 0.333333, 0.388889,
   0.444444, 0.5, 0.555556, 0.611111, 0.666667, 0.722222, 0.777778, 0.833333,
   0.888889, 0.944444, 1, 1.05556, 1.11111, 1.16667, 1.22222, 1.27778, 1.33333,
   1.38889, 1.44444, 1.5, 1.55556, 1.61111, 1.66667, 1.72222, 1.77778, 1.83333, 
   1.88889};

   const vector<double> y = {
   0.333034, 0.422844, 1.48633, 0.729215, 1.00652, 1.18649, 1.43956, 0.263414,
   1.34644, 1.82042, 0.805843, 1.68948, 1.57449, 0.140161, 0.419817, 1.18763,
   0.630543, 1.13713, 0.132873, 1.20804, 0.222821, 1.53384, 0.542998, 0.467257,
   1.03713, 1.53649, 1.01036, 1.03909, 0.589998, 0.55461, 0.996414, 1.21134, 
   1.27194, 1.20289}; */

   //DNA iniziale ordinato
   auto DNA_iniziale = DNA(x, y);
   DNA_iniziale.geni_iniziali();

   //dati per inizializzare la classe algoritmo_genetico
   const size_t m = 200;
   const double percentuale_candidati = 50.0;
   const double prob_ricombinazione = 0.9;  
   const double prob_mutazione = 0.02;       

   //inizializzazione algoritmo_genetico
   auto ga = unique_ptr<algoritmo_genetico<DNA>>(new algoritmo_genetico<DNA>(
            &DNA_iniziale, m, percentuale_candidati, prob_ricombinazione,
            prob_mutazione));

   cout << "Generazione: " << ga->get_generazione();
   cout << "  Fitness migliore: " << ga->get_fitness_migliore_DNA() << endl;
   ga->print();

   //generazioni successive
   for(int i = 1; i < 15; i++){
      ga->evoluzione();
      ga->print();
      cout << "Generazione: " << ga->get_generazione();
      cout << "  Fitness migliore: " << ga->get_fitness_migliore_DNA() << endl;
   }

   //migliore sequenza
   string str = ga->get_migliore_DNA();
   cout << "Migliore sequenza: " << str << endl;

   return 0;
}
