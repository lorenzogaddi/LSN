#include "algoritmo_genetico.h"
#include <mpi.h>

//funzione per inizializzare il generatore di numeri casuali
void setup(int rank){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   for(int i = 0; i < rank+1; i++){
      if (Primes.is_open()){
         Primes >> p1 >> p2 ;
      } else cerr << "PROBLEM: Unable to open Primes" << endl;
   }
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

//funzione che evita di avere duplicati in un vettore
vector<int> no_duplicati(vector<int> vec) {
    vector<int> risultato(vec.size());
    unordered_set<int> set;
    for (int i = 0; i < vec.size(); i++) {
        int numero = vec[i];
        if (set.count(numero) == 0) {
            // Il numero non è un duplicato, lo inserisco nel set e nel vettore risultante
            set.insert(numero);
            risultato[i] = numero;
        } else {
            // Il numero è un duplicato, cerco un numero non presente nel set
            int prossimo_numero = 0;
            while (set.count(prossimo_numero) > 0) {
                prossimo_numero++;
            }
            // Inserisco il numero nel set e nel vettore risultante
            set.insert(prossimo_numero);
            risultato[i] = prossimo_numero;
        }
    }
    return risultato;
}

int main(int argc, char* argv[]){
   MPI_Init(&argc, &argv);
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   //coordinate x città
   const vector<double> x = 
      {-86.300568, -112.096962, -92.288986, -121.493629, -104.984856, -72.682198,
      -75.519722, -84.281296, -84.388229, -116.199722, -89.654961, -86.162643, 
      -93.603729, -95.677956, -84.875374, -91.187393, -69.781693, -76.490936, 
      -71.063698, -84.555328, -93.102211, -90.182106, -92.172935, -112.018417, 
      -96.699654, -119.766121, -71.537994, -74.769913, -105.939728, -78.639099, 
      -100.783318, -73.757874, -82.999069, -97.503342, -123.030403, -76.883598, 
      -71.414963, -81.033211, -100.346405, -86.784241, -97.740349, -111.888237, 
      -72.580536, -77.43364, -122.905014, -81.612328, -89.384445, -104.820236,
      -157.857376, -134.420212};
   //coordinate y città
   const vector<double> y = {32.377716, 33.448143, 34.746613, 38.576668, 39.739227,
   41.764046, 39.157307, 30.438118, 33.749027, 43.617775, 39.798363, 39.768623, 
   41.591087, 39.048191, 38.186722, 30.457069, 44.307167, 38.978764, 42.358162, 
   42.733635, 44.955097, 32.303848, 38.579201, 46.585709, 40.808075, 39.163914, 
   43.206898, 40.220596, 35.68224, 35.78043, 46.82085, 42.652843, 39.961346, 
   35.492207, 44.938461, 40.264378, 41.830914, 34.000343, 44.367031, 36.16581, 
   30.27467, 40.777477, 44.262436, 37.538857, 47.035805, 38.336246, 43.074684, 
   41.140259, 21.307442, 58.301598};

   //p1 e p2 del generatore di numeri casuali dipende dal rank: in questo
   //modo per ciascun rank viene utilizzato un generatore di numeri casuali
   //indipendente dagli altri rank
   setup(rank);
   
   //DNA iniziale ordinato
   auto DNA_iniziale = DNA(x, y);
   DNA_iniziale.geni_iniziali();

   //dati per inizializzare la classe algoritmo_genetico
   const size_t m = 200;
   const double percentuale_candidati = 50.0;
   const double prob_ricombinazione = 0.9;  
   const double prob_mutazione = 0.02;
   int N_migr = 2;
   //se p=0, le ricerche per la soluzione sono indipendenti, se p=1, ogni N_migr
   //i continenti si scambiano le sequenze migliori
   int p = 0;       

   //inizializzazione algoritmo_genetico
   auto ga = unique_ptr<algoritmo_genetico<DNA>>(new algoritmo_genetico<DNA>(
            &DNA_iniziale, m, percentuale_candidati, prob_ricombinazione,
            prob_mutazione));
   cout << "Rank: " << rank << "   Generazione: " << ga->get_generazione() << "  Fitness migliore: " << ga->get_fitness_migliore_DNA() << endl;

   //generazioni successive
   for(int i = 1; i < 10; i++){
         ga->evoluzione();
         //scambio delle migliori sequenze tra diversi ranks
         if(ga->get_generazione()%N_migr == 0){
            //DNA migliore
            DNA* dna = ga->get_DNA(ga->get_indice_migliore());
            //sequenza migliore
            vector<size_t> vec = dna->get_geni();
            //inserimento di un elemento casuale che determina lo scambio delle
            //sequenze migliori
            vec.push_back(floor(rnd.Rannyu()*size));
            vector<size_t> dati_condivisi((x.size()+1)*size);
            int MPI_Barrier(MPI_Comm);
            //condivide le migliori sequenze di tutti i processi tra i ranks
            MPI_Allgather(vec.data(), x.size()+1, MPI_UNSIGNED_LONG, dati_condivisi.data(), 
                        x.size()+1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
            //vettore che determina lo scambio delle migliori sequenze
            vector<int> sequenza(size);
            for(int j = 0; j < size; j++){
               sequenza[j] = dati_condivisi[(j+1)*x.size() + j];
            }
            //evito che ci siano duplicati nel vettore sequenza
            sequenza = no_duplicati(sequenza);
            //n determina da quale rank è copiata la nuova migliore sequenza
            int n = 0;
            for(int j = 0; j < size; j++){
               if(rank==sequenza[j]){
                  if(j==(size-1)){
                     n = sequenza[0];
                  } else{
                     n = sequenza[j+1];
                  }
               }
            }
            //sostituzione sequenza migliore
            vector<size_t> nuovo(x.size());
            for(int j = 0; j < x.size(); j++){
               nuovo[j] = dati_condivisi[n*x.size() + j + n];
            }
            //impostazione delle nuove migliori sequenze
            if(p==1){
               auto dna2 = ga->get_DNA(ga->get_indice_migliore());
               dna2->set_sequenza(nuovo);
               dna2->calcola_fitness();
               ga->fitness_popolazione();
               int MPI_Barrier(MPI_Comm);
            }
         }
         cout << "Rank: " << rank << "   Generazione: " << ga->get_generazione() << "  Fitness migliore: " << ga->get_fitness_migliore_DNA() << endl;
   }
   //migliore sequenza
   string str = ga->get_migliore_DNA();
   ofstream sequenza;
   sequenza.open("best_sequence" + to_string(rank) + ".out");
   sequenza << str << endl;
   sequenza << ga->get_fitness_migliore_DNA();
   sequenza.close();


   MPI_Finalize();
   return 0;
}
