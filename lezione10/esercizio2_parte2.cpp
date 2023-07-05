#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <numeric>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TH1F.h"

using namespace std;

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

//calcolo della funzione costo per una determinata sequenza di città
double funzione_costo(vector<int> vec, int n) {
  double somma = 0;
  for (int i = 0; i < n - 1; i++) {
    int indice1 = vec[i];
    int indice2 = vec[i+1];
    double dx = x[indice1] - x[indice2];
    double dy = y[indice1] - y[indice2];
    double d = sqrt(dx*dx + dy*dy);
    somma += d;
  }
  double dx = x[0] - x[n-1];
  double dy = y[0] - y[n-1];
  somma += sqrt(dx*dx + dy*dy);
  return somma;
}


//questa funzione trova in quale file è presente la migliore sequenza
int trova(int rank, int n){
    //riempio un vettore con le migliori distanze per ciascun rank
    int a = 0;
    vector<double> distanza(rank);
    for(int j = 0; j < rank; j++){
        ifstream sequence;
        sequence.open("best_sequence" + to_string(j) + ".out");
        for(int i = 0; i < n; i++){
            sequence >> a;
        }
        sequence >> distanza[j];
        sequence.close();
    }
    //trova il file la cui sequenza minimizza la distanza
    double minVal = distanza[0];
    int minIndex = 0;
    for (int i = 1; i < distanza.size(); ++i) {
        if (distanza[i] < minVal) {
            minVal = distanza[i];
            minIndex = i;
        }
    }
    return minIndex;
}

//Questo programma mostra il miglior percorso trovato dal programma
//esercizio2_parte2.exe
int main(){
    TApplication app("app", 0, 0);
    ifstream path, sequence;
    int n = 50;
    int size = 6; //numero di rank utilizzati
    int a = 0;
    double L = 0;
    int best_rank = trova(size, n);
    cout << "Rank migliore: " << best_rank << endl;
    sequence.open("best_sequence" + to_string(best_rank) + ".out");

    //grafico del sentiero percorso per la miglior sequenza
    vector<int> sequ;
    int best_generazione = 0;
    double b = 0;
    for(int i = 0; i < n; i++){
        sequence >> a;
        sequ.push_back(a);
    }
    TGraph* graph = new TGraph();

    for(int i = 0; i < n; i++){
        a = sequ[i];
        graph->SetPoint(i, x[a], y[a]);
    }
    graph->SetPoint(n, x[sequ[0]], y[sequ[0]]);

    TCanvas percorso("percorso", "percorso");

    graph->SetTitle("best_path");
    graph->GetXaxis()->SetTitle("x");
    graph->GetYaxis()->SetTitle("y");
    graph->Draw("ALP*");
    sequence.close();

    app.Run();
    return 0;
}