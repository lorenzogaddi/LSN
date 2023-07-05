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

//Questo programma fa i grafici del fitness medio della miglior metà della 
//popolazione e del fitness migliore in funzione del numero di generazioni e del 
//percorso della sequenza migliore. Sostituendo square con circle nei nomi dei 
//file, si ottengono i grafici per 34 città disposte all'interno di un quadrato
int main(){
    TApplication app("app", 0, 0);
    ifstream fitness_medio, fitness_migliore;
    fitness_medio.open("average_fitness.dat");
    fitness_migliore.open("best_fitness.dat");
    int q = 0;
    int generazioni = 15;
    int n = 34;
    int a = 0;
    double fitness = 0; 
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
    -0.18375}; 

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

    //grafico fitness medio della miglior metà della popolazione in funzione del
    //numero di generazioni
    TGraph* graph1 = new TGraph();

    for(int i = 0; i < generazioni; i++){
        fitness_medio >> fitness;
        graph1->SetPoint(i, (i+1), fitness);
    }

    TCanvas average("<L>", "<L>");

    graph1->SetTitle("<L>");
    graph1->GetXaxis()->SetTitle("Generazione");
    graph1->GetYaxis()->SetTitle("<L>");
    graph1->Draw("ALP");

    //grafico del fitness migliore in funzione del numero di generazioni
    TGraph* graph2 = new TGraph();
    for(int i = 0; i < generazioni; i++){
        fitness_migliore >> fitness;
        graph2->SetPoint(i, (i+1), fitness);
    }

    TCanvas best("L_best", "L_best");

    graph2->SetTitle("L_best");
    graph2->GetXaxis()->SetTitle("Generazione");
    graph2->GetYaxis()->SetTitle("L_best");
    graph2->Draw("ALP");

    //grafico del sentiero percorso per la miglior sequenza
    //sequenza per città disposte su una circonferenza
    vector<int> sequ = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 0, 1, 2, 3, 4, 5};
    //sequenza per città disposte all'interno di un quadrato
    //vector<int> sequ = {2, 9, 11, 12, 15, 17, 19, 21, 25, 32, 33, 31, 30, 27,
    //26, 24, 28, 29, 23, 22, 20, 18, 13, 14, 16, 10, 7, 0, 1, 3, 4, 5, 8, 6};

    TGraph* graph3 = new TGraph();

    for(int i = 0; i < n; i++){
        a = sequ[i];
        graph3->SetPoint(i, x[a], y[a]);
    }
    a = sequ[0];
    graph3->SetPoint(n, x[a], y[a]);

    TCanvas percorso("percorso", "percorso");

    graph3->SetTitle("best_path");
    graph3->GetXaxis()->SetTitle("x");
    graph3->GetYaxis()->SetTitle("y");
    graph3->Draw("ALP*");

    fitness_medio.close();
    fitness_migliore.close();

    app.Run();
    return 0;
}