#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <numeric>
#include "random.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TH1F.h"

using namespace std;

//funzione per calcolare la deviazione standard della media a partire da 
//media, somma dei quadrati e numero di elementi di una serie di elementi


double error(double mean, double squares, unsigned int n) {
    double error = 0;
    error = ((double) squares / n) - (mean*mean);
    error = sqrt((double) error / (n-1));
    return error;
}

//funzione per generare theta

bool angolo(double x, double y){
    if((x*x + y*y) <= 1){
        return true;
    } else{
        return false;
    }
    
}


int main(){
    double L = 1.5;
    double d = 2;
    //Di seguito il ragionamento utilizzato per risolvere l'esercizio
    /*Per risolvere questo problema, possiamo generare un numero casuale s
    compreso tra 0 e d che rappresenta la distanza di un'estremità dell'ago 
    dalla linea orizzontale inferiore. Successivamente, possiamo 
    generare un numero compreso tra zero e due pi greco con il metodo
    accept-reject. L'ago interseca una linea se la quantità (s+Lsin(theta))
    t è maggiore di d o minore di zero.
    */

   //inizializzazione generatore
    Random rnd;
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

   //fine inizializzazione gen

    unsigned int N = 50; //numero di blocchi
    unsigned int N_throws = 10000; //numero di lanci dell'ago per calcolare ogni volta pi greco
    unsigned int M = 10000; //Numero di lanci totali
    unsigned int steps = (int) M / N;
    unsigned int N_hit = 0;
    double pi_greco = 0;
    double x = 0;
    double y = 0;
    double theta = 0;
    double s = 0;
    double somma = 0;
    double somma_totale = 0;
    double somma_quadrati = 0;
    vector<double> vec; //vettore per contenere step
    vector<double> media; //vettore che contiene i valori medi di pi greco per ogni valore di N
    vector<double> stddev; //vettore che contiene le deviazioni standard della media di pi greco per ogni valore di N

    //Il seguente ciclo fa ciò che è stato spiegato nel ragionamento iniziale
    for(int i = 0; i < N; i++){
        for(int n = 0; n < steps; n++){  
            for(int m = 0; m < N_throws; m++){
                s = rnd.Rannyu(0, d);
                bool logic = false;
                while(logic == false){
                    x = rnd.Rannyu();
                    y = rnd.Rannyu();
                    theta = 2*acos((double) x / (sqrt(x*x + y*y)));
                    logic = angolo(x, y);
                }
                if((s+L*sin(theta)) < 0 or (s+L*sin(theta)) > d){
                    N_hit = N_hit + 1;
                }
            }
        pi_greco = (double) (2*L*N_throws) / (d*N_hit);
        vec.push_back(pi_greco);
        N_hit = 0;
        }
        somma = accumulate(vec.begin(), vec.end(), 0.0);
        somma_totale += somma;
        media.push_back((double) somma_totale / (steps*(i+1)));
        somma_quadrati += (double) (somma*somma) / (steps*steps);
        if(i > 0){
            stddev.push_back(error(media[i], somma_quadrati, i+1));
        } else{
            stddev.push_back(0);
        }
        vec.clear();
        cout << "Block number " << i + 1 << endl << endl;
        cout << "----------------------------" << endl << endl;

    }
    cout << "Valore finale di pi greco: " << media[N-1] << "±" << stddev[N-1] << endl;

    //Disegno del grafico dei valori medi di pi greco col relativo errore in funzione di N

    TApplication app("app", 0, 0);
    TGraphErrors *graph = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph->SetPoint(i, (i+1), media[i]);
        graph->SetPointError(i, 0.0, stddev[i]);
    }

    TCanvas Pi_greco("Pi_greco", "Pi_greco");

    graph->SetTitle("Pi_greco");
    graph->GetXaxis()->SetTitle("N");
    graph->GetYaxis()->SetTitle("Pi greco");
    graph->Draw("ALP");

    app.Run();

    return 0;
}