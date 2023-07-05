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


int main(){
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

   //calcolo integrale con campionamento di x uniforme su zero e uno

   unsigned int N = 100; //numero di blocchi
    unsigned int M = 100000; //Numero di lanci
    unsigned int steps = (unsigned int) M / N;
    double s = 0;
    double integrale = 0;
    double somma_totale = 0;
    double somma_quadrati = 0;
    vector<double> media; //vettore che contiene i valori medi dell'integrale per ogni valore di N
    vector<double> stddev; //vettore che contiene le deviazioni standard della media dell'integrale per ogni valore di N

    //Riempimento dei vettori media e stddev per il calcolo dell'integrale

    for(int i = 0; i < N; i++){
        for(int m = 0; m < steps; m++){
            s = rnd.Rannyu(0, 1);
            integrale += (double) (M_PI*cos((double) (M_PI*s) / 2)) / 2;
            }
        integrale = (double) integrale / steps;
        somma_totale += integrale;
        somma_quadrati += integrale*integrale;
        media.push_back((double) somma_totale / (i+1));
        if(i > 0){
            stddev.push_back(error(media[i], somma_quadrati, i+1));
        } else{
            stddev.push_back(0);
        }
        integrale = 0;
    }


    //Grafico del valore dell'integrale con la relativa deviazione standard
    //della media in funzione di N

    TApplication app("app", 0, 0);
    TGraphErrors *graph1 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph1->SetPoint(i, (i+1), media[i]);
        graph1->SetPointError(i, 0.0, stddev[i]);
    }

    TCanvas Distrib_unif("Distrib_unif", "Distrib_unif");

    graph1->SetTitle("Distrib_unif");
    graph1->GetXaxis()->SetTitle("N");
    graph1->GetYaxis()->SetTitle("Valore integrale");
    graph1->Draw("ALP");

    cout << "Valore finale dell'integrale per distribuzione uniforme: " << media[N-1] << " ± " << stddev[N-1] << endl;


    //Riempimento dei vettori media2 e stddev2 per il calcolo dell'integrale
    //campionando x con distribuzione di probabilità p(x) = 2*(1-x)
    double r = 0;
    double integrale2 = 0;
    double somma_totale2 = 0;
    double somma_quadrati2 = 0;
    vector<double> media2; //vettore che contiene i valori medi dell'integrale per ogni valore di N
    vector<double> stddev2; //vettore che contiene le deviazioni standard della media dell'integrale per ogni valore di N
    for(int i = 0; i < N; i++){
        for(int m = 0; m < steps; m++){
            r = rnd.Polinomio();
            integrale2 += (double) ((M_PI*cos((double) M_PI*r / 2)) / (4*(1-r)));
            }
        integrale2 = (double) integrale2 / steps;
        somma_totale2 += integrale2;
        somma_quadrati2 += integrale2*integrale2;
        media2.push_back((double) somma_totale2 / (i+1));
        if(i > 0){
            stddev2.push_back(error(media2[i], somma_quadrati2, i+1));
        } else{
            stddev2.push_back(0);
        }
        integrale2 = 0;
    }

    //Grafico del valore dell'integrale con la relativa deviazione standard
    //della media in funzione di N

    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph2->SetPoint(i, (i+1), media2[i]);
        graph2->SetPointError(i, 0.0, stddev2[i]);
    }

    TCanvas distrib_polin("Distrib_polin", "Distrib_polin");

    graph2->SetTitle("Distrib_polin");
    graph2->GetXaxis()->SetTitle("N");
    graph2->GetYaxis()->SetTitle("valore integrale");
    graph2->Draw("ALP");

    cout << "Valore finale dell'integrale per distribuzione polinomiale: " << media2[N-1] << " ± " << stddev2[N-1] << endl;


    app.Run();
   return 0;
}