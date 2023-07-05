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
    //dati per metodo dei blocchi
    unsigned int N = 100; //numero di blocchi
    unsigned int M = 100000; //Numero di lanci
    unsigned int steps = (unsigned int) M / N;
    double s = 0;
    double integrale = 0;
    double somma_totale = 0;
    double somma_quadrati = 0;
    vector<double> media; //vettore che contiene i valori medi dell'integrale per ogni valore di N
    vector<double> stddev; //vettore che contiene le deviazioni standard della media dell'integrale per ogni valore di N

    //Riempimento dei vettori media e stddev per il calcolo del
    //valore medio di r

    for(int i = 0; i < N; i++){
        for(int m = 0; m < steps; m++){
            s = rnd.Rannyu(0, 1);
            integrale += s;
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


    //Grafico del valore medio di r con la relativa deviazione standard
    //della media in funzione di N

    TApplication app("app", 0, 0);
    TGraphErrors *graph1 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph1->SetPoint(i, (i+1), media[i]);
        graph1->SetPointError(i, 0.0, stddev[i]);
    }

    TCanvas valor_medio_r("valor_medio_r", "valor_medio_r");

    graph1->SetTitle("Valore medio di r");
    graph1->GetXaxis()->SetTitle("N");
    graph1->GetYaxis()->SetTitle("Valore medio di r");
    graph1->Draw("ALP");

    cout << "Per N = 100, <r> = " << media[N-1] << " ± " << stddev[N-1] << endl;


    //Riempimento dei vettori media2 e stddev2 per il calcolo della
    //varianza di r
    double r = 0;
    double integrale2 = 0;
    double somma_totale2 = 0;
    double somma_quadrati2 = 0;
    vector<double> media2; //vettore che contiene i valori medi dell'integrale per ogni valore di N
    vector<double> stddev2; //vettore che contiene le deviazioni standard della media dell'integrale per ogni valore di N
    for(int i = 0; i < N; i++){
        for(int m = 0; m < steps; m++){
            r = rnd.Rannyu();
            integrale2 += (r - 0.5)*(r - 0.5);
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

    //Grafico della varianza di r con la relativa deviazione standard
    //della media in funzione di N

    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph2->SetPoint(i, (i+1), media2[i]);
        graph2->SetPointError(i, 0.0, stddev2[i]);
    }

    TCanvas varianza_r("varianza_r", "varianza_r");

    graph2->SetTitle("Varianza r");
    graph2->GetXaxis()->SetTitle("N");
    graph2->GetYaxis()->SetTitle("Varianza r");
    graph2->Draw("ALP");

    cout << "Per N = 100, varianza = " << media2[N-1] << " ± " << stddev2[N-1] << endl;


    //istogramma valori chi-quadro ottenuti
    const int n = 10000;
    const int k = 100;
    const double bin_size = 1.0 / k;

    TH1F histo("chi_quadrato", "chi_quadrato", 20, 65, 130);

    for(int j = 0; j < 100; j++){  
        double observed[k] = {0};
        for (int i = 0; i < n; ++i) {
        double x = rnd.Rannyu();
        // Determinazione dell'indice del bin in cui cade il numero generato
        int bin = floor(x / bin_size);
        // Incremento del conteggio degli elementi del bin
        observed[bin] += 1;
       }
       // Calcolo del chi quadro con la formula di Pearson
       double expected = n / k;
       double chi_square = 0;

       for (int i = 0; i < k; ++i) {
           chi_square += pow(observed[i] - expected, 2) / expected;
       }
       //Riempimento istogramma con i chi quadro calcolati
       histo.Fill(chi_square);
    }

    cout << "Il valore medio del chi quadrato è: " << histo.GetMean() << endl;


    TCanvas mycanvas("chi_quadrato", "chi_quadrato");
    histo.StatOverflows(kTRUE);
    histo.GetXaxis()->SetTitle("Chi_quadrato");
    histo.GetYaxis()->SetTitle("conteggi");
    histo.Draw();


    app.Run();

    return 0;
}