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


    //C e P nel caso in cui campiono direttamente S(T)

    double S_zero = 100;
    double T = 1;
    double K = 100;
    double r = 0.1;
    double volatility = 0.25;
    double C = 0;
    double P = 0;
    unsigned int N = 100; //numero di blocchi
    unsigned int M = 100000; //Numero di lanci
    unsigned int steps = (unsigned int) M / N;
    double s = 0;
    double somma_totale = 0;
    double somma_totale2 = 0;
    double somma_quadrati = 0;
    double somma_quadrati2 = 0;
    vector<double> media; //vettore che contiene i valori medi di C per ogni valore di N
    vector<double> stddev; //vettore che contiene le deviazioni standard della media di C per ogni valore di N
    vector<double> media2; //vettore che contiene i valori medi di P per ogni valore di N
    vector<double> stddev2; //vettore che contiene le deviazioni standard della media di P per ogni valore di N

    //Riempimento dei vettori media, media2, stddev2 e stddev per il calcolo 
    //di C e P

    for(int i = 0; i < N; i++){
        for(int m = 0; m < steps; m++){
            s = rnd.Gauss(0, 1);
            double S_T = S_zero*exp(T*(r - (double) (volatility*volatility) / 2) + volatility*sqrt(T)*s);
            if(S_T >= K){
            C += exp(-r*T)*(S_T - K);
            } else{
                P += exp(-r*T)*(K - S_T);
            }
            }
        C = (double) C / steps;
        P = (double) P / steps;
        somma_totale += C;
        somma_totale2 += P;
        somma_quadrati += C*C;
        somma_quadrati2 += P*P;
        media.push_back((double) somma_totale / (i+1));
        media2.push_back((double) somma_totale2 / (i+1));
        if(i > 0){
            stddev.push_back(error(media[i], somma_quadrati, i+1));
            stddev2.push_back(error(media2[i], somma_quadrati2, i+1));
        } else{
            stddev.push_back(0);
            stddev2.push_back(0);
        }
        C = 0;
        P = 0;
        }
        cout << "Call-option price per uno step: " << media[N-1] << " ± " << stddev[N-1] << endl;
        cout << "Put-option price per uno step: " << media2[N-1] << " ± " << stddev2[N-1] << endl;


    //Grafico del valore medio di C e P con la relativa deviazione standard
    //della media in funzione di N

    TApplication app("app", 0, 0);

    TGraphErrors *graph1 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph1->SetPoint(i, (i+1), media[i]);
        graph1->SetPointError(i, 0.0, stddev[i]);
    }

    TCanvas valor_medio_C1("C", "C");

    graph1->SetTitle("C_uno_step");
    graph1->GetXaxis()->SetTitle("N");
    graph1->GetYaxis()->SetTitle("C");
    graph1->Draw("ALP");

    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph2->SetPoint(i, (i+1), media2[i]);
        graph2->SetPointError(i, 0.0, stddev2[i]);
    }

    TCanvas valor_medio_P1("P", "P");

    graph2->SetTitle("P_uno_step");
    graph2->GetXaxis()->SetTitle("N");
    graph2->GetYaxis()->SetTitle("P");
    graph2->Draw("ALP");

    //azzeramento dei precedenti vettori e loro riempimento campionando
    //S(T) a passi
    media.clear();
    media2.clear();
    stddev.clear();
    stddev2.clear();
    somma_totale = 0;
    somma_totale2 = 0;
    somma_quadrati = 0;
    somma_quadrati2 = 0;
    double S = 0;
        for(int i = 0; i < N; i++){
        for(int m = 0; m < steps; m++){
            s = rnd.Gauss(0, 1);
            double S_T = S_zero*exp(0.01*(r - (double) (volatility*volatility) / 2) + volatility*sqrt(0.01)*s);
            for(int q = 0; q < 99; q++){
                s = rnd.Gauss(0, 1);
                S_T = S_T*exp(0.01*(r - (double) (volatility*volatility) / 2) + volatility*sqrt(0.01)*s);
            }

            if(S_T >= K){
                C += exp(-r*T)*(S_T - K);
            } else{
                P += exp(-r*T)*(K - S_T);
            }
            }
        C = (double) C / steps;
        P = (double) P / steps;
        somma_totale += C;
        somma_totale2 += P;
        somma_quadrati += C*C;
        somma_quadrati2 += P*P;
        media.push_back((double) somma_totale / (i+1));
        media2.push_back((double) somma_totale2 / (i+1));
        if(i > 0){
            stddev.push_back(error(media[i], somma_quadrati, i+1));
            stddev2.push_back(error(media2[i], somma_quadrati2, i+1));
        } else{
            stddev.push_back(0);
            stddev2.push_back(0);
        }
        C = 0;
        P = 0;
        }

        cout << "Call-option price per 100 step: " << media[N-1] << " ± " << stddev[N-1] << endl;
        cout << "Put-option price per 100 step: " << media2[N-1] << " ± " << stddev2[N-1] << endl;


    //Grafico del valore medio di C e P con la relativa deviazione standard
    //della media in funzione di N campionando S(T) per passi


    TGraphErrors *graph3 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph3->SetPoint(i, (i+1), media[i]);
        graph3->SetPointError(i, 0.0, stddev[i]);
    }

    TCanvas valor_medio_C2("C", "C");

    graph3->SetTitle("C_100_step");
    graph3->GetXaxis()->SetTitle("N");
    graph3->GetYaxis()->SetTitle("C");
    graph3->Draw("ALP");

    TGraphErrors *graph4 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        graph4->SetPoint(i, (i+1), media2[i]);
        graph4->SetPointError(i, 0.0, stddev2[i]);
    }

    TCanvas valor_medio_P2("P", "P");

    graph4->SetTitle("P_100_step");
    graph4->GetXaxis()->SetTitle("N");
    graph4->GetYaxis()->SetTitle("P");
    graph4->Draw("ALP");


    app.Run();

    return 0;
}