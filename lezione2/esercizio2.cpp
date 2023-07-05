#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <numeric>
#include <sstream>
#include <iomanip>
#include "random.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TH1F.h"

using namespace std;


// Funzione per calcolare la media degli elementi di un vettore
double mean(vector<double> arr) {
    double sum = accumulate(arr.begin(), arr.end(), 0.0);
    double media = sum / arr.size();
    return media;
}

// Funzione per calcolare la deviazione standard della media degli elementi di un vettore
double deviazioneStandardm(vector<double> arr) {
    double m = mean(arr);
    double acc = 0.0;
    for (int i = 0; i < arr.size(); i++) {
        acc += pow(arr[i], 2);
    }
    double error = sqrt( (( acc / arr.size())- m*m) / (arr.size() - 1));
    return error;
}


int main(){

    TApplication app("app", 0, 0);

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

    //Random walk for a cubic lattice
    unsigned int M = 10000; //simulazioni totali
    unsigned int N = 100; //numero blocchi
    unsigned int step = 100; //numero di salti fatti per la simulazione di un singolo random walk
    unsigned int steps = (int) M / N; //numero di elementi contenuti in un singolo blocco
    vector<int> xyz(3*M); //vettore che contiene le coordinate x, y e z
                            //durante ogni step per le M simulazioni.
                            //I primi tre elementi del vettore corrispondono
                            //rispettivamente alle coordinate x, y e z
                            //per la prima simulazione, i
                            //successivi tre per la seconda e così via
    vector<double> A_n; //vettore che contiene gli elementi A_i ricavati
                           //da ciascun blocco
    vector<double> media;
    vector<double> stddev;
    for(int i = 0; i < step; i++){
        //riempio il vettore xyz
        for(int j = 0; j < M; j++){
            double x = rnd.Rannyu();
            int bin = floor(x*3);
            int indice = (3*j) + bin;
            double r = rnd.Rannyu();
            if(r >= 0.5){
            xyz[indice] = xyz[indice] + 1;
            } else{
                xyz[indice] = xyz[indice] - 1;
            }
        }
        //riempio il vettore A_n
        for(int j = 0; j < N; j++){
            int u = 0;
            for(int m = 0; m < steps; m++){
                int x = xyz[(j*steps) + 3*m];
                int y = xyz[(j*steps) + 3*m + 1];
                int z = xyz[(j*steps) + 3*m + 2];
                u = u + x*x + y*y + z*z;
            }
            double A = (double) u / steps;
            A_n.push_back(A);
        }
        double a = round(mean(A_n)*100) / 100;
        media.push_back(a);
        stddev.push_back(deviazioneStandardm(A_n));
        A_n.clear();
    }

    //Grafico della radice quadrata del valore medio di r^2 in funzione
    //degli step
    TGraphErrors *graph1 = new TGraphErrors();

    graph1->SetPoint(0, 0, 0);
    graph1->SetPointError(0, 0, 0);
    for(int i = 0; i < step; i++){
        graph1->SetPoint(i+1, (i+1), sqrt(media[i]));
        //L'errore di sqrt(<r^2>) si può trovare utilizzando la propagazione
        //degli errori. Il nuovo errore è stddev[i] / (2*sqrt(<r^2>))
        double a = round((stddev[i] / (2*sqrt(media[i])))*100) / 100;
        graph1->SetPointError(i+1, 0.0, a);
    }

    TCanvas Cubic_lattice("Cubic lattice", "Cubic lattice");

    graph1->SetTitle("Cubic lattice");
    graph1->GetXaxis()->SetTitle("Step");
    graph1->GetYaxis()->SetTitle("sqrt(<r^2>)");
    graph1->Draw("ALP");
    //fit del grafico con una funzione del tipo f(x) = k*sqrt(x)
    TF1 *fitfunction = new TF1("fitfunction", "[0]*sqrt(x)", 0, 100);
    graph1->Fit("fitfunction");
    cout << "Chi quadro: " << fitfunction->GetChisquare() << endl;


    //Per risolvere la seconda parte dell'esercizio usiamo lo stesso ragionamento 
    //e anche procedimento, però cambiando il modo di riempire il vettore xyz2: 
    //generiamo theta e phi campionando l'angolo solido e poi troviamo x, y e z 
    //con le coordinate sferiche: x = sin(theta)*cos(phi), y = sin(theta)*sin(phi) 
    //e z = cos(theta). Il nuovo vettore xyz2 ha 3*M componenti, i cui primi tre 
    //elementi sono le coordinate x, y e z della prima simulazione, i successivi 
    //tre elementi sono le coordinate della seconda simulazione e così via, come
    //nel punto 1) dell'esercizio.

    vector<double> xyz2(3*M);
    vector<double> A_n2; 
    vector<double> media2;
    vector<double> stddev2;

    for(int i = 0; i < step; i++){
        //riempio il vettore xyz2 uno step alla volta per poi ricavare
        //media e deviazione standard della media di sqrt(<r^2>) con il
        //metodo dei blocchi
        for(int j = 0; j < M; j++){
            double theta = acos(2*rnd.Rannyu() - 1);
            double phi = rnd.Rannyu(0, 2*M_PI);
            xyz2[3*j] = xyz2[3*j] + sin(theta)*cos(phi);
            xyz2[(3*j) + 1] = xyz2[(3*j) + 1] + sin(theta)*sin(phi);
            xyz2[(3*j) + 2] = xyz2[(3*j) + 2] + cos(theta);
        }
        //riempio il vettore A_i
        for(int j = 0; j < N; j++){
            double A = 0;
            for(int m = 0; m < steps; m++){
                double x = xyz2[(j*steps) + 3*m];
                double y = xyz2[(j*steps) + 3*m + 1];
                double z = xyz2[(j*steps) + 3*m + 2];
                A = A + pow(x, 2) + pow(y, 2) + pow(z, 2);
            }
            A = (double) A / steps;
            A_n2.push_back(A);
        }
        double a = round(mean(A_n2)*100) / 100;
        media2.push_back(mean(A_n2));
        stddev2.push_back(deviazioneStandardm(A_n2));
        A_n2.clear();
    }

    //riempimento grafico per random walk nel continuo
    TGraphErrors *graph2 = new TGraphErrors();
    TF1 *fitfunction2 = new TF1("fitfunction2", "[0]*sqrt(x)", 0, 100);

    graph2->SetPoint(0, 0, 0);
    graph2->SetPointError(0, 0, 0);
    for(int i = 0; i < step; i++){
        graph2->SetPoint(i+1, (i+1), sqrt(media2[i]));
        //L'errore di sqrt(<r^2>) si può trovare utilizzando la propagazione
        //degli errori. Il nuovo errore è stddev2[i] / (2*sqrt(<r^2>))
        double a = round((stddev2[i] / (2*sqrt(media2[i])))*100) / 100;
        graph2->SetPointError(i+1, 0.0, a);
    }

    TCanvas Continuum("Continuum", "Continuum");

    graph2->SetTitle("Continuum");
    graph2->GetXaxis()->SetTitle("Step");
    graph2->GetYaxis()->SetTitle("sqrt(<r^2>)");
    graph2->Draw("ALP");

    //fit del grafico con una funzione del tipo f(x) = k*sqrt(x)
    graph2->Fit("fitfunction2");
    cout << "Chi quadro: " << fitfunction2->GetChisquare() << endl;

    app.Run();
    return 0;
}