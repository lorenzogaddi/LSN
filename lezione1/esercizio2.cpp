#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TH1F.h"

using namespace std;

int main(){
    unsigned int N = 100;
    unsigned int M = 10000;
    double Sn = 0;

    //Inizializzazione generatore
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

   //Fine inizializzazione gen

  
    TApplication app("app", 0, 0);

    //Generazione e riempimento istogramma con M input di Sn campionati 
    //tramite lanci di dado classico

    TH1F histo1("standard_dice", "standard_dice", 100, 1.5, 5.5);

   
    for(int o = 0; o < M; o++){   
        for(int i = 0; i < N; i++){
            double r = rnd.Rannyu();
            int x = floor(r*6);
            Sn = Sn + (x + 1);
        }
        Sn = (double) Sn / N;
        histo1.Fill(Sn);
        Sn = 0;
    }

    TCanvas mycanvas1("standard_dice", "standard_dice");
    histo1.StatOverflows(kTRUE);
    histo1.GetXaxis()->SetTitle("Sn");
    histo1.GetYaxis()->SetTitle("conteggi");
    histo1.Fit("gaus");  //fit con una gaussiana dell'istogramma
    histo1.Draw();
    TF1 *fitFunc1 = histo1.GetFunction("gaus");
    cout << "Chi quadro per dado standard: " << fitFunc1->GetChisquare() << endl;


    //Generazione e riempimento istogramma con M input di Sn campionati 
    //tramite lanci di dado esponenziali
    TH1F histo2("exponential_dice", "exponential_dice", 200, 0.0, 2.0);

    for(int o = 0; o < M; o++){   
        for(int i = 0; i < N; i++){
            Sn = Sn + (rnd.Exponential(1.0));
        }
        Sn = (double) Sn / N;
        histo2.Fill(Sn);
        Sn = 0;
    }

    TCanvas mycanvas2("exponential_dice", "exponential_dice");
    histo2.StatOverflows(kTRUE);
    histo2.GetXaxis()->SetTitle("Sn");
    histo2.GetYaxis()->SetTitle("conteggi");
    histo2.Fit("gaus");   //fit con una gaussiana dell'istogramma
    histo2.Draw();
    TF1 *fitFunc2 = histo2.GetFunction("gaus");
    cout << "Chi quadro per dado esponenziale: " << fitFunc2->GetChisquare() << endl;

    //Generazione e riempimento istogramma con M input di Sn campionati 
    //tramite lanci di dado lorentziani

    TH1F histo3("lorentzian_dice", "lorentzian_dice", 200, -20.0, 20.0);

   
    for(int o = 0; o < M; o++){   
        for(int i = 0; i < N; i++){
            Sn = Sn + (rnd.Lorentzian(0.0, 1.0));
        }
        Sn = (double) Sn / N;
        histo3.Fill(Sn);
        Sn = 0;
    }

    TCanvas mycanvas3("lorentzian_dice", "lorentzian_dice");
    histo3.StatOverflows(kTRUE);
    histo3.GetXaxis()->SetTitle("Sn");
    histo3.GetYaxis()->SetTitle("conteggi");
    histo3.Draw();

    app.Run();
    return 0;
}