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

//La temperatura iniziale per la fase liquida tale per cui negli istanti
//temporali successivi si stabilizza attorno a 1.1 è 1.981

//La temperatura iniziale per la fase solida tale per cui negli istanti
//temporali successivi si stabilizza attorno a 0.8 è 1.543

//La temperatura iniziale per la fase gassosa tale per cui negli istanti
//temporali successivi si stabilizza attorno a 1.2 è 0.958

//questo programma fa i grafici delle quantità richieste dall'esercizio
//usando i dati dei file creati dal programma NVE_NVT.exe. Per visualizzare
//i grafici di ciascuna fase, è sufficiente scrivere alle righe 39, 
//40, 41, 42, 43 _S  per i grafici della fase solida, _L per quella liquida 
//e _G per quella gassosa.


int main(){
    TApplication app("app", 0, 0);
    ifstream Epot, Ekin, Etot, Temp, Pres;
    unsigned int N = 20;
    double media = 0;
    double stddev = 0;
    double a = 0;
    Epot.open("output_epot_G.out",ios::app);
    Ekin.open("output_ekin_G.out",ios::app);
    Temp.open("output_temp_G.out",ios::app);
    Etot.open("output_etot_G.out",ios::app);
    Pres.open("output_pres_G.out",ios::app);

    //grafico energia interna per particella

    TGraphErrors *graph1 = new TGraphErrors();

    //Faccio il fit dell'energia interna per particella con una funzione del
    //tipo f(x) = cost in quanto l'energia interna del sistema si deve
    //conservare. Ricavo successivamente il chi-quadro per verificare la
    //compatibilità del grafico di E_tot/N con f(x)
    TF1 *fitfunction = new TF1("fitfunction", "[0]", 0, 100);

    for(int i = 0; i < N; i++){
        Etot >> a >> a >> media >> stddev;
        graph1->SetPoint(i, (i+1), media);
        graph1->SetPointError(i, 0.0, stddev);
    }

    TCanvas Etotal("E/Npart", "E/Npart");

    graph1->SetTitle("E/Npart");
    graph1->GetXaxis()->SetTitle("N");
    graph1->GetYaxis()->SetTitle("E/Npart");
    graph1->Draw("ALP");
    graph1->Fit("fitfunction");
    cout << "Chi quadro ricavato dal fit del grafico di E_tot/N con una";
    cout << " funzione costante: " << fitfunction->GetChisquare() << endl;

    //grafico energia potenziale per particella

    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Epot >> a >> a >> media >> stddev;
        graph2->SetPoint(i, (i+1), media);
        graph2->SetPointError(i, 0.0, stddev);
        if(i==N-1){
            cout << "Valore finale di U/N: " << media << " ± " << stddev << endl;
        }
    }

    TCanvas Epotential("U/Npart", "U/Npart");

    graph2->SetTitle("U/Npart");
    graph2->GetXaxis()->SetTitle("N");
    graph2->GetYaxis()->SetTitle("U/Npart");
    graph2->Draw("ALP");

    //grafico energia cinetica per particella

    TGraphErrors *graph3 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Ekin >> a >> a >> media >> stddev;
        graph3->SetPoint(i, (i+1), media);
        graph3->SetPointError(i, 0.0, stddev);
        if(i==N-1){
            cout << "Valore finale di K/N: " << media << " ± " << stddev << endl;
        }
    }

    TCanvas Ekinetic("K/Npart", "K/Npart");

    graph3->SetTitle("K/Npart");
    graph3->GetXaxis()->SetTitle("N");
    graph3->GetYaxis()->SetTitle("K/Npart");
    graph3->Draw("ALP");

    //grafico temperatura

    TGraphErrors *graph4 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Temp >> a >> a >> media >> stddev;
        graph4->SetPoint(i, (i+1), media);
        graph4->SetPointError(i, 0.0, stddev);
        if(i==N-1){
            cout << "Valore finale di T: " << media << " ± " << stddev << endl;
        }
    }

    TCanvas Temperature("T", "T");

    graph4->SetTitle("T");
    graph4->GetXaxis()->SetTitle("N");
    graph4->GetYaxis()->SetTitle("T");
    graph4->Draw("ALP");

    //grafico pressione

    TGraphErrors *graph5 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Pres >> a >> a >> media >> stddev;
        graph5->SetPoint(i, (i+1), media);
        graph5->SetPointError(i, 0.0, stddev);
        if(i==N-1){
            cout << "Valore finale di P: " << media << " ± " << stddev << endl;
        }
    }

    TCanvas Pressure("P", "P");

    graph5->SetTitle("P");
    graph5->GetXaxis()->SetTitle("N");
    graph5->GetYaxis()->SetTitle("P");
    graph5->Draw("ALP");


    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Pres.close();

    app.Run();
    return 0;
}