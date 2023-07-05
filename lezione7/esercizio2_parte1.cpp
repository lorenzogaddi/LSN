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

//Lo step per avere circa un'acceptance rate del 50% col Metropolis per
//la fase liquida è 0.2

//Lo step per avere circa un'acceptance rate del 50% col Metropolis per
//la fase solida è 0.11

//Per la fase gassosa non è possibile avere un'acceptance rate del 50%.
//L'acceptance rate minore ottenibile è circa del 60% ed è ottenuta
//per uno step pari a L/2, ossia per 6.4633

//questo programma fa i grafici delle quantità dell'esercitazione
//usando i dati dei file creati dal programma NVE_NVT.exe. 

int main(){
    TApplication app("app", 0, 0);
    ifstream Epot, Ekin, Etot, Temp, Pres;
    unsigned int N = 20;
    double media = 0;
    double stddev = 0;
    double a = 0;
    Epot.open("output_epot.dat",ios::app);
    Ekin.open("output_ekin.dat",ios::app);
    Temp.open("output_temp.dat",ios::app);
    Etot.open("output_etot.dat",ios::app);
    Pres.open("output_pres.dat",ios::app);

    //grafico energia interna per particella

    TGraphErrors *graph1 = new TGraphErrors();

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

    //grafico energia potenziale per particella

    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Epot >> a >> a >> media >> stddev;
        graph2->SetPoint(i, (i+1), media);
        graph2->SetPointError(i, 0.0, stddev);
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