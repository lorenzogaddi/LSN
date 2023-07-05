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
//Questo programma fa i grafici dell'energia interna, della capacità termica,
//della magnetizzazione e della sucettività in funzione del numero
//di blocchi utilizzando i dati dei file output.ene.O, output.heat.O, 
//output.magn.O e output.susc.O
int main(){
    TApplication app("app", 0, 0);
    ifstream Ene, Heat, Mag, Chi;
    unsigned int N = 20;
    double media = 0;
    double stddev = 0;
    double a = 0;
    Ene.open("output.ene.0",ios::app);
    Heat.open("output.heat.0",ios::app);
    Mag.open("output.magn.0",ios::app);
    Chi.open("output.susc.0",ios::app);

    //grafico energia

    TGraphErrors *graph1 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Ene >> a >> a >> media >> stddev;
        graph1->SetPoint(i, (i+1), media);
        graph1->SetPointError(i, 0.0, stddev);
    }

    TCanvas Energy("Energy", "Energy");

    graph1->SetTitle("Ising 1D, Internal energy");
    graph1->GetXaxis()->SetTitle("N");
    graph1->GetYaxis()->SetTitle("U/N");
    graph1->Draw("ALP");

    //grafico capacità termica

    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Heat >> a >> a >> media >> stddev;
        graph2->SetPoint(i, (i+1), media);
        graph2->SetPointError(i, 0.0, stddev);
    }

    TCanvas Heat_capacity("C", "C");

    graph2->SetTitle("Ising 1D, Heat capacity");
    graph2->GetXaxis()->SetTitle("N");
    graph2->GetYaxis()->SetTitle("C");
    graph2->Draw("ALP");

    //grafico magnetizzazione

    TGraphErrors *graph3 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Mag >> a >> a >> media >> stddev;
        graph3->SetPoint(i, (i+1), media);
        graph3->SetPointError(i, 0.0, stddev);
    }

    TCanvas Magnetization("M", "M");

    graph3->SetTitle("Ising 1D, Magnetizaztion");
    graph3->GetXaxis()->SetTitle("N");
    graph3->GetYaxis()->SetTitle("M");
    graph3->Draw("ALP");

    //grafico magnetic susceptibility

    TGraphErrors *graph4 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Chi >> a >> a >> media >> stddev;
        graph4->SetPoint(i, (i+1), media);
        graph4->SetPointError(i, 0.0, stddev);
    }

    TCanvas Susceptibility("X", "X");

    graph4->SetTitle("Ising 1D, susceptibility");
    graph4->GetXaxis()->SetTitle("N");
    graph4->GetYaxis()->SetTitle("X");
    graph4->Draw("ALP");

    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();

    app.Run();
    return 0;
}