#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <numeric>
#include "random.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TH1F.h"

using namespace std;

//questo programma fa i grafici delle diverse grandezze in funzione della
//temperatura prendendo i dati dai file energy_Gibbs.dat, energy_metro.dat
//e così via. Sul grafico di ciascuna grandezza è possibile vedere sia
//la versione fatta con il campionamento di Gibbs sia la versione fatta
//con l'algoritmo di Metropolis

int main(){
    TApplication app("app", 0, 0);
    ifstream Ene1, Heat1, Mag1, Chi1, Ene2, Heat2, Mag2, Chi2;
    unsigned int N = 16;
    double media = 0;
    double stddev = 0;
    double T = 0;
    Ene1.open("energy_metro.dat",ios::app);
    Heat1.open("heat_metro.dat",ios::app);
    Mag1.open("magn_metro.dat",ios::app);
    Chi1.open("susc_metro.dat",ios::app);
    Ene2.open("energy_Gibbs.dat",ios::app);
    Heat2.open("heat_Gibbs.dat",ios::app);
    Mag2.open("magn_Gibbs.dat",ios::app);
    Chi2.open("susc_Gibbs.dat",ios::app);

    //grafico energia

    TGraphErrors *graph1 = new TGraphErrors();
    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Ene1 >> T >> media >> stddev;
        stddev *= 100; //Gli errori sono moltiplicati per 100 in modo tale che
                       //possono essere visualizzati ad occhio nudo sul grafico
        graph1->SetPoint(i, T, media);
        graph1->SetPointError(i, 0.0, stddev);
        Ene2 >> T >> media >> stddev;
        stddev *= 100;
        media *= 1.05; //i valori medi di U/N ottenuti con il campionamento di
                       //Gibbs sono moltiplicati per 1.05 in modo tale che
                       //sulla canvas si possano distinguere le due curve
                       //ottenute con Metropolis e Gibbs
        graph2->SetPoint(i, T, media);
        graph2->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Energy = new TCanvas("Energy", "Energy", 800, 600);

    graph1->SetLineColor(kBlack);
    graph2->SetLineColor(kRed);
    graph1->SetTitle("Ising 1D, Internal energy");
    graph1->GetXaxis()->SetTitle("T");
    graph1->GetYaxis()->SetTitle("U/N");
    graph1->Draw("ALP");
    graph2->Draw("LP SAME");
    TLegend* legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend->AddEntry(graph1, "Metropolis", "l");
    legend->AddEntry(graph2, "Gibbs", "l");
    legend->Draw();
    Energy->Update();

    //grafico capacità termica

    TGraphErrors *graph3 = new TGraphErrors();
    TGraphErrors *graph4 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Heat1 >> T >> media >> stddev;
        stddev *= 10; //Per la stessa considerazione fatta per l'energia interna
                      //gli errori di C sono moltiplicati per 10
        graph3->SetPoint(i, T, media);
        graph3->SetPointError(i, 0.0, stddev);
        Heat2 >> T >> media >> stddev;
        stddev *= 10;
        media *= 1.05; //i valori medi di C ottenuti con il campionamento di
                       //Gibbs sono moltiplicati per 1.05 in modo tale che
                       //sulla canvas si possano distinguere le due curve
                       //ottenute con Metropolis e Gibbs
        graph4->SetPoint(i, T, media);
        graph4->SetPointError(i, 0.0, stddev);
    }

    TCanvas Heat_capacity("Heat_capacity", "Heat_capacity", 800, 600);
    
    graph3->SetLineColor(kBlack);
    graph4->SetLineColor(kRed);
    graph3->SetTitle("Ising 1D, Heat capacity");
    graph3->GetXaxis()->SetTitle("T");
    graph3->GetYaxis()->SetTitle("C");
    graph3->Draw("ALP");
    graph4->Draw("LP SAME");
    TLegend* legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend2->AddEntry(graph3, "Metropolis", "l");
    legend2->AddEntry(graph4, "Gibbs", "l");
    legend2->Draw();
    Heat_capacity.Update();


    //grafico magnetizzazione

    TGraphErrors *graph5 = new TGraphErrors();
    TGraphErrors *graph6 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Mag1 >> T >> media >> stddev;
        stddev *= 10; //Per la stessa considerazione fatta per l'energia interna
                      //gli errori di M sono moltiplicati per 10
        graph5->SetPoint(i, T, media);
        graph5->SetPointError(i, 0.0, stddev);
        Mag2 >> T >> media >> stddev;
        stddev *= 10;
        media *= 1.15; //i valori medi di M ottenuti con il campionamento di
                       //Gibbs sono moltiplicati per 1.15 in modo tale che
                       //sulla canvas si possano distinguere le due curve
                       //ottenute con Metropolis e Gibbs
        graph6->SetPoint(i, T, media);
        graph6->SetPointError(i, 0.0, stddev);
    }

    TCanvas Magnetization("Magnetization", "Magnetization", 800, 600);

    graph5->SetLineColor(kBlack);
    graph6->SetLineColor(kRed);
    graph5->SetTitle("Ising 1D, Magnetization");
    graph5->GetXaxis()->SetTitle("T");
    graph5->GetYaxis()->SetTitle("M");
    graph5->Draw("ALP");
    graph6->Draw("LP SAME");
    TLegend* legend3 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend3->AddEntry(graph3, "Metropolis", "l");
    legend3->AddEntry(graph4, "Gibbs", "l");
    legend3->Draw();
    Magnetization.Update();
    
    
    //grafico suscettività

    TGraphErrors *graph7 = new TGraphErrors();
    TGraphErrors *graph8 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Chi1 >> T >> media >> stddev;
        stddev *= 20; //Per la stessa considerazione fatta per l'energia interna
                      //gli errori di X sono moltiplicati per 20
        graph7->SetPoint(i, T, media);
        graph7->SetPointError(i, 0.0, stddev);
        Chi2 >> T >> media >> stddev;
        stddev *= 20;
        media *= 1.15; //i valori medi di X ottenuti con il campionamento di
                       //Gibbs sono moltiplicati per 1.15 in modo tale che
                       //sulla canvas si possano distinguere le due curve
                       //ottenute con Metropolis e Gibbs
        graph8->SetPoint(i, T, media);
        graph8->SetPointError(i, 0.0, stddev);
    }

    TCanvas Susceptibility("Susceptibility", "Susceptibility", 800, 600);

    graph7->SetLineColor(kBlack);
    graph8->SetLineColor(kRed);
    graph7->SetTitle("Ising 1D, Susceptibility");
    graph7->GetXaxis()->SetTitle("T");
    graph7->GetYaxis()->SetTitle("X");
    graph7->Draw("ALP");
    graph8->Draw("LP SAME");
    TLegend* legend4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend4->AddEntry(graph3, "Metropolis", "l");
    legend4->AddEntry(graph4, "Gibbs", "l");
    legend4->Draw();
    Susceptibility.Update();
    

    Ene1.close();
    Heat1.close();
    Mag1.close();
    Chi1.close();
    Ene2.close();
    Heat2.close();
    Mag2.close();
    Chi2.close();

    app.Run();
    return 0;
}