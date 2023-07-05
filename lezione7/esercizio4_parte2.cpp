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

//Questo programma fa i grafici di U/N, P, g(r) per NVE e NVT. Tutte
//le grandezze sono espresse in unità SI
int main(){
    TApplication app("app", 0, 0);
    ifstream Epot, Pres, Gdir1, Gdir2;
    unsigned int N = 20;
    double media = 0;
    double stddev = 0;
    double a = 0;
    double sigma = 0.34E-9;
    double k_B = 1.380649E-23;
    double temp = 120;
    double epsilon = temp*k_B;
    vector<double> gdir1;
    vector<double> gdir2;
    vector<double> err_gdir1;
    vector<double> err_gdir2;
    Epot.open("output_epot_G.out",ios::app);
    Pres.open("output_pres_G.out",ios::app);
    Gdir1.open("gdir_NVT_G.out",ios::app);
    Gdir2.open("gdir_NVE_G.out",ios::app);

    //grafico energia potenziale per particella

    TGraphErrors *graph1 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Epot >> a >> a >> media >> stddev;
        media *= epsilon;
        stddev *= epsilon;
        graph1->SetPoint(i, (i+1), media);
        graph1->SetPointError(i, 0.0, stddev);
    }

    cout << "U/N: " << media << " ± " << stddev << endl;

    TCanvas Epotential("U/Npart", "U/Npart");

    graph1->SetTitle("U/Npart");
    graph1->GetXaxis()->SetTitle("N");
    graph1->GetYaxis()->SetTitle("U/Npart (J)");
    graph1->Draw("ALP");

    //grafico pressione

    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Pres >> a >> a >> media >> stddev;
        media *= epsilon/pow(sigma,3);
        stddev *= epsilon/pow(sigma,3);
        graph2->SetPoint(i, (i+1), media);
        graph2->SetPointError(i, 0.0, stddev);
    }

    cout << "P: " << media << " ± " << stddev << endl;

    TCanvas Pressure("P", "P");

    graph2->SetTitle("P");
    graph2->GetXaxis()->SetTitle("N");
    graph2->GetYaxis()->SetTitle("P (Pa)");
    graph2->Draw("ALP");

    //grafico g(r) NVT
    TGraphErrors *graph3 = new TGraphErrors();

    for(int i=0; i < 100; i++){
        Gdir1 >> a >> a >> media >> stddev;
        gdir1.push_back(media);
        err_gdir1.push_back(stddev);
        graph3->SetPoint(i, a*sigma, media);
        graph3->SetPointError(i, 0.0, stddev);
    }

    TCanvas radial_distribution1("g(r)", "g(r)");

    graph3->SetTitle("g(r) NVT");
    graph3->GetXaxis()->SetTitle("r (m)");
    graph3->GetYaxis()->SetTitle("g(r)");
    graph3->Draw("ALP");

    //grafico g(r) NVE
    TGraphErrors *graph4 = new TGraphErrors();

    for(int i=0; i < 100; i++){
        Gdir2 >> a >> a >> media >> stddev;
        gdir2.push_back(media);
        err_gdir2.push_back(stddev);
        graph4->SetPoint(i, a*sigma, media);
        graph4->SetPointError(i, 0.0, stddev);
    }

    TCanvas radial_distribution2("g(r)", "g(r)");

    graph4->SetTitle("g(r) NVE");
    graph4->GetXaxis()->SetTitle("r (m)");
    graph4->GetYaxis()->SetTitle("g(r)");
    graph4->Draw("ALP");

    //calcolo del chi_quadro per (g(r)_NVT - g(r)_NVE) con la funzione f(r)=0
    double chi = 0;
    for(int i = 0; i < 100; i++){
        media = gdir1[i] - gdir2[i];
        stddev = sqrt(pow(err_gdir1[i],2)+pow(err_gdir2[i],2));
        if(stddev==0){
            stddev=1;
        }
        chi += pow(media,2)/pow(stddev,2);
    }

    cout << "Chi quadrato di g1(r)-g2(r) con la funzione f(r)=0: " << chi << endl;

    Epot.close();
    Pres.close();
    Gdir1.close();
    Gdir2.close();

    app.Run();
    return 0;
}