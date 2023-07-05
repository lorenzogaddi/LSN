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

//questo programma fa i fit dei diversi grafici ottenuti dai dati contenuti
//nei file energy_metro.dat, energy_Gibbs.dat e così via con le funzioni
//ricavate risolvendo esattamente il modello di Ising 1D

// Dichiarazione delle funzioni di fit
Double_t energy(Double_t *x, Double_t *par) {
    int Ns = 50;
    double T = x[0];
    double J = par[0];
    double th = tanh(J / T);
    double thN = pow(th, Ns);
    double ch = 1.0 / th;
    double e = -J * (th + ch * thN) / (1.0 + thN);
    return e;
}

Double_t heat(Double_t *x, Double_t *par) {
    int Ns = 50;
    double T = x[0];
    double J = par[0];
    double th = tanh(J / T);
    double thN = pow(th, Ns);
    double ch = 1.0 / th;
    double heat=(pow((J / T), 2))*(((1+thN+(Ns-1)*(pow(th, 2))+(Ns-1)*(pow(ch, 2))*thN)/(1+thN))-Ns*pow(((th+ch*thN)/(1+thN)), 2));
    return heat;
}

double_t magn(Double_t *x, Double_t *par){
    double J = par[0];
    double h = par[1];
    double T = x[0];
    int Ns = 50;
    double b = 1.0 / T;
    double l1 = exp(b * J) * cosh(b * h) + sqrt(exp(2 * b * J) * cosh(b * h) * cosh(b * h) - 2 * sinh(2 * b * J));
    double l2 = exp(b * J) * cosh(b * h) - sqrt(exp(2 * b * J) * cosh(b * h) * cosh(b * h) - 2 * sinh(2 * b * J));
    double Z = pow(l1, Ns) + pow(l2, Ns);
    double M = (exp(b * J) * sinh(b * h) * ((pow(l1, Ns - 1)) * (1 + exp(b * J) * cosh(b * h) / sqrt(exp(2 * b * J) * cosh(b * h) * cosh(b * h) - 2 * sinh(2 * b * J)))
    + (pow(l2, Ns - 1)) * (1 - exp(b * J) * cosh(b * h) / sqrt(exp(2 * b * J) * cosh(b * h) * cosh(b * h) - 2 * sinh(2 * b * J))))) / (Z);
    return M;
}

Double_t susc(Double_t *x, Double_t *par) {
    int Ns = 50;
    double T = x[0];
    double beta = 1 / T;
    double J = par[0];
    double th = tanh(J / T);
    double thN = pow(th, Ns);
    double ch = 1.0 / th;
    double X = beta * exp(2.0 * beta * J) * (1.0 - thN) / (1.0 + thN);
    return X;
}

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

    //grafico energia Metropolis

    TGraphErrors *graph1 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Ene1 >> T >> media >> stddev;
        graph1->SetPoint(i, T, media);
        graph1->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Energy_Metro = new TCanvas("Energy_Metropolis", "Energy_Metropolis", 800, 600);
    TF1 *energy_fit = new TF1("energy_fit", energy, 0.5, 2.0, 1);
    energy_fit->SetParameter(0, 1);

    graph1->SetTitle("Ising 1D, Internal energy Metropolis");
    graph1->GetXaxis()->SetTitle("T");
    graph1->GetYaxis()->SetTitle("U/N");
    graph1->Draw("ALP");
    graph1->Fit("energy_fit", "R");
    cout << "Chi quadro per U/N per Metropolis: " << energy_fit->GetChisquare() << endl;

    //grafico energia Gibbs
    
    TGraphErrors *graph2 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Ene2 >> T >> media >> stddev;
        graph2->SetPoint(i, T, media);
        graph2->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Energy_Gibbs = new TCanvas("Energy_Gibbs", "Energy_Gibbs", 800, 600);
    TF1 *energy_fit2 = new TF1("energy_fit2", energy, 0.5, 2.0, 1);
    energy_fit2->SetParameter(0, 1);

    graph2->SetTitle("Ising 1D, Internal energy Gibbs");
    graph2->GetXaxis()->SetTitle("T");
    graph2->GetYaxis()->SetTitle("U/N");
    graph2->Draw("ALP");
    graph2->Fit("energy_fit2");
    cout << "Chi quadro per U/N per Gibbs: " << energy_fit2->GetChisquare() << endl;



    //grafico capacità termica Metropolis

    TGraphErrors *graph3 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Heat1 >> T >> media >> stddev;
        graph3->SetPoint(i, T, media);
        graph3->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Heat_Metro = new TCanvas("Heat_Metropolis", "Heat_Metropolis", 800, 600);
    TF1 *heat_fit = new TF1("heat_fit", heat, 0.5, 2.0, 1);
    heat_fit->SetParameter(0, 1);

    graph3->SetTitle("Ising 1D, Heat Capacity Metropolis");
    graph3->GetXaxis()->SetTitle("T");
    graph3->GetYaxis()->SetTitle("C");
    graph3->Draw("ALP");
    graph3->Fit("heat_fit", "R");
    cout << "Chi quadro per C per Metropolis: " << heat_fit->GetChisquare() << endl;
    
    //grafico capacità termica Gibbs

    TGraphErrors *graph4 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Heat2 >> T >> media >> stddev;
        graph4->SetPoint(i, T, media);
        graph4->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Heat_Gibbs = new TCanvas("Heat_Gibbs", "Heat_Gibbs", 800, 600);
    TF1 *heat_fit2 = new TF1("heat_fit2", heat, 0.5, 2.0, 1);
    heat_fit2->SetParameter(0, 1);

    graph4->SetTitle("Ising 1D, Heat Capacity Gibbs");
    graph4->GetXaxis()->SetTitle("T");
    graph4->GetYaxis()->SetTitle("C");
    graph4->Draw("ALP");
    graph4->Fit("heat_fit2", "R");
    cout << "Chi quadro per C per Gibbs: " << heat_fit2->GetChisquare() << endl;

    //grafico magnetizzazione Metropolis

    TGraphErrors *graph5 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Mag1 >> T >> media >> stddev;
        graph5->SetPoint(i, T, media);
        graph5->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Magn_Metro = new TCanvas("Magn_Metropolis", "Magn_Metropolis", 800, 600);
    TF1 *magn_fit = new TF1("magn_fit", magn, 0.5, 2.0, 2);
    magn_fit->SetParameter(0, 1);
    magn_fit->SetParameter(1, 0.02);

    graph5->SetTitle("Ising 1D, Magnetization Metropolis");
    graph5->GetXaxis()->SetTitle("T");
    graph5->GetYaxis()->SetTitle("M");
    graph5->Draw("ALP");
    graph5->Fit("magn_fit", "R");
    cout << "Chi quadro per M per Metropolis: " << magn_fit->GetChisquare() << endl;

    //grafico magnetizzazione Gibbs

    TGraphErrors *graph6 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Mag2 >> T >> media >> stddev;
        graph6->SetPoint(i, T, media);
        graph6->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Magn_Gibbs = new TCanvas("Magn_Gibbs", "Magn_Gibbs", 800, 600);
    TF1 *magn_fit2 = new TF1("magn_fit2", magn, 0.5, 2.0, 2);
    magn_fit2->SetParameter(0, 1);
    magn_fit2->SetParameter(1, 0.02);

    graph6->SetTitle("Ising 1D, Magnetization Gibbs");
    graph6->GetXaxis()->SetTitle("T");
    graph6->GetYaxis()->SetTitle("M");
    graph6->Draw("ALP");
    graph6->Fit("magn_fit2", "R");
    cout << "Chi quadro per M per Gibbs: " << magn_fit2->GetChisquare() << endl;

    //grafico suscettività Metropolis

    TGraphErrors *graph7 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Chi1 >> T >> media >> stddev;
        graph7->SetPoint(i, T, media);
        graph7->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Susc_Metro = new TCanvas("Susc_Metropolis", "Susc_Metropolis", 800, 600);
    TF1 *susc_fit = new TF1("susc_fit", susc, 0.5, 2.0, 1);
    susc_fit->SetParameter(0, 1);

    graph7->SetTitle("Ising 1D, Susceptivity Metropolis");
    graph7->GetXaxis()->SetTitle("T");
    graph7->GetYaxis()->SetTitle("X");
    graph7->Draw("ALP");
    graph7->Fit("susc_fit", "R");
    cout << "Chi quadro per X per Metropolis: " << susc_fit->GetChisquare() << endl;


    //grafico suscettività Gibbs

    TGraphErrors *graph8 = new TGraphErrors();

    for(int i = 0; i < N; i++){
        Chi2 >> T >> media >> stddev;
        graph8->SetPoint(i, T, media);
        graph8->SetPointError(i, 0.0, stddev);
    }

    TCanvas* Susc_Gibbs = new TCanvas("Susc_Gibbs", "Susc_Gibbs", 800, 600);
    TF1 *susc_fit2 = new TF1("susc_fit2", susc, 0.5, 2.0, 1);
    susc_fit2->SetParameter(0, 1);

    graph8->SetTitle("Ising 1D, Susceptivity Gibbs");
    graph8->GetXaxis()->SetTitle("T");
    graph8->GetYaxis()->SetTitle("X");
    graph8->Draw("ALP");
    graph8->Fit("susc_fit2", "R");
    cout << "Chi quadro per X per Gibbs: " << susc_fit2->GetChisquare() << endl;
    
    

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