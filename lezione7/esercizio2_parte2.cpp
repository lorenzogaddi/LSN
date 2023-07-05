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


//sono importanti solo i primi 1000 lags

//questo programma fa il grafico dell'autocorrelazione in funzione del tempo
//prendendo i dati dell'energia potenziale dai file instantaneous_v_S.out,
//instantaneous_v_L.out e instantaneous_v_G.out


//calcola varianza per i 500000 dati
double varianza(vector<double> vec){
    double mean = 0;
    double squares = 0;
    for(int i = 0; i<500000; i++){
        squares += vec[i]*vec[i];
        mean += vec[i];
    }
    mean = mean /(double) 500000;
    double varianza = (squares /(double) 500000) - pow(mean,2);
    return varianza;
}

//calcola una delle quantità richieste per l'autocorrelazione
double par1(vector<double> vec, int t){
    double p1 = 0;
    for(int i = 0; i < (500000-t); i++){
        p1 += vec[i]*vec[i+t];
    }
    return p1;
}

int main(){
    TApplication app("app", 0, 0);
    //riempio vettore con i valori dell'energia potenziale
    ifstream Epot;
    Epot.open("instantaneous_v_G.out",ios::app);
    vector<double> energia;
    double v = 0;
    for(int i = 0; i < 500000; i++){
        Epot >> v;
        energia.push_back(v);
    }
    Epot.close();

    //grafico autocorrelazione in funzione di t
    double var = varianza(energia);
    int a = 1;
    double autocor = 0;
    double p1 = 0;
    double p2 = 0;
    double p3 = 0;
    for(int i = 0; i < 500000; i++){
        p2 += energia[i];
    }
    p3 = p2;
    TGraph *graph = new TGraph();
    graph->SetPoint(0, 0, 1);

    for(int i = 1; i < 1000; i++){
        p1 = (double) par1(energia, i) / (500000-i);
        p2 -= energia[500000-i];
        p3 -= energia[i-1];
        double den = pow((500000-i),2);
        autocor = (double) (p1-((p2*p3)/den))/ var;
        graph->SetPoint(a, i, autocor);
        a += 1;
        autocor = 0;
        p1 = 0;
    }

    TCanvas Autocorrelation("Autocorrelazione", "Autocorrelazione");

    graph->SetTitle("Autocorrelazione");
    graph->GetXaxis()->SetTitle("t");
    graph->GetYaxis()->SetTitle("A(t)");
    graph->Draw("ALP");
    app.Run();
    return 0;
}