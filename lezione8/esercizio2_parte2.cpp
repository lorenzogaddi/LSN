#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TH1F.h"

using namespace std;

//Questo programma disegna il grafico di <H> in funzione del numero di blocchi
//per i parametri mu e sigma che minimizzano l'energia. Inoltre, riempie un
//istogramma con i valori delle x ottenute campionando il modulo al quadrato
//della funzione d'onda. Le funzioni presenti in questo programma sono le
//stesse di quelle del programma esercizio2_parte1.cpp

//Random numbers
Random rnd;

//variables
int accepted, attempted;
int nstep, nblk;
double mu, sigma, x, delta;
TH1F histo("p(x)", "p(x)", 60, -3, 3);
vector<double> sum;
vector<double> sum2;

void Input(){
  ifstream Primes, Seed, ReadInput;
  //Read seed for random numbers
  int seed[4];
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  mu = 0.789132;
  sigma = 0.600949;
  nblk = 100;
  nstep = 1000;
  delta = 5.25;
  return;
}

double psi(double z){
  double esp1 = pow((double) (z-mu)/(sqrt(2)*sigma),2);
  double esp2 = pow((double) (z+mu)/(sqrt(2)*sigma),2);
  double y = exp(-esp1)+exp(-esp2);
  return y;
}

//sample the squared module of psi with Metropolis
void Move()
{
  double x_new = x + delta*rnd.Rannyu(-0.5, 0.5);
  double p = pow(psi(x_new),2) / pow(psi(x),2);
  attempted += 1;
  if(p>=1){
    x = x_new;
    accepted += 1;
  } else{
    double z = rnd.Rannyu();
    if(z<p){
      x = x_new;
      accepted += 1;
    }
  }
  histo.Fill(x);
}

//funzione che ha lo scopo di scartare le prime mosse del Metropolis per
//campionare x
void equilibrazione(){
  for(int i = 0; i < 1000; i++){
    Move();
  }
}

double V(){
  return (pow(x,2)-2.5)*pow(x,2);
}

//Calculate H*psi/psi
double E_loc(){
  double esp1 = pow((double) (x-mu)/(sqrt(2)*sigma),2);
  double esp2 = pow((double) (x+mu)/(sqrt(2)*sigma),2);
  double p1 = -exp(-esp2)/pow(sigma,2);
  double p2 = exp(-esp2)*pow((double) (x+mu)/(sigma*sigma),2);
  double p3 = -exp(-esp1)/pow(sigma,2);
  double p4 = exp(-esp1)*pow((double) (x-mu)/(sigma*sigma),2);
  double E = ((-0.5*(p1+p2+p3+p4))/psi(x)) + V();
  return E;
}

vector<double> Energia(){
    vector<double> vec(2);
    double E = 0;
    double somma = 0;
    double stddev = 0;
    double somma2 = 0;
    for(int i = 0; i < nblk; i++){
        for(int j = 0; j < nstep; j++){
            Move();
            E += E_loc();
        }
        E = (double) E / nstep;
        somma += E;
        sum.push_back(somma);
        somma2 += E*E;
        sum2.push_back(somma2);
        E = 0;
    }
    stddev = sqrt(fabs(somma2/(double)nblk - pow(somma/(double)nblk,2))/(double)nblk);
    vec[0] = (double) somma / nblk;
    vec[1] = stddev;
    return vec;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

int main(){
    TApplication app("app", 0, 0);
    Input();
    equilibrazione();
    vector<double> vec = Energia();
    cout << "Energia: " << vec[0] << " ± " << vec[1] << endl;
    cout << "Rate d'accettazione: " << (double) accepted/attempted << endl;
    //istogramma contenente i valori di x campionati
    TCanvas mycanvas("psi", "psi");
    histo.StatOverflows(kTRUE);
    histo.GetXaxis()->SetTitle("x");
    histo.GetYaxis()->SetTitle("conteggi");
    histo.Draw();
    
    //grafico di <H> in funzione del numero di blocchi
    TGraphErrors *graph = new TGraphErrors();
    for(int i = 0; i < nblk; i++){
      double media = (double) sum[i]/(i+1);
      double stddev = Error(sum[i], sum2[i], i+1);
      graph->SetPoint(i, (i+1), media);
      graph->SetPointError(i, 0.0, stddev);
    }
    TCanvas Energia("<H>", "<H>");
    graph->SetTitle("<H>");
    graph->GetXaxis()->SetTitle("N");
    graph->GetYaxis()->SetTitle("<H>");
    graph->Draw("ALP");

    app.Run();
    return 0;
}
