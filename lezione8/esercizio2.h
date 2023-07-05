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

//Random numbers
Random rnd;

//variables
int accepted, attempted;
double best_E = INFINITY; //stima energia più bassa
double best_mu = 0; //stima mu per la quale si ottiene l'energia più bassa
double best_sigma = 0; //stima sigma per la quale si ottiene l'energia più bassa
double best_err = 0; //errore per l'energia più bassa ottenuta
double energy = 0; //energia ad ogni step del simulated annealing
double errore = 0; //errore energia ad ogni step del simulated annealing
int nstep, nblk; //numero di lanci per ciascun blocco e numero di blocchi
double mu, sigma, x, delta, delta2;

//pigreco
const double pi=3.1415927;

//Funzione che inizializza il generatore e le variabili mu, sigma, x, delta,
//delta2, nblk e nstep prendendole dal file input.dat
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

  //read data from file input.dat
  ReadInput.open("input.dat");

  ReadInput >> mu >> sigma >> x >> delta >> nblk >> nstep >> delta2;
  ReadInput.close();
  return;
}

//funzione che calcola il valore di psi per una determinata coordinata
double psi(double z){
  double esp1 = pow((double) (z-mu)/(sqrt(2)*sigma),2);
  double esp2 = pow((double) (z+mu)/(sqrt(2)*sigma),2);
  double y = exp(-esp1)+exp(-esp2);
  return y;
}

//campionamento del modulo al quadrato di psi attraverso metropolis usando
//una probabilità di transizione uniforme
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
}

//funzione che ha lo scopo di scartare le prime mosse del Metropolis per
//campionare x
void equilibrazione(){
  for(int i = 0; i < 1000; i++){
    Move();
  }
}

//funzione che restituisce il valore del potenziale per la coordinata x
double V(){
  return (pow(x,2)-2.5)*pow(x,2);
}

//funzione che calcola un singolo valore dell'energia locale, 
//ossia (H applicato a psi)/psi
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

//funzione che calcola e restituisce la media e la deviazione standard della 
//media dell'energia locale campionando x con la funzione Move() per determinati
//valori di sigma e mu
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
        somma2 += E*E;
        E = 0;
    }
    stddev = sqrt(fabs(somma2/(double)nblk - pow(somma/(double)nblk,2))/(double)nblk);
    vec[0] = (double) somma / nblk;
    vec[1] = stddev;
    //salva il valor migliore trovato dell'energia e i migliori valori di
    //mu e sigma
    if(vec[0] <= best_E){
      best_E = vec[0];
      best_mu = mu;
      best_sigma = sigma;
      best_err = vec[1];
    }
    accepted = 0;
    attempted = 0;
    return vec;
}

//funzione che compie uno step del simulated annealing per una certa temperatura
//T utilizzando Metropolis
void Move2(double T){
  double beta = 1.0/T;
  double old_mu = mu;
  double old_sigma = sigma;
  mu = mu + delta2*rnd.Rannyu(-0.05,0.05);
  sigma = sigma + delta2*rnd.Rannyu(-0.05,0.05);
  equilibrazione();
  vector<double> vec = Energia();
  double new_energy = vec[0];
  double p = exp(-beta*(new_energy-energy));
  if(rnd.Rannyu()<p){
    energy=new_energy;
    errore = vec[1];
  } else{
    mu = old_mu;
    sigma = old_sigma;
  }
}