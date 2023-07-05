#ifndef __DNA__
#define __DNA__
#include <iostream>
#include <numeric>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <fstream>
#include "random.h"

using namespace std;

Random rnd;

class DNA{
public:
//costruttori
DNA(const vector<double>& x, const vector<double>& y);
DNA(const vector<double>& x, const vector<double>& y, vector<int> sequenza);
DNA(const DNA& copia);
//definizione operatore =
DNA& operator = (const DNA& copia);
//distruttore
~DNA();

//imposta i geni iniziali casualmente
void geni_iniziali();
//copia i geni a partire da un altro DNA
void copia_geni(const DNA& copia);
//crossover
void ricombina_geni(const DNA& Genitore1, const DNA& Genitore2);
//mutazione
void muta_geni(double prob);
//calcola il fitness
double calcola_fitness();
//restituisce il fitness
double get_fitness() const;
//restituisce i geni
vector<size_t> get_geni();
//metodo che permette di definire la sequenza di città
void set_sequenza(vector<size_t>);
//trasforma il DNA in una stringa
string to_stringa() const;
//calcola la lunghezza del percorso data una sequenza di DNA
double lunghezza_percorso(const vector<size_t>& sequenza);
//calcola la distanza tra due città
double dist(double x1, double y1, double x2, double y2);
//crossover
void crossover(const DNA& dna1, const DNA& dna2);
//cerca un gene in una sezione del DNA
bool trova_gene(size_t gene, size_t inizio, size_t fine) const;
//cambia due geni casualmente
void cambio_geni();
//funzione che ottimizza il DNA
void ricerca_locale();
//inverte una sequenza di geni
void scambio_ricerca_locale(const vector<size_t>& geni_input, vector<size_t>& geni_output,
                size_t gene1, size_t gene2);

private:
const vector<double>* X; //puntatore alle x delle città
const vector<double>* Y; //puntatore alle y delle città
vector<size_t> geni; //sequenza città
double fitness; //fitness sequenza
};

//costruttore
DNA::DNA(const vector<double>& x, const vector<double>& y){
    X = &x;
    Y = &y;
    geni.resize(x.size());
    //riempitura geni in ordine
    iota(geni.begin(), geni.end(), 0);
    //fitness iniziale impostato ad infinito
    fitness = INFINITY;
}

//costruttore
DNA::DNA(const vector<double>& x, const vector<double>& y, vector<int> sequenza){
    X = &x;
    Y = &y;
    geni.resize(x.size());
    for(int i = 0; i < sequenza.size(); i++){
        geni[i] = sequenza[i];
    }
    calcola_fitness();
}

//costruttore
DNA::DNA(const DNA& copia){
    X = copia.X;
    Y = copia.Y;
    fitness = copia.fitness;
    copia_geni(copia);
}

//operatore =
DNA& DNA::operator = (const DNA& copia){
    X = copia.X;
    Y = copia.Y;
    fitness = copia.fitness;
    copia_geni(copia);
}

//distruttore
DNA::~DNA(){}

// mischia casualmente i geni
void DNA::geni_iniziali(){
    size_t a = 0;
    size_t b = 0;
    for (size_t i = 0; i < geni.size(); i++){
        a = floor(rnd.Rannyu()*geni.size());       
        b = geni[i];
        geni[i] = geni[a];
        geni[a] = b;
    }
}

//copia i geni di un altro DNA
void DNA::copia_geni(const DNA& copia){
    geni.resize(X->size());
    for (size_t i = 0; i < geni.size(); i++)
        geni[i] = copia.geni[i];
}

//crossover
void DNA::ricombina_geni(const DNA& Genitore1, const DNA& Genitore2){
    crossover(Genitore1, Genitore2);
}

//mutazione geni e miglioramento DNA
void DNA::muta_geni(double prob_mutazione){
    if (rnd.Rannyu() <= prob_mutazione){
        cambio_geni();
        }
    ricerca_locale();
}

//calcola il fitness del DNA
double DNA::calcola_fitness(){
    fitness = lunghezza_percorso(geni);
    return fitness;
}

//restituisce il fitness
double DNA::get_fitness() const{
    return fitness;
}

//restituisce la sequenza del DNA
vector<size_t> DNA::get_geni(){
    return geni;
}

//permette di settare la sequenza delle città
void DNA::set_sequenza(vector<size_t> vec){
    for(int i = 0; i < vec.size(); i++){
        geni[i] = vec[i];
    }
    calcola_fitness();
}

//restituisce la sequenza del DNA sotto forma di stringa
string DNA::to_stringa() const{
    string str = " ";
    for (size_t i = 0; i < geni.size(); i++){
        str += to_string(geni[i]) + " ";
    }
    return str;
}

//calcola il fitness per la sequenza del DNA
double DNA::lunghezza_percorso(const vector<size_t>& sequenza){
    double lunghezza = 0.0;
    for (size_t i = 0; i < sequenza.size() - 1; i++){
        lunghezza += dist(X->at(sequenza[i]), Y->at(sequenza[i]), X->at(sequenza[i + 1]), Y->at(sequenza[i + 1]));
    }
    lunghezza += dist(X->at(sequenza[sequenza.size() - 1]), Y->at(sequenza[sequenza.size() - 1]), X->at(sequenza[0]),
        Y->at(sequenza[0]));
    return lunghezza;
}

//calcola la distanza tra 2 città
double DNA::dist(double x1, double y1, double x2, double y2){
    return sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));
}

//operatore di crossover: una sezione dei geni del dna1 viene copiato a questo
//dna e i rimanenti geni sono copiati a questo dna nell'ordine in cui appaiono
//nel dna2
void DNA::crossover(const DNA& dna1, const DNA& dna2){
    // inizio e fine casuale per copiare una sezione del dna1
    size_t inizio = 0;
    size_t fine = 0;
    while (inizio >= fine){
        inizio = floor(rnd.Rannyu()*dna1.geni.size());
        fine = floor(rnd.Rannyu()*dna1.geni.size());
    }
    //copia una sequenza dal dna1 a questo dna
    for (size_t i = inizio; i <= fine; i++){
        geni[i] = dna1.geni[i];
    }

    const size_t differenza = inizio - fine;

    //copia i geni del dna2 che non sono stati presi dal dna1. Inizia a copiare
    //dopo l'indice fine fino alla fine del dna2 e poi copia dall'inizio del dna2
    //fino all'indice fine del dna2
    vector<size_t> differenza_DNA;
    differenza_DNA.reserve(dna2.geni.size() - differenza);

    if (fine + 1 <= dna2.geni.size() - 1){
        for (size_t i = fine + 1; i < dna2.geni.size(); i++){
            if (!trova_gene(dna2.geni[i], inizio, fine)){
                differenza_DNA.push_back(dna2.geni[i]);
            }
        }
    }

    for (size_t i = 0; i <= fine; i++){
        if (!trova_gene(dna2.geni[i], inizio, fine)){
            differenza_DNA.push_back(dna2.geni[i]);
        }
    }

    //copia i geni rimanenti dal vettore differenza_DNA a questo dna iniziando ad
    //aggiungere i nuovi geni dall'indice (fine+1)
    size_t i = 0;
    if (fine + 1 <= dna2.geni.size() - 1){
        i = fine + 1;
    }

    for (size_t j = 0; j < differenza_DNA.size(); j++){
        geni[i] = differenza_DNA[j];
        i++;
        if (i > geni.size() - 1){
            i = 0;
        }
    }
}

//cerca gene in una sezione del DNA
bool DNA::trova_gene(size_t gene, size_t inizio, size_t fine) const{
    for (size_t i = inizio; i <= fine; i++){
        if (gene == geni[i]){
            return true;
        }
    }
    return false;
}

//scambia due geni del DNA
void DNA::cambio_geni(){
    const auto gene1 = floor(rnd.Rannyu()*geni.size());
    const auto gene2 = floor(rnd.Rannyu()*geni.size());
    const auto a = geni[gene1];
    geni[gene1] = geni[gene2];
    geni[gene2] = a;
}

//funzione che ottimizza il DNA scambiando la sequenza di geni se tale scambio
//migliora il fitness del DNA
void DNA::ricerca_locale(){
    bool migliorato = true;
    double fitness_scambio_geni = 0.0;
    vector<size_t> geni_scambiati(geni.size());

    while (migliorato){
        for (size_t i = 1; i < geni.size() - 1; i++){
            for (size_t j = i + 1; j < geni.size(); j++){
                scambio_ricerca_locale(geni, geni_scambiati, i, j);
                fitness_scambio_geni = lunghezza_percorso(geni_scambiati);

                if (fitness_scambio_geni < fitness){
                    geni.swap(geni_scambiati);
                    fitness = fitness_scambio_geni;
                    migliorato = true;
                    }else{
                        migliorato = false; 
                    }
            }
        }
    }
}

//funzione che inverte la sequenza di geni di un DNA dall'indice gene1 all'indice
//gene2
void DNA::scambio_ricerca_locale(const vector<size_t>& geni_input, vector<size_t>& geni_output,
                    size_t gene1, size_t gene2){
    //copia i geni dall'indice 0 all'indice (gene1-1)
    for (size_t i = 0; i <= gene1 - 1; i++){
        geni_output[i] = geni_input[i];
        }

    //inversione della sequenza da gene1 a gene2
    size_t contatore = 0;
    for (size_t i = gene1; i <= gene2; i++){
        geni_output[i] = geni_input[gene2 - contatore];
        contatore++;
    }

    //copia dei geni dall'indice (gene2+1) fino alla fine
    for (size_t i = gene2 + 1; i < geni_input.size(); i++){
        geni_output[i] = geni_input[i];
    }
}

#endif
