#ifndef __algoritmo_genetico__
#define __algoritmo_genetico__
#include "dna.h"

using namespace std;

//Random rnd;


template<typename DNA> class algoritmo_genetico{
    public:
    //costruttore
    algoritmo_genetico(const DNA* DNAiniziale, size_t dim_popolazione, 
                    double perc_candidati, double probricombinazione, 
                    double probmutazione);
    //metodo per passare alla prossima generazione
    void evoluzione();
    //metodo che restituisce la dimensione della popolazione
    size_t get_m() const;
    //metodo che restituisce il numero della generazione
    size_t get_generazione() const;
    //metodo che restituisce il fitness della migliore sequenza
    double get_fitness_migliore_DNA() const;
    //metodo che restituisce il fitness medio della metà migliore
    double get_fitness_medio() const;
    //metodo che restituisce la migliore sequenza sotto forma di stringa
    string get_migliore_DNA() const;
    //metodo che crea la popolazione iniziale
    void crea_popolazione_iniziale();
    //metodo che calcola il fitness di tutta la popolazione e salva il fitness 
    //migliore
    void fitness_popolazione();
    //funzione che calcola e salva il fitness medio della metà migliore della
    //popolazione
    void calcola_fitness_medio();
    //metodo che riempie il vettore scelta_DNA con gli indici dei DNA che
    //possono riprodursi
    void creazione_scelta_DNA();
    //selezione dei genitori per produrre i figli
    void selezione();
    //funzione che salva in un file i dati richiesti dal problema
    void print();

    private:
    const DNA* DNA_iniziale; //DNA iniziale
    size_t m; //dimensione popolazione
    size_t meta_m; //m diviso due
    double percentuale_candidati; //percentuale di candidati: solo una parte dei
                                  //DNA è selezionata per riprodursi
    double prob_ricombinazione; //probabilità crossover
    double prob_mutazione; //probabilità mutazione
    double fitness_medio; //fitness della metà migliore della popolazione
    vector<unique_ptr<DNA>> popolazione; //vettore che contiene la popolazione
    vector<unique_ptr<DNA>> nuova_generazione; //vettore della nuova generazione
    vector<size_t> scelta_DNA; //vettore che contiene i DNA che possono ricombinarsi
    DNA* genitore1; //DNA1 per la ricombinazione
    DNA* genitore2; //DNA2 per la ricombinazione
    DNA* migliore_DNA; //DNA con il migliore fitness
    size_t numero_evoluzioni; //numero di generazioni
};


//costruttore che inizializza tutte le variabili private e la popolazione
template<typename DNA> algoritmo_genetico <DNA>::algoritmo_genetico(const DNA* DNAiniziale, 
size_t dim_popolazione, double perc_candidati, double probricombinazione, 
double probmutazione){
    //inizializzazioni variabili
    DNA_iniziale = DNAiniziale;
    m = dim_popolazione;
    meta_m = m / 2;
    percentuale_candidati = perc_candidati;
    prob_ricombinazione = probricombinazione;
    prob_mutazione = probmutazione;
    fitness_medio = 0;
    scelta_DNA.reserve(m);
    genitore1 = nullptr;
    genitore2 = nullptr;
    migliore_DNA = nullptr;
    numero_evoluzioni = 0;
    //inizializzazione popolazione
    crea_popolazione_iniziale();
}

//funzione che restituisce la dimensione della popolazione
template<typename DNA> size_t algoritmo_genetico<DNA>::get_m() const{
    return m;
}

//funzione che restituisce il numero della generazione
template<typename DNA> size_t algoritmo_genetico<DNA>::get_generazione() const{
    return numero_evoluzioni;
}

//funzione che restituisce il fitness del migliore DNA
template<typename DNA> double algoritmo_genetico<DNA>::get_fitness_migliore_DNA() const{
    return migliore_DNA->get_fitness();
}

//
template<typename DNA> double algoritmo_genetico<DNA>::get_fitness_medio() const{
    return fitness_medio;
}

//funzione che restituisce il migliore DNA
template<typename DNA> string algoritmo_genetico<DNA>::get_migliore_DNA() const{
    return migliore_DNA->to_stringa();
}

//funzione che genera la popolazione iniziale
template<typename DNA> void algoritmo_genetico<DNA>::crea_popolazione_iniziale(){
    popolazione.resize(m);
    nuova_generazione.resize(m);
    for (size_t i = 0; i < m; i++){
        popolazione[i] = unique_ptr<DNA>( new DNA(*DNA_iniziale) );
        popolazione[i]->geni_iniziali();
        nuova_generazione[i] = unique_ptr<DNA>( new DNA(*DNA_iniziale) );
        nuova_generazione[i]->geni_iniziali();
    }
    fitness_popolazione();
    calcola_fitness_medio();
    numero_evoluzioni = 1;
}

//funzione che calcola il fitness della popolazione e salva il migliore
template<typename DNA> void algoritmo_genetico<DNA>::fitness_popolazione(){
    double fitness_minimo = INFINITY;
    double fitness = 0.0;
    for (size_t i = 0; i < m; i++){
            fitness = popolazione[i]->calcola_fitness();
            if (fitness < fitness_minimo){
                fitness_minimo = fitness;
                migliore_DNA = popolazione[i].get();
            }
        }
}

//funzione che calcola il fitness medio della metà migliore della popolazione
template<typename DNA> void algoritmo_genetico<DNA>::calcola_fitness_medio(){
    vector<double> vec(m);
    for(int i = 0; i < m; i++){
        vec[i] = popolazione[i]->get_fitness();
    }
    sort(vec.begin(), vec.end());
    double media = 0;
    for(int i = 0; i < meta_m; i++){
        media += vec[i];
    }
    fitness_medio = media /(double) meta_m;
}

//funzione che fa passare alla nuova generazione
template<typename DNA> void algoritmo_genetico<DNA>::evoluzione(){
    //creazione del vettore con gli indici dei DNA che possono essere genitori
    creazione_scelta_DNA();
    //ciclo che crea la nuova generazione
    for (size_t i = 0; i < meta_m; i++){
        //selezione dei genitori
        selezione();
        //creazione figli
        auto figlio1 = nuova_generazione[i].get();
        auto figlio2 = nuova_generazione[i + meta_m].get();
        if (rnd.Rannyu() <= prob_ricombinazione){
            figlio1->ricombina_geni(*genitore1, *genitore2);
            figlio2->ricombina_geni(*genitore2, *genitore1);
        }
        else{
            figlio1->copia_geni(*genitore1);
            figlio2->copia_geni(*genitore2);
        }
        //mutazione figli
        figlio1->muta_geni(prob_mutazione);
        figlio2->muta_geni(prob_mutazione);
    }
    popolazione.swap(nuova_generazione);
    fitness_popolazione();
    calcola_fitness_medio();
    numero_evoluzioni++;
}


template<typename DNA> void algoritmo_genetico<DNA>::creazione_scelta_DNA(){
    //vettore di coppie: una copia rappresenta rispettivamente l'indice del DNA
    //nella popolazione e il suo fitness 
    vector<pair<size_t, double>> coppia(m);
    for (size_t i = 0; i < m; i++){
        coppia[i] = make_pair(i, popolazione[i]->get_fitness());
    }

    //funzione lambda che ordina il vettore di coppie in base al fitness
    sort(coppia.begin(), coppia.end(),
            [](pair<size_t, double>& a, pair<size_t, double>& b){
                return a.second < b.second;
            });
    
    //reset vecchio vettore
    scelta_DNA.clear();
    
    const size_t numero_candidati = static_cast<size_t>(coppia.size()*(percentuale_candidati*0.01));
    const double q = 1.0 / numero_candidati;
    double prob_inserimento = 0.0;

    for (size_t i = 0; i < numero_candidati; i++){
        //la probabilità di inserire un DNA in scelta_DNA è più alta per
        //DNA con fitness migliori
        prob_inserimento = 1.0 - i*q;
        prob_inserimento = pow(prob_inserimento, 2.0);
        if (rnd.Rannyu() <= prob_inserimento){
            scelta_DNA.push_back(coppia[i].first);
        }
    }
}

//funzione che seleziona i due genitori
template<typename DNA> void algoritmo_genetico<DNA>::selezione(){
    const auto indice1 = floor(rnd.Rannyu()*scelta_DNA.size());
    const auto indice2 = floor(rnd.Rannyu()*scelta_DNA.size());
    const auto ind1 = scelta_DNA[indice1];
    const auto ind2 = scelta_DNA[indice2];
    genitore1 = popolazione[ind1].get();
    genitore2 = popolazione[ind2].get();
}

template<typename DNA> void algoritmo_genetico<DNA>::print(){
    ofstream fit_migliore;
    ofstream fit_medio;
    fit_migliore.open("best_fitness.dat", ios::app);
    fit_medio.open("average_fitness.dat", ios::app);
    fit_migliore << get_fitness_migliore_DNA() << endl;
    fit_medio << fitness_medio << endl;
    fit_migliore.close();
    fit_medio.close();
}

#endif
