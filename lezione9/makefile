esercizio1_parte1.o: esercizio1_parte1.cpp algoritmo_genetico.h random.h dna.h
	g++ -c esercizio1_parte1.cpp
random.o: random.h
	g++ -c random.cpp
esercizio1_parte1.exe: esercizio1_parte1.o algoritmo_genetico.h random.o
	g++ -o esercizio1_parte1.exe esercizio1_parte1.o random.o
esercizio1_parte2.exe: esercizio1_parte2.cpp
	g++ -o esercizio1_parte2.exe esercizio1_parte2.cpp `root-config --cflags` `root-config --libs`

compila: esercizio1_parte1.exe esercizio1_parte2.exe