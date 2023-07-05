esercizio2_parte1.o: esercizio2_parte1.cpp algoritmo_genetico.h random.h dna.h
	mpic++ -c esercizio2_parte1.cpp
random.o: random.h
	g++ -c random.cpp
esercizio2_parte1.exe: esercizio2_parte1.o algoritmo_genetico.h random.o
	mpic++ -o esercizio2_parte1.exe esercizio2_parte1.o random.o
esercizio2_parte2.exe: esercizio2_parte2.cpp
	g++ -o esercizio2_parte2.exe esercizio2_parte2.cpp `root-config --cflags` `root-config --libs`

compila: esercizio2_parte1.exe esercizio2_parte2.exe

esegui:
	mpirun --mca btl_vader_single_copy_mechanism none -np 6 esercizio2_parte1.exe