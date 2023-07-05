CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
esercizio1_parte1.exe: esercizio1_parte1.cpp
	$(CC) -o esercizio1_parte1.exe esercizio1_parte1.cpp `root-config --cflags` `root-config --libs`
esercizio1_parte2.exe: esercizio1_parte2.cpp
	$(CC) -o esercizio1_parte2.exe esercizio1_parte2.cpp `root-config --cflags` `root-config --libs`
esercizio1_parte3.exe: esercizio1_parte3.cpp
	$(CC) -o esercizio1_parte3.exe esercizio1_parte3.cpp `root-config --cflags` `root-config --libs`
compila: Monte_Carlo_ISING_1D.exe esercizio1_parte1.exe esercizio1_parte2.exe esercizio1_parte3.exe
clean :
	rm *.o Monte_Carlo_ISING_1D.exe
