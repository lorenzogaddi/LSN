random.o: random.h random.cpp
	g++ -c random.cpp

esercizio1.o: random.h esercizio1.cpp
	g++ -c esercizio1.cpp `root-config --cflags` `root-config --libs`

esercizio1.exe: random.o esercizio1.o random.h
	g++ -o esercizio1.exe esercizio1.o random.o `root-config --cflags` `root-config --libs`

compila: esercizio1.exe