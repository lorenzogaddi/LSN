random.o: random.h random.cpp
	g++ -c random.cpp

esercizio2_parte1.o: random.h esercizio2.h esercizio2_parte1.cpp
	g++ -c esercizio2_parte1.cpp `root-config --cflags` `root-config --libs`

esercizio2_parte1.exe: random.o esercizio2_parte1.o esercizio2.h random.h
	g++ -o esercizio2_parte1.exe esercizio2_parte1.o random.o `root-config --cflags` `root-config --libs`

esercizio2_parte2.o: random.h esercizio2_parte2.cpp
	g++ -c esercizio2_parte2.cpp `root-config --cflags` `root-config --libs`

esercizio2_parte2.exe: random.o esercizio2_parte2.o random.h
	g++ -o esercizio2_parte2.exe esercizio2_parte2.o random.o `root-config --cflags` `root-config --libs`

compila: esercizio2_parte1.exe esercizio2_parte2.exe