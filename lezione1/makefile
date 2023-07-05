random.o: random.h random.cpp
	g++ -c random.cpp

esercizio2.o: random.h esercizio2.cpp
	g++ -c esercizio2.cpp `root-config --cflags` `root-config --libs`

esercizio2.exe: random.o esercizio2.o random.h
	g++ -o esercizio2.exe esercizio2.o random.o `root-config --cflags` `root-config --libs`

esercizio3.o: random.h esercizio3.cpp
	g++ -c esercizio3.cpp `root-config --cflags` `root-config --libs`

esercizio3.exe: random.o esercizio3.o random.h
	g++ -o esercizio3.exe esercizio3.o random.o `root-config --cflags` `root-config --libs`

esercizio1.o: random.h esercizio1.cpp
	g++ -c esercizio1.cpp `root-config --cflags` `root-config --libs`

esercizio1.exe: random.o esercizio1.o random.h
	g++ -o esercizio1.exe esercizio1.o random.o `root-config --cflags` `root-config --libs`

compila: esercizio2.exe esercizio3.exe esercizio1.exe