CC = g++
CFLAGS = -Wall -O3 --std=c++11


esercizio2.exe : esercizio2.cpp
	$(CC) esercizio2.cpp -o esercizio2.exe `root-config --cflags` `root-config --libs`
NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

compila: NVE_NVT.exe esercizio2.exe
clean :
	rm *.o NVE_NVT.exe
