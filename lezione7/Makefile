CC = g++
CFLAGS = -Wall -O3 --std=c++11


esercizio2_parte1.exe : esercizio2_parte1.cpp
	$(CC) esercizio2_parte1.cpp -o esercizio2_parte1.exe `root-config --cflags` `root-config --libs`
esercizio2_parte2.exe : esercizio2_parte2.cpp
	$(CC) esercizio2_parte2.cpp -o esercizio2_parte2.exe `root-config --cflags` `root-config --libs`
NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
esercizio4_parte1.exe : esercizio4_parte1.cpp
	$(CC) esercizio4_parte1.cpp -o esercizio4_parte1.exe `root-config --cflags` `root-config --libs`
esercizio4_parte2.exe : esercizio4_parte2.cpp
	$(CC) esercizio4_parte2.cpp -o esercizio4_parte2.exe `root-config --cflags` `root-config --libs`

compila: NVE_NVT.exe esercizio2_parte1.exe esercizio2_parte2.exe esercizio4_parte1.exe esercizio4_parte2.exe
clean :
	rm *.o NVE_NVT.exe
