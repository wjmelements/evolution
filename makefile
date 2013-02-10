CC=g++
CFLAGS=-O3

all: depend evolution microscope
depend: dna env org
env:
	$(CC) $(CFLAGS) -c src/Environment.cpp -o lib/Environment.o
dna:
	$(CC) $(CFLAGS) -c src/DNA.cpp -o lib/DNA.o
org:
	$(CC) $(CFLAGS) -c src/Organism.cpp -o lib/Organism.o
evo: depend evolution
evolution:
	$(CC) $(CFLAGS) lib/*.o main.cpp -o bin/evolution
mic: depend microscope
microscope:
	$(CC) $(CFLAGS) lib/*.o Microscope.cpp -o bin/microscope
run:
	bin/evolution
data:
	bin/microscope
clean:
	rm lib/*.o
	rm bin/*
