# MAX. K-MER LENGTH, NUMBER OF FILES
 CC = g++ -O3 -march=native -DmaxK=32 -DmaxN=5000

## IF DEBUG
# CC = g++ -g -march=native -DmaxK=33 -DmaxN=1024

## IF BIFROST LIBRARY SHOULD BE USED
# CC = g++ -O3 -march=native -DmaxK=32 -DmaxN=1024 -DuseBF
# BF = -lbifrost -lpthread -lz

# Wrap Windows / Unix commands
ifeq ($(OS), Windows_NT)
	TD = obj
	MK = mkdir obj
	RM = rmdir /s /q obj 
	MV = cmd /C move *.o obj

else
	TD = obj/
	MK = mkdir obj/
	RM = rm -rf obj/
	MV = mv *.o obj/
endif

ifeq ("$(wildcard $(TD))", "")
    RM = @echo ""
endif

SANS: main.o
	$(CC) -o SANS main.o graph.o kmer32.o kmerXX.o kmerAminoXX.o kmerAmino12.o color64.o colorXX.o util.o translator.o cleanliness.o $(BF)
	$(RM)
	$(MK)
	$(MV)

main.o: src/main.cpp src/main.h translator.o graph.o util.o cleanliness.o
	$(CC) -c src/main.cpp

graph.o: src/graph.cpp src/graph.h kmer32.o kmerXX.o kmerAmino12.o kmerAminoXX.o color64.o colorXX.o
	$(CC) -c src/graph.cpp

kmer32.o: src/kmer32.cpp src/kmer32.h
	$(CC) -c src/kmer32.cpp

kmerXX.o: src/kmerXX.cpp src/kmerXX.h
	$(CC) -c src/kmerXX.cpp

kmerAmino12.o: src/kmerAmino12.cpp src/kmerAmino12.h util.o
	$(CC) -c src/kmerAmino12.cpp

kmerAminoXX.o: src/kmerAminoXX.cpp src/kmerAminoXX.h util.o
	$(CC) -c src/kmerAminoXX.cpp

color64.o: src/color64.cpp src/color64.h
	$(CC) -c src/color64.cpp

colorXX.o: src/colorXX.cpp src/colorXX.h
	$(CC) -c src/colorXX.cpp

util.o: src/util.cpp src/util.h
	$(CC) -c src/util.cpp

translator.o: src/translator.cpp src/translator.h src/gc.h
	$(CC) -c src/translator.cpp

cleanliness.o: src/cleanliness.cpp src/cleanliness.h
	$(CC) -c src/cleanliness.cpp

.PHONY: clean
clean:
	$(RM)
