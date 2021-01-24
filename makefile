# MAX. K-MER LENGTH, NUMBER OF FILES
CC = g++ -O3 -march=native -DmaxK=32 -DmaxN=64

## IF BIFROST LIBRARY SHOULD BE USED
# CC = g++ -O3 -march=native -DmaxK=32 -DmaxN=64 -DuseBF
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
	$(RM)
	$(CC) -o SANS main.o graph.o kmer32.o kmerXX.o color64.o colorXX.o util.o $(BF)
	$(MK)
	$(MV)

main.o: src/main.cpp src/main.h graph.o util.o
	$(CC) -c src/main.cpp

graph.o: src/graph.cpp src/graph.h kmer32.o kmerXX.o color64.o colorXX.o
	$(CC) -c src/graph.cpp

kmer32.o: src/kmer32.cpp src/kmer32.h
	$(CC) -c src/kmer32.cpp

kmerXX.o: src/kmerXX.cpp src/kmerXX.h
	$(CC) -c src/kmerXX.cpp

color64.o: src/color64.cpp src/color64.h
	$(CC) -c src/color64.cpp

colorXX.o: src/colorXX.cpp src/colorXX.h
	$(CC) -c src/colorXX.cpp

util.o: src/util.cpp src/util.h
	$(CC) -c src/util.cpp

.PHONY: clean
clean:
	$(RM)