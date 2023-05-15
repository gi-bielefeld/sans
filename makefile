# MAX. K-MER LENGTH, NUMBER OF FILES
CC = g++ -O3 -march=native -DmaxK=32 -DmaxN=64 -std=c++14
BF = -lpthread

## IF DEBUG
# CC = g++ -g -march=native -DmaxK=33 -DmaxN=64 -std=c++14

## IF BIFROST LIBRARY SHOULD BE USED
# CC = g++ -O3 -march=native -DmaxK=64 -DmaxN=64 -DuseBF -std=c++14
# BF = -lbifrost -lpthread

# GZ STREAM LIB
CFLAGS = gcc -O3 -march=native

# Wrap Windows / Unix commands
ifeq ($(OS), Windows_NT)
	TD = obj
	MK = rmdir /s /q obj && mkdir obj
	RM = rmdir /s /q obj
	MV = cmd /C move *.o obj
	CP = cp makefile obj
else
	TD = obj/
	MK = mkdir -p obj/
	RM = rm -rf obj/
	MV = mv *.o obj/
	CP = cp makefile obj/makefile
endif


ifeq ("$(wildcard $(TD))", "")
    RM = @echo ""
endif

SANS: start makefile obj/main.o done
	$(CC) -o SANS obj/main.o obj/graph.o obj/index.o obj/kmer.o obj/kmerAmino.o obj/color.o obj/util.o obj/translator.o obj/cleanliness.o obj/gzstream.o -lz $(BF)


obj/main.o: makefile src/main.cpp src/main.h obj/translator.o obj/graph.o obj/util.o obj/cleanliness.o obj/gzstream.o
	$(CC) -c src/main.cpp

obj/graph.o: makefile src/graph.cpp src/graph.h obj/index.o obj/kmer.o obj/kmerAmino.o obj/color.o
	$(CC) -c src/graph.cpp

obj/index.o: makefile src/index.cpp src/index.h obj/kmer.o obj/color.o
	$(CC) -c src/index.cpp

obj/kmer.o: makefile src/kmer.cpp src/kmer.h
	$(CC) -c src/kmer.cpp

obj/kmerAmino.o: makefile src/kmerAmino.cpp src/kmerAmino.h obj/util.o
	$(CC) -c src/kmerAmino.cpp

obj/color.o: makefile src/color.cpp src/color.h
	$(CC) -c src/color.cpp

obj/util.o: makefile src/util.cpp src/util.h
	$(CC) -c src/util.cpp

obj/translator.o:  makefile src/translator.cpp src/translator.h src/gc.h
	$(CC) -c src/translator.cpp

obj/cleanliness.o: makefile src/cleanliness.cpp src/cleanliness.h
	$(CC) -c src/cleanliness.cpp

obj/gzstream.o: makefile src/gz/gzstream.C src/gz/gzstream.h	
	$(CFLAGS) -c src/gz/gzstream.C

# [Internal rules]

# Print info at compile start
start:
	@echo "";
	@echo "   ________________________________";
	@echo "     <<< BUILDING SANS SERIF >>>   ";
	@echo "   ________________________________";
	@echo "";
	$(MK)

# Print info when done
done:
	$(MV)
	@echo "";
	@echo "   _______________";
	@echo "    <<< Done! >>> ";
	@echo "   _______________";
	@echo "";


.PHONY: clean

# Remove build files
clean:
	$(RM)



