# MAX. K-MER LENGTH, NUMBER OF FILES
CC = g++ -O3 -flto=auto -march=native -DmaxK=32 -DmaxN=64
XX = -lpthread

## IF BIFROST LIBRARY SHOULD BE USED
# CC = g++ -O3 -flto=auto -march=native -DmaxK=32 -DmaxN=64 -DuseBF
# XX = -lbifrost -lpthread -lz

SANS: obj/ obj/main.o
	$(CC) obj/main.o obj/index.o obj/graph.o obj/tree.o obj/kmer.o obj/color.o obj/util.o -o SANS $(XX)

obj/main.o: src/main.cpp src/main.h obj/index.o obj/graph.o obj/util.o
	$(CC) -c src/main.cpp -o obj/main.o

obj/graph.o: src/graph.cpp src/graph.h obj/tree.o obj/kmer.o
	$(CC) -c src/graph.cpp -o obj/graph.o

obj/tree.o: src/tree.cpp src/tree.h obj/color.o
	$(CC) -c src/tree.cpp -o obj/tree.o

obj/index.o: src/index.cpp src/index.h
	$(CC) -c src/index.cpp -o obj/index.o

obj/kmer.o: src/kmer.cpp src/kmer.h src/byte.h
	$(CC) -c src/kmer.cpp -o obj/kmer.o

obj/color.o: src/color.cpp src/color.h src/byte.h
	$(CC) -c src/color.cpp -o obj/color.o

obj/util.o: src/util.cpp src/util.h
	$(CC) -c src/util.cpp -o obj/util.o

obj/: makefile src/ansi.h src/*/*.h
	@rm -rf obj/ && mkdir obj/
