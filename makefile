# MAX. K-MER LENGTH, NUMBER OF FILES
CC = g++ -O3 -march=native -DmaxK=64 -DmaxN=64 -std=c++17
XX = -lpthread -lz

## IF DEBUG
# CC = g++ -g -march=native -DmaxK=32 -DmaxN=64 -std=c++14

## IF BIFROST LIBRARY SHOULD BE USED
# CC = g++ -O3 -march=native -DmaxK=32 -DmaxN=64 -DuseBF -std=c++14
# XX = -lbifrost -lpthread -lz

# GZ STREAM LIB
CFLAGS = gcc -O3 -march=native

# Directories
SRCDIR		:= src
BUILDDIR 	:= obj


# Wrap Windows / Unix commands
ifeq ($(OS), Windows_NT)
	TD = $(BUILDDIR)
	MK = rmdir /s /q $(BUILDDIR) && mkdir $(BUILDDIR)
	RM = rmdir /s /q $(BUILDDIR)
	MV = cmd /C move *.o $(BUILDDIR)
	CP = cp makefile $(BUILDDIR)
else
	TD = $(BUILDDIR)/
	MK = mkdir -p $(BUILDDIR)/
	RM = rm -rf $(BUILDDIR)/
	MV = mv *.o $(BUILDDIR)/
	CP = cp makefile $(BUILDDIR)/makefile
endif


ifeq ("$(wildcard $(TD))", "")
    RM = @echo ""
endif

all: makefile start SANS done

SANS: makefile $(BUILDDIR)/main.o
	$(CC) -o SANS $(BUILDDIR)/nexus_color.o $(BUILDDIR)/main.o $(BUILDDIR)/graph.o $(BUILDDIR)/kmer.o $(BUILDDIR)/kmerAmino.o $(BUILDDIR)/color.o $(BUILDDIR)/util.o $(BUILDDIR)/translator.o $(BUILDDIR)/cleanliness.o $(BUILDDIR)/gzstream.o  $(BUILDDIR)/PCTree_basic.o $(BUILDDIR)/PCTree_construction.o $(BUILDDIR)/PCTreeForest.o $(BUILDDIR)/PCTree_restriction.o $(BUILDDIR)/PCTree_intersect.o $(BUILDDIR)/PCNode.o $(XX)

$(BUILDDIR)/main.o: makefile $(SRCDIR)/main.cpp $(SRCDIR)/main.h $(BUILDDIR)/color.o $(BUILDDIR)/translator.o $(BUILDDIR)/graph.o $(BUILDDIR)/util.o $(BUILDDIR)/cleanliness.o $(BUILDDIR)/gzstream.o $(BUILDDIR)/nexus_color.o $(BUILDDIR)/PCTree_construction.o $(BUILDDIR)/PCTree_basic.o $(BUILDDIR)/PCTreeForest.o $(BUILDDIR)/PCTree_restriction.o $(BUILDDIR)/PCTree_intersect.o $(BUILDDIR)/PCNode.o
	$(CC) -c $(SRCDIR)/main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/graph.o: makefile $(SRCDIR)/graph.cpp $(SRCDIR)/graph.h $(BUILDDIR)/kmer.o $(BUILDDIR)/kmerAmino.o $(BUILDDIR)/color.o
	$(CC) -c $(SRCDIR)/graph.cpp -o $(BUILDDIR)/graph.o

$(BUILDDIR)/kmer.o: makefile $(SRCDIR)/kmer.cpp $(SRCDIR)/kmer.h $(BUILDDIR)/util.o
	$(CC) -c $(SRCDIR)/kmer.cpp -o $(BUILDDIR)/kmer.o

$(BUILDDIR)/kmerAmino.o: makefile $(SRCDIR)/kmerAmino.cpp $(SRCDIR)/kmerAmino.h $(BUILDDIR)/util.o
	$(CC) -c $(SRCDIR)/kmerAmino.cpp -o $(BUILDDIR)/kmerAmino.o

$(BUILDDIR)/color.o: makefile $(SRCDIR)/color.cpp $(SRCDIR)/color.h
	$(CC) -c $(SRCDIR)/color.cpp -o $(BUILDDIR)/color.o
	
$(BUILDDIR)/nexus_color.o: $(SRCDIR)/nexus_color.cpp $(SRCDIR)/nexus_color.h
	$(CC) -c $(SRCDIR)/nexus_color.cpp -o $(BUILDDIR)/nexus_color.o

$(BUILDDIR)/util.o: $(SRCDIR)/util.cpp $(SRCDIR)/util.h
	$(CC) -c $(SRCDIR)/util.cpp -o $(BUILDDIR)/util.o

$(BUILDDIR)/translator.o: $(SRCDIR)/translator.cpp $(SRCDIR)/translator.h $(SRCDIR)/gc.h
	$(CC) -c $(SRCDIR)/translator.cpp -o $(BUILDDIR)/translator.o

$(BUILDDIR)/cleanliness.o: $(SRCDIR)/cleanliness.cpp $(SRCDIR)/cleanliness.h
	$(CC) -c $(SRCDIR)/cleanliness.cpp -o $(BUILDDIR)/cleanliness.o

$(BUILDDIR)/gzstream.o: $(SRCDIR)/gz/gzstream.C $(SRCDIR)/gz/gzstream.h	
	$(CFLAGS) -c $(SRCDIR)/gz/gzstream.C  -o $(BUILDDIR)/gzstream.o


# PC-Tree

	
$(BUILDDIR)/PCNode.o: $(SRCDIR)/pctree/PCNode.cpp $(SRCDIR)/pctree/PCNode.h $(BUILDDIR)/PCTreeForest.o
	$(CC) -c $(SRCDIR)/pctree/PCNode.cpp -o $(BUILDDIR)/PCNode.o
	
$(BUILDDIR)/PCTree_basic.o: $(SRCDIR)/pctree/PCTree_basic.cpp $(BUILDDIR)/PCNode.o 
	$(CC) -c $(SRCDIR)/pctree/PCTree_basic.cpp -o $(BUILDDIR)/PCTree_basic.o
	
$(BUILDDIR)/PCTree_construction.o: $(SRCDIR)/pctree/PCTree_construction.cpp $(BUILDDIR)/PCNode.o 
	$(CC) -c $(SRCDIR)/pctree/PCTree_construction.cpp -o $(BUILDDIR)/PCTree_construction.o
	
$(BUILDDIR)/PCTreeForest.o: $(SRCDIR)/pctree/PCTreeForest.cpp $(SRCDIR)/pctree/PCTreeForest.h
	$(CC) -c $(SRCDIR)/pctree/PCTreeForest.cpp -o $(BUILDDIR)/PCTreeForest.o
	
$(BUILDDIR)/PCTree_intersect.o: $(SRCDIR)/pctree/PCTree_intersect.cpp
	$(CC) -c $(SRCDIR)/pctree/PCTree_intersect.cpp -o $(BUILDDIR)/PCTree_intersect.o
		
$(BUILDDIR)/PCTree_restriction.o: $(SRCDIR)/pctree/PCTree_restriction.cpp $(BUILDDIR)/PCNode.o
	$(CC) -c $(SRCDIR)/pctree/PCTree_restriction.cpp -o $(BUILDDIR)/PCTree_restriction.o
	


# [Internal rules]

# Print info at compile start
start:
	@echo "";
	@echo "   ________________________________";
	@echo "     <<< BUILDING SANS AMBAGES >>>   ";
	@echo "   ________________________________";
	@echo "";
	@$(MK)

# Print info when done
done:
	@echo "";
	@echo "   _______________";
	@echo "    <<< Done! >>> ";
	@echo "   _______________";
	@echo "";


# [Remove current build files]

.PHONY: clean

# Remove build files
clean:
	$(RM)
