.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		=
SOURCE		= main.cpp Mapping.cpp AlignmentCandidates.cpp GetData.cpp tools.cpp PairwiseAlignment.cpp bwt_search.cpp bwt_index.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -lpthread -o kart-rna $(LIB)
%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<
clean:
		rm -f *.o *~
eva:
		$(Compiler) $(FLAGS) SamEvaulation.cpp -o eva
eva2:
		$(Compiler) $(FLAGS) SamEvaulation2.cpp -o eva2
