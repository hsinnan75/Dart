.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread
SOURCE		= main.cpp Mapping.cpp AlignmentCandidates.cpp GetData.cpp tools.cpp nw_alignment.cpp bwt_search.cpp bwt_index.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main index

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o dart $(LIB)

static:		$(OBJECT)
			$(Compiler) $(FLAGS) -static $(OBJECT) -o dart_static $(LIB)
                        
%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .

clean:
		rm -f *.o *~ && make clean -C BWT_Index
