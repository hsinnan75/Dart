.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lstdc++
SOURCE		= main.cpp Mapping.cpp AlignmentCandidates.cpp GetData.cpp tools.cpp nw_alignment.cpp bwt_search.cpp bwt_index.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o dart $(LIB) htslib/libhts.a -lbz2 -llzma

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<
clean:
		rm -f *.o *~
