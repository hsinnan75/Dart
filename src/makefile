.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -I${PREFIX}/include -L${PREFIX}/lib -Wall -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lstdc++
SOURCE		= main.cpp Mapping.cpp AlignmentCandidates.cpp GetData.cpp tools.cpp nw_alignment.cpp bwt_search.cpp bwt_index.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o dart $(LIB) BWT_Index/libbwa.a htslib/libhts.a -lbz2 -llzma

clean:
		rm -f *.o *~
