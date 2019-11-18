.KEEP_STAT:

all: main

CXX		= g++
CFLAGS		= -I${PREFIX}/include -L${PREFIX}/lib -Wall -D NDEBUG -O3 -m64
LDFLAGS		= -L${PREFIX}/lib -lz -lm -lpthread -lstdc++
SOURCE		= main.cpp Mapping.cpp AlignmentCandidates.cpp GetData.cpp tools.cpp nw_alignment.cpp bwt_search.cpp bwt_index.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

main:		$(OBJECT)
			$(CXX) $(CFLAGS) $(OBJECT) -o dart $(LDFLAGS) htslib/libhts.a -lbz2 -llzma

%.o:		%.cpp $(HEADER)
			$(CXX) $(CFLAGS) -c $<
clean:
		rm -f *.o *~
