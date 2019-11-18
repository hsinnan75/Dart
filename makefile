.KEEP_STAT:

all: dart bwt_index

htslib:
		mkdir -p bin
		$(MAKE) -C src/htslib libhts.a

dart: htslib
		$(MAKE) -C src && mv -f src/$@ bin

bwt_index:
		$(MAKE) -C src/BWT_Index && mv -f src/BWT_Index/$@ bin/

clean:
		rm -f bin/dart bin/bwt_index
		$(MAKE) clean -C src
		$(MAKE) clean -C src/htslib
		$(MAKE) clean -C src/BWT_Index
