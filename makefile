.KEEP_STAT:

all: bwt_index dart

htslib:
		mkdir -p bin/
		$(MAKE) -C src/htslib libhts.a

bwt_index:
		$(MAKE) -C src/BWT_Index && mv -f src/BWT_Index/$@ bin/

dart: htslib bwt_index
		$(MAKE) -C src && mv -f src/$@ bin

clean:
		rm -f bin/dart bin/bwt_index
		$(MAKE) clean -C src
		$(MAKE) clean -C src/htslib
		$(MAKE) clean -C src/BWT_Index
