CC= g++ --std=c++11 -O3 -Wall
INCLUDE= -I ${SEQLIB}/ -I ${SEQLIB}/htslib/ ${SEQLIB}/lib/libseqlib.a ${SEQLIB}/lib/libhts.a ${SEQLIB}/lib/libbwa.a ${SEQLIB}/lib/libfml.a -I ${ALGLIB}/src ${ALGLIB}/src/*.o -I ${LDPATH}/include  -L ${LDPATH}/lib
LIBS= -llzma -lbz2 -L. -lgzstream -lz -lpthread

bamDCS:  main.cc main.h Makefile
	$(CC) $< -o $@ $(INCLUDE) $(LIBS)
compass_search.o:

.PHONY: clean

clean:
	-rm bamDCS
