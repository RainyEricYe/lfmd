CC= g++ --std=c++11 -O3
INCLUDE= -I ${SEQLIB}/ -I ${SEQLIB}/htslib/ ${SEQLIB}/lib/libseqlib.a ${SEQLIB}/lib/libhts.a ${SEQLIB}/lib/libbwa.a ${SEQLIB}/lib/libfml.a
LIBS= -llzma -lbz2 -L. -lz -lpthread

bamRdf: bamRdf.cc Makefile
	$(CC) $< $(INCLUDE) $(LIBS) -o $@

.PHONY: clean

clean:
	-rm bamRdf
