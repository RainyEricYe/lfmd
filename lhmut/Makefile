CC= g++ -std=c++11 -Wall -O3

INCLUDE= -I ${LDPATH}/include ${ALGLIB}/src/*.o -I ${ALGLIB}/src
LIBS= -L ${LDPATH}/lib -L.
OBJS= main.o likelihood.o option.o

lhmut: lhmut.cc option.h $(OBJS) Makefile
	$(CC) $< -o $@ $(OBJS) $(INCLUDE) $(LIBS)

main.o: main.cc main.h option.h
	$(CC) -c $< -o $@

likelihood.o: likelihood.cc likelihood.h main.h option.h
	$(CC) -c $< $(INCLUDE) -o $@

option.o: option.cc option.h
	$(CC) -c $< -o $@

.PHONY: clean

clean:
	-rm lhmut $(OBJS)
