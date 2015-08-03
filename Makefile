CC=			gcc
CFLAGS=		-g -Wall -Wno-unused-function -O2
PROG=		naivepca
INCLUDES=	
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

naivepca:eigen.o naivepca.o
		$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

naivepca.o: ksort.h kseq.h
