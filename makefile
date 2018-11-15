CC      =  gcc
CFLAGS  =  -g -Wall -D VERBOSE -D PRINT_COVERAGE #-O2  
LDFLAGS = -lz

#OBJS = gfa.o opt.o paf.o sdict.o eg.o 
PROG = scaff_10x  scaff_hic

.SUFFIXS:.c .o

all:$(PROG)

scaff_hic: bamlite.o bed.o cdict.o graph.o scaff_hic.o sdict.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

scaff_10x: bamlite.o bed.o cdict.o graph.o scaff_10x.o sdict.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	


