CC      =  gcc
CFLAGS  =  -g -Wall -D VERBOSE -D PRINT_COVERAGE #-O2  
LDFLAGS = -lz -lm

#OBJS = gfa.o opt.o paf.o sdict.o eg.o 
PROG = pin_10x  pin_hic pin_ld

.SUFFIXS:.c .o

all:$(PROG)

pin_hic: bamlite.o bed.o cdict.o graph.o pin_hic.o sdict.o col_hic_lnks.o build_graph.o get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

pin_10x: bamlite.o bed.o cdict.o graph.o pin_10x.o sdict.o build_graph.o col_10x_lnks.o  get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

pin_ld: ld_io.o bed.o cdict.o graph.o pin_ld.o sdict.o build_graph.o col_ld_lnks.o get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	


