CC      =  gcc
CFLAGS  =  -g -Wall -D VERBOSE -D PRINT_COVERAGE #-O2  
LDFLAGS = -lz -lm

#OBJS = gfa.o opt.o paf.o sdict.o eg.o 
PROG = scaff_10x  scaff_hic scaff_gm

.SUFFIXS:.c .o

all:$(PROG)

scaff_hic: bamlite.o bed.o cdict.o graph.o scaff_hic.o sdict.o col_hic_lnks.o build_graph.o get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

scaff_10x: bamlite.o bed.o cdict.o graph.o scaff_10x.o sdict.o build_graph.o col_10x_lnks.o  get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

scaff_gm: ld_io.o bed.o cdict.o graph.o scaff_gm.o sdict.o build_graph.o col_gm_lnks.o get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	


