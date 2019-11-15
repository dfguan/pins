CC      =  gcc
CFLAGS  =  -g -Wall -D VERBOSE -D PRINT_COVERAGE #-O2  
LDFLAGS = -lz -lm

#OBJS = gfa.o opt.o paf.o sdict.o eg.o 
PROG = pin_10x  pin_hic  pin_hic_it  pin_ld pin_ld_it #pin_hic2

.SUFFIXS:.c .o

all:$(PROG)
	mkdir -p ../bin && cp $(PROG) ../bin

pin_hic_it: bamlite.o bed.o cdict.o graph.o pin_hic_it.o sdict.o col_hic_lnks.o build_graph.o get_seq.o make_brk.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

pin_hic: bamlite.o bed.o cdict.o graph.o pin_hic.o sdict.o col_hic_lnks.o build_graph.o get_seq.o make_brk.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#pin_hic2: bamlite.o bed.o cdict.o graph.o pin_hic.o sdict.o col_hic_lnks.o build_graph2.o get_seq.o
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

pin_10x: bamlite.o bed.o cdict.o graph.o pin_10x.o sdict.o build_graph.o col_10x_lnks.o  get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

pin_ld: ld_io.o bed.o cdict.o graph.o pin_ld.o sdict.o build_graph.o col_ld_lnks.o get_seq.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lhts -lpthread -lz -lm -lbz2 -llzma -lcurl

pin_ld_it: ld_io.o bed.o cdict.o graph.o sdict.o build_graph.o col_ld_lnks.o get_seq.o pin_ld_it.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lhts -lpthread -lz -lm -lbz2 -llzma -lcurl

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(PROG) ../bin

	


