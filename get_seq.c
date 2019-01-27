/*
 * =====================================================================================
 *
 *       Filename:  get_seq.c
 *
 *    Description:  get path sequence from gfa
 *
 *        Version:  1.0
 *        Created:  20/11/2018 13:12:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "graph.h"


int main_get_seq(int argc, char *argv[])
{

	char *seq_fn = 0, *out_fn = 0;
	uint32_t min_l = 0;	
	
	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	
	--argc, ++argv;
	int c;
	while (~(c = getopt(argc, argv, "s:o:l:h"))) {
		switch (c) {
			case 's':
				seq_fn = optarg;
				break;
			case 'o':
				out_fn = optarg;
				break;
			case 'l':
				min_l = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
		help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <GFA> ...\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -s    STR      sequence source file\n");
				fprintf(stderr, "         -o    STR      output file [stdout]\n");
				fprintf(stderr, "         -l    STR      minimum output sequence length [stdout]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "[E::%s] Require GFA file", __func__);
		goto help;
	}
	char *gfa_fn = argv[optind];
	fprintf(stderr, "[M::%s] program starts\n", __func__);
	graph_t *g = load_gfa(gfa_fn);
	if (seq_fn) read_seq(g, seq_fn);
	get_path(g, min_l);

	graph_destroy(g);	
	fprintf(stderr, "[M::%s] program ends\n", __func__);
	return 0;


}

