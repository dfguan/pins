/*
 * =====================================================================================
 *
 *       Filename:  pin_ld_it.c
 *
 *    Description:  scaffolding with LD information iteratively
 *
 *        Version:  1.0
 *        Created:  13/04/2019 09:25:09
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
int main(int argc, char *argv[]) 
{
	double pv = 1e-6;
	int min_wt = 5, iter= 3;
	int c, use_sat = 0;
	uint32_t ws = 50000;
	char *program, *sat_fn = 0, *faidx_fn= 0, *seq_fn = 0;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	while (~(c=getopt(argc, argv, "p:l:i:s:x:r:w:h"))) {
		switch (c) {
			case 'p': 
				pv = atof(optarg);
				break;
			case 'w': 
				min_wt = atoi(optarg);
				break;
			case 'l': 
				ws = atoi(optarg);
				break;
			case 'i': 
				iter = atoi(optarg);
				break;
			case 's': 
				sat_fn = optarg;
				use_sat = 1;
				break;
			case 'x': 
				faidx_fn = optarg;
				break;
			case 'r':
				seq_fn = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <BCF> <LD>\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -p    FLOAT    maximum p value [1e-6]\n");
				fprintf(stderr, "         -i    INT      iteration times [3]\n");
				fprintf(stderr, "         -l    INT      windows size [50K]\n");
				fprintf(stderr, "         -w    INT      minimum linkage weight [5]\n");
				fprintf(stderr, "         -s    STR      sat file [nul]\n");
				fprintf(stderr, "         -x    STR      reference fa index file [nul]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] required files missing!\n", __func__); goto help;
	}
	
	char *bcf_fn = argv[optind++];
	char **ld_fn = argv + optind;
	int n_ld = argc - optind;
	
	int sat_fnl = sat_fn ? strlen(sat_fn) : 0;
	char mat_fn[] = "links.01.mat";
	char *sat_ofn = malloc(sizeof(char) * (max(sat_fnl, 13) + 1));
	char *sat_nfn = malloc(sizeof(char) * 14);
	if (sat_fn) 
		strcpy(sat_ofn, sat_fn);
	else
		sat_ofn[0] = 0;

	int i;
	uint32_t min_l = 0;
	for ( i = 1; i <= iter; ++i) {
		sprintf(sat_nfn, "scaffs.%02d.sat", i);
		sprintf(mat_fn, "links.%02d.mat", i);
		col_ld_lnks(bcf_fn, sat_ofn, ld_fn, n_ld, pv, ws, mat_fn);
		buildg(use_sat ? sat_ofn : faidx_fn, mat_fn, min_wt, use_sat, sat_nfn);
		if (i == iter) get_seq(sat_nfn, seq_fn, min_l, "scaffolds_final.fa");
		strcpy(sat_ofn, sat_nfn);
		use_sat = 1;
	}
	
	free(sat_ofn);
	free(sat_nfn);
	fprintf(stderr, "Program ends\n");	
	return 0;	
}
