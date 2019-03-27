/*
 * =====================================================================================
 *
 *       Filename:  pin_hic_it.c
 *
 *    Description:  pin hic iteratively
 *
 *        Version:  1.0
 *        Created:  26/03/2019 11:12:31
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
#include <stdint.h>
#include <string.h>
#include <getopt.h>

#define max(a, b) ((a) > (b) ? (a) : (b))

int main(int argc, char *argv[])
{
	int c;
	int  min_mq = 10;
	int  iter = 3;
	char *program;
	char *sat_fn = 0, *ctg_fn = 0, *seq_fn = 0;
	int use_sat = 0;
	int min_wt = 5;
	uint32_t min_l  = 0;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	--argc, ++argv;
	while (~(c=getopt(argc, argv, "q:i:o:s:h"))) {
		switch (c) {
			case 'q':
				min_mq = atoi(optarg);
				break;
			case 'w': 
				min_wt = atoi(optarg);
				break;
			case 's': 
				sat_fn = optarg;
				use_sat = 1;
				break;
			case 'c': 
				ctg_fn = optarg;
				break;
			case 'r':
				seq_fn = optarg;
				break;
			case 'l':
				min_l = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [options] <BAM_FILE>\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -q    INT      minimum alignment quality [10]\n");
				fprintf(stderr, "         -w    INT      minimum linkage weight [5]\n");
				fprintf(stderr, "         -s    STR      sat file\n");
				fprintf(stderr, "         -c    STR      reference fa index file \n");
				fprintf(stderr, "         -r    STR      reference file \n");
				fprintf(stderr, "         -l    STR      minimum scaffold length [0]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require at least one bam file!\n", __func__); goto help;
	}
	char **bam_fn = &argv[optind];
	int n_bam = argc - optind;
	//check parameters	
	if (!sat_fn && !ctg_fn) {
		fprintf(stderr,"[E::%s] require a sat file or a reference index file!\n", __func__); goto help;
		
	}
	if (sat_fn && !seq_fn) {
		fprintf(stderr,"[W::%s] reference file not supplied, contigs should be all in sat file!\n", __func__); goto help;
	}
	fprintf(stderr, "Program starts\n");
	//start iteration
	int sat_fnl = sat_fn ? strlen(sat_fn) : 0;
	char mat_fn[] = "links.01.mat";
	/*char mat_fn[] = "links.01.mat";*/
	//scaffs.01.sat 14
	char *sat_ofn = malloc(sizeof(char) * (max(sat_fnl, 13) + 1));
	char *sat_nfn = malloc(sizeof(char) * 14);

	if (sat_fn) 
		strcpy(sat_ofn, sat_fn);
	else
		sat_ofn[0] = 0;

	//scaffs.01.fa
	
	int i;
	for ( i = 1; i <= iter; ++i) {
		//input sat_fn output mat_fn
		sprintf(sat_nfn, "scaffs.%02d.sat", i);
		sprintf(mat_fn, "links.%02d.mat", i);
		col_hic_lnks(sat_ofn, bam_fn, n_bam, min_mq, 5000, mat_fn);
		buildg(sat_ofn, mat_fn, min_wt, use_sat, sat_nfn);
		//get seq at the final round
		if (i == iter) get_seq(sat_nfn, seq_fn, min_l, "scaffolds_final.fa");
		strcpy(sat_ofn, sat_nfn);
	}
	free(sat_ofn);
	free(sat_nfn);
	fprintf(stderr, "Program ends\n");	
	return 0;	
}



