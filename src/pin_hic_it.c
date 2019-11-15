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

#include "col_hic_lnks.h"
#include "build_graph.h"
#include "get_seq.h"
#include "make_brk.h"

#define max(a, b) ((a) > (b) ? (a) : (b))

int main(int argc, char *argv[])
{
	int c;
	int  min_mq = 10;
	int  iter = 3;
	int norm = 1, brk = 1, limn = 4;
	char *program;
	char *sat_fn = 0, *faidx_fn = 0, *seq_fn = 0;
	int use_sat = 0;
	char *outdir = ".";
	int use_min_dist = 1;
	int cann = 5;
	int min_wt = 100;
	uint32_t min_l  = 0;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	while (~(c=getopt(argc, argv, "O:q:nc:ds:br:a:w:i:l:x:h"))) {
		switch (c) {
			case 'q':
				min_mq = atoi(optarg);
				break;
			case 'b':
				brk = 0;
				break;
			case 'a':
				limn = atoi(optarg);
				break;
			case 'O':
				outdir = optarg;
				break;
			case 'd':
				use_min_dist = 0;
				break;
			case 'c': 
				cann = atoi(optarg);
				break;
			case 's': 
				sat_fn = optarg;
				use_sat = 1;
				break;
			case 'w': 
				min_wt = atoi(optarg);
				break;
			case 'x': 
				faidx_fn = optarg;
				break;
			case 'r':
				seq_fn = optarg;
				break;
			case 'i':
				iter = atoi(optarg);
				break;
			case 'n':
				norm = 0;
				break;
			case 'l':
				min_l = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s [options] <BAM_FILE>\n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -i    INT      iteration times\n");
				fprintf(stderr, "         -a    INT      allowed top N candidates [4]\n");
				fprintf(stderr, "         -O    STR      output directory [.]\n");
				fprintf(stderr, "         -q    INT      minimum alignment quality [10]\n");
				fprintf(stderr, "         -n    BOOL     do not use normalized weight [TRUE]\n");
				fprintf(stderr, "         -d    BOOL     do not use minimum distance to normalize weight [FALSE]\n");
				fprintf(stderr, "         -b    BOOL     do not break at the final step [FALSE]\n");
				fprintf(stderr, "         -c    INT      candidate number [5]\n");
				fprintf(stderr, "         -s    STR      sat file [nul]\n");
				fprintf(stderr, "         -x    STR      reference fa index file [nul]\n");
				fprintf(stderr, "         -r    STR      reference file [nul]\n");
				fprintf(stderr, "         -l    INT      minimum scaffold length [0]\n");
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
	if (!sat_fn && !faidx_fn) {
		fprintf(stderr,"[E::%s] require a sat file or a reference index file!\n", __func__); goto help;
	}
	if (sat_fn && !seq_fn) {
		fprintf(stderr,"[W::%s] reference file not supplied, contigs should be all in sat file!\n", __func__); goto help;
	}
	fprintf(stderr, "Program starts\n");
	//start iteration
	int loutdir = strlen(outdir);
	int sat_nfnl = loutdir + 1 + strlen("scaffs.01.sat");
	int sat_fnl = sat_fn ? max(strlen(sat_fn), sat_nfnl) : sat_nfnl;
	int scf_fnl = loutdir + 1 + strlen("scaffolds_final.fa");
	int mat_fnl = loutdir + 1 + strlen("links.01.mat");
	/*char mat_fn[] = "links.01.mat";*/
	//scaffs.01.sat 14
	char *sat_ofn = malloc(sizeof(char) * (sat_fnl + 1));
	char *sat_nfn = malloc(sizeof(char) * (sat_nfnl + 1));
	char *mat_fn = malloc(sizeof(char) * (mat_fnl + 1));
	char *scf_fn = malloc(sizeof(char) * (scf_fnl + 1));
	if (sat_fn) 
		strcpy(sat_ofn, sat_fn);
	else
		sat_ofn[0] = 0;
		
	//scaffs.01.fa
	
	int i;
	if (use_min_dist) min_wt = -1, norm = 0;
	for ( i = 1; i <= iter; ++i) {
		//input sat_fn output mat_fn
		sprintf(sat_nfn, "%s/scaffs.%02d.sat", outdir, i);
		sprintf(mat_fn, "%s/links.%02d.mat", outdir, i);
		col_hic_lnks(sat_ofn, bam_fn, n_bam, min_mq, 5000, use_min_dist, mat_fn);
		/*fprintf(stderr, "%p\n", sat_nfn);*/
		buildg_hic(use_sat ? sat_ofn : faidx_fn, mat_fn, min_wt, use_sat, norm, 0, cann, sat_nfn);
		//get seq at the final round
		
		strcpy(sat_ofn, sat_nfn);
		if (i == iter) {
			sprintf(mat_fn, "%s/links.%02d.mat", outdir, 1);
			sprintf(scf_fn, "%s/scaffolds_final.fa", outdir);
			if (brk) sprintf(sat_nfn, "%s/scaffs.bk.sat", outdir), mk_brks(sat_ofn, mat_fn, limn, sat_nfn);
			get_seq(sat_nfn, seq_fn, min_l, scf_fn);
			
		}
		use_sat = 1;
		/*use_min_dist = 0;*/
		/*min_wt = 100;*/
		/*norm = 1;*/
	}
	free(sat_ofn);
	free(sat_nfn);
	free(scf_fn);
	free(mat_fn);
	fprintf(stderr, "Program ends\n");	
	return 0;	
}



