/*
 * =====================================================================================
 *
 *       Filename:  ld_scaffold.c
 *
 *    Description:  scaffolding ld information
 *
 *        Version:  1.0
 *        Created:  21/10/2018 10:04:55
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

#include "sdict.h"
#include "cdict.h"
#include "bed.h"
#include "ld_io.h"
#include "graph.h"

sdict_t *col_snp_distri(char *fn)
{
	bed_file_t *bf = bed_open(fn);
	if (!bf) return 0;
	sdict_t *ctgs = sd_init();		
	bed_rec_t r;
	while (bed_read(bf, &r) >= 0) 
		sd_put2(ctgs, r.ctgn, r.len, r.le, r.rs, r.l_snp_n, r.r_snp_n);
	bed_close(bf);
	return ctgs;
}

uint32_t check_left_half(uint32_t le, uint32_t rs, uint32_t p) // 1 for left half 0 for right half 2 for middle
{
	if (p > le && p < rs) return 2;
	else if (p <= le) return 1;
	else return 0;	
}

int add_record(ld_t *r, sdict_t *ctgs, cdict_t *cds)
{
	uint32_t ind1 = sd_get(ctgs, r->ctgn1);
	uint32_t ind2 = sd_get(ctgs, r->ctgn2);
	/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
	uint32_t is_l1 = check_left_half(ctgs->seq[ind1].le, ctgs->seq[ind1].rs, r->pos1);
	if (is_l1 > 1) return 0; //middle won't be added
	uint32_t is_l2 = check_left_half(ctgs->seq[ind2].le, ctgs->seq[ind2].rs, r->pos2);
	if (is_l2 > 1) return 0; //middle won't be added
	
	cd_add(&cds[ind1<<1|is_l1], r->ctgn2, is_l2, is_l2?ctgs->seq[ind2].l_snp_n:ctgs->seq[ind2].r_snp_n);		
	cd_add(&cds[ind2<<1|is_l2], r->ctgn1, is_l1, is_l1?ctgs->seq[ind1].l_snp_n:ctgs->seq[ind1].r_snp_n);		
	return 0;
}

void out_record(cdict_t *cds, sdict_t *ctgs, uint32_t n)
{
	uint32_t i;
	cdict_t *c;
	for ( i = 0; i < n; ++i) {
		c = cds + i;
		uint32_t j;
		for ( j = 0; j < c->n_cnt; ++j) {
			fprintf(stdout, "%s\t%c\t%s\t%c\t%u\t%u\t%u\n", ctgs->seq[i>>1].name, i&1?'+':'-', c->cnts[j].name, j&1?'+':'-',c->cnts[j].cnt, c->cnts[j].snp_n, i&1?ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n);				
		}	
	}

}

int proc_ld(char *ld_fn, double pv, sdict_t *ctgs, cdict_t *cds)
{
	ld_file_t *fp = open_ld(ld_fn);	
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open ld file: %s\n", __func__, ld_fn);
		return 1;
	}
	ld_t r;
	int ret;
	long n_records = 0;	
	while ((ret = read_ld1(fp, &r) >= 0)) { // we should add the snps no matter whether it  
		++n_records;
		if (r.p < pv &&  r.pos1 < r.pos2  && strcmp(r.ctgn1, r.ctgn2) < 0) 
			add_record(&r, ctgs, cds);	
		if (n_records % 1000000 == 0) fprintf(stderr, "Process %lu records\n", n_records);
	}
	fprintf(stderr, "Finish processing %lu records\n", n_records);
	return 0;
}


int col_gm_lnks(char *snps_fn, char **ld_fn, int n_ld, double pv) //will be changed into two_fn
{
	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting snps information\n", __func__);
#endif
	sdict_t *ctgs = col_snp_distri(snps_fn);	
	if (!ctgs) return 1;
	
	uint32_t i, n_cds = ctgs->n_seq<<1;
	cdict_t* cds = calloc(n_cds, sizeof(cdict_t)); 
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing linkage disequilibrium file\n", __func__);
#endif
	
	for ( i = 0; i < n_ld; ++i) {
		if (proc_ld(ld_fn[i], pv, ctgs, cds))
			return 1;
	}	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] outputing linkage information\n", __func__);
#endif
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
	out_record(cds, ctgs, n_cds);
	for (i = 0; i < n_cds; ++i)  cd_destroy(cds +i);	
	if (cds) free(cds);
	return 0;
}


int main_gm_lnks(int argc, char *argv[]) 
{
	double pv = 1e-6;
	int c;
	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	--argc, ++argv;	
	while (~(c=getopt(argc, argv, "p:h"))) {
		switch (c) {
			case 'p': 
				pv = atof(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <SNP_DISTRI> <LD>\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -p    FLOAT      maximum p value [1e-6]\n");
				fprintf(stderr, "         -h               help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require at least one bam file!\n", __func__); goto help;
	}
	
	char *snps_fn = argv[optind++];
	char **ld_fn = argv + optind;
	int n_ld = argc - optind;
	col_gm_lnks(snps_fn, ld_fn, n_ld, pv);
	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	return 0;
}
