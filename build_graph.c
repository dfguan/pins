/*
 * =====================================================================================
 *
 *       Filename:  build_graph.c
 *
 *    Description:  build graph with links information
 *
 *        Version:  1.0
 *        Created:  19/11/2018 19:39:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "bed.h"
#include "graph.h"
#include "sdict.h"
#include "cdict.h"

int get_links(char *links_fn, cdict_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1;
		if (r.is_l) 
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, 0);		
		else
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, 0, r.llen);		
		if (r.is_l2)
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, r.rlen, 0);		
		else
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, 0, r.rlen);		
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd_add2(&cds[ind1<<1|r.is_l], r.ctgn2, r.is_l2, r.wt, r.llen);	//this has been normalized	
	} 
	/*fprintf(stderr, "%u links\n", line_n);*/
	bed_close(bf);
	return 0;	
}
int anothernorm(cdict_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq << 1;
	uint32_t i;
	cdict_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;
		uint32_t j;
		c = cds + i;
		uint32_t icnt;
		for (j = 0; j < c->n_cnt; ++j) {
			char *name2 = c->cnts[j].name; 
			icnt = c->cnts[j].cnt ; 
			uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;
			fprintf(stderr, "%s\t%c\t%s\t%c\t%u\t%u\t%u\t%lf\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, snpn, snp2, 100000.0*(double)icnt/(snp2*snpn));
		}
	}
	return 0;
}
graph_t *build_graph(cdict_t *cds, sdict_t *ctgs)
{
	graph_t *g = graph_init();
	uint32_t n = ctgs->n_seq << 1;
	uint32_t i, j;
	for ( i = 0; i < n; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint8_t is_l = i & 1;	
		cdict_t *c = cds + i;	
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2) << 1 | c->cnts[j].is_l;
			uint32_t k;
			uint8_t hand_shaking = 0;
			uint32_t ocnt = 0;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0 && cds[ind].cnts[k].is_l == is_l) {
						hand_shaking = 1;
						ocnt = cds[ind].cnts[k].cnt;
						break;
					}
			}	
			/*if (hand_shaking) fprintf(stderr, "I hand shaking\n");*/
			if (hand_shaking) add_edge2(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	
		}		
	}	
	
	return g;
}



int run(int n_ctg, char *edge_fn, int min_wt)
{
	sdict_t *ctgs = sd_init();	
	if (!ctgs) return 1;
	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict_t* cds = calloc(n_ctg<<1, sizeof(cdict_t)); 
	uint32_t n_cds = n_ctg<<1;
	
	uint32_t i;
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
	get_links(edge_fn, cds, ctgs);
	/*anothernorm(cds, ctgs);*/
	/*return 0;*/
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i); 
	cd_set_lim(cds, n_cds, min_wt); 
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
	graph_t *g = build_graph(cds, ctgs);
	process_graph(g);
	

	for (i = 0; i < n_cds; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd_destroy(cds +i);	

	} 
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
	graph_destroy(g);
	return 0;

}

int main_bldg(int argc, char *argv[])
{
	int c;
	uint32_t min_wt = 5;
	while (~(c=getopt(argc, argv, "b:B:c:C:q:S:a:L:l:h"))) {
		switch (c) {
			case 'w': 
				min_wt = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <CONTIG_NUM> <LINKS_MATRIX> \n", argv[0], argv[1]);
				fprintf(stderr, "Options:\n");
				/*fprintf(stderr, "         -L    INT      maximum insertion length [10000]\n");*/
				fprintf(stderr, "         -w    INT      minimum weight for links [5]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require number of contig and links file!\n", __func__); goto help;
	}
	int n_ctg = atoi(argv[optind++]);
	char *lnk_fn = argv[optind];
	fprintf(stderr, "Program starts\n");	
	run(n_ctg, lnk_fn, min_wt);
	return 0;	

}

