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

sdict_t *col_ctgs(char *fn)
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

sdict_t *col_ctgs_from_graph(graph_t *g)
{
	sdict_t *ctgs = sd_init();
	vertex_t *vs = g->vtx.vertices;
	uint32_t n_vs = g->vtx.n;
	uint32_t i;
	for ( i = 0; i < n_vs; ++i) 
		sd_put2(ctgs, vs[i].name, vs[i].len, 0, 0, 0, 0);
	return ctgs;
}

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
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) 
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len);
	//create edges	
	for ( i = 0; i < n; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint8_t is_l = i & 1;	
		cdict_t *c = cds + i;	
		if (!c->n_cnt) continue;	
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
			ocnt = 1;
			if (hand_shaking) add_dedge(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.
		}		
	}	
	
	return g;
}



int buildg(char *fn, char *edge_fn, int min_wt, int use_sat)
{
	graph_t *og; 
	sdict_t *ctgs;
	if (use_sat) {
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contigs from sat file\n", __func__);
#endif
		og = load_sat(fn);
				
	} else {
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contigs from faidx file\n", __func__);
#endif
		ctgs = col_ctgs(fn);	
		og = graph_init();	
	}
	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict_t* cds = calloc(n_ctg<<1, sizeof(cdict_t)); 
	uint32_t n_cds = n_ctg<<1;
	uint32_t i;
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	get_links(edge_fn, cds, ctgs);
	/*anothernorm(cds, ctgs);*/
	/*return 0;*/
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i); 
	cd_set_lim(cds, n_cds, min_wt); 
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] building graph\n", __func__);
#endif
	graph_t *g = build_graph(cds, ctgs);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing graph\n", __func__);
#endif
	process_graph(g);
	

	for (i = 0; i < n_cds; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd_destroy(cds +i);	

	} 
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	graph_destroy(g);
	return 0;

}

int main_bldg(int argc, char *argv[])
{
	int c;
	uint32_t min_wt = 5; char *program = argv[0];
	char *sat_fn = 0, *ctg_fn = 0;
	int use_sat = 0;
	--argc, ++argv;
	while (~(c=getopt(argc, argv, "w:h"))) {
		switch (c) {
			case 'w': 
				min_wt = atoi(optarg);
				break;
			case 's': 
				sat_fn = optarg;
				use_sat =  1;
				break;
			case 'c': 
				ctg_fn = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <LINKS_MATRIX> \n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -w    INT      minimum weight for links [5]\n");
				fprintf(stderr, "         -c    FILE     reference index file\n");
				fprintf(stderr, "         -s    FILE     sat file\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require number of contig and links file!\n", __func__); goto help;
	}
	char *lnk_fn = argv[optind];
	fprintf(stderr, "Program starts\n");	
	if (use_sat) 
		buildg(sat_fn, lnk_fn, min_wt, 1);
	else 
		buildg(ctg_fn, lnk_fn, min_wt, 0);

	fprintf(stderr, "Program ends\n");	
	return 0;	

}

