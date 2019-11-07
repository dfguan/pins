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
#include "utls.h"

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
	asm_t *ca = &g->as.asms[g->as.casm];
	path_t *pt = g->pt.paths;
	uint32_t n = ca->n;
	uint32_t i, j;
	for ( i = 0; i < n; ++i) {
		//get path length
		/*uint32_t path_len = 0;*/
		/*for (j = 0; j < pt[ca->pn[i]>>1].n; ++j) */
			/*path_len += ctgs->seq[sd_get(ctgs, pt[ca->pn[i]>>1].name)].len;	*/
		/*path_len += (pt[ca->pn[i]>>1].n - 1) * 200;	 //gap size = 200	*/
		sd_put2(ctgs, pt[ca->pn[i]>>1].name, pt[ca->pn[i]>>1].len, 0, 0, 0, 0);
	} 
	return ctgs;
}

int get_links_hic(char *links_fn, cdict2_t *cds, sdict_t *ctgs)
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
		cd2_add(&cds[ind1], r.is_l, r.ctgn2, r.is_l2, r.fwt);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
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
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.rlen, 0); //wrong length add to contig
		else
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, 0, r.rlen);		
		if (r.is_l2)
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, r.llen, 0);		
		else
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, 0, r.llen);		
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd_add2(&cds[ind1<<1|r.is_l], r.ctgn2, r.is_l2, r.fwt, r.llen);	//this has been normalized	
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

int print_cdict2(cdict2_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		char *name1 = ctgs->seq[i].name;
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j;
		c = cds + i;
		uint32_t icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			/*fprintf(stderr, "ctg_idx: %u\n", ctg2_idx);*/
			icnt = c->cnts[j].cnt[0] + c->cnts[j].cnt[1] + c->cnts[j].cnt[2] + c->cnts[j].cnt[3]; 
            c->cnts[j].ncnt = (float) icnt / ctgs->seq[ctg2_idx].len;
            /*c->cnts[j].ncnt = (float) icnt / (ctgs->seq[ctg2_idx].l_snp_n + ctgs->seq[ctg2_idx].r_snp_n);*/
            /*uint32_t z;*/
            /*for ( z = 0; z < 4; ++z) c->cnts[j].fcnt[z] = (float) c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[i].l_snp_n : ctgs->seq[i].r_snp_n) / ( z & 0x1 ? ctgs->seq[ctg2_idx].l_snp_n : ctgs->seq[ctg2_idx].r_snp_n);  */
			/*uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;*/
			fprintf(stderr, "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%u\n", name1, name2, c->cnts[j].ncnt, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3], ctgs->seq[ctg2_idx].len);
		}
	}
	return 0;
}

int print_cdict(cdict_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq << 1;
	uint32_t i;
	cdict_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;
		uint32_t j;
		c = cds + i;
		float icnt;
		uint32_t intcnt;
		for (j = 0; j < c->n_cnt; ++j) {
			char *name2 = c->cnts[j].name; 
			icnt = c->cnts[j].cnt ; 
			intcnt = 1; //c->cnts[j].intcnt;
			uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;
			fprintf(stderr, "%s\t%c\t%s\t%c\t%f\t%u\t%u\t%u\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, intcnt, snpn, snp2);
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
			float ocnt = 0;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0 && cds[ind].cnts[k].is_l == is_l) {
						hand_shaking = 1;
						ocnt = cds[ind].cnts[k].cnt;
						break;
					}
			}	
			/*if (hand_shaking) fprintf(stderr, "I hand shaking\n");*/
			/*ocnt = 1;*/
			if (hand_shaking) add_dedge(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.
		}		
	}	
	return g;
}

int norm_links(cdict2_t *cds, sdict_t *ctgs, int norm)
{
	uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		/*char *name1 = ctgs->seq[i>>1].name;*/
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j;
		c = cds + i;
		uint32_t icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			icnt = c->cnts[j].cnt[0] + c->cnts[j].cnt[1] + c->cnts[j].cnt[2] + c->cnts[j].cnt[3]; 
            c->cnts[j].ncnt = norm ? (float) icnt / ctgs->seq[ctg2_idx].len : (float) icnt;
            /*c->cnts[j].ncnt = (float) icnt / (ctgs->seq[ctg2_idx].l_snp_n + ctgs->seq[ctg2_idx].r_snp_n);*/
            /*uint32_t z;*/
            /*for ( z = 0; z < 4; ++z) c->cnts[j].fcnt[z] = (float) c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[i].l_snp_n : ctgs->seq[i].r_snp_n) / ( z & 0x1 ? ctgs->seq[ctg2_idx].l_snp_n : ctgs->seq[ctg2_idx].r_snp_n);  */
			/*uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;*/
			/*fprintf(stderr, "%s\t%c\t%s\t%c\t%u\t%u\t%u\t%lf\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, snpn, snp2, 100000.0*(double)icnt/(snp2*snpn));*/
		}
	}
	return 0;
}
int det_ori(float *ws)
{
	float l, sl;
	l = ws[0];
	sl = 0;
	int i, maxi = 0;
	for (i = 1; i < 4; ++i ) 
		if (ws[i] >= l) 
				maxi = i, sl = l, l = ws[i];  
		else if (ws[i] > sl) sl = ws[i];
	if (sl != l)
		return maxi;
	else
		return -1;
	/*if (norm_cdf(l, 0.5, sl + l) <= 0.95) */
			/*return -1; */
	/*else*/
		   /*return maxi;	*/
}
graph_t *build_graph_hic(cdict2_t *cds, sdict_t *ctgs)
{
	graph_t *g = graph_init();
	
	uint32_t i, j;
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) 
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len);
	//create edges
		
	for ( i = 0; i < ctgs->n_seq; ++i) {
		char *name1 = ctgs->seq[i].name;
		cdict2_t *c = cds + i;	

		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2);
			uint32_t k;
			uint8_t hand_shaking = 0;
			uint32_t is_l, is_l2;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0) {
						hand_shaking = 1;
						break;
					}
			}	
			/*if (hand_shaking) fprintf(stderr, "I hand shaking\n");*/
			//determine joint direction, here we use a very easy model
			if (!hand_shaking) continue;
			int idx;
			if (~(idx = det_ori(c->cnts[j].cnt)))
			/*uint32_t hh = c->cnts[j].cnt[0];*/
			/*uint32_t ht = c->cnts[j].cnt[1];*/
			/*uint32_t th = c->cnts[j].cnt[2];*/
			/*uint32_t tt = c->cnts[j].cnt[3];*/
			/*uint32_t tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];*/
			/*uint32_t hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];	*/
			/*uint32_t hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];*/
			/*uint32_t tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];*/
			/*if (tl != hd) */
					/*is_l = tl > hd ? 0 : 1; */
			/*else */
					/*continue;*/
			/*if (hd2 != tl2) */
					/*is_l2 = tl2 > hd ? 0 : 1; 	*/
			/*else */
					/*continue;*/
			is_l = idx >> 1, is_l2 = idx & 1, add_dedge(g, name1, is_l, name2, is_l2, c->cnts[j].cnt[is_l << 1 | is_l2]);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.
		}		
	}	
	return g;
}
int buildg_hic(char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn)
{
	/*fprintf(stderr, "%s %s\n", fn, edge_fn);*/
	/*fprintf(stderr, "%d %d\n", min_wt, use_sat);*/
	/*fprintf(stderr, "%d %f\n", norm, min_mdw);*/
	/*fprintf(stderr, "%d\n", mlc);*/
	/*fprintf(stderr, "%p %s\n", out_fn, out_fn);*/
	graph_t *og; 
	sdict_t *ctgs = 0;
	if (use_sat) {
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] collecting contigs from sat file\n", __func__);
#endif
		og = load_sat(fn);
		ctgs = col_ctgs_from_graph(og);
	} else {
		if (!fn) {
			fprintf(stderr, "[E::%s] please set reference index file with -c\n", __func__);
			return 1;
		}
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contigs from faidx file\n", __func__);
#endif
		ctgs = col_ctgs(fn);	
		og = graph_init();	
	}
	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict2_t* cds = calloc(n_ctg, sizeof(cdict2_t)); 
	uint32_t i;
	for ( i = 0; i < n_ctg; ++i) cd2_init(cds+i); 
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	get_links_hic(edge_fn, cds, ctgs);
	/*anothernorm(cds, ctgs);*/
	/*return 0;*/
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*print_cdict2(cds, ctgs);	*/
	norm_links(cds, ctgs, norm);
	for ( i = 0; i < n_ctg; ++i) cd2_sort(cds+i); 
	print_cdict2(cds, ctgs);	
	cd2_set_lim(cds, n_ctg, mlc); 
	/*if (norm) */
	/*if (norm) cd_filt(cds, n_cds, min_rat); */
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] building graph\n", __func__);
#endif
	graph_t *g = build_graph_hic(cds, ctgs);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing graph\n", __func__);
#endif
	process_graph(g);
			
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] merging graph\n", __func__);
#endif
	merge_graph(og, g, 1);

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] output graph\n", __func__);
#endif

	dump_sat(og, out_fn);
	for (i = 0; i < n_ctg; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd2_destroy(cds +i);	

	} 
	/*fprintf(stderr, "CheckPoint2\n");*/
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	graph_destroy(og);
	return 0;

}

int buildg(char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn)
{
	fprintf(stderr, "%s %s\n", fn, edge_fn);
	fprintf(stderr, "%d %d\n", min_wt, use_sat);
	fprintf(stderr, "%d %f\n", norm, min_mdw);
	fprintf(stderr, "%d\n", mlc);
	fprintf(stderr, "%p\n", out_fn);
	graph_t *og; 
	sdict_t *ctgs = 0;
	if (use_sat) {
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] collecting contigs from sat file\n", __func__);
#endif
		og = load_sat(fn);
		ctgs = col_ctgs_from_graph(og);
	} else {
		if (!fn) {
			fprintf(stderr, "[E::%s] please set reference index file with -c\n", __func__);
			return 1;
		}
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
	if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i); 
	print_cdict(cds, ctgs);	
	cd_set_lim(cds, n_cds, min_wt, min_mdw, mlc, norm); 
	/*if (norm) */
	/*if (norm) cd_filt(cds, n_cds, min_rat); */
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] building graph\n", __func__);
#endif
	graph_t *g = build_graph(cds, ctgs);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing graph\n", __func__);
#endif
	process_graph(g);
			
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] merging graph\n", __func__);
#endif
	merge_graph(og, g, 1);

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] output graph\n", __func__);
#endif

	dump_sat(og, out_fn);

	for (i = 0; i < n_cds; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd_destroy(cds +i);	

	} 
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	graph_destroy(og);
	return 0;

}

int main_bldg(int argc, char *argv[], int ishic)
{
	int c;
	uint32_t min_wt = 5; char *program = argv[0];
	char *sat_fn = 0, *ctg_fn = 0, *out_fn = 0;
	int use_sat = 0, mlc = 1;
	int norm = 0;
	float msn = .7, mdw = 0.95;
	--argc, ++argv;
	while (~(c=getopt(argc, argv, "w:o:s:c:nm:f:k:h"))) {
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
			case 'f': 
				mdw = atof(optarg);
				break;
			case 'k': 
				mlc = atoi(optarg);
				break;
			case 'n': 
				norm = 1;
				break;
			/*case 'm': */
				/*msn = atof(optarg);*/
				/*break;*/
			case 'o': 
				out_fn = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <LINKS_MATRIX> \n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -w    INT      minimum weight for links [10]\n");
				fprintf(stderr, "         -k    INT      maximum linking candiates [1]\n");
				fprintf(stderr, "         -c    FILE     reference index file [nul]\n");
				fprintf(stderr, "         -n    BOOL     normalize weight [false]\n");
				/*fprintf(stderr, "         -m    FLOAT    minimum barcode ratio [.7]\n");*/
				fprintf(stderr, "         -f    FLOAT    minimum weight difference [0.95]\n");
				fprintf(stderr, "         -s    FILE     sat file [nul]\n");
				fprintf(stderr, "         -o    FILE     output file [stdout]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require number of contig and links file!\n", __func__); goto help;
	}
	char *lnk_fn = argv[optind];
	fprintf(stderr, "Program starts\n");	
	if (!sat_fn) sat_fn = ctg_fn;
	int ret;
	if (!ishic)
		ret = buildg(sat_fn, lnk_fn, min_wt, use_sat, norm, mdw, mlc, out_fn);
	else 
		ret = buildg_hic(sat_fn, lnk_fn, min_wt, use_sat, norm, mdw, mlc, out_fn);

	fprintf(stderr, "Program ends\n");	
	return ret;	

}
