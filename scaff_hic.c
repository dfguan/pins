/*
 * =====================================================================================
 * *       Filename:  ld_scaffold.c
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
#include <getopt.h>
#include <math.h>


#include "sdict.h"
#include "cdict.h"
#include "bamlite.h"
#include "bed.h"
#include "kvec.h"
#include "graph.h"


typedef struct {
	/*int mq:15, rev:1, as:16;*/
	uint32_t s, ns, tid, ntid;
}aln_inf_t;

uint32_t check_left_half(uint32_t le, uint32_t rs, uint32_t p) // 1 for left half 0 for right half 2 for middle
{
	if (p > le && p < rs) return 2;
	else if (p <= le) return 1;
	else return 0;	
}

graph_t *gen_graph(cdict_t *cds, sdict_t *ctgs)
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
			if (hand_shaking) fprintf(stderr, "I hand shaking\n");
			if (hand_shaking) add_edge2(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	
		}		
	}	
	

	return g;
}


int get_edge_txt(char *links_fn, cdict_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1 = sd_get(ctgs, r.ctgn);		
		line_n += 1;
		cd_add2(&cds[ind1<<1|r.is_l], r.ctgn2, r.is_l2, r.wt);	//this has been normalized	
	} 
	/*fprintf(stderr, "%u links\n", line_n);*/
	bed_close(bf);
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
			fprintf(stdout, "%s\t%c\t%s\t%c\t%u\n", ctgs->seq[i>>1].name, i&1?'+':'-', c->cnts[j].name, j&1?'+':'-', c->cnts[j].cnt);				
		}	
	}

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

}
/*  
int core(char *snps_fn, char *edge_fn)
{
	sdict_t *ctgs = init_snps(snps_fn);	
	if (!ctgs) return 1;
	
	fprintf(stderr, "%u\n", ctgs->n_seq);
	cdict_t* cds = calloc(ctgs->n_seq<<1, sizeof(cdict_t)); 
	uint32_t n_cds = ctgs->n_seq<<1;
	
	uint32_t i;
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
	get_edge_txt(edge_fn, cds, ctgs);
	anothernorm(cds, ctgs);
	return 0;
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i); 
	cd_set_lim(cds, n_cds); 
	graph_t *g = gen_graph(cds, ctgs);
	process_graph(g);
	
	fprintf(stderr, "enter\n");

	for (i = 0; i < n_cds; ++i) {
		fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);
		cd_destroy(cds +i);	

	} 
	fprintf(stderr, "leave\n");
	if (cds) free(cds);
	graph_destroy(g);
	return 0;

}
*/
int col_joints(aln_inf_t *a, int a_cnt, aln_inf_t *f, int f_cnt, sdict_t *ctgs, cdict_t *cs)
{
	if (a_cnt == 2) {
		uint32_t ind1 = a[0].tid;
		uint32_t ind2 = a[1].tid;
		/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
		uint32_t is_l1 = check_left_half(ctgs->seq[ind1].le, ctgs->seq[ind1].rs, a[0].s);
		if (is_l1 > 1) return 1; //middle won't be added
		uint32_t is_l2 = check_left_half(ctgs->seq[ind2].le, ctgs->seq[ind2].rs, a[1].s);
		if (is_l2 > 1) return 1; //middle won't be added
		
		cd_add(&cs[ind1<<1|is_l1], ctgs->seq[ind2].name, is_l2, is_l2?ctgs->seq[ind2].l_snp_n:ctgs->seq[ind2].r_snp_n);		
		cd_add(&cs[ind2<<1|is_l2], ctgs->seq[ind1].name, is_l1, is_l1?ctgs->seq[ind1].l_snp_n:ctgs->seq[ind1].r_snp_n);		
		return 0;
	} else if (f_cnt == 2){
		uint32_t ind1 = f[0].tid;
		uint32_t ind2 = f[1].tid;
		/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
		uint32_t is_l1 = check_left_half(ctgs->seq[ind1].le, ctgs->seq[ind1].rs, f[0].s);
		if (is_l1 > 1) return 1; //middle won't be added
		uint32_t is_l2 = check_left_half(ctgs->seq[ind2].le, ctgs->seq[ind2].rs, f[1].s);
		if (is_l2 > 1) return 1; //middle won't be added
		
		cd_add(&cs[ind1<<1|is_l1], ctgs->seq[ind2].name, is_l2, is_l2?ctgs->seq[ind2].l_snp_n:ctgs->seq[ind2].r_snp_n);		
		cd_add(&cs[ind2<<1|is_l2], ctgs->seq[ind1].name, is_l1, is_l1?ctgs->seq[ind1].l_snp_n:ctgs->seq[ind1].r_snp_n);		
		
		/*cd_add(&cs[ind1<<1|is_l1], ctgs->seq[ind1].name, is_l2, is_l2?ctgs->seq[ind2].l_snp_n:ctgs->seq[ind2].r_snp_n);		*/
		/*cd_add(&cs[ind2<<1|is_l2], ctgs->seq[ind2].name, is_l1, is_l1?ctgs->seq[ind1].l_snp_n:ctgs->seq[ind1].r_snp_n);		*/
		return 0;	
	}
	return 1;
}

int proc_bam(char *bam_fn, int min_mq, uint32_t ws, sdict_t *ctgs, cdict_t **cs)
{
	
	bamFile fp;
	bam_header_t *h;
	bam1_t *b;
	fp = bam_open(bam_fn, "r"); //should check if bam is sorted
	if (fp == 0) {
		fprintf(stderr, "[E::%s] fail to open %s\n", __func__, bam_fn);
		return -1;
	}
	
	h = bam_header_read(fp);
	b = bam_init1();
	
	int i;
	for ( i = 0; i < h->n_targets; ++i) {
		char *name = h->target_name[i];
		uint32_t len = h->target_len[i];
		if (ws > len) ws = len;
		uint32_t le = (len - ws) >> 1;
		uint32_t rs = (len + ws) >> 1;
		uint32_t lenl, lenr;
		lenl = lenr = (len - ws) >> 1;
		sd_put2(ctgs, name, len, le, rs, lenl, lenr);
	}
	if (!*cs) {
		*cs = calloc(ctgs->n_seq << 1, sizeof(cdict_t));
		for ( i = 0; i < ctgs->n_seq << 1; ++i) cd_init(cs[i]);
	}		
	/*if (!ns->ct) { //not initiate yet*/
		/*init_gaps(gap_fn, ns, ctgs, max_ins_len);*/
	/*}*/
	

	char *cur_qn = NULL;
	long bam_cnt = 0;
	int is_set = 0;
	/*aln_inf_t aln[2];*/
	/*int aln_cnt;*/
	
	kvec_t(aln_inf_t) all;
	kv_init(all);
	kvec_t(aln_inf_t) five;
	kv_init(five);

	uint8_t rev;
	uint64_t rdp_counter  = 0;
	uint64_t used_rdp_counter = 0;
	while (1) {
		//segment were mapped 
		if (bam_read1(fp, b) >= 0 ) {
			if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {
				if (!col_joints(all.a, all.n, five.a, five.n, ctgs, *cs)) ++used_rdp_counter;
				/*aln_cnt = 0;	*/
				/*rev = 0;*/
				/*is_set = 0;*/
				kv_reset(all);
				kv_reset(five);
				if (cur_qn) ++rdp_counter, free(cur_qn); 
				cur_qn = strdup(bam1_qname(b));
			}
			if (b->core.flag & 0x4 || b->core.flag & 0x8 || b->core.qual < min_mq || b->core.tid == b->core.mtid) continue; //not aligned
			rev = !!(b->core.flag & 0x10);
			uint32_t *cigar = bam1_cigar(b);
			//only collects five prime
			aln_inf_t tmp;
			tmp.tid = b->core.tid;
			tmp.ntid = b->core.mtid;
			tmp.s = b->core.pos + 1;
			tmp.ns = b->core.mpos + 1;
			kv_push(aln_inf_t, all, tmp);
			
			if ((rev && bam_cigar_op(cigar[0]) == BAM_CMATCH) || (!rev && bam_cigar_op(cigar[b->core.n_cigar-1]) == BAM_CMATCH)) 
				kv_push(aln_inf_t, five, tmp);
			
			/*aln_cnt = (aln_cnt + 1 ) & 1;*/
			/*if ((++bam_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld bams\n", __func__, bam_cnt); */
		} else {
			if (!col_joints(all.a, all.n, five.a, five.n, ctgs, *cs)) ++used_rdp_counter;
			break;	
		}
	}
	fprintf(stderr, "[M::%s] finish processing %ld read pairs %ld (%.2f) passed\n", __func__, rdp_counter, used_rdp_counter, (double)used_rdp_counter/rdp_counter); 
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);

	return 0;
}

/*int aa_10x_hic(char *bam_fn, int min_as, int min_mq, int min_cov, float min_cov_rat, int max_cov, float max_cov_rat)*/
/*int aa_hic(char *bam_fn, int min_as, int min_mq, int min_cov, int max_cov, uint32_t max_ins_len)*/
int scaff_hic(char **bam_fn, int n_bam, int min_mq)
{

	/*uint32_t n_cds = ctgs->n_seq<<1;*/
	/*cdict_t* cds = calloc(ctgs->n_seq<<1, sizeof(cdict_t)); */
	sdict_t *ctgs = sd_init();	
	cdict_t *cds = 0;

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing bam file\n", __func__);
#endif
	int i;	
	uint32_t ws = 100;
	for ( i = 0; i < n_bam; ++i) {
		if (proc_bam(bam_fn[i], min_mq, ws, ctgs, &cds)) {
			return -1;	
		}	
	}	
	uint32_t n_cds = ctgs->n_seq << 1;
	for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);
	out_record(cds, ctgs, n_cds);
	for (i = 0; i < n_cds; ++i)  cd_destroy(cds +i);	
	if (cds) free(cds);
	return 0;

}

int main(int argc, char *argv[])
{
	int c;
	int max_cov = 10000000, min_cov = 7, min_mq = 30;
	int min_as = 0;
	uint32_t max_ins_len = 10000;
	/*int max_cov = 100, min_cov = 0, min_mq = 0;*/
	/*int min_as = 0;*/
	/*uint32_t max_ins_len = 10000;*/
	char *r;
	while (~(c=getopt(argc, argv, "c:C:q:s:L:h"))) {
		switch (c) {
			case 'q':
				min_mq = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: aa_hic [options] <BAM_FILE>\n");
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -q    INT      minimum alignment quality [0]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require at least one bam file!\n", __func__); goto help;
	}
	char **bam_fn = &argv[optind];
	int n_bam = argc - optind;
	fprintf(stderr, "Program starts\n");	
	scaff_hic(bam_fn, n_bam, min_mq);
	return 0;	
}

