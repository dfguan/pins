/*
 * =====================================================================================
 *
 *       Filename:  ld_scaffold.c
 *
 *    Description:  scaffolding based on ld information
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

#include "htslib/synced_bcf_reader.h"
#include "sdict.h"
#include "cdict.h"
#include "bed.h"
#include "ld_io.h"
#include "graph.h"

int col_snp_distr(char *fn, uint32_t ws, sdict_t *ctgs)
{
	bcf_srs_t *sr = bcf_sr_init();
	if (!bcf_sr_add_reader(sr, fn)) {
		fprintf(stderr, "failed to open bcf file");
		return 0;
	}
	bcf_hdr_t *hr = bcf_sr_get_header(sr, 0);
	bcf_idpair_t *ctg = hr->id[BCF_DT_CTG];
	uint32_t i;
	uint32_t n_ctg = hr->n[BCF_DT_CTG];
	for ( i = 0; i < n_ctg; ++i) {
		const char *ctg_name = ctg[i].key;
		uint32_t ctg_len = ctg[i].val->info[0];
		uint32_t cur_ws = ws;
		//exclude middle part (now) or ends 
		
		if (ctg_len < (cur_ws << 1)) cur_ws = ctg_len >> 1;
		uint32_t le = ctg_len > cur_ws ? (ctg_len - cur_ws) >> 1 : 0;
		uint32_t rs = ctg_len > cur_ws ? (ctg_len + cur_ws) >> 1: ctg_len;
		sd_put3(ctgs, ctg_name, ctg_len, le, rs, 0,  0, 0);
	}
	while (bcf_sr_next_line(sr)) {
		bcf1_t *rec = bcf_sr_get_line(sr, 0);
		if (bcf_get_variant_types(rec) == VCF_SNP && rec->n_allele == 2) {
			sd_seq_t *seq = &ctgs->seq[rec->rid];
			++seq->snp_n;
			if (rec->pos < seq->le) ++seq->l_snp_n;	
			else if (rec->pos > seq->rs) ++seq->r_snp_n;
		}
	}
	bcf_sr_destroy(sr);
	return 0;
}

uint32_t check_left_half(uint32_t le, uint32_t rs, uint32_t p) // 1 for left half 0 for right half 2 for middle
{
	if (p > le && p < rs) return 2;
	else if (p <= le) return 1;
	else return 0;	
}

int add_lnks(ld_t *r, sdict_t *ctgs, sdict_t *scfs, cdict_t *cds)
{
	if (scfs->n_seq) {
		uint32_t ind1 = sd_get(ctgs, r->ctgn1);
		uint32_t ind2 = sd_get(ctgs, r->ctgn2);
		sd_seq_t *sq1 = &ctgs->seq[ind1];
		sd_seq_t *sq2 = &ctgs->seq[ind2];	
		uint32_t scf_ind1 = sq1->le;
		uint32_t scf_ind2 = sq2->le;
		if (scf_ind1 == scf_ind2) return 1;	
	
		uint32_t a0s = sq1->l_snp_n == 0 ? sq1->rs + r->pos1 : sq1->rs + sq1->len - r->pos1; 
		uint32_t a1s = sq2->l_snp_n == 0 ? sq2->rs + r->pos2: sq2->rs + sq2->len - r->pos2; 
		uint32_t is_l1 = check_left_half(scfs->seq[ind1].le, scfs->seq[ind1].rs, a0s);
		if (is_l1 > 1) return 1; //middle won't be added
		uint32_t is_l2 = check_left_half(scfs->seq[ind2].le, scfs->seq[ind2].rs, a1s);
		if (is_l2 > 1) return 1; //middle won't be added
		
		cd_add(&cds[ind1<<1|is_l1], scfs->seq[ind2].name, is_l2, is_l2?scfs->seq[ind2].l_snp_n:scfs->seq[ind2].r_snp_n);		
		cd_add(&cds[ind2<<1|is_l2], scfs->seq[ind1].name, is_l1, is_l1?scfs->seq[ind1].l_snp_n:scfs->seq[ind1].r_snp_n);		
		return 0;
	} else {
	
		uint32_t ind1 = sd_get(ctgs, r->ctgn1);
		uint32_t ind2 = sd_get(ctgs, r->ctgn2);
		/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
		uint32_t is_l1 = check_left_half(ctgs->seq[ind1].le, ctgs->seq[ind1].rs, r->pos1);
		if (is_l1 > 1) return 1; //middle won't be added
		uint32_t is_l2 = check_left_half(ctgs->seq[ind2].le, ctgs->seq[ind2].rs, r->pos2);
		if (is_l2 > 1) return 1; //middle won't be added
		
		cd_add(&cds[ind1<<1|is_l1], r->ctgn2, is_l2, is_l2?ctgs->seq[ind2].l_snp_n:ctgs->seq[ind2].r_snp_n);		
		cd_add(&cds[ind2<<1|is_l2], r->ctgn1, is_l1, is_l1?ctgs->seq[ind1].l_snp_n:ctgs->seq[ind1].r_snp_n);		
		return 0;
	}
	
}

void out_mat(cdict_t *cds, sdict_t *ctgs, uint32_t n, char *out_fn)
{
	uint32_t i;
	cdict_t *c;
	FILE *fout = out_fn ? fopen(out_fn, "w") : stdout;
	for ( i = 0; i < n; ++i) {
		c = cds + i;
		uint32_t j;
		for ( j = 0; j < c->n_cnt; ++j) {
			fprintf(fout, "%s\t%c\t%s\t%c\t%.0f\t%u\t%u\n", ctgs->seq[i>>1].name, i&1?'+':'-', c->cnts[j].name, j&1?'+':'-',c->cnts[j].cnt, c->cnts[j].snp_n, i&1?ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n);				
		}	
	}
	if (out_fn) fclose(fout);
}

int col_ld(char *ld_fn, double pv, sdict_t *ctgs, sdict_t *scfs, cdict_t *cds)
{
	ld_file_t *fp = open_ld(ld_fn);	
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open ld file: %s\n", __func__, ld_fn);
		return 1;
	}
	ld_t r;
	int ret;
	long n_records = 0, passed = 0;	
	while ((ret = read_ld1(fp, &r) >= 0)) { // we should add the snps no matter whether it  
		++n_records;
		if (r.p < pv &&  r.pos1 < r.pos2  && strcmp(r.ctgn1, r.ctgn2) < 0) 
			if (!add_lnks(&r, ctgs, scfs, cds)) ++passed;	
		if (n_records % 1000000 == 0) fprintf(stderr, "Process %lu records\n", n_records);
	}
	fprintf(stderr, "Finish processing %lu records %lu (%lf) passed\n", n_records, passed, (double)passed/n_records);
	return 0;
}
//don't merge with the one in col_hic not the same
int init_scaffs(graph_t *g, char *bcf_fn, sdict_t *ctgs, sdict_t *scfs, uint32_t ws)
{
	/*dump_sat(g);*/
	//get contigs 
	bcf_srs_t *sr = bcf_sr_init();
	if (!bcf_sr_add_reader(sr, bcf_fn)) {
		fprintf(stderr, "failed to open bcf file");
		return 0;
	}
	bcf_hdr_t *hr = bcf_sr_get_header(sr, 0);
	bcf_idpair_t *ctg = hr->id[BCF_DT_CTG];
	uint32_t i;
	uint32_t n_ctg = hr->n[BCF_DT_CTG];
	for ( i = 0; i < n_ctg; ++i) {
		const char *ctg_name = ctg[i].key;
		uint32_t ctg_len = ctg[i].val->info[0];
		uint32_t cur_ws = ws;
		//exclude middle part (now) or ends 
		
		if (ctg_len < (cur_ws << 1)) cur_ws = ctg_len >> 1;
		uint32_t le = ctg_len > cur_ws ? (ctg_len - cur_ws) >> 1 : 0;
		uint32_t rs = ctg_len > cur_ws ? (ctg_len + cur_ws) >> 1: ctg_len;
		sd_put3(ctgs, ctg_name, ctg_len, le, rs, 0,  0, 0);
	}
	
	//parse_path 
	asm_t *as = &g->as.asms[g->as.casm];
	vertex_t *vt = g->vtx.vertices;
	path_t *pt = g->pt.paths;
	uint32_t n = as->n;
	for ( i = 0; i < n; ++i) {
		uint32_t m;
		/*fprintf(stderr, "%d\n", as->pn[i]); */
		uint32_t *p = parse_path(g, as->pn[i], &m);
		uint32_t j, len = 0, len_ctg;
		//push scaffold name to scfs 
		int32_t scf_id =sd_put2(scfs, pt[as->pn[i]>>1].name, 0, 0, 0, 0, 0);
		for ( j = 0; j < m; ++j ) { // pid, length,   
			len_ctg = vt[p[j]>>2].len;	
			uint32_t sid = sd_get(ctgs, vt[p[j]>>2].name);
			//contig points to scaffold and set its start and direction
			/*fprintf(stderr, "sid: %u\t%s\t%s\n", sid, pt[as->pn[i]>>1].name, vt[p[j]>>2].name);*/
			ctgs->seq[sid].le = scf_id; //difference
			ctgs->seq[sid].rs = len;
			ctgs->seq[sid].l_snp_n = p[j] & 1;
			if (j == m - 1) 
				len += len_ctg;
			else 
				len += len_ctg + 200;
		}
		uint32_t cur_ws = ws;
		//exclude middle part (now) or ends 
		
		if (len < (cur_ws << 1)) cur_ws = len >> 1;
		uint32_t le = len > cur_ws ? (len - cur_ws) >> 1 : 0;
		uint32_t rs = len > cur_ws ? (len + cur_ws) >> 1: len;
		//reset scaffold length le rs l_snp_n, r_snp_n
		sd_put3(scfs, pt[as->pn[i]>>1].name, len, le, rs, 0, 0, 0);
		free(p);
	}
	//update snp number for both sides
	while (bcf_sr_next_line(sr)) {
		bcf1_t *rec = bcf_sr_get_line(sr, 0);
		if (bcf_get_variant_types(rec) == VCF_SNP && rec->n_allele == 2) {
			sd_seq_t *seq = &ctgs->seq[rec->rid];
			sd_seq_t *scfseq = &scfs->seq[seq->le];
			uint32_t p = seq->l_snp_n ? seq->rs + seq->len - rec->pos : seq->rs + rec->pos; 
			++seq->snp_n;
			if (p < scfseq->le) ++scfseq->l_snp_n;	
			else if (p > scfseq->rs) ++scfseq->r_snp_n;
		}
	}
	bcf_sr_destroy(sr);
	return 0;
}

int init_seqs(char *bcf_fn, char *sat_fn, uint32_t ws, sdict_t *ctgs, sdict_t *scfs)
{
	int use_sat = 0;
	if (sat_fn && *sat_fn) use_sat = 1; 
	if (!use_sat) {
		col_snp_distr(bcf_fn, ws, ctgs);	
		if (!ctgs->n_seq) return 1;
	}
	else {
		graph_t *g = load_sat(sat_fn);
		init_scaffs(g, bcf_fn, ctgs, scfs, ws);
		graph_destroy(g);
		if (!ctgs->n_seq || !scfs->n_seq) return 1;
	}	
	return 0;
}

int col_ld_lnks(char *bcf_fn, char *sat_fn, char **ld_fn, int n_ld, double pv, uint32_t ws, char *out_fn) //will be changed into two_fn
{
	
/*#ifdef VERBOSE*/
	/*fprintf(stderr, "[M::%s] collecting snps information\n", __func__);*/
/*#endif*/
	/*sdict_t *ctgs = col_snp_distr(bcf_fn, ws);	*/
	/*if (!ctgs) return 1;*/
	sdict_t *ctgs = sd_init();	
	sdict_t *scfs = sd_init();
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collect scaffolds\n", __func__);
#endif
	init_seqs(bcf_fn, sat_fn, ws, ctgs, scfs);	
	if (!scfs->n_seq && !ctgs->n_seq) {
		fprintf(stderr, "[E::%s] fail to collect sequences\n", __func__);	
		return 1;
	}
	int32_t i;
	
   	uint32_t n_cds = scfs->n_seq ? scfs->n_seq << 1 : ctgs->n_seq<<1;
	
	cdict_t* cds = calloc(n_cds, sizeof(cdict_t)); 
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing linkage disequilibrium file\n", __func__);
#endif
	
	for ( i = 0; i < n_ld; ++i) {
		if (col_ld(ld_fn[i], pv, ctgs, scfs, cds))
			return 1;
	}	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] outputing linkage information\n", __func__);
#endif
	sdict_t *_sd = scfs->n_seq ? scfs : ctgs;
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
	out_mat(cds, _sd, n_cds, out_fn);
	for (i = 0; i < n_cds; ++i)  cd_destroy(cds +i);	
	if (cds) free(cds);
	sd_destroy(ctgs);
	sd_destroy(scfs);
	return 0;
}


int main_ld_lnks(int argc, char *argv[]) 
{
	double pv = 1e-6;
	int c;
	uint32_t ws = 50000;
	
	char *program, *out_fn = 0, *sat_fn = 0;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	--argc, ++argv;	
	while (~(c=getopt(argc, argv, "p:s:w:o:h"))) {
		switch (c) {
			case 'p': 
				pv = atof(optarg);
				break;
			case 's': 
				sat_fn = optarg;
				break;
			case 'w': 
				ws = atoi(optarg);
				break;
			case 'o': 
				out_fn = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <> <LD>\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -p    FLOAT      maximum p value [1e-6]\n");
				fprintf(stderr, "         -s    FILE       sat file [NULL]\n");
				fprintf(stderr, "         -w    FLOAT      window size [50K]\n");
				fprintf(stderr, "         -o    FILE       output to a file [stdout]\n");
				fprintf(stderr, "         -h               help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] required files missing!\n", __func__); goto help;
	}
	
	char *bcf_fn = argv[optind++];
	char **ld_fn = argv + optind;
	int n_ld = argc - optind;
	//bcf file is not necessary when you can read two
	col_ld_lnks(bcf_fn, sat_fn, ld_fn, n_ld, pv, ws, out_fn); 
	
	return 0;
}
