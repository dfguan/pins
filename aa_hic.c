/*
 * =====================================================================================
 *
 *       Filename:  aa_hic.c
 *
 *    Description:  assembly assessment with hic 
 *
 *        Version:  1.0
 *        Created:  15/09/2018 15:42:59
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
#include <string.h>
#include <zlib.h>


#include "aa.h"
#include "bamlite.h"
#include "sdict.h"
#include "kvec.h"
#include "bed.h"

typedef struct {
	int mq:15, rev:1, as:16;
	uint32_t s, e, tid;
}aln_inf_t;

typedef struct {
	cord_t *ct;
	size_t n, m;	
}ns_t;


/*void col_pos(aln_inf_t  *fal, int min_as, int min_mq, uint32_t max_ins_len, ctg_pos_t *d)*/
void col_pos(aln_inf_t  *fal, int min_mq, uint32_t max_ins_len, ns_t *ns, ctg_pos_t *d)
{
	if (fal->mq > min_mq ) {
		uint32_t s = fal->s;
		uint32_t e = fal->e;
		/*fprintf(stderr, "%u\t%u\t%u\t%d\t%u\n", s, e, e - s, min_mq, max_ins_len);*/
		if (e - s < max_ins_len || span_gaps(s, e, &ns->ct[fal->tid])) {
	/*if (e-s >= max_ins_len) fprintf(stderr, "%u\t%u\t%u\n", s, e, e - s);*/
			s = s << 1;
			e = e << 1 | 1; //
			pos_push(&d->ctg_pos[fal->tid], s);
			pos_push(&d->ctg_pos[fal->tid], e); // we init 
		}
	}
}

uint32_t get_target_end(uint32_t *cigar, int n_cigar, uint32_t s)
{
	int i = 0;
	for ( i = 0; i < n_cigar; ++i) {
		uint8_t c  = bam_cigar_opchr(cigar[i]);
		if (c == 'M' || c == 'D') 
			s += cigar[i] >> BAM_CIGAR_SHIFT;
	}	
	return s;
}

void init_gaps(char *gap_fn, ns_t *ns, sdict_t *ctgs, uint32_t max_ins_len)
{
	bed_file_t* bf = bed_open(gap_fn);
	bed_rec_t r;
	ns->ct = calloc(ctgs->n_seq, sizeof(cord_t));
	while (bed_read(bf, &r) >= 0) {
		uint32_t ind = sd_put(ctgs, r.ctgn, 10);	
		if (r.e - r.s >= max_ins_len) {
			cors tmp = (cors){r.s, r.e};
			cord_push(&ns->ct[ind], &tmp);				
		}
	}
	bed_close(bf);
}




int proc_bam(char *bam_fn, int min_mq, uint32_t max_ins_len, sdict_t *ctgs, ctg_pos_t *d, char *gap_fn, ns_t *ns)
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
		sd_put(ctgs, h->target_name[i], h->target_len[i]);
		ctg_pos_push(d, i);
	}
	
	if (!ns->ct) { //not initiate yet
		init_gaps(gap_fn, ns, ctgs, max_ins_len);
	}


	char *cur_qn = NULL;
       /**cur_bc = NULL;*/
	/*int32_t cur_l = 0;*/
	/*kvec_t(aln_inf_t) fsa = {0, 0, 0}; // fisrt segment*/
	/*kvec_t(aln_inf_t) lsa = {0, 0, 0}; //last segemnt*/
	/*sdict_t* bc_n = sd_init();*/
	long bam_cnt = 0;
	int is_set = 0;
	aln_inf_t aln;
	int aln_cnt;
	uint8_t rev;
	while (1) {
		//segment were mapped 
		if (bam_read1(fp, b) >= 0 ) {
			if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {
				if (aln_cnt == 2 && (rev == 1 || rev == 2)) col_pos(&aln, min_mq, max_ins_len, ns, d);
				aln_cnt = 0;	
				rev = 0;
				is_set = 0;
				if (cur_qn) free(cur_qn); 
				cur_qn = strdup(bam1_qname(b));
			}
			if (b->core.flag & 0x4) continue; //not aligned
			++aln_cnt;
			rev = (rev << 1) | !!(b->core.flag & 0x10);
			if (b->core.isize > 0 && !is_set) {
				uint8_t *s = bam_aux_get(b, "AS");
				if (s) aln.as = *(int32_t *)(s+1); else aln.as = -1;	
				aln.s = b->core.pos + 1;
				aln.mq = b->core.qual;
				aln.tid = b->core.tid;
				/*uint32_t e = get_target_end(bam1_cigar(b), b->core.n_cigar, aln.s);			*/
				aln.e = aln.s + b->core.isize - 1;	
			} 
			if ((++bam_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld bams\n", __func__, bam_cnt); 
		} else {
			if (aln_cnt == 2 && (rev == 1 || rev == 2)) col_pos(&aln, min_mq, max_ins_len, ns, d);
			break;	
		}
	}
	fprintf(stderr, "[M::%s] finish processing %ld bams\n", __func__, bam_cnt); 
	/*while (1) {*/
		//segment were mapped 
		/*if (bam_read1(fp, b) >= 0 ) {*/
			/*fprintf(stderr, "%s\t%s\n",cur_qn, bam1_qname(b));*/
			/*if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {*/
				/*col_pos(fsa.a, fsa.n, lsa.a, lsa.n, min_as, min_mq,max_ins_len, d);*/
				/*if (cur_qn) free(cur_qn); */
				/*cur_qn = strdup(bam1_qname(b));*/
				/*cur_bc = cur_qn + b->core.l_qname - BC_LEN;*/
				/*cur_l = b->core.l_qseq;*/
				/*lsa.n = fsa.n = 0;*/
			/*}*/
			/*if (b->core.flag & 0x4) continue; //not aligned*/
			
			/*aln_inf_t tmp;*/
			/*tmp.rev = !!(b->core.flag & 0x10);*/
			/*uint8_t *s = bam_aux_get(b, "AS");*/
			/*if (s) tmp.as = *(int32_t *)(s+1); else tmp.as = -1;	*/
			/*tmp.s = b->core.pos;*/
			/*tmp.mq = *(uint16_t *)bam1_qual(b);*/
			/*tmp.tid = b->core.tid;*/
			/*uint32_t e = get_target_end(bam1_cigar(b), b->core.n_cigar, tmp.s);			*/
			/*tmp.e = e;	*/
			/*if (tmp.rev) {*/
				/*uint32_t tmp_pos = tmp.s;*/
				/*tmp.s = ctgs->seq[b->core.tid].len - tmp.e;*/
				/*tmp.e = ctgs->seq[b->core.tid].len - tmp_pos;	*/
			/*}*/
			/*if ((b->core.flag & 0x40) && !(b->core.flag & 0x80)) {*/
				/*add reverse or forward infor*/
				/*tmp.tn = (uint32_t)gfa_name2id(g, h->target_name[b->core.tid]);*/
					
				/*uint8_t *s = bam_aux_get(b, "NM");*/
				/*if (s) tmp.nm = *(int32_t *)(s+1); else tmp.nm = -1;*/
				/*kv_push(aln_inf_t,fsa,tmp);				*/
			/*} else if ((b->core.flag & 0x80) && !(b->core.flag & 0x40)) {*/
				/*aln_inf_t tmp;*/
				/*tmp.tn = (uint32_t)gfa_name2id(g, h->target_name[b->core.tid]);*/
				/*tmp.rc = !!(b->core.flag & 0x10);*/
				/*uint8_t *s = bam_aux_get(b, "NM");*/
				/*if (s) tmp.nm = *(int32_t *)(s+1); else tmp.nm = -1;*/
				/*uint8_t *s = bam_aux_get(b, "AS");*/
				/*if (s) tmp.as = *(int32_t *)(s+1); else tmp.as = -1;	*/
				/*tmp.s = b->core.pos;*/
				/*kv_push(aln_inf_t, lsa, tmp); 				*/
			/*}*/
			/*if ((++bam_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld bams\n", __func__, bam_cnt); */
		/*} else {*/
			/*col_pos(fsa.a, fsa.n, lsa.a, lsa.n, min_as, min_mq, max_ins_len, d);*/
			/*calc_alns(g, fsa.a, fsa.n, lsa.a, lsa.n, cur_l, fas, sd_put(bc_n,cur_bc, BC_LEN), seg_bc);*/
			/*break;	*/
		/*}*/
	/*}*/
	
	/*kv_destroy(fsa);*/
	/*kv_destroy(lsa);	*/
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	/*sd_destroy(bc_n);	*/

	return 0;
}

/*int aa_10x_hic(char *bam_fn, int min_as, int min_mq, int min_cov, float min_cov_rat, int max_cov, float max_cov_rat)*/
/*int aa_hic(char *bam_fn, int min_as, int min_mq, int min_cov, int max_cov, uint32_t max_ins_len)*/
int aa_hic(char **bam_fn, int n_bam, char *gap_fn, int min_mq, int min_cov, int max_cov, uint32_t max_ins_len)
{
	sdict_t *ctgs = sd_init();
		
	ctg_pos_t *d = ctg_pos_init();
	int i;	
	ns_t *ns = calloc(1, sizeof(ns_t));

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing bam file\n", __func__);
#endif

	for ( i = 0; i < n_bam; ++i) {
		if (proc_bam(bam_fn[i], min_mq, max_ins_len, ctgs, d, gap_fn, ns)) {
			return -1;	
		}	
	}	
	
	ns_destroy(ns);	

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] calculating coverage\n", __func__);
#endif
	cov_ary_t *ca = cal_cov(d, ctgs);

	ctg_pos_destroy(d);
	if (!ca) {
		fprintf(stderr, "[W::%s] low quality alignments\n", __func__);
		return 0;	
	}
		
	/*sel_sup_reg(ca, min_cov_rat, min_cov, max_cov_rat, max_cov, ctgs);*/
	char *type = "HC";
	char *desc = "hic data";
#ifdef PRINT_COVERAGE
	print_coverage(ca, ctgs, "hic");
#endif	
	sel_sup_reg(ca, min_cov, max_cov, ctgs, type, desc);
		
	cov_ary_destroy(ca, ctgs->n_seq); //a little bit messy
	sd_destroy(ctgs);

	fprintf(stderr, "Program finished successfully\n");
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
			case 'c':
				min_cov = atoi(optarg); 
				break;
			case 'C':
				max_cov = atoi(optarg); 
				break;
			case 'q':
				min_mq = atoi(optarg);
				break;
			case 'L':
				max_ins_len = strtol(optarg, &r, 10);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: aa_hic [options] <BAM_FILE>\n");
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c    INT      minimum coverage [3]\n");
				fprintf(stderr, "         -C    INT      maximum coverage [100]\n");
				fprintf(stderr, "         -q    INT      minimum alignment quality [0]\n");
				/*fprintf(stderr, "         -s    INT      minimum alignment score [0]\n");*/
				fprintf(stderr, "         -L    INT      maximum insertion length [10000]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require at least one bam file!\n", __func__); goto help;
	}
	char *gap_fn = argv[optind++];
	char **bam_fn = &argv[optind];
	int n_bam = argc - optind;
	fprintf(stderr, "Program starts\n");	
	aa_hic(bam_fn, n_bam,  gap_fn, min_mq, min_cov, max_cov, max_ins_len);
	return 0;	
}




