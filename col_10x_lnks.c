/*
 * =====================================================================================
 *
 *       Filename:  col_10x_lnks.c 
 *
 *    Description:  collect links from bam files 
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
#include <math.h>

#include "bamlite.h"
#include "sdict.h"
#include "kvec.h"
#include "ksort.h"
#include "cdict.h"
#include "utls.h"

#define BC_LEN 16

typedef struct {
	int mq:15, rev:1, as:16;
	uint32_t s, e, tid;
	int iden;
}aln_inf_t;



typedef struct {
	uint32_t s, e; //no more than 
	uint64_t bctn;
}bc_t;


#define bc_key(a) ((a).bctn)
KRADIX_SORT_INIT(bct, bc_t, bc_key, 8)

/*#define cord_key(a) ((a).s)*/
/*KRADIX_SORT_INIT(cord, cors, cord_key, 4);*/
	
typedef struct {
	uint32_t n, m;
	bc_t *ary;
}bc_ary_t;

void bc_ary_push(bc_ary_t *l, bc_t *z)
{
	uint32_t max = -1;
	if (l->n >= l->m) {
		if (l->m > (max >> 1)) {
			fprintf(stderr, "Too many values here\n");
			exit(1);
		} else 
			l->m = l->m ? l->m << 1 : 16;
		l->ary = realloc(l->ary, sizeof(bc_t) * l->m);//beware of overflow
	}
	l->ary[l->n++] = *z;
}

int col_bcnt(aln_inf_t  *fal, uint32_t bc, int min_mq, uint32_t max_is, bc_ary_t *l)
{
	/*if (fal->iden < 98)*/
		/*return 1;*/
	if (fal->mq > min_mq) {
		uint32_t s = fal->s;
		uint32_t e = fal->e;
		/*s = s << 1;*/
		/*e = e << 1 | 1; */
		/*if (e - s < max_ins_len) {*/
			/*if (opt) {*/
				/*uint32_t tmp = s;*/
				/*s = e;*/
				/*e = tmp;*/
			/*}*/
		if (e - s < max_is) {
			bc_t t = (bc_t) {s, e, (uint64_t)bc << 32 | fal->tid};
			
			bc_ary_push(l, &t);	
			return 0;
		}
		/*}*/
		/*pos_push(&d->ctg_pos[fal[0].tid], s);*/
		/*pos_push(&d->ctg_pos[fal[0].tid], e); // we init */
		/*ctg_pos_push(&d[fal[0].tid], s);k*/
		/*ctg_pos_push(&d[fal[0].tid], e);*/
	}
	return 1;
}

void out_matrix(cdict_t *cds, sdict_t *ctgs, uint32_t n)
{
	uint32_t i;
	cdict_t *c;
	for ( i = 0; i < n; ++i) {
		c = cds + i;
		uint32_t j;
		for ( j = 0; j < c->n_cnt; ++j) {
			if (c->cnts[j].cnt) fprintf(stdout, "%s\t%c\t%s\t%c\t%u\t%u\t%u\n", ctgs->seq[i>>1].name, i&1?'+':'-', c->cnts[j].name, j&1?'+':'-', c->cnts[j].cnt, c->cnts[j].snp_n, ctgs->seq[i>>1].l_snp_n);				
		}	
	}
}

uint32_t check_left_half(uint32_t le, uint32_t rs, uint32_t p) // 1 for left half 0 for right half 2 for middle
{
	if (p > le && p < rs) return 2;
	else if (p <= le) return 1;
	else return 0;	
}

int get_identity(uint32_t *cigar, int n_cigar, int len, int ed)
{
	int i, s = 0;
	for ( i = 0; i < n_cigar; ++i) {
		uint8_t c  = bam_cigar_opchr(cigar[i]);
		if (c == 'M' || c == '=' || c == 'X' || c == 'I') 
			s += bam_cigar_oplen(cigar[i]);
	}	
	/*fprintf(stderr, "cigar_v: %d\n", s);*/
	return (int)((float)(s - ed)/len * 100);
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


/*sdict_t *init_snps(char *fn)*/
/*{*/
	/*bed_file_t *bf = bed_open(fn);*/
	/*if (!bf) return 0;*/
	/*sdict_t *ctgs = sd_init();		*/
	/*bed_rec_t r;*/
	/*while (bed_read(bf, &r) >= 0) */
		/*sd_put2(ctgs, r.ctgn, r.len, r.le, r.rs, r.l_snp_n, r.r_snp_n);*/
	/*bed_close(bf);*/
	/*return ctgs;*/
/*}*/


int proc_bam(char *bam_fn, int min_mq, uint32_t max_is, uint32_t ws, sdict_t *ctgs, int opt, bc_ary_t *bc_l)
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
	
	/*ctg_pos_t *d = ctg_pos_init();*/
	/*for ( i = 0; i < h->n_targets; ++i) */
		/*sd_put(ctgs, h->target_name[i], h->target_len[i]);*/
		/*ctg_pos_push(d, i);*/
	int i;
	//50k
	uint32_t cur_ws;
	for ( i = 0; i < h->n_targets; ++i) {
		char *name = h->target_name[i];
		uint32_t len = h->target_len[i];
		cur_ws = ws;
		if (len < (cur_ws << 1)) cur_ws = len >> 1;
		uint32_t le = cur_ws;
		uint32_t rs = len - cur_ws + 1;
		uint32_t lenl, lenr;
		lenl = lenr = cur_ws;
		sd_put2(ctgs, name, len, le, rs, lenl, lenr);
	}
	/*for ( i = 0; i < h->n_targets; ++i) {*/
		/*char *name = h->target_name[i];*/
		/*uint32_t len = h->target_len[i];*/
		/*if (ws > len) ws = len;*/
		/*uint32_t le = (len - ws) >> 1;*/
		/*uint32_t rs = (len + ws) >> 1;*/
		/*uint32_t lenl, lenr;*/
		/*lenl = lenr = (len - ws) >> 1;*/
		/*sd_put2(ctgs, name, len, le, rs, lenl, lenr);*/
	/*}*/
	char *cur_qn = NULL, *cur_bc = NULL;
	/*int32_t cur_l = 0;*/
	sdict_t* bc_n = sd_init();
	long rdp_cnt = 0;
	long used_rdp_cnt = 0;
	int is_set = 0;
	aln_inf_t aln;
	int aln_cnt = 0;
	uint8_t rev;
	while (1) {
		//segment were mapped 
		if (bam_read1(fp, b) >= 0 ) {
			if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {
				/*fprintf(stderr, "%d\t%d\t%d\n", aln_cnt, rev, aln.mq);*/
				if (aln_cnt == 2 && (rev == 1 || rev == 2) && !col_bcnt(&aln, sd_put(bc_n, cur_bc), min_mq, max_is, bc_l)) ++used_rdp_cnt;
				aln_cnt = 0;	
				rev = 0;
				is_set = 0;
				if (cur_qn) {
					++rdp_cnt, free(cur_qn);
					if ((rdp_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld read pairs\n", __func__, rdp_cnt); 
				} 
				cur_qn = strdup(bam1_qname(b));
				cur_bc = cur_qn + b->core.l_qname - BC_LEN;
			}
			if (b->core.flag & 0x4) continue; //not aligned
			++aln_cnt;
			rev = (rev << 1) | !!(b->core.flag & 0x10);
		
			if (b->core.isize > 0 && !is_set) {
				uint8_t *s; //= bam_aux_get(b, "AS");
				/*if (s) aln.as = bam_aux2i(s); else aln.as = -1;	*/
				s = bam_aux_get(b, "NM");
				if (s) aln.iden = bam_aux2i(s); else aln.iden = 0;	
				/*fprintf(stderr, "%s %d %d\n", cur_qn,aln.as, aln.iden);*/
				/*aln.iden = get_identity(bam1_cigar(b), b->core.n_cigar, b->core.l_qseq, aln.iden);*/
				/*fprintf(stderr, "%s %d\n", cur_qn, aln.iden);*/
				aln.s = opt ? get_target_end(bam1_cigar(b), b->core.n_cigar, b->core.pos) :  b->core.pos + 1;
				aln.mq = b->core.qual;
				aln.tid = b->core.tid;
				/*uint32_t e = get_target_end(bam1_cigar(b), b->core.n_cigar, aln.s);			*/
				aln.e = aln.s + b->core.isize - 1;	
				is_set = 1;
			} 
			if (opt && b->core.isize < 0) 
				aln.e = b->core.pos + 1;
		} else {
			/*fprintf(stderr, "%d\t%d\t%d\n", aln_cnt, rev, aln.mq);*/
			if (aln_cnt == 2 && (rev == 1 || rev == 2) && !col_bcnt(&aln, sd_put(bc_n, cur_bc), min_mq, max_is, bc_l)) ++used_rdp_cnt;
			break;	
		}
	}
	fprintf(stderr, "[M::%s] finish processing %ld read pairs, %ld (%.2f) passed\n", __func__, rdp_cnt, used_rdp_cnt, (double)used_rdp_cnt/rdp_cnt); 
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	sd_destroy(bc_n);	
	return 0;
}


void insert_sort(bc_t *b, bc_t *e)
{
	bc_t *i, *j, swap_tmp;
		for (i = b + 1; i < e; ++i)										
			for (j = i; j > b && ((j-1)->s > j->s || ((j-1)->s == j->s && (j-1)->e > j->e)); --j) {			
				swap_tmp = *j; *j = *(j-1); *(j-1) = swap_tmp;			
			}															
}
void srt_by_nm_loc(bc_t *s, bc_t *e)
{
	bc_t *i = s;
	while (i < e) {
		bc_t *j;
		for (j = i + 1; j < e && j->bctn == i->bctn; ++j);
		insert_sort(i, j);
		i = j;
	}	
}

int col_joints(uint32_t ind1l, uint32_t ind2l, sdict_t *ctgs, cdict_t *cs)
{
	/*fprintf(stderr, "%u\t%u\t%u\n", ctgs->n_seq << 1, ind1l, ind2l);*/
	cd_add(&cs[ind1l], ctgs->seq[ind2l>>1].name, ind2l & 1, ind2l & 1 ? ctgs->seq[ind2l>>1].l_snp_n:ctgs->seq[ind2l>>1].r_snp_n);		
	/*fprintf(stderr, "leave2\n");*/
	cd_add(&cs[ind2l], ctgs->seq[ind1l>>1].name, ind1l & 1, ind1l & 1 ? ctgs->seq[ind1l>>1].l_snp_n:ctgs->seq[ind1l>>1].r_snp_n);		
	return 0;
}

//from  https://github.com/bcgsc/arcs
//
/*float norm_cdf(int x, float p, int n) {*/
    /*float mean = n * p;*/
    /*float sd = sqrt(n * p * (1 - p));*/
    /*return 0.5 * (1 + erf((x - mean)/(sd * sqrt(2))));*/
/*}*/


/*ctg_pos_t *col_pos(bc_ary_t *bc_l, uint32_t min_bc, uint32_t max_bc, uint32_t min_inner_bcn, uint32_t min_mol_len, int n_targets)*/
	/*cord_t *cc = col_cords(bc_l, min_bc, max_bc, min_inner_bcn, max_span, min_mol_len, ctgs->n_seq, ctgs);	*/
cdict_t *col_cds(bc_ary_t *bc_l, uint32_t min_bc, uint32_t max_bc, uint32_t min_inner_bcn, int n_targets, sdict_t *ctgs)
{
	/*cord_t *cc = calloc(n_targets, sizeof(cord_t));*/
	
	radix_sort_bct(bc_l->ary, bc_l->ary + bc_l->n);	
	uint32_t n = bc_l->n;
	bc_t *p = bc_l->ary;
	
	kvec_t(uint32_t) ctgl;	
	kv_init(ctgl);
	uint32_t i;
	cdict_t *cs = calloc(ctgs->n_seq << 1, sizeof(cdict_t));	
	for ( i = 0; i < ctgs->n_seq << 1; ++i) cd_init(&cs[i]);
	i = 0;
	while (i < n) {
		uint32_t z;
		for ( z = i; z < n && (p[z].bctn >> 32) == (p[i].bctn >> 32); ++z);
		uint32_t n_bc = z - i;
		/*fprintf(stderr, "enter a new barcode\n");*/
		if (n_bc > min_bc && n_bc < max_bc) {
			srt_by_nm_loc(p+i, p+z); //sort by ctg name and locus
			kv_reset(ctgl);
			uint32_t lt_cnt, re_cnt, mid_cnt;	
			lt_cnt = re_cnt = mid_cnt = 0;
			uint32_t is_l = check_left_half(ctgs->seq[p[i].bctn & 0xFFFF].le, ctgs->seq[p[i].bctn & 0xFFFF].rs, p[i].s); 
			if (is_l == 1) ++lt_cnt;
			else if (is_l == 0) ++re_cnt;
			else ++mid_cnt;
			uint32_t j = i + 1;
			/*fprintf(stderr, "%u\t%u\n",j,z);*/
			while (j <= z) {
				if (j == z || p[j].bctn != p[i].bctn) {
					uint32_t is_hd = lt_cnt > re_cnt ? 1 : 0;//is_head

					if (j - i > min_inner_bcn && lt_cnt != re_cnt && 1 - norm_cdf(is_hd ? lt_cnt : re_cnt, 0.5, lt_cnt + re_cnt) < 0.05)  {
					
					fprintf(stderr, "%s\t%u\t%u\t%u\t%u\n",ctgs->seq[p[i].bctn&0xFFFF].name, lt_cnt,re_cnt, mid_cnt, j - i);
						kv_push(uint32_t, ctgl, (p[i].bctn & 0xFFFF) << 1 | is_hd);
					}
					 					
					if (j == z) break; //without this will lead to infinate loop
					i = j;
					lt_cnt = re_cnt = mid_cnt = 0;
					is_l = check_left_half(ctgs->seq[p[i].bctn & 0xFFFF].le, ctgs->seq[p[i].bctn & 0xFFFF].rs, p[i].s); 
					if (is_l == 1) ++lt_cnt;
					else if (is_l == 0) ++re_cnt;
					else ++mid_cnt;
					j = i + 1;
				} else {
					is_l = check_left_half(ctgs->seq[p[j].bctn & 0xFFFF].le, ctgs->seq[p[j].bctn & 0xFFFF].rs, p[j].s);
					if (is_l == 1) ++lt_cnt;
					else if (is_l == 0) ++re_cnt;
					else ++mid_cnt;
					++j;
				}
			}
			uint32_t ctgl_s = kv_size(ctgl);
		/*fprintf(stderr, "enter\n");*/
			/*if (ctgl_s > 1) */
				/*for ( j = 0; j < ctgl_s; ++j) fprintf(stderr, "%s%c,",ctgs->seq[ctgl.a[j]>>1].name, ctgl.a[j]&1 ? '+':'-');*/
			/*fprintf(stderr, "\n");	*/
			for (j = 0; j < ctgl_s; ++j) {
					uint32_t w;
					for ( w = j + 1; w < ctgl_s; ++w) col_joints(ctgl.a[j], ctgl.a[w], ctgs, cs); 
			}
		/*fprintf(stderr, "leave\n");*/
		}	
		i = z;	
	}
	kv_destroy(ctgl);	
	return cs;
}

/*int aa_10x(char *srt_bam_fn, int min_as, int min_mq, int min_cov, float min_cov_rat, int max_cov, float max_cov_rat)*/
int col_10x_lnks(char *bam_fn[], int n_bam, int min_mq,  uint32_t win_s, uint32_t max_is, int min_bc, int max_bc, uint32_t min_inner_bcn, int opt)
{
	sdict_t *ctgs = sd_init();

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing bam file\n", __func__);
#endif
	int i;
	bc_ary_t *bc_l = calloc(1, sizeof(bc_ary_t));
	for ( i = 0; i < n_bam; ++i) {
		if (proc_bam(bam_fn[i], min_mq, max_is, win_s, ctgs, opt, bc_l)) {
			return -1;	
		}	
	}	

	if (!(bc_l&&bc_l->n)) {
		fprintf(stderr, "[W::%s] none useful information, quit\n", __func__);
		return -1;
	}
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting coordinates\n", __func__);
#endif
	cdict_t *cds = col_cds(bc_l, min_bc, max_bc, min_inner_bcn, ctgs->n_seq, ctgs);	
	free(bc_l->ary); free(bc_l);	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] normalizing joints\n", __func__);
#endif
	uint32_t n_cds = ctgs->n_seq << 1;
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
	out_matrix(cds, ctgs, n_cds);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing joints\n", __func__);
#endif
	for (i = 0; i < n_cds; ++i)  cd_destroy(cds +i);	
	if (cds) free(cds);
	sd_destroy(ctgs);
	fprintf(stderr, "Program finished successfully\n");
	return 0;

}



int main_10x_lnks(int argc, char *argv[])
{
	int c;
	int  min_mq = 20;
	uint32_t win_s = 50000, min_inner_bcn = 8; 
	uint32_t min_bc = 20, max_bc = 100000, max_is=1000;
	char *r;
	
	int option = 0; //the way to calculate molecule length //internal parameters not allowed to adjust by users
	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	
	--argc, ++argv;
	//optind points at argv[1]
	while (~(c=getopt(argc, argv, "b:B:q:L:w:a:h"))) {
		switch (c) {
			case 'b': 
				min_bc = strtol(optarg, &r, 10);
				break;
			case 'B':
				max_bc = strtol(optarg, &r, 10);
				break;
			case 'q':
				min_mq = atoi(optarg);
				break;
			case 'L':
				max_is = strtol(optarg, &r, 10);
				break;
			case 'w':
				win_s = strtol(optarg, &r, 10);
				break;
			case 'a':
				min_inner_bcn = strtol(optarg, &r, 10);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <BAM_FILE> ...\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -b    INT      minimum barcode number for each molecule [20]\n");	
				fprintf(stderr, "         -B    INT      maximum barcode number for each molecule [100000]\n");
				fprintf(stderr, "         -q    INT      minimum mapping quality [20]\n");
				fprintf(stderr, "         -w    INT      edge length [50000]\n");
				fprintf(stderr, "         -L    INT      maximum insertion length [1000]\n");
				fprintf(stderr, "         -a    INT      minimum barcode for contig [8]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require at least one bam file!\n", __func__); goto help;
	}
	char **bam_fn = argv+optind;
	int n_bam = argc - optind;
	fprintf(stderr, "Program starts\n");	
	col_10x_lnks(bam_fn, n_bam, min_mq,  win_s,  max_is, min_bc, max_bc, min_inner_bcn, option);
	fprintf(stderr, "Program ends\n");	
	return 0;	
}



