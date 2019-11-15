/*
 * =====================================================================================
 *
 *       Filename:  make_brk.c
 *
 *    Description:	generate breaks for a scaffold
 *
 *        Version:  1.0
 *        Created:  14/11/2019 09:39:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>	

#include "make_brk.h"
#include "graph.h" 
#include "sdict.h"
#include "cdict.h"
#include "bed.h"
#include "kvec.h"
#include "utls.h"

int mb_get_links_hic(char *links_fn, cdict2_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1;
		ind1 = sd_get(ctgs, r.ctgn);		
		sd_put(ctgs, r.ctgn2);		
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd2_add(&cds[ind1], r.is_l, r.ctgn2, r.is_l2, r.wt);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
}

int mb_norm_links(cdict2_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		/*char *name1 = ctgs->seq[i>>1].name;*/
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j;
		c = cds + i;
		float icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			icnt = c->cnts[j].cnt[0] + c->cnts[j].cnt[1] + c->cnts[j].cnt[2] + c->cnts[j].cnt[3]; 
            c->cnts[j].ncnt = icnt / ctgs->seq[ctg2_idx].len;
            /*c->cnts[j].ncnt = (float) icnt / (ctgs->seq[ctg2_idx].l_snp_n + ctgs->seq[ctg2_idx].r_snp_n);*/
            /*uint32_t z;*/
            /*for ( z = 0; z < 4; ++z) c->cnts[j].fcnt[z] = (float) c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[i].l_snp_n : ctgs->seq[i].r_snp_n) / ( z & 0x1 ? ctgs->seq[ctg2_idx].l_snp_n : ctgs->seq[ctg2_idx].r_snp_n);  */
			/*uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;*/
			/*fprintf(stderr, "%s\t%c\t%s\t%c\t%u\t%u\t%u\t%lf\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, snpn, snp2, 100000.0*(double)icnt/(snp2*snpn));*/
		}
	}
	return 0;
}
uint32_t *find_breaks2(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks, int lim)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        char *name = ctgs->seq[i].name;
        uint32_t scf_idx = ctgs->seq[i].le;
        uint32_t ctg_idx = i; 
		uint32_t j;
		c = cds + i;
        uint32_t fgt, flt, s_ok, p_ok;
        fgt = flt = s_ok = p_ok = 0;
        int ctgn = 0;
		for (j = 0; j < c->n_cnt; ++j) {
		   uint32_t ctg_idx2 =  sd_get(ctgs, c->cnts[j].name);	
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
			/*if (ctgn <= 50) fprintf(stderr, "%u\t%s\t%u\t%s\t%f\t%f\t%f\t%f\n", ctg_idx, name, ctg_idx2, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
				++ctgn;
				if (ctg_idx2 > ctg_idx) {
					++fgt;
					if (ctg_idx2 == ctg_idx + 1 && fgt < lim) s_ok = 1;
				} else {
					++flt;	
					if (ctg_idx2 + 1 == ctg_idx && flt < lim) p_ok = 1;
				}
				/*if (ctg_idx2 > ctg_idx + 1) */
					/*fgt = 1;*/
				/*else if (ctg_idx2 + 1 == ctg_idx) */
					/*flt = 1;	*/
			}
			if ((p_ok && s_ok) || (fgt + 1 > lim  && flt + 1 > lim)) 
				break;
        }
		/*fprintf(stderr, "%s\t%d\n", name, ctgs->seq[ctg_idx].rs >> 1 & 1);*/
		/*if (!s_ok) fprintf(stderr, "miss successor add a break\n");	*/
		/*if (!p_ok) fprintf(stderr, "miss predecessor add a break\n");*/
		//1 tail 2 head 3 head + tail 0 middle
        if (!s_ok && (ctgs->seq[ctg_idx].rs & 1) == 0)  
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (!p_ok && (ctgs->seq[ctg_idx].rs >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        /*if (ctgs->seq[ctg_idx].rs & 1) */
            /*kv_push(uint32_t, v, ctg_idx << 1 | 1);  */
	}
    *n_brks = v.n;
	return v.a;
}
uint32_t *find_breaks(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        char *name = ctgs->seq[i].name;
		/*uint32_t seq_idx = i;*/
        uint32_t scf_idx = ctgs->seq[i].le;
        uint32_t ctg_idx = i; 
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j, z;
		c = cds + i;
        uint32_t fgt, flt;
        fgt = flt = 0;
        int ctgn = 0;
        /*float density[4];*/
        uint32_t susp_hd = -1, susp_tl = -1;
		for (j = 0; j < c->n_cnt; ++j) {
            uint32_t sum = 0;
            for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z];
		   uint32_t ctg_idx2 =  sd_get(ctgs, c->cnts[j].name);	
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
                //is_suc 
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                ctgn += 1;
                /*for (z = 0; z < 4; ++z) density[z] = (float)c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[ctg_idx].r_snp_n : ctgs->seq[ctg_idx].l_snp_n) / (z&0x1 ? ctgs->seq[ctg_idx2].r_snp_n : ctgs->seq[ctg_idx2].l_snp_n);*/
                /*for (z = 0; z < 4; ++z) density[z] = (float)c->cnts[j].cnt[z]/((z >> 1 ? ctgs->seq[ctg_idx].r_snp_n : ctgs->seq[ctg_idx].l_snp_n) + (z&0x1 ? ctgs->seq[ctg_idx2].r_snp_n : ctgs->seq[ctg_idx2].l_snp_n));*/
                /*for (z = 0; z < 4; ++z) density[z] = c->cnts[j].cnt[z]/(ctgs->seq[ctg_idx].r_snp_n + ctgs->seq[ctg_idx2].r_snp_n);*/
				if (ctgn <= 50) fprintf(stderr, "%u\t%s\t%u\t%s\t%f\t%f\t%f\t%f\n", ctg_idx, name, ctg_idx2, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);

                    /*if (ctgn >= 10) break;*/
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                {
                        if (ctg_idx2 > ctg_idx) {
                            /*if (get_max(c->cnts[j].cnt, 4) > 1) {*/
                            /*if (get_fmax(density, 4) > 1) {*/
                            /*if (tl > hd || get_max(c->cnts[j].cnt, 4) > 1) {*/
                            if (tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95) { //a very loose condition for successor otherwise cause many false positives
                                if (fgt) continue;
                                ++fgt;
                                if (ctg_idx2 != ctg_idx + 1) {
                                    fprintf(stderr, "n1: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                    if (norm_cdf(tl, 0.5, tl + hd) > 0.95) kv_push(uint32_t, v, ctg_idx << 1 | 1);     
                                    else susp_tl = j;
                                } 
                            } // a successor
                            else {
                                if (flt) continue;
                                ++flt;
                                fprintf(stderr, "n0: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                if (!(ctgs->seq[ctg_idx].rs & 0x2)) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                            }  //precessor
                        } else { //potential to be a precessor 
                            /*if (get_max(c->cnts[j].cnt, 4) < 2) {*/
                            /*if (get_fmax(density, 4) < 2) {*/
                            /*if (tl < hd || get_max(c->cnts[j].cnt, 4) < 2) {*/
                            if (tl < hd || norm_cdf(tl, 0.5, tl + hd) <= 0.95) {
                                if (flt) continue;
                                ++flt;
                                if (ctg_idx2 + 1 != ctg_idx) {
                                    fprintf(stderr, "m1: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                    if (!(ctgs->seq[ctg_idx].rs & 0x2)) {
                                        if (norm_cdf(hd, 0.5, tl + hd) > .95) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                                        else susp_hd = j;
                                    }
                                } 
                            } // a predecessor
                            else {
                                if (fgt) continue;
                                ++fgt;
                                fprintf(stderr, "m0: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                            } // a successor 
                        }
                        //validate 
                }
            }
            /*if (fgt && flt) break;*/
        }
		
        if (susp_hd != 0xFFFFFFFF) {
            int insert = 1;
            for (j = susp_hd + 1; j < susp_hd + 2 && j < c->n_cnt; ++j) {
                uint32_t sum = 0;
                for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
				uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
                uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                if (scf_idx2 == scf_idx && ctg_idx2 + 1 == ctg_idx) {
                fprintf(stderr, "val:%s\t%s\t%f\t%f\t%f\t%f\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                    if (hd > tl && norm_cdf(hd, 0.5, tl + hd) > 0.95) 
                        insert = 0;
                    break;   
                } 
            }
            if (insert) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);      
        }
        if (susp_tl != 0xFFFFFFFF) {
            int insert = 1;
            for (j = susp_tl + 1; j < susp_tl + 2 && j < c->n_cnt; ++j) {
                uint32_t sum = 0;
                for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
				uint32_t seq_idx2 = sd_get(ctgs, c->cnts[j].name);
                uint32_t ctg_idx2 = ctgs->seq[seq_idx2].r_snp_n;
                uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                if (scf_idx2 == scf_idx && ctg_idx2 == ctg_idx + 1) {
                fprintf(stderr, "val:%s\t%s\t%f\t%f\t%f\t%f\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                    if (hd < tl && norm_cdf(tl, 0.5, tl + hd) >= 0.95) 
                        insert = 0;
                    break;   
                } 
            }
            if (insert) kv_push(uint32_t, v, (ctg_idx) << 1 | 1);      
        }
        if (fgt == 0 && (ctgs->seq[ctg_idx].rs & 1) == 0) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (flt == 0 && (ctgs->seq[ctg_idx].rs >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        /*if (ctgs->seq[ctg_idx].rs & 1) */
            /*kv_push(uint32_t, v, ctg_idx << 1 | 1);  */
	}
    *n_brks = v.n;
	return v.a;
}
int mb_init_scaffs(graph_t *g, sdict_t *ctgs, sdict_t *scfs)
{
	/*dump_sat(g);*/
	asm_t *as = &g->as.asms[g->as.casm];
	vertex_t *vt = g->vtx.vertices;
	path_t *pt = g->pt.paths;
	uint32_t n = as->n;
	uint32_t i;	
	/*for (i = 0; i < g->vtx.n; ++i) sd_put(ctgs, vt[i].name);	*/
	for ( i = 0; i < n; ++i) {
		uint32_t m = pt[as->pn[i] >> 1].n;
		/*fprintf(stderr, "%d\n", as->pn[i]); */
		uint32_t *p = pt[as->pn[i] >> 1].ns;
		uint32_t j, len = 0, len_ctg;
		//push scaffold name to scfs 
		int32_t scf_id = as->pn[i] >> 1;
		uint32_t sid;
		for ( j = 0; j < m; ++j ) { // pid, length,   
			len_ctg = vt[p[j]>>2].len;	
			sid = sd_put(ctgs, vt[p[j]>>2].name);
			//contig points to scaffold and set its start and direction
			/*fprintf(stderr, "sid: %u\t%s\t%s\n", sid, pt[as->pn[i]>>1].name, vt[p[j]>>2].name, p[j]);*/
			/*fprintf(stdout, "%s\t%s\t%u\n", pt[as->pn[i]>>1].name, vt[p[j]>>2].name, len);*/
			ctgs->seq[sid].len = len_ctg;
			ctgs->seq[sid].le = scf_id;
			ctgs->seq[sid].l_snp_n = p[j] & 1;
			ctgs->seq[sid].r_snp_n = j;
			ctgs->seq[sid].rs = 0;	
			if (!j) 
				ctgs->seq[sid].rs = 2;
			if (j == m - 1) 
				ctgs->seq[sid].rs |= 1, len += len_ctg;
			else 
				len += len_ctg + 200;
		}
		//reset scaffold length le rs l_snp_n, r_snp_n
		/*free(p);*/
	}
	return 0;
}
int cmp_brks(const void *a, const void *b)
{
    uint32_t p = *(uint32_t *)a;
    uint32_t q = *(uint32_t *)b;
    if (p < q) return -1;
    else if (p == q) return 0;
    else return 1;
}

int cut_paths(graph_t *g, uint32_t *brks, uint32_t n_brks, sdict_t *ctgs) 
{
    qsort(brks, n_brks, sizeof(uint32_t), cmp_brks);  
    
    uint32_t i, j;
    /*for (i = 0; i < n_brks; ++i) {*/
        /*fprintf(stdout, "%s\t%d\t%d\n", ctgs->seq[brks[i]>>1].name, brks[i]>>1, brks[i]&1);*/
    /*}*/
    /*name_t nt, nt1;*/
	char *name, *name2;
    /*for ( i = 1, j = 0; i <= n_brks; ++i) {*/
        /*if (i == n_brks || brks[i] != brks[j]) {*/
            /*uint32_t ctg_idx = brks[j];*/
			/*fprintf(stdout, "%s\t%s\n", ctgs->seq[ctg_idx >> 1].name, ctgs->seq[(ctg_idx >> 1)+1].name);*/
			/*j = i;        */
        /*} */
    /*}*/
	uint32_t scf_id = -1;	
	kvec_t(uint32_t) pos;
	kv_init(pos);
	int brkn = 0;
	int adpn = 0;
	for ( i = 1, j = 0; i <= n_brks; ++i) {
		if (i == n_brks || brks[i] != brks[j]) {
			uint32_t ctg_idx = brks[j] >> 1;
			fprintf(stderr, "[M::%s] break: %s\t%s\n", __func__, ctgs->seq[ctg_idx].name, ctgs->seq[(ctg_idx)+1].name);
			if (ctgs->seq[ctg_idx].le != scf_id) {
				if (~scf_id) brkn += pos.n, adpn += break_path(g, scf_id, pos.a, pos.n);
				kv_reset(pos);
				scf_id = ctgs->seq[ctg_idx].le;	
			}	
			kv_push(uint32_t, pos, ctgs->seq[ctg_idx].r_snp_n);
			j = i;	
		} 			
		/*uint32_t ctg_idx = brks[i];*/
			/*fprintf(stdout, "%u\t%s\t%s\n", ctgs->seq[ctg_idx>>1].le, ctgs->seq[ctg_idx >> 1].name, ctgs->seq[(ctg_idx >> 1)+1].name);*/
	}
	fprintf(stderr, "[M::%s] find %d breaks, add %d paths\n", __func__, brkn, adpn);
	kv_destroy(pos);
    return 0;
}
int mk_brks(char *sat_fn, char *links_fn, int limn, char *out_fn)
{
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] loading SAT graph\n", __func__);
#endif
	graph_t *og = load_sat(sat_fn);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] simplify SAT graph\n", __func__);
#endif
	simp_graph(og);
	sdict_t *ctgs = sd_init(); 
	sdict_t *scfs = sd_init();
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] initiate scaffolds\n", __func__);
#endif
	mb_init_scaffs(og, ctgs, scfs);
	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	uint32_t i = 0;
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict2_t* cds = calloc(n_ctg, sizeof(cdict2_t)); 
	for ( i = 0; i < n_ctg; ++i) cd2_init(cds+i); 
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	mb_get_links_hic(links_fn, cds, ctgs);

	/*for ( i = 0 ; i < n_ctg; ++i) {*/
		/*fprintf(stderr, "%s\t%u\t%u\n", ctgs->seq[i].name, ctgs->seq[i].le, ctgs->seq[i].r_snp_n);*/
	/*}*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] normalize contact matrix\n", __func__);
#endif
    mb_norm_links(cds, ctgs);
	/*return 0;*/
     //sort cd2
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] sort contact matrix\n", __func__);
#endif
	for ( i = 0; i < n_ctg; ++i) cd2_sort(cds+i); 

	uint32_t n_brks;
	uint32_t *brks = find_breaks2(cds, ctgs, &n_brks, limn);
		
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] cut paths\n", __func__);
#endif
	cut_paths(og, brks, n_brks, ctgs);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] dump sat graph\n", __func__);
#endif
	dump_sat(og, out_fn);	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] release memory resource\n", __func__);
#endif
	graph_destroy(og);
	if (brks) free(brks);
	return 0;
}

int main_brks(int argc, char *argv[])
{
	int c;
	char *sat_fn = 0, *links_fn = 0, *out_fn = 0;
	int limn = 4;
	char *program = argv[0];
	--argc, ++argv;
	while (~(c=getopt(argc, argv, "n:o:h"))) {
		switch (c) {
			case 'n':
				limn = atoi(optarg);
				break;	
			case 'o':
				out_fn = (optarg);
				break;	
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <LINKS_MATRIX> <SAT>\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -n    INT      allow successor or predecessor be top N candidates [4]\n");
				fprintf(stderr, "         -o    FILE     output file [stdout]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require a SAT and link matrix file!\n", __func__); goto help;
	}
	links_fn = argv[optind++];
	sat_fn = argv[optind];
	fprintf(stderr, "Program starts\n");	
	mk_brks(sat_fn, links_fn, limn, out_fn);
	fprintf(stderr, "Program ends\n");	
	return 0;	
}

