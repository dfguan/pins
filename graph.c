/*
 * =====================================================================================
 *
 *       Filename:  graph.c
 *
 *    Description:  realization of graph functions  
 *
 *        Version:  1.0
 *        Created:  21/10/2018 12:55:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <zlib.h>

#include "graph.h"
#include "khash.h"
#include "ksort.h"
#include "kdq.h"
#include "kseq.h"
#include "kvec.h"



KDQ_INIT(uint32_t)


KSEQ_INIT(gzFile, gzread, gzseek)	
/*KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000)*/



#define edge_key(a) ((a).v)
KRADIX_SORT_INIT(edge, edge_t, edge_key, 4)

KHASH_MAP_INIT_STR(str, uint32_t)

typedef khash_t(str) shash_t;

uint8_t rc_table[128]={
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,
	4,4,4,4,65,4,4,4,4,4,4,4,4,4,4,4,4,
	84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,4,
	4,4,4,65,4,4,4,4,4,4,4,4,4,4,4
};

graph_t *graph_init(void) 
{	
	graph_t *g = (graph_t *)calloc(1, sizeof(graph_t));	
	g->vtx.h = kh_init(str);
	g->pt.h = kh_init(str);
	return g;
}

void free_node(vertex_t *n)
{
	if (!n) return;
	if (n->name) free(n->name);
	if (n->seq) free(n->seq);
}


void graph_destroy(graph_t *g)
{
	/*fprintf(stderr, "enter\n");*/
	if (g) {
		//free nodes
		uint32_t n_nd = g->vtx.n;
		uint32_t i;
		for ( i = 0 ; i < n_nd; ++i) free_node(&g->vtx.vertices[i]);
		free(g->vtx.vertices);
		if (g->vtx.h) kh_destroy(str, (shash_t *)g->vtx.h);		
		//free edges
		if (g->eg.edges) free(g->eg.edges);
		if (g->eg.edge_idx) free(g->eg.edge_idx);
		//free paths
		if (g->pt.paths) {
			for (i = 0; i < g->pt.n; ++i) free(g->pt.paths[i].ns);
			free(g->pt.paths);
		} 
		free(g);
	}

	/*fprintf(stderr, "leave\n");*/
}


int out_vetices(graph_t *g)
{

	vertex_t *vs = g->vtx.vertices;
	uint32_t n_vs = g->vtx.n;
	uint32_t i;
	for ( i = 0; i < n_vs; ++i) 
		fprintf(stdout, "S\t%s\t%s\n",  vs[i].name, vs[i].seq ? vs[i].seq : "*");
	
	return 0;
}

int out_edges(graph_t *g, int all)
{
	vertex_t *vs = g->vtx.vertices;
	edge_t *edg = g->eg.edges;
	uint32_t n_edges = g->eg.n;		
	if (all) n_edges += g->eg.n_del;	
	uint32_t i;
	for ( i = 0; i < n_edges; ++i) {
		uint32_t v = edg[i].v, w = edg[i].w;
		if (v>>1 != w >> 1) fprintf(stdout, "L\t%s\t%c\t%s\t%c\t%s\twt:%u\n", vs[v>>1].name, v&1?'+':'-', vs[w>>1].name, w&1?'+':'-', "*", edg[i].wt); // + head of sequence - tail of sequqnce	
	}
	return 0;
}

int out_paths(graph_t *g)
{
	vertex_t *vs = g->vtx.vertices;
	path_t *p = g->pt.paths;
	uint32_t n_p = g->pt.n;
	uint32_t i;	
	for ( i = 0; i < n_p; ++i) {
		uint32_t j;
		if (!p[i].name) fprintf(stdout, "P\t%c%06u\t", p[i].is_circ?'c':'u',i);
		else fprintf(stdout, "P\t%s\t", p[i].name);
		uint32_t v;	
		for ( j = 0; j + 1 < p[i].n; ++j)  {
			v = p[i].ns[j];
			fprintf(stdout, "%s%c,", vs[v>>1].name, v&1?'+':'-'); // + head of sequence - tail of sequqnce	
		}
		v = p[i].ns[j];
		fprintf(stdout, "%s%c\n", vs[v>>1].name, v&1?'+':'-');
	}
	return 0;
}

int out_graph(graph_t *g)
{
	fprintf(stdout, "H\tVN:Z:1.0\n");	
	out_vetices(g);
	out_edges(g, 0);
	out_paths(g);
	return 0;
}

uint32_t get_name2id(graph_t *g, char *nm)
{
	shash_t *h = (shash_t *)g->vtx.h;
	khint_t k = kh_get(str, h, nm);
	return k == kh_end(h) ? -1 : kh_val(h, k);
}
uint32_t add_node(graph_t *g, char* name, char *seq, uint32_t len)
{
	shash_t *h = (shash_t *)g->vtx.h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		vertices_t * ns = &g->vtx;
		if (ns->n == ns->m) {
			ns->m = ns->m ? ns->m << 1 : 16;
			ns->vertices = realloc(ns->vertices, sizeof(vertex_t) * ns->m);	
		}
		vertex_t *n = &ns->vertices[ns->n]; 
		if (seq) n->seq = strdup(seq);
		else n->seq = 0;
		n->len = len;
		kh_key(h, k) =n->name = strdup(name);	
		kh_val(h, k) = ns->n++;	
	} else {
		uint32_t ind = kh_val(h, k);
		if (seq && !g->vtx.vertices[ind].seq) g->vtx.vertices[ind].seq = strdup(seq);
		if (len && !g->vtx.vertices[ind].len) g->vtx.vertices[ind].len = len; 	
	}	
	return kh_val(h, k);
}

int idx_edge(graph_t *g)
{
	uint32_t n_vtx = g->vtx.n;
	if (!g->eg.edge_idx) g->eg.edge_idx = calloc(n_vtx<<1, sizeof(uint64_t));
	
	edges_t *es = &g->eg;
	uint32_t n_edges = es->n;
	if (!es->is_srt) {
		radix_sort_edge(es->edges, es->edges+n_edges);
		es->is_srt = 1;
	} 
	uint64_t *idx = g->eg.edge_idx;
	uint32_t i, j;
	/*fprintf(stderr, "enter %u edges\n",n_edges);*/
	if (n_edges) 	
		for ( i = 1, j = 0; i <= n_edges; ++i) {
				/*fprintf(stderr, "%u\n", i);*/
			if (i == n_edges || es->edges[i].v != es->edges[j].v) {
				/*fprintf(stderr, "%u\n", es->edges[j].v);*/
				idx[es->edges[j].v] = (uint64_t) j << 32 | (i - j);	
				j = i;	
			}
		}
	/*fprintf(stderr, "leave\n");*/
}

int add_edge1(graph_t *g, edge_t *e)
{
	edges_t *es = &g->eg;
	if (es->n == es->m) {
		es->m = es->m ? es->m << 1 : 16;
		es->edges = realloc(es->edges, sizeof(edge_t) * es->m);
	}
	es->edges[es->n++] = *e;
}

int add_udedge(graph_t *g, char *sname, uint32_t sl, char *ename, uint32_t er, uint32_t wt)
{
	uint32_t sind = add_node(g, sname, 0, 0);
	uint32_t eind = add_node(g, ename, 0, 0);
	
	edge_t e = (edge_t) {sind << 1 | sl, eind << 1 | er, wt, 0, 0};
	edge_t re = (edge_t) {eind << 1 | er, sind << 1 | sl, wt, 0, 0}; //undirect graph
	
	add_edge1(g, &e);
	add_edge1(g, &re);
	return 0;
}


int add_dedge(graph_t *g, char *sname, uint32_t sl, char *ename, uint32_t er, uint32_t wt)
{
	uint32_t sind = add_node(g, sname, 0, 0);
	uint32_t eind = add_node(g, ename, 0, 0);
	
	edge_t e = (edge_t) {sind << 1 | sl, eind << 1 | er, wt, 0, 0};
	/*edge_t re = (edge_t) {eind << 1 | er, sind << 1 | sl, wt, 0, 0}; //undirect graph*/
	
	add_edge1(g, &e);
	/*add_edge1(g, &re);*/
	return 0;
}


int add_path(graph_t *g, char *name,  uint32_t *nodes, uint32_t n) 
{
	paths_t *ps = &g->pt;
	shash_t *h = (shash_t *)ps->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		if (ps->n == ps->m) {
			ps->m = ps->m ? ps->m << 1 : 16;
			ps->paths = realloc(ps->paths, sizeof(path_t) * ps->m);
		}	
		path_t *p = &ps->paths[ps->n];
		p->ns = malloc(n*sizeof(uint32_t));
		memcpy(p->ns, nodes, sizeof(uint32_t) * n);
		p->n = n;
		kh_key(h, k) =p->name = strdup(name);	
		kh_val(h, k) = ps->n++;	
	} else 
		fprintf(stderr, "[W::%s] edge has been added\n", __func__);
	return 0;
}

int del_r(graph_t *g, edge_t *a) // remove reverse edge
{
	edge_t *e = edges(g, a->w);
	uint32_t ne = edge_n(g, a->w);
	edge_t *i;
	for ( i = e; i < e + ne; ++i) {
		if (i->w == a->v) {
			i->is_del = 1;
			break;
		}
	}
	return 0;
}

//only keeps maximum weighted edges
int clean_edges(graph_t *g)
{	
	uint32_t n_vtx = g->vtx.n; 
	uint32_t i;
	/*out_edges(g);*/
	uint32_t n_del = 0;
	for (i = 0; i < n_vtx<<1; ++i) {
		edge_t *e = edges(g, i);
		uint32_t en = edge_n(g, i);	
		edge_t *a;
		//find maximum weight
		//if there are mulitple maximum weights break them all we can't decide which way to choose
		uint32_t mwt = 0;
		uint32_t n_mwt = 0;	
		for ( a = e; a < e + en; ++a) {
			if (a->is_del) continue;
			if (a->wt > mwt) {
				mwt = a->wt;
				n_mwt = 1;	
			} else if (a->wt == mwt) {
				++n_mwt;
			} 	
		} 

		for (a = e; a < e + en; ++a) {
			if (a->is_del) continue;
			if (n_mwt > 1 || a->wt != mwt) {
				++n_del;
				a->is_del = 1;
				del_r(g, a);			
			} 		
		}
	}	
	g->eg.n_del = n_del;
	/*fprintf(stderr, "%u nodes %u del\n", n_vtx, n_del);	*/
	return 0;
}

int join_ends(graph_t *g)
{
	uint32_t n_vtx = g->vtx.n;
	
	uint32_t i;
	for ( i = 0; i < n_vtx << 1; i+=2) {
		edge_t e = (edge_t) {i, i+1, 1000, 0, 0};	
		add_edge1(g, &e);
		edge_t re = (edge_t) {i+1, i, 1000, 0, 0};
		add_edge1(g, &re);
	}
	fprintf(stderr, "[M::%s] %u edges\n", __func__, g->eg.n);
	g->eg.is_srt = 0; //set unsorted
	return 0;
}

int update_graph(graph_t *g) // update index
{
	uint32_t n_edges = g->eg.n;		
	edge_t *es = g->eg.edges;
	edge_t *end = es + n_edges; 
	uint32_t n_del = 0;
	
	for (; es < end; ++es) {
		if (es->is_del) {
			edge_t tmp = *es;
			*es = *(end - 1);
			*(end - 1) = tmp;
			--end;
			--es;	
			++n_del;
		}	
	}
	fprintf(stderr, "%u edges, %u del\n", n_edges, n_del);	
	g->eg.n -= n_del;
	idx_edge(g);
	return 0;
}


int vis_r(graph_t *g, edge_t *a)
{
	edge_t *e = edges(g, a->w);
	uint32_t ne = edge_n(g, a->w);
	edge_t *i;
	for ( i = e; i < e + ne; ++i) {
		if (i->w == a->v) {
			i->is_vis = 1;
			break;
		}
	}
	return 0;
}


int srch_path(graph_t *g)
{
	
	uint32_t n_vtx = g->vtx.n;
	
	uint8_t *mark = calloc(n_vtx << 1, sizeof(uint8_t));
	
	uint32_t i;
	
	kdq_t(uint32_t) *q;
	q = kdq_init(uint32_t);
	uint8_t is_circ;
	for ( i = 0; i < n_vtx<<1; ++i) {
		if (mark[i]) continue;
		uint32_t p, s;
		p = s = i; q->count = 0; is_circ = 0;
		while (1) {
			if (mark[p]) break;
			else {
				mark[p] = 1;
				kdq_push(uint32_t, q, p);				
			}
		 //traverse forwardly
			edge_t *es = edges(g, p);
			uint32_t es_n = edge_n(g, p);
			uint32_t j;
			/*fprintf(stderr, " enter %u\n", es_n);*/
			for (j = 0; j < es_n; ++j) if (!(es + j)->is_vis) break;
			/*fprintf(stderr, " leave \n");*/
			if (j == es_n) break;
			else {
				es = es + j;
				es->is_vis = 1;
				vis_r(g, es);
				p = es->w;
			}
		}
		/*fprintf(stderr, "2 enter \n");*/
		if (p != s || kdq_size(q) == 1) {
			p = s;
			while (1) { // traverse the other way
				edge_t *es = edges(g, p);
				uint32_t es_n = edge_n(g, p);
				uint32_t j;
				for (j = 0; j < es_n; ++j) if (!(es + j)->is_vis) break;
				/*while (j < es_n) if (!(es + j)->is_vis) break;*/
				if (j == es_n) break;
				else {
					es = es + j;
					es->is_vis = 1;
					vis_r(g, es);
					p = es->w;
				}
				if (mark[p]) break;
				else {
					mark[p] = 1;
					kdq_unshift(uint32_t, q, p);
				}	
			}	
		} else is_circ = 1;
		//add path
		if (kdq_size(q) >= 2) {
			paths_t *pths = &g->pt;
			if (pths->n == pths->m) {
				pths->m = pths->m ? pths->m << 1 : 16;
				pths->paths = realloc(pths->paths, sizeof(path_t) * pths->m);
			}
			path_t *pth = &pths->paths[pths->n++];
			pth->name = 0;
			pth->ns = malloc(sizeof(uint32_t) * (kdq_size(q) >> 1));	
			uint32_t k, m;
			for ( k = 0, m = 0; k < kdq_size(q); k += 2, ++m) pth->ns[m] = kdq_at(q, k);
			pth->n = m;
			pth->is_circ = is_circ;
		}
	}
	return 0;
}


int process_graph(graph_t *g)
{
	idx_edge(g);
	/*out_edges(g);*/
	clean_edges(g);
	join_ends(g);
	update_graph(g);
	srch_path(g);
	out_graph(g);
	return 0;
}

// GFA IO
int add_s(graph_t *g, char *s)
{	
	/*fprintf(stderr, "enters\n");*/
	char *name = 0;
	char *seq = 0;
	char *p, *q;
	int i;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) name = q;
			else if (i == 1) {
				seq = q[0] == '*' ? 0 : strdup(q);
				break;
			}	
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	
	uint32_t len = seq == 0 ? 0 : strlen(seq);
	add_node(g, name, seq, len);
	/*fprintf(stderr, "leaves\n");*/
	return 0;

}

int add_e(graph_t *g, char *s)
{
	char *n1, *n2;
	char d1, d2;
	char *p, *q;
	int i;
	uint32_t wt;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) n1 = q;
			else if (i == 1) d1 = q[0];
		   	else if (i == 2) n2 = q;
			else if (i == 3) d2 = q[0];
			else if (i == 5) wt = strtoul(q, NULL, 10);	
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	

	add_dedge(g, n1, d1 == '+', n2, d2 == '+', wt);
	return 0;
}

int add_p(graph_t *g, char *s)
{
	char *name;	
	char *p, *q;
	int i;
	kvec_t(uint32_t) ns;
	kv_init(ns);
	char *nodes_str;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) name = q;
			else if (i == 1) nodes_str = q;
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	
	for (p = q = nodes_str;; ++p) {
		int e = *p; 
		if (*p == 0 || *p == ',') {
			*p = 0;
				/*int c = q[ql - 1];*/
			int c =*(p-1);
			*(p-1) = 0;
			/*q[ql-1] = 0;*/
			uint32_t n_id = add_node(g, q, 0, 0);
			/*fprintf(stderr, "node: %s  %u\n",q, n_id);*/
			n_id = n_id << 1 | (c == '-');
			kv_push(uint32_t, ns, n_id);
			q = p + 1;
			
		}
	  	if (!e)	break;
	}
	/*fprintf(stderr, "enterp %d\n", ns.n);*/
	
	if (ns.n) add_path(g, name, ns.a, ns.n);
	kv_destroy(ns);
	/*fprintf(stderr, "leavep\n");*/
	return 0;
}


graph_t  *load_gfa(char *fn)
{
	graph_t *g = graph_init();
	
	kstream_t *ks;
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	kstring_t buf = {0, 0, 0};
	int dret;
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		if (buf.s[0] == 'S') add_s(g, buf.s);
		else if (buf.s[0] == 'L') add_e(g, buf.s);	
		else if (buf.s[0] == 'P') add_p(g, buf.s);
	}
	return g;
}

int read_seq(graph_t *g, char *seqfn)
{
	gzFile fp = seqfn && strcmp(seqfn, "-") ? gzopen(seqfn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 1;
	kseq_t *seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) 
		add_node(g, seq->name.s, seq->seq.s, seq->seq.l);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}


int cp_seq(char *s, char *t, uint32_t len, int is_rc)
{
	if (is_rc) {
		uint32_t i;
		for ( i = 0; i < len; ++i) s[i] = rc_table[t[len -i - 1]];	
	} else 
		memcpy(s, t, len * sizeof(char));
	return 0;
}

int get_path(graph_t *g, uint32_t min_l)
{	
	paths_t *ps = &g->pt;
	vertex_t *vs = g->vtx.vertices;
	uint32_t i, j;
	for ( i = 0; i < ps->n; ++i) {
		path_t *p = &ps->paths[i];
	   	char *ref_nm = p->name;
		uint32_t ref_len = 0; 	
		for ( j = 0; j < p->n; ++j) {
			uint32_t seq_len = vs[p->ns[j] >> 1].len;
			if (seq_len) {
				ref_len += seq_len;//200 'N' s	
				if (j != p->n - 1) ref_len += 200;
			} else {
				ref_len = 0;
				break;
			}
		}
		char *ref_seq = NULL;
		if (ref_len && ref_len >= min_l) {
			ref_seq = malloc(sizeof(char) * (ref_len+1));
			char *s = ref_seq;
		
			for (j = 0; j < p->n; ++j) {
				cp_seq(s, vs[p->ns[j]>> 1].seq, vs[p->ns[j]>>1].len, p->ns[j] & 1) , s += vs[p->ns[j]>>1].len;
				if (j != p->n - 1) memset(s,'N', 200), s += 200;
			}
			*s = 0;
		} 	
		fprintf(stdout, ">%s_%u\n",ref_nm, ref_len);
		if (ref_seq) {
			fprintf(stdout, "%s\n",ref_seq);
			free(ref_seq);		
		}
	}
	return 0;
}

