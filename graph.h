/*
 * =====================================================================================
 *
 *       Filename:  graph.h
 *
 *    Description:  graph	 
 *
 *        Version:  1.0
 *        Created:  21/10/2018 12:04:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _GRAPH_H
#define _GRAPH_H

#include <stdio.h>
#include <stdint.h>

typedef struct {
	char *name;
	char *seq;
	uint32_t len;
}vertex_t;

typedef struct {
	vertex_t *vertices;
	void	*h;
	uint32_t n, m;	
}vertices_t;

typedef struct {
	uint32_t v; //
	uint32_t w; //
	uint32_t wt:30, is_vis:1, is_del:1;
}edge_t;
typedef struct {
	uint32_t n, m;
	uint32_t is_srt;
	uint32_t n_del;
	edge_t *edges;
	uint64_t *edge_idx; // idx + num 
}edges_t;

typedef struct {
	uint32_t *ns;
	uint32_t n:31, is_circ:1; //be careful with these uint32_t 
	char *name;
}path_t;

typedef struct {
	path_t *paths;
	void	*h;
	uint32_t n,m; //be careful with these uint32_t
}paths_t;


typedef struct {
	vertices_t vtx;
	edges_t eg;
	paths_t pt;
} graph_t;

#define edges(g, v) ((g)->eg.edges + (g->eg.edge_idx[(v)] >> 32))
#define edge_n(g, v) ((uint32_t)g->eg.edge_idx[(v)])
#ifdef __cplusplus 
extern "C" {
#endif
	graph_t *graph_init(void);
	void graph_destroy(graph_t *g);
	int add_udedge(graph_t *g, char *sname, uint32_t sl, char *ename, uint32_t er, uint32_t wt);
	int add_dedge(graph_t *g, char *sname, uint32_t sl, char *ename, uint32_t er, uint32_t wt);
	uint32_t add_node(graph_t *g, char* name, char *seq, uint32_t len);
	int srch_path(graph_t *g);
	int out_graph(graph_t *g); // print path
	int process_graph(graph_t *g);
	uint32_t get_name2id(graph_t *g, char *nm);
	int chk_edge(uint32_t id1, uint32_t id2);
//gfa 
	graph_t  *load_gfa(char *fn);
	int get_path(graph_t *g, uint32_t minl);
	int read_seq(graph_t *g, char *fn);
#ifdef __cplusplus
}
#endif

#endif
