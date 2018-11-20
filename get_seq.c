/*
 * =====================================================================================
 *
 *       Filename:  get_seq.c
 *
 *    Description:  get path sequence from gfa
 *
 *        Version:  1.0
 *        Created:  20/11/2018 13:12:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include "graph.h"


int main_get_seq(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "%s %s <GFA>", argv[0], argv[1]);
		return 1;
	}
	char *gfa_fn = argv[2];
	fprintf(stderr, "[M::%s] program starts\n", __func__);
	graph_t *g = load_gfa(gfa_fn);
	get_path(g);

	graph_destroy(g);	
	fprintf(stderr, "[M::%s] program ends\n", __func__);
	return 0;


}

