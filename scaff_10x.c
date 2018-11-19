/*
 * =====================================================================================
 *
 *       Filename:  scaff_10x.c
 *
 *    Description:  scaffolding with 10x data
 *
 *        Version:  1.0
 *        Created:  19/11/2018 19:36:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <string.h>

#include "col_10x_lnks.h"
#include "build_graph.h"


int main(int argc, char *argv[])
{
	if (argc < 2) {
help:
		fprintf(stderr, "\n  scaff_10x [-v] [-h] <command> [<args>]\n");
		fprintf(stderr, "  commands:\n");
		fprintf(stderr, "           link        generate links for contigs\n");
		fprintf(stderr, "           build       generate graph with links\n");
		return 1;
	} else {
		if (!strcmp(argv[1], "link")) main_10x_lnks(argc - 1, argv + 1);
	   	else if (!strcmp(argv[1], "build")) main_bldg(argc - 1, argv + 1);	
	   	else if (!strcmp(argv[1], "-h")) goto help;	
	   	else if (!strcmp(argv[1], "-v")) fprintf(stderr, "version: 0.0.0\n");	
		else {
			fprintf(stderr, "  [E::%s] unrecognized command %s\n", __func__, argv[1]);
			goto help;	
		}	
	}
	return 0;
}


