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
#include "get_seq.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
help:
		fprintf(stderr, "\n  pin_10x [-v] [-h] <command> [<args>]\n");
		fprintf(stderr, "  commands:\n");
		fprintf(stderr, "           link        generate link matrix for pairs of contigs\n");
		fprintf(stderr, "           build       generate a scaffolding graph with links\n");
		fprintf(stderr, "           gets        get scaffolds from a scaffolding graph\n");
		return 1;
	} else {
		if (!strcmp(argv[1], "link")) main_10x_lnks(argc , argv);
	   	else if (!strcmp(argv[1], "build")) main_bldg(argc , argv, 0);	
		else if (!strcmp(argv[1], "gets")) main_get_seq(argc, argv);
	   	else if (!strcmp(argv[1], "-h")) goto help;	
	   	else if (!strcmp(argv[1], "-v")) fprintf(stderr, "version: 0.0.0\n");	
		else {
			fprintf(stderr, "  [E::%s] unrecognized command %s\n", __func__, argv[1]);
			goto help;	
		}	
	}
	return 0;
}


