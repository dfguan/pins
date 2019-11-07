/*
 * =====================================================================================
 *
 *       Filename:  two.h
 *
 *    Description:  two header
 *
 *        Version:  1.0
 *        Created:  17/05/2018 09:10:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dg539@cam.ac.uk
 *   Organization:  Department of Genetics, Cambridge University
 *
 * =====================================================================================
 */
#ifndef _LD_IO_H
#define _LD_IO_H


#include <stdio.h>

#include "ld.h"

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif


typedef struct {
	void *fp;
	kstring_t buf;	
} ld_file_t;

#ifdef __cplusplus
extern "C" {
#endif
ld_file_t * open_ld(const char *ld_fn);
int close_ld(ld_file_t *fp);
/**
 * @func	parse_ld
 * @brief	parse ld file 
 *
 */
int read_ld1(ld_file_t *ld_fp, ld_t *t);
#ifdef __cplusplus
}
#endif
#endif




