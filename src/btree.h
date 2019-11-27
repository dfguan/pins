/*
 * =====================================================================================
 *
 *       Filename:  btree.h
 *
 *    Description:  binary tree 
 *
 *        Version:  1.0
 *        Created:  25/11/2019 09:43:14
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

// this is not a good abstract of a btree, consider to modify it 

#ifndef _MBTREE_H
#define _MBTREE_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


typedef struct {
	uint32_t tid, nhd, ntl, nall;
	uint32_t left, right;	
}ele_t; 

typedef struct {
	uint32_t rot, n_bc, n_ele, m_ele;
	ele_t *data;			
} btree_t;

#ifdef __cplusplus
extern "C" {
#endif

uint32_t srch_bt(btree_t *bt, uint32_t tid, int *fnd);
uint32_t app_ele(btree_t *bt, uint32_t tid); //append element to data arrayhh

uint32_t insert_bt(btree_t *bt, uint32_t pidx, uint32_t tid);

int update_bt(btree_t *bt, uint32_t idx, int which);
#ifdef __cplusplus
}
#endif


#endif
