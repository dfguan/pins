/*
 * =====================================================================================
 *
 *       Filename:  btree.c
 *
 *    Description:  btree implementation
 *
 *        Version:  1.0
 *        Created:  25/11/2019 10:57:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include "btree.h"

uint32_t srch_bt(btree_t *bt, uint32_t tid, int *fnd)
{	
	/*fprintf(stderr, "ele: %lu\t%p\n",bt->n_ele, &bt->data[bt->rot]);*/
	*fnd = 0;
	if (!bt->n_ele) return -1;
	uint32_t i;
	ele_t *z;
	/*ele_t * z = &bt->data[bt->rot];*/
	for (i = bt->rot;;) {
		/*fprintf(stderr, "z: %p %u\t%u\t%u\n",z, z->tid, z->left, z->right);	*/
		z = &bt->data[i];
		if (z->tid  == tid) {
			*fnd = 1;
			return i;
		} else if (z->tid > tid) {
			if (~z->left) 
				i = z->left;
			else 
				return i;
		} else {
			if (~z->right) 
				i = z->right;
			else
				return i;
		}
	} 
}

uint32_t app_ele(btree_t *bt, uint32_t tid) //append element to data array
{
	if (bt->n_ele >= bt->m_ele) {
		bt->m_ele = bt->m_ele ? bt->m_ele << 1 : 16;
		bt->data =(ele_t *)realloc(bt->data, bt->m_ele * sizeof(ele_t));
		/*fprintf(stderr, "zco\n");*/
	}
	bt->data[bt->n_ele]= (ele_t) {tid, 0, 0, 0, 0xFFFFFFFF, 0xFFFFFFFF};	
	uint32_t rtnv = bt->n_ele ++;
	return rtnv;
}

uint32_t insert_bt(btree_t *bt, uint32_t pidx, uint32_t tid)
{
	//either insert or update a vertex
	uint32_t idx = app_ele(bt, tid);
	/*fprintf(stderr, "to:%p %d %lu\n", pe, !!pe, pe-bt->data);*/
	if (~pidx) {//need to think about how to add root
		ele_t *pe = &bt->data[pidx];
		if (pe->tid > tid) pe->left = idx;
		else pe->right = idx;
	}	
	/*pe->tid > tid ? pe>left = idx : pe->right = idx;*/
	return idx;	
}

int update_bt(btree_t *bt, uint32_t idx, int which) 
{
	ele_t *d = &bt->data[idx];	
	++bt->n_bc;
	if (which == 1) ++d->nhd;
	else if (which == 0) ++d->ntl;
	++d->nall;	
	return 0;
}
