/*
 * =====================================================================================
 *
 *       Filename:  ld_io.c
 *
 *    Description:	parse two text file 
 *
 *        Version:  1.0
 *        Created:  17/05/2018 09:10:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dg539@cam.ac.uk
 *   Organization:  Department of Genetics, Cambridge University
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "ld.h"
#include "ld_io.h"
#include "kseq.h"


KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000)

ld_file_t * open_ld(const char *ld_fn)
{
	gzFile fp = *ld_fn == '-' ? gzdopen(fileno(stdin), "r") : gzopen(ld_fn, "r");
	if (!fp) return NULL;
	kstream_t *ks = ks_init(fp);
	ld_file_t *f = (ld_file_t *)calloc(1, sizeof(ld_file_t));
	f->fp = ks;
	return f;		
}


int close_ld(ld_file_t *fp)
{
	if (fp) {
		kstream_t *ks = (kstream_t *)fp->fp;
		if (ks) gzclose(ks->f),ks_destroy(ks);
		if (fp->buf.s) free(fp->buf.s);	
		return 0;
	} else 
		return 0;
}
/**
 * @func	parse_ld_t
 * @brief	parse ld records
 *
 */


int parse_ld_t(kstring_t *buf, ld_t *t)
{
	char *s = buf->s, *r, *q;
	int l = buf->l;
	int i, j;
	/*fprintf(stderr, "%s\n", buf->s);*/
	for( i = j = 0, q = s; i < l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		switch (j) {
			case 0:
				t->flags = strtoul(q, &r, 10);
				break;
			case 1:
				/*t->ctgn1 = strtoul(q, &r, 10);*/
				t->ctgn1 = q;
				break;
			case 2:
				t->pos1 = strtoul(q, &r, 10);
				break;
			case 3:
				/*t->ctgn2 = strtoul(q, &r, 10);*/
				t->ctgn2 = q;
				break;
			case 4:
				t->pos2 = strtoul(q, &r, 10);
				break;
			case 5:
				t->p1 = strtof(q, &r);
				break;
			case 6:
				t->p2 = strtof(q, &r);
				break;
			case 7:
				t->q1 = strtof(q, &r);
				break;
			case 8:
				t->q2 = strtof(q, &r);
				break;
			case 9:
				t->r = strtof(q, &r);
				break;
			case 10:
				t->r2 = strtof(q, &r);
				break;
				//d dprime
			case 11:
				t->d = strtof(q, &r);
				break;
			case 12:
				t->dprime = strtof(q, &r);
			case 13:
				t->p = strtod(q, &r);
				break;
			case 14:
				t->chiSqFisher = strtod(q, &r);
				break;
			case 15:
				t->chiSqModel = strtod(q, &r);
				break;
		}	
		++j, q = s + i + 1;	
	}
	if (j < 13) return -1;	
	else return 0;
}
int read_ld1(ld_file_t *ld_fp, ld_t *t)
{
	int ret, dret;
read_more:	
	ret = ks_getuntil((kstream_t *)ld_fp->fp, KS_SEP_LINE, &ld_fp->buf, &dret);
		if (ret < 0) return ret;
		if (ld_fp->buf.s[0] == '#') goto read_more;
		ret = parse_ld_t(&ld_fp->buf, t);
		if (ret < 0) goto read_more;	
		return ret;
}	
