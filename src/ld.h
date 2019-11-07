#ifndef _LD_H
#define _LD_H

#include <stdint.h>
typedef struct {
	//uint32_t    ctgn1;
	//uint32_t    Amissing: 1, Aphased: 1, Aposition: 30;
	//uint32_t    ctgn2;
	//uint32_t    Bmissing: 1, Bphased: 1, Bposition: 30;
	
	uint16_t    flags;
	char		*ctgn1;
	uint32_t	pos1;
	char		*ctgn2;
	uint32_t	pos2;
	float		p1, p2, q1, q2;
	float		d, dprime; // D and D'
	float		r, r2;     // Correlation coefficient
	double		p;         // P-value
	double		chiSqFisher;
	double		chiSqModel;		
} ld_t;

#endif
