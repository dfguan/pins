/*
 * =====================================================================================
 *
 *       Filename:  make_brk.h
 *
 *    Description:  make_brk header
 *
 *        Version:  1.0
 *        Created:  14/11/2019 15:28:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef MAKE_BRK_H
#define MAKE_BRK_H

#ifdef __cplusplus
extern "C" {
#endif
int main_brks(int argc, char *argv[]);
int mk_brks(char *sat_fn, char *links_fn, int limn, char *out_fn);
#ifdef __cplusplus
}
#endif

#endif
