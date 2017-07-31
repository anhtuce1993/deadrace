/*******************************************************************************
*
* File net.h
*
* Copyright (C) 2008 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#include <stdlib.h>

#define MIN_PACKET_SIZE 128
#define MAX_PACKET_SIZE 1024

extern int safe_mod(int x,int y);
extern void error(int test,int no,char *name,char *format,...);
extern void error_root(int test,int no,char *name,char *format,...);
extern char *alloc_buffer(size_t n);
extern void init_buffer(int n,double *buf);

