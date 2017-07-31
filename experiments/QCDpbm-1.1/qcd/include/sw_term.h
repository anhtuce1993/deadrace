
/*******************************************************************************
*
* File sw_term.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef SW_TERM_H
#define SW_TERM_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef START_H
#include "start.h"
#endif

#ifndef PAULI_C
extern void mul_pauli(pauli *m,weyl *s,weyl *r);
extern void assign_pauli(int vol,pauli_dble *md,pauli *m);
extern void apply_sw(int vol,pauli *m,spinor *s,spinor *r);
#endif

#ifndef PAULI_DBLE_C
extern void mul_pauli_dble(pauli_dble *m,weyl_dble *s,weyl_dble *r);
extern int inv_pauli_dble(pauli_dble *m,pauli_dble *im);
extern double det_pauli_dble(pauli_dble *m);
extern void apply_sw_dble(int vol,pauli_dble *m,spinor_dble *s,spinor_dble *r);
#endif

#ifndef SWINIT_C
extern void alloc_sw(void);
extern void alloc_swd(void);
extern void free_sw(void);
extern void free_swd(void);
extern void assign_swd2sw(void);
extern int invert_swd(ptset_t set);
#endif

#ifndef SW_TERM_C
extern void sw_term(void);
extern void free_swbufs(void);
#endif

#endif
