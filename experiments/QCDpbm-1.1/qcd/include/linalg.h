
/*******************************************************************************
*
* File linalg.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef LINALG_H
#define LINALG_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef LINALG_C
extern complex spinor_prod(int vol,int icom,spinor *pk,spinor *pl);
extern float norm_square(int vol,int icom,spinor *pk);
extern void mulc_spinor_add(int vol,spinor *pk,spinor *pl,complex z);
extern float normalize(int vol,int icom,spinor *pk);
extern void mulg5(int vol,spinor *pk);
extern void mulmg5(int vol,spinor *pk);
#endif

#ifndef LINALG_DBLE_C
extern complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *pk,
                                     spinor_dble *pl);
extern double norm_square_dble(int vol,int icom,spinor_dble *pk);
extern void mulc_spinor_add_dble(int vol,spinor_dble *pk,spinor_dble *pl,
                                 complex_dble z);
extern double normalize_dble(int vol,int icom,spinor_dble *pk);
extern void mulg5_dble(int vol,spinor_dble *pk);
extern void mulmg5_dble(int vol,spinor_dble *pk);
#endif

#endif
