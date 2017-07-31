
/*******************************************************************************
*
* File dirac.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef DIRAC_H
#define DIRAC_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef PBND_C
extern void (*assign_s2w[])(int imb[],int vol,spinor *pk,weyl *pl);
extern void (*add_assign_w2s[])(int imb[],int vol,weyl *pk,spinor *pl);
extern void (*sub_assign_w2s[])(int imb[],int vol,weyl *pk,spinor *pl);
extern void (*mulg5_sub_assign_w2s[])(int imb[],int vol,weyl *pk,spinor *pl);
#endif

#ifndef PBND_DBLE_C
extern void (*assign_sd2wd[])(int imb[],int vol,spinor_dble *pk,
                              weyl_dble *pl);
extern void (*add_assign_wd2sd[])(int imb[],int vol,weyl_dble *pk,
                                  spinor_dble *pl);
extern void (*sub_assign_wd2sd[])(int imb[],int vol,weyl_dble *pk,
                                  spinor_dble *pl);
extern void (*mulg5_sub_assign_wd2sd[])(int imb[],int vol,weyl_dble *pk,
                                        spinor_dble *pl);
#endif

#ifndef QHAT_C
extern void Qhat(int k,int l);
extern void Qnohat(int k,int l);
extern void Qoe(int k,int l);
extern void Qeo(int k,int l);
#endif

#ifndef QHAT_DBLE_C
extern void Qhat_dble(int k,int l);
extern void Qnohat_dble(int k,int l);
extern void Qoe_dble(int k,int l);
extern void Qeo_dble(int k,int l);
#endif

#ifndef SCOM_C
extern void cps_int_bnd(int k);
extern void cps_ext_bnd(int k);
extern void free_sbufs(void);
#endif

#ifndef SDCOM_C
extern void cpsd_int_bnd(int k);
extern void cpsd_ext_bnd(int k);
extern void free_sdbufs(void);
#endif

#endif
