
/*******************************************************************************
*
* File parms.c
*
* Copyright (C) 2005, 2007 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Lattice parameters
*
* The externally accessible functions are
*
*   lat_parms_t set_lat_parms(double beta,double kappa,double csw)
*     Sets the basic lattice parameters. The return value is a structure
*     that contains the lattice parameters
*
*   lat_parms_t lat_parms(void)
*     Returns the current lattice parameters in a structure with elements
*     *.beta,*.kappa,*.m0,*.csw where m0 is bare quark mass determined by
*     kappa
*
* Notes:
*
* The lattice parameters must be set simultaneously on all processes. The
* type lat_parms_t is defined in the file start.h
*
*******************************************************************************/

#define PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "start.h"
#include "global.h"

static lat_parms_t lat={0.0,0.0,0.0,0.0};


lat_parms_t set_lat_parms(double beta,double kappa,double csw)
{
   double dprms[3];

   if (NPROC>1)
   {
      dprms[0]=beta;
      dprms[1]=kappa;
      dprms[2]=csw;

      MPI_Bcast(dprms,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((dprms[0]!=beta)||(dprms[1]!=kappa)||(dprms[2]!=csw),1,
             "set_lat_parms [parms.c]","Parameters are not global");
   }

   error_root(kappa<=0.0,1,
              "set_lat_parms [parms.c]","kappa must be positive");

   lat.beta=beta;
   lat.kappa=kappa;
   lat.m0=1.0/(2.0*kappa)-4.0;
   lat.csw=csw;

   return lat;
}


lat_parms_t lat_parms(void)
{
   return lat;
}

