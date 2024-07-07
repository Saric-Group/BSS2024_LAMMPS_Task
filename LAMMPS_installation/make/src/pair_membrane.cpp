/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "pair_membrane.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include <string.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMembrane::PairMembrane(LAMMPS *lmp) : Pair(lmp)
{

  single_enable = 0;  // 1 if single() routine exists
}

/* ---------------------------------------------------------------------- */

PairMembrane::~PairMembrane()
{
  if (allocated) {

    // pw359: NOTE that you will get a segfault, if you don't allocate cutsq.
    memory->destroy(cutsq);
    memory->destroy(setflag);
    memory->destroy(rc);
    memory->destroy(rc_sq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(rmin);
    memory->destroy(zeta);
    memory->destroy(nu);
    memory->destroy(sintheta0);
  }
}

/* ---------------------------------------------------------------------- */

void PairMembrane::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul, fx,fy,fz;
  double rsq,rinv,r2inv,r6inv,r3inv,r5inv,r7inv;
  double forcex,forcey,forcez,crossx,crossy,crossz;
  double forcelj_uR, forcelj_uA,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double gradax, graday, gradaz, gradPhix, gradPhiy, gradPhiz;
  double nirhat, njrhat, nijrhat, rhatx, rhaty, rhatz;
  double normmui, normmuj, nix, niy, niz, njx, njy, njz, nijx, nijy, nijz;
  double a, Phi, r, ftni, cosarg, powcosarg, arg, rminsqorsq, uA, uR;
  double nicrossrhatx, nicrossrhaty, nicrossrhatz, njcrossrhatx, njcrossrhaty, njcrossrhatz;
  double nicrossanix, nicrossaniy, nicrossaniz, njcrossanjx, njcrossanjy, njcrossanjz,
         gradanix, gradaniy, gradaniz, gradanjx, gradanjy, gradanjz;


  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i     = ilist[ii];
    xtmp  = x[i][0];
    ytmp  = x[i][1];
    ztmp  = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum  = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j         = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j        &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        r2inv = 1.0/rsq;
        rinv  = sqrt(r2inv);
        r     = 1.0/rinv;

        if (r < rc[itype][jtype]) {

          // calculate Phi and dPhi/dri

          forcex = forcey = forcez = 0.0;
          a = uA  = uR  = Phi = 0;
          gradPhix    = gradPhiy    = gradPhiz = 0.0;
          nicrossanix = nicrossaniy = nicrossaniz = 0;
          njcrossanjx = njcrossanjy = njcrossanjz = 0;
          

          rhatx = delx * rinv;
          rhaty = dely * rinv;
          rhatz = delz * rinv;

          // NOTE: no contribution if magnitude of orientation vector is 0
          if (mu[i][3] > 0.0 && mu[j][3] > 0.0) {


            nix     = mu[i][0] / mu[i][3];
            niy     = mu[i][1] / mu[i][3];
            niz     = mu[i][2] / mu[i][3];

            njx     = mu[j][0] / mu[j][3];
            njy     = mu[j][1] / mu[j][3];
            njz     = mu[j][2] / mu[j][3];

  
            nijx    = nix - njx;
            nijy    = niy - njy;
            nijz    = niz - njz;

            nirhat   = nix *rhatx  + niy *rhaty  + niz *rhatz;
            njrhat   = njx *rhatx  + njy *rhaty  + njz *rhatz;
            nijrhat  = nijx*rhatx  + nijy*rhaty  + nijz*rhatz;

            gradax    = 2.0 * nirhat * njrhat * rhatx - nix * njrhat  - njx * nirhat + sintheta0[itype][jtype] * (nijrhat * rhatx - nijx);
            graday    = 2.0 * nirhat * njrhat * rhaty - niy * njrhat  - njy * nirhat + sintheta0[itype][jtype] * (nijrhat * rhaty - nijy);
            gradaz    = 2.0 * nirhat * njrhat * rhatz - niz * njrhat  - njz * nirhat + sintheta0[itype][jtype] * (nijrhat * rhatz - nijz);

            gradax   *= rinv;
            graday   *= rinv;
            gradaz   *= rinv;

            gradPhix  = nu[itype][jtype] * gradax;
            gradPhiy  = nu[itype][jtype] * graday;
            gradPhiz  = nu[itype][jtype] * gradaz;


            nicrossrhatx = niy * rhatz - niz * rhaty;
            nicrossrhaty = niz * rhatx - nix * rhatz;
            nicrossrhatz = nix * rhaty - niy * rhatx;

            njcrossrhatx = njy * rhatz - njz * rhaty;
            njcrossrhaty = njz * rhatx - njx * rhatz;
            njcrossrhatz = njx * rhaty - njy * rhatx;

            a   =   nicrossrhatx * njcrossrhatx + nicrossrhaty * njcrossrhaty +  nicrossrhatz * njcrossrhatz 
                  + sintheta0[itype][jtype] * (-nijrhat -sintheta0[itype][jtype]);

            Phi = 1.0 + nu[itype][jtype] * (a - 1.0);


            // torque

            gradanix = njx - rhatx * (njrhat + sintheta0[itype][jtype]);
            gradaniy = njy - rhaty * (njrhat + sintheta0[itype][jtype]);
            gradaniz = njz - rhatz * (njrhat + sintheta0[itype][jtype]);
            
            nicrossanix = niy * gradaniz - niz * gradaniy;
            nicrossaniy = niz * gradanix - nix * gradaniz;
            nicrossaniz = nix * gradaniy - niy * gradanix;

      
            if (newton_pair || j < nlocal) {

              gradanjx    = nix + rhatx * (-nirhat + sintheta0[itype][jtype]);
              gradanjy    = niy + rhaty * (-nirhat + sintheta0[itype][jtype]);
              gradanjz    = niz + rhatz * (-nirhat + sintheta0[itype][jtype]);

              njcrossanjx = njy * gradanjz - njz * gradanjy;
              njcrossanjy = njz * gradanjx - njx * gradanjz;
              njcrossanjz = njx * gradanjy - njy * gradanjx;

            }
          } 

          rminsqorsq = rmin[itype][jtype] * rmin[itype][jtype] * r2inv;

          if (rminsqorsq > 1.0) {

            // Eq (7): grad(u_R) 
            forcelj_uR = 4.*epsilon[itype][jtype] * rminsqorsq * ( 1. - rminsqorsq) * rinv;

            // Eq (6): -grad(u_R) + epsilon * grad(Phi)
            forcex -=   forcelj_uR * rhatx  - epsilon[itype][jtype] * gradPhix;
            forcey -=   forcelj_uR * rhaty  - epsilon[itype][jtype] * gradPhiy;
            forcez -=   forcelj_uR * rhatz  - epsilon[itype][jtype] * gradPhiz;

            // factor for torque
            ftni       = -epsilon[itype][jtype];
            uR         = epsilon[itype][jtype] * rminsqorsq * (rminsqorsq - 2.0); 

          }

          else {

            arg        =  M_PI/2.0 * (r - rmin[itype][jtype])/(rc[itype][jtype]-rmin[itype][jtype]);
            cosarg     =  cos(arg);
            powcosarg  =  pow(cosarg, 2.0*zeta[itype][jtype]-1.0);
            forcelj_uA =  zeta[itype][jtype]*M_PI*epsilon[itype][jtype]/(rc[itype][jtype]-rmin[itype][jtype]) 
                        * powcosarg           
                        * sin(arg);
            uA         =  -epsilon[itype][jtype] * powcosarg*cosarg;
            ftni       = uA;

            forcex -=  forcelj_uA * rhatx * Phi +    uA * gradPhix;
            forcey -=  forcelj_uA * rhaty * Phi +    uA * gradPhiy;
            forcez -=  forcelj_uA * rhatz * Phi +    uA * gradPhiz;

          }  

          fx = forcex * factor_lj;
          fy = forcey * factor_lj;
          fz = forcez * factor_lj;
          
          f[i][0]  += fx; 
          f[i][1]  += fy;
          f[i][2]  += fz;

          torque[i][0] += -nu[itype][jtype] * nicrossanix * ftni;
          torque[i][1] += -nu[itype][jtype] * nicrossaniy * ftni;
          torque[i][2] += -nu[itype][jtype] * nicrossaniz * ftni;

          if (newton_pair || j < nlocal) {

            f[j][0] -= fx;
            f[j][1] -= fy;
            f[j][2] -= fz;

            torque[j][0] += -nu[itype][jtype] * njcrossanjx * ftni;
            torque[j][1] += -nu[itype][jtype] * njcrossanjy * ftni;
            torque[j][2] += -nu[itype][jtype] * njcrossanjz * ftni;

          }

          evdwl = ecoul = 0;
          if (eflag) {
            if (rminsqorsq > 1.0)   
               evdwl = uR + epsilon[itype][jtype] *(1.0 - Phi);
            else 
              evdwl = uA * Phi;

            evdwl  *= factor_lj;
          }

          if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,ecoul,fx,fy,fz,delx,dely,delz);
        }
      } 
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMembrane::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,  n+1,  n+1,  "pair:cutsq");
  memory->create(rc,     n+1,  n+1,  "pair:rc");
  memory->create(rc_sq,  n+1,  n+1,  "pair:rc_sq");
  memory->create(epsilon,n+1,  n+1,  "pair:epsilon");
  memory->create(sigma,  n+1,  n+1,  "pair:sigma");
  memory->create(rmin,   n+1,  n+1,  "pair:rmin");
  memory->create(zeta,   n+1,  n+1,  "pair:zeta");
  memory->create(nu,     n+1,  n+1,  "pair:nu");
  memory->create(sintheta0, n+1,  n+1,  "pair:sintheta0");

}




/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMembrane::settings(int narg, char **arg)
{

  if (narg !=  1)
    error->all(FLERR,"Incorrect args in pair_style membrane command");

  // global cutoff 

  rc_global  = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = i+1; j <= atom->ntypes; j++) {
        if (setflag[i][j]) {
          rc[i][j] = rc_global;
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMembrane::coeff(int narg, char **arg)
{

  if (narg != 9)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one   = force->numeric(FLERR,arg[2]);
  double sigma_one     = force->numeric(FLERR,arg[3]);
  double rmin_one      = force->numeric(FLERR,arg[4]);
  double rc_one        = force->numeric(FLERR,arg[5]);
  double zeta_one      = force->numeric(FLERR,arg[6]);
  double nu_one        = force->numeric(FLERR,arg[7]);
  double theta0_one    = force->numeric(FLERR,arg[8]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j]   = epsilon_one;
      sigma[i][j]     = sigma_one;
      rmin[i][j]      = rmin_one;
      rc[i][j]        = rc_one;
      zeta[i][j]      = zeta_one;
      nu[i][j]        = nu_one;
      sintheta0[i][j] = sin(theta0_one);
      setflag[i][j]   = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMembrane::init_style()
{
  // q could potentially be useful to encode any particle specific information. 
  // mu = n
  // NOTE: At the moment I just require these quantities to be set.
  if (!atom->q_flag || !atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair membrane requires atom attributes q, n, torque");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMembrane::init_one(int i, int j)
{
  //pw359:  This function will be called for i<=j. 
  //        The lower diagonal part of the matrix needs to be inialised. 

  epsilon[j][i]   = epsilon[i][j];
  sigma[j][i]     = sigma[i][j];
  rmin[j][i]      = rmin[i][j];
  rc[j][i]        = rc[i][j];
  zeta[j][i]      = zeta[i][j];
  nu[j][i]        = nu[i][j];
  sintheta0[j][i] = sintheta0[i][j];
  rc_sq[i][j]     = rc[i][j] * rc[i][j];
  rc_sq[j][i]     = rc_sq[i][j];

  return rc[i][j];
}


