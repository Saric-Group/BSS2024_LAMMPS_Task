/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.



   This class implements a pair style for the membrane model described in
   Yuan et al., Phys. Rev. E 82, 011905 (2010). 
   DOI: 10.1103/PhysRevE.82.011905

------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(membrane,PairMembrane)

#else

#ifndef LMP_PAIR_MEMBRANE_H
#define LMP_PAIR_MEMBRANE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMembrane : public Pair {
 public:
  PairMembrane(class LAMMPS *);
  virtual ~PairMembrane();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:

  void   allocate();

  double rc_global;
  double **rc, **rc_sq;
  double **epsilon;
  double **sigma;
  double **rmin;
  double **sintheta0;
  double **zeta;
  double **nu;

};

}

#endif
#endif

