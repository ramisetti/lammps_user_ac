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

#include <stdlib.h>
#include <string.h>
#include "fix_momentum_new.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

FixMomentumNew::FixMomentumNew(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix momentum command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix momentum command");

  dynamic = linear = angular = rescale = 0;
  vcm_factor=1.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"linear") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix momentum command");
      linear = 1;
      xflag = force->inumeric(FLERR,arg[iarg+1]);
      yflag = force->inumeric(FLERR,arg[iarg+2]);
      zflag = force->inumeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"angular") == 0) {
      angular = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"rescale") == 0) {
      rescale = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      vcm_factor = 1.0/force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix momentum command");
  }

  if (linear == 0 && angular == 0)
    error->all(FLERR,"Illegal fix momentum command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1)
      error->all(FLERR,"Illegal fix momentum command");

  dynamic_group_allow = 1;
/*  if (group->dynamic[igroup]) {
    dynamic = 1;
  } else {
   if (group->count(igroup) == 0)
     error->all(FLERR,"Fix momentum group has no atoms");
  }
*/
  // vcm_factor=1.0;
  // if(narg == 5)
  // 	vcm_factor = 1.0/force->numeric(FLERR,arg[5]);
}

/* ---------------------------------------------------------------------- */

int FixMomentumNew::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMomentumNew::init()
{
  if (group->dynamic[igroup]) {
    dynamic = 1;
  } else {
   if (group->count(igroup) == 0)
     error->all(FLERR,"Fix momentum group has no atoms");
  }

  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixMomentumNew::end_of_step()
{
  double **v = atom->v;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  double ekin_old,ekin_new;
  ekin_old = ekin_new = 0.0;

  if (dynamic)
    masstotal = group->mass(igroup);

  // do nothing is group is empty, i.e. mass is zero;

  if (masstotal == 0.0) return;

  // compute kinetic energy before momentum removal, if needed

  if (rescale) {

    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    double ke=0.0;

    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          ke += rmass[i] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          ke +=  mass[type[i]] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    }
    MPI_Allreduce(&ke,&ekin_old,1,MPI_DOUBLE,MPI_SUM,world);
  }

  if (linear) {
    double vcm[3];
    group->vcm(igroup,masstotal,vcm);

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xflag) v[i][0] -= vcm[0]*vcm_factor;
        if (yflag) v[i][1] -= vcm[1]*vcm_factor;
        if (zflag) v[i][2] -= vcm[2]*vcm_factor;
      }
  }

  if (angular) {
    double xcm[3],angmom[3],inertia[3][3],omega[3];
    group->xcm(igroup,masstotal,xcm);
    group->angmom(igroup,xcm,angmom);
    group->inertia(igroup,xcm,inertia);
    group->omega(angmom,inertia,omega);

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    double **x = atom->x;
    double **v = atom->v;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    double dx,dy,dz;
    double unwrap[3];

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        domain->unmap(x[i],image[i],unwrap);
        dx = unwrap[0] - xcm[0];
        dy = unwrap[1] - xcm[1];
        dz = unwrap[2] - xcm[2];
        v[i][0] -= omega[1]*dz - omega[2]*dy;
        v[i][1] -= omega[2]*dx - omega[0]*dz;
        v[i][2] -= omega[0]*dy - omega[1]*dx;
      }
  }

  // compute kinetic energy after momentum removal, if needed

  if (rescale) {

    double ke=0.0, factor=1.0;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;

    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          ke += rmass[i] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          ke +=  mass[type[i]] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    }
    MPI_Allreduce(&ke,&ekin_new,1,MPI_DOUBLE,MPI_SUM,world);

    if (ekin_new != 0.0) factor = sqrt(ekin_old/ekin_new);
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] *= factor;
        v[i][1] *= factor;
        v[i][2] *= factor;
      }
    }
  }
}
