#include <omp.h>
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "data_structures.h"
#include "sfunction.h"
#include "macro.h"


int cube_inside(dvec n, dvec cell, dvec x, double hblock)
{
  dvec n2;
  for(int d = 0; d < DIM; d++)
	  n2[d] = n[d] + cell[d];

  ivec f;

  for(int d = 0; d < DIM; d++) {
  
//    if (c[d]<=x[d]-hblock && x[d]+hblock<c2[d])
      f[d] = 1;
//    else
      if (x[d]-hblock<n[d])
        f[d] = 0;
      else
        if (x[d]+hblock>n2[d])
          f[d] = 0;
  }
  
  return f[0] && f[1] && f[2];
}

double cubes_intersect_vol(dvec n, dvec cell, dvec x, double hblock)
{
  dvec n2;
  for(int d = 0; d < DIM; d++)
	  n2[d] = n[d] + cell[d];

  dvec f;

  for(int d = 0; d < DIM; d++) {
  
//    if (c[d]<=x[d]-hblock && x[d]+hblock<c2[d])
      f[d] = hblock + hblock;
//    else
      if (x[d]-hblock<n[d])
        f[d] = (x[d] - n[d]) + hblock;  // /block
      else
        if (x[d]+hblock>n2[d])
          f[d] = (n2[d] - x[d]) + hblock;  // /block

    if (f[d]<0 || f[d]>hblock+hblock)
      f[d] = 0;
  }
  
  return f[0]*f[1]*f[2];
}

// double ro_interp(FHMD *fh, int k) {

// 	double x = (fh->step_MD % fh->FH_step) / (double) fh->FH_step;
//     return (1.0-x)*fh->arr[k].ro_fh + x*fh->arr_next[k].ro_fh;
// }

/* moved to mdlib/update.cpp for OpenMP parallelization
// maryna rho-interpretation

static double DRNOR() {
  double R1 = (double)rand() / ((double)(RAND_MAX) + 1.0);
  double R2 = (double)rand() / ((double)(RAND_MAX) + 1.0);
  double R11 = sqrt(2.0 * (-log(1.0 - R1)));
  double R22 = 2.0 * (3.1415926535897932384626433832795) * R2;

  return R11 * cos(R22);
}

double ro_interp(FHMD *fh, int k) {
    double x = (fh->step_MD % fh->FH_step) / (double) fh->FH_step;

    // Deterministic linear interpolation
    double ro_mean = (1.0 - x) * fh->arr[k].ro_fh + x * fh->arr_next[k].ro_fh;

    // Add thermal fluctuation to restore correct variance
    double fluct_rho = sqrt(1.0 / 3.0) * fh->std_rho * DRNOR();  // DRNOR() should return N(0,1)

    // return ro_mean + fluct_rho;
    return fh->arr_next[k].ro_fh;
}
*/
void fhmd_update_MD_in_FH(rvec x[], rvec v[], float mass[], rvec f[], int N_atoms, FHMD *fh)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    int        ind;

    double   block = 0.24 ; //0.24
    double   hblock = block*0.5;
    double   iblock3 = 1.0/(block*block*block);
	dvec     cell;
    for(int d = 0; d < DIM; d++)
		cell[d] = fh->box[d]/(double)(fh->N[d]);

    /* Collect statistics */
    for(int n = 0; n < fh->Ntot; n++) {
        arr[n].natoms = 0;
		arr[n].ro_md = 0;
		arr[n].u_md[0] = 0;
		arr[n].u_md[1] = 0;
		arr[n].u_md[2] = 0 ;
		arr[n].uro_md[0] = 0;
		arr[n].uro_md[1] = 0;
		arr[n].uro_md[2] = 0;
    }

    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        for(int d = 0; d < DIM; d++)
        {
            fh->indv[n][d] = (int)(xn[d]/cell[d]) + fh->N_shift[d];
        }

        ind = I(fh->indv[n], fh->N);
if (ind<0 || ind>=fh->Ntot) {
	printf("outofbounds %d %d %f %f %f  ",n,ind,xn[0],xn[1],xn[2]);
	ind=0;
}

        fh->ind[n] = ind;

		if (cube_inside(fh->grid.n[ind],cell,xn,hblock))
		{
		    arr[ind].natoms += 1.0;
			arr[ind].ro_md += fh->m_hybr[n];
			arr[ind].uro_md[0] += v[n][0]*mass[n];
			arr[ind].uro_md[1] += v[n][1]*mass[n];
			arr[ind].uro_md[2] += v[n][2]*mass[n];
            // fprintf(stdout, "cube_inside");
            // fprintf(stdout, "volume in mass_star: K: %d, volume: %.5f\n", n, arr[n].vol_star);
		}
		else
		{
			for (int di=-1; di<=1; di++)
				for (int dj=-1; dj<=1; dj++)
					for (int dk=-1; dk<=1; dk++)
					{
			        ivec indv2;
		            indv2[0] = fh->indv[n][0] + di;
		            indv2[1] = fh->indv[n][1] + dj;
		            indv2[2] = fh->indv[n][2] + dk;
/*			        for(int d = 0; d < DIM; d++)
			            if (indv2[d]<0)
			            	indv2[d] += fh->N_md[d];
			            else
				            if (indv2[d]>=fh->N_md[d])
				            	indv2[d] -= fh->N_md[d];
no need, this is included in I() */
			        int ind2 = I(indv2, fh->N);

// we need negative bottom-left corner n[] for cells like (-1,-1,-1), not PBC'ed!!
                    dvec n2;
                    n2[0] = fh->grid.n[ind][0] + di*cell[0];
                    n2[1] = fh->grid.n[ind][1] + dj*cell[1];
                    n2[2] = fh->grid.n[ind][2] + dk*cell[2];
					double fr = cubes_intersect_vol(n2/*fh->grid.n[ind2]*/,cell,xn,hblock) * iblock3;
				    arr[ind2].natoms += fr;
                    // fprintf(stdout, "cube_inside");


				    fr *= fh->m_hybr[n];
					arr[ind2].ro_md += fr;
					arr[ind2].uro_md[0] += v[n][0]*fr;
					arr[ind2].uro_md[1] += v[n][1]*fr;
					arr[ind2].uro_md[2] += v[n][2]*fr;
					}
		}

	}

    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md *= fh->grid.ivol[i];
        for(int d = 0; d < DIM; d++)
            arr[i].uro_md[d] *= fh->grid.ivol[i];
    }
    
}

void fhmd_sum_arrays(t_commrec *cr, FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    /* Pack FHMD arrays to linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->mpi_linear[i]              = arr[i].natoms;
        fh->mpi_linear[i + fh->Ntot]   = arr[i].ro_md;
        fh->mpi_linear[i + fh->Ntot*2] = arr[i].uro_md[0];
        fh->mpi_linear[i + fh->Ntot*3] = arr[i].uro_md[1];
        fh->mpi_linear[i + fh->Ntot*4] = arr[i].uro_md[2];
    }
    /* Broadcast linear array */
    gmx_sumd(fh->Ntot*5, fh->mpi_linear, cr);

    /* Unpack linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].natoms    = fh->mpi_linear[i];
        arr[i].ro_md     = fh->mpi_linear[i + fh->Ntot];
        arr[i].uro_md[0] = fh->mpi_linear[i + fh->Ntot*2];
        arr[i].uro_md[1] = fh->mpi_linear[i + fh->Ntot*3];
        arr[i].uro_md[2] = fh->mpi_linear[i + fh->Ntot*4];
    }

//	double nn=0;    
//    for(int i = 0; i < fh->Ntot; i++)
//        nn += arr[i].natoms;
//    for(int i = 0; i < 10; i++)
//        printf("%.3f ", arr[i].uro_md[0]);
//    if (MASTER(cr) && (fh->Noutput > 0))
//      if(!(fh->step_MD % fh->Noutput))
//        fprintf(stderr,"%f ", nn);
//        fprintf(stderr,"y%.3f ", arr[0].u_md[0]);
//        fprintf("y%d ", arr[2].natoms);
}

void fhmd_update_MD_in_FH2(rvec x[], rvec v[], float mass[], rvec f[], int N_atoms, FHMD *fh)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    int        ind;

return;

    double   block = 0.1; //0.24;
    double   hblock = block*0.5;
    double   iblock3 = 1.0/(block*block*block);
	dvec     cell;
    for(int d = 0; d < DIM; d++)
		cell[d] = fh->box[d]/(double)(fh->N_md[d]);

    /* Collect statistics */
    for(int n = 0; n < fh->Ntot; n++)
        arr[n].ro_mixed = 0;

    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        for(int d = 0; d < DIM; d++)
        {
            fh->indv[n][d] = (int)(xn[d]/cell[d]) + fh->N_shift[d];
        }

        ind = I(fh->indv[n], fh->N);

if (ind<0 || ind>=fh->Ntot) {
	printf("outofbounds %d %d %f %f %f  ",n,ind,xn[0],xn[1],xn[2]);
	ind=0;
}
        fh->ind[n] = ind;

		if (cube_inside(fh->grid.n[ind],cell,xn,hblock))
		{
			arr[ind].ro_mixed += fh->ro_hybr[n];
		}
		else
		{
			for (int di=-1; di<=1; di++)
				for (int dj=-1; dj<=1; dj++)
					for (int dk=-1; dk<=1; dk++)
					{
			        ivec indv2;
		            indv2[0] = fh->indv[n][0] + di;
		            indv2[1] = fh->indv[n][1] + dj;
		            indv2[2] = fh->indv[n][2] + dk;
/*			        for(int d = 0; d < DIM; d++)
			            if (indv2[d]<0)
			            	indv2[d] += fh->N_md[d];
			            else
				            if (indv2[d]>=fh->N_md[d])
				            	indv2[d] -= fh->N_md[d];
no need, this is included in I() */
			        int ind2 = I(indv2, fh->N);

// we need negative bottom-left corner n[] for cells like (-1,-1,-1), not PBC'ed!!
                    dvec n2;
                    n2[0] = fh->grid.n[ind][0] + di*cell[0];
                    n2[1] = fh->grid.n[ind][1] + dj*cell[1];
                    n2[2] = fh->grid.n[ind][2] + dk*cell[2];
					double fr = cubes_intersect_vol(n2/*fh->grid.n[ind2]*/,cell,xn,hblock) * iblock3;

					arr[ind2].ro_mixed += fh->ro_hybr[n];
					}
		}

	}
}

void fhmd_sum_arrays2(t_commrec *cr, FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    /* Pack FHMD arrays to linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->mpi_linear[i]   =   arr[i].ro_mixed;
    }
    /* Broadcast linear array */
    gmx_sumd(fh->Ntot*5, fh->mpi_linear, cr);

    /* Unpack linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_mixed   =   fh->mpi_linear[i];
    }
}

void fhmd_calculate_MDFH_terms(rvec x[], rvec v[], float mass[], rvec f[], int N_atoms, FHMD *fh, int th)
{
// moved to mdlib/update.cpp for OpenMP parallelization
}

