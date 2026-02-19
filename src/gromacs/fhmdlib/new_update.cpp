#include "data_structures.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "interpolation.h"
#include "macro.h"
#include "sfunction.h"

// include "gromacs/hhmd_G/libtessel.h"

int start_step = 0;  // VF: was 500

double norm(double vect[3]) {
  return sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
}

static double DRNOR() {
  double R1 = (double)rand() / ((double)(RAND_MAX) + 1.0);
  double R2 = (double)rand() / ((double)(RAND_MAX) + 1.0);
  double R11 = sqrt(2.0 * (-log(1.0 - R1)));
  double R22 = 2.0 * (3.1415926535897932384626433832795) * R2;

  return R11 * cos(R22);
}

// double u_interp(FHMD* fh, int k, int d) {
//   double x = (fh->step_MD % fh->FH_step) / (double)fh->FH_step;
//   return (1.0 - x) * fh->arr[k].u_fh[d] + x * fh->arr_next[k].u_fh[d];
// }

double u_interp(FHMD *fh, int k, int d)
{
    double x = (fh->step_MD % fh->FH_step) / (double)fh->FH_step;

    // Mean (linear interpolation)
    double u_mean =
        (1.0 - x) * fh->arr[k].u_fh[d] +
         x        * fh->arr_next[k].u_fh[d];

    // Restore missing variance
    double sigma_missing =
        sqrt(2.0 * x * (1.0 - x)) * fh->std_u;

    double fluct_u = sigma_missing * DRNOR();

    return u_mean + 1.0 * fluct_u;  //  empirical factor to prevent divergence, can be adjusted or removed
}

double s_k_t(FHMD* fh, double t, double k, double S) {
  double a = log(2.0) / log(t);  // 0.5 - so starts from 0.5
  return 1.0 / (1.0 + pow(pow(S, a) - 1.0, k));
}

/*double u_deriv(FHMD *fh, int k, int d) {

    return (fh->arr[k].u_fh[d] - fh->arr_prev[k].u_fh[d]) / fh->dt_FH;
}*/

void fhmd_do_update_md(int start, int nrend, double dt, int nstpcouple,
                       t_grp_tcstat* tcstat, double nh_vxi[], gmx_bool bNEMD,
                       t_grp_acc* gstat, rvec accel[], ivec nFreeze[],
                       float invmass[], unsigned short ptype[],
                       unsigned short cFREEZE[], unsigned short cACC[],
                       unsigned short cTC[], rvec x[], rvec xprime[], rvec v[],
                       rvec f[], matrix M, gmx_bool bNH, gmx_bool bPR,
                       t_commrec* cr, FHMD* fh) {
  double w_dt;
  int gf = 0, ga = 0, gt = 0;
  rvec vrel;
  float /*vn,*/ vv, va, vb, vnrel;
  float lg, vxi = 0, u;
  int n, d;

  /* FHMD variables */
  FH_arrays* arr = fh->arr;
  int ind;
  double invro_dt;
  double S;
  double k = 0.5, t = 7.0;  // for lambda
  double gamma_u, gamma_x;
  int nbr[8];
  dvec xi;
  const double g_eps = 1e-10;
  double lambda = 1.0;
  // double tau = 1; //#for s=0.5 was 0.1 - artefact, but temp is ok
  // double T_0 = -1;
  double tau = fh->tau_m;

  // double lambda_s = s_k_t(fh, 0.5, 7, !S!);
  double lambda_s = 0;

  // Declare T_0 outside of the function to keep its value across multiple calls
  static double T_0 = 0.730094;  // T_0 is initialized only once

  FILE* th;
  char fname[256];

  // Build filename with parameters (no rounding)
  if (start == 0)  //  the first thread
  {
    snprintf(fname, sizeof(fname), "output/thermo_S%f_tau%f_fluct%f.dat", fh->S,
             fh->tau_m, fh->fluct_m);

    // Open the file for appending
    th = fopen(fname, "a");
    fprintf(th, "step = %d\n", fh->step_MD);

    // th2 = fopen("thermo_vel.dat", "a");
  }

  // Ensure that step_MD > 0 for processing
  if (fh->step_MD > 0) {
    // Print the current step to the file
    // fprintf(th, "step = %d\n", fh->step_MD);

    // Initialize variables for momentum and energy calculations
    dvec total_mol_p = {0.0, 0.0, 0.0};  // Initialize total momentum
    double total_kinetic_energy = 0.0;
    // double lambda;

    // Calculate the total momentum and kinetic energy
    for (n = start; n < nrend; n++) {
      for (d = 0; d < DIM; d++) {
        total_mol_p[d] +=
            fh->p_mol_m[n][d];  // Summing momentum for each dimension
        total_kinetic_energy += (fh->p_mol_m[n][d] * fh->p_mol_m[n][d]) *
                                invmass[n] / 2.0;  // p^2 / 2m
      }
    }

    // Calculate temperature based on kinetic energy
    double temperature = (2.0 / 3.0) * (total_kinetic_energy / (nrend - start));

    // Set T_0 only during the first step (step 1)
    if (fh->step_MD == 1) {
      T_0 = 0.730094;  // Set reference temperature
                       // fprintf(th, "T_0 in if: %f\n", T_0);
    }
    lambda = sqrt(1.0 + (dt / tau) * (T_0 / temperature - 1.0));
    // lambda = 1.0;
    // Calculate lambda using the reference temperature T_0
    // if (T_0 > 0) {
    //     lambda = sqrt(1.0 + (dt / tau) * (T_0 / temperature - 1.0));  //
    //     Scaling factor for velocity
    // } else {
    //     lambda = 1.0;  // Default to 1 if T_0 hasn't been set yet
    // }

    if (start == 0)  //  the first thread
    {
      // fprintf(th, "step = %d\n", fh->step_MD);
      fprintf(th, "Total momentum: x = %f, y = %f, z = %f\n", total_mol_p[0],
              total_mol_p[1], total_mol_p[2]);

      fprintf(th, "Total Temperature: %f\n", temperature);
      fprintf(th, "T_0: %f\n", T_0);
      fprintf(th, "Total kinetic energy: %f\n", total_kinetic_energy);
      fprintf(th, "lambda: %f\n", lambda);
      fprintf(th, "lambda_s, %f\n", lambda_s);
      fprintf(th, "dt: %f\n", dt);
      // fprintf(th, "tau: %f\n", tau);//check
      // fprintf(th, "fluct_mag: %f\n", fh->fluct_m);//check
    }
  }

  if (bNH || bPR) {
    /* Update with coupling to extended ensembles, used for
     * Nose-Hoover and Parrinello-Rahman coupling
     * Nose-Hoover uses the reversible leap-frog integrator from
     * Holian et al. Phys Rev E 52(3) : 2338, 1995
     */

    /* FHMD Error */
    printf(MAKE_RED
           "\nFHMD: ERROR: FH-MD coupling doesn't support Nose-Hoover and "
           "Parrinello-Rahman\n" RESET_COLOR "\n");
    exit(11);
  } else if (cFREEZE != NULL || nFreeze[0][XX] || nFreeze[0][YY] ||
             nFreeze[0][ZZ] || bNEMD) {
    /* Update with Berendsen/v-rescale coupling and freeze or NEMD */

    /* FHMD Error */
    printf(MAKE_RED
           "\nFHMD: ERROR: FH-MD coupling doesn't support freeze or "
           "NEMD\n" RESET_COLOR "\n");
    exit(12);
  } else {
    /* Plain update with Berendsen/v-rescale coupling */

    for (n = start; n < nrend; n++) {
      if ((ptype[n] != eptVSite) && (ptype[n] != eptShell)) {
        w_dt = invmass[n] * dt;
        ind = fh->ind[n];

        S = fh->S;
        if (fh->S_function == moving_sphere)
          S = fhmd_Sxyz_r(x[n], fh->protein_com,
                          fh);  // MD/FH sphere follows protein
        else if (fh->S_function == fixed_sphere)
          S = fhmd_Sxyz_r(x[n], fh->box05, fh);  // Fixed MD/FH sphere
        if (fh->step_MD < fh->nsteps_temper)     // VF:
          S = S * (double)fh->step_MD / (double)fh->nsteps_temper;

        if (cTC) {
          gt = cTC[n];
        }
        lg = tcstat[gt].lambda;  // Thermostat

        /* Local thermostat */
        if (fh->S_berendsen >= 0) {
          if (S > fh->S_berendsen) lg = 1;
        } else {
          lg = lg * (1 - pow(S, -fh->S_berendsen)) + pow(S, -fh->S_berendsen);
        }

        // // // //Maryna
        // FILE *th;
        // th = fopen("thermo.dat", "a");

        // if (fh->step_MD > 0){

        // 	fprintf(th,"step = %d\n", fh->step_MD);
        // 	dvec total_mol_p;
        // 	total_mol_p[0] = 0;
        // 	total_mol_p[1] = 0;
        // 	total_mol_p[2] = 0;

        // 	fprintf(th, "Initial total momentum: x = %f, y = %f, z =
        // %f\n", total_mol_p[0], total_mol_p[1], total_mol_p[2]);

        // // for (n = start; n < nrend; n++){
        // // 	for (d = 0; d < DIM; d++){
        // // 		// total_mol_p[d] += fh->p_mol_m[n][d];
        // // 		fprintf(th, "d: d = %d\n, n = %d\n", d, n);
        // // 	}
        // // 	}
        // 	// fprintf(th, "Final total momentum: x = %f, y = %f, z = %f\n",
        // total_mol_p[0], total_mol_p[1], total_mol_p[2]); 	fprintf(th, "c:
        // c = %d\n",  c);
        // }
        // FILE *th;
        // th = fopen("thermo.dat", "a");
        // if (fh->step_MD > 0) {
        //   // Open file for printing

        //   // Print the current step to the file
        //   fprintf(th, "step = %d\n", fh->step_MD);

        //   // Initialize counter for iterations
        //   double c = 0;

        //   // Debug: Print the start and nrend values
        //   fprintf(th, "start = %d, nrend = %d\n", start, nrend);

        //   // Iterate over the range from start to nrend
        //   for (n = start; n < nrend; n++) {
        //     // Debug: print current iteration
        //     fprintf(th, "In loop: n = %d\n", n);

        //     // Increment the counter c for each iteration
        //     c += 1;
        //   }

        //   // After the loop, print the final value of c
        //   fprintf(th, "Final value of c: %d\n", (int)c);
        //   fclose(th);
        // }

        //  hhMD addtions:
        /*
        int cells[27] =
        //{94,187,210,315,354,477,599,739,741,988,1018,1156,1571,1825,2240,2259,2275,2414,2477,2530,2725,2780,2782,3231,3313,4134,4439};
        {1780,4275,3079,276,878,4244,1198,1887,1457,2821,1374,2190,1738,2347,2768,4055,48,146,3516,2610,3617,3008,2166,658,221,2840,3591};

        int outp=0;
        int outq=0;
        for (int q=0; q<27; q++)
          if (n==cells[q]-1) {
                outp=1;
                outq=q;
                break;
                }

        if (fh->step_MD%100!=0)
          outp=0;
        */

        /*if (outp)
                printf("%d:",n+1);
        */

        dvec r, vn, vn1;

        PBC(r, x[n], fh->box);
        if (r[0] <= 0) r[0] = 0.0000000000000000001;
        if (r[1] <= 0) r[1] = 0.0000000000000000001;
        if (r[2] <= 0) r[2] = 0.0000000000000000001;

        /*
                if (outp)
for (int d=0; d<3; d++)
  {
  r[d] = fh->grid.c[outq][d];//+0.00001;
//  if (r[d]-0.0<=0.5)  r[d] = 10.0;  else
//  if (r[d]-20.0<=0.5)  r[d] = 30.0;  else
//  if (r[d]-40.0<=0.5)  r[d] = 50.0;
  }

                if (outp)
                        printf("%.3f  %.3f  %.3f     ",r[0],r[1],r[2]);
*/

        /* this is for velocity update - do not use it at present
                        vn[0] = lg*v[n][0];
                        vn[1] = lg*v[n][1];
                        vn[2] = lg*v[n][2];

                        if (dnorm2(vn)>1e20) {
                                fprintf (stderr,"\nvel inf   : %d\n",n);
                                switch (fpclassify(dnorm2(vn))) {
                                        case FP_INFINITE:  fprintf
           (stderr,"infinite");  break; case FP_NAN:       fprintf
           (stderr,"NaN");       break; case FP_ZERO:      fprintf
           (stderr,"zero");      break; case FP_SUBNORMAL: fprintf
           (stderr,"subnormal"); break; case FP_NORMAL:    fprintf
           (stderr,"normal");    break;
                                }
                        }
        */
        /*		double ro_fh;
                trilinear_find_neighbours(x[n], n, xi, nbr, fh);
                trilinear_interpolation_scalar(ro_fh, xi, INTERPOLATE(ro_fh));
        //        if (n%1000==0)
        //            printf("%.3f %.3f %.3f  ", xi[0], xi[2], ro_fh);

                    double invden = 1.0/(1.0 + fh->S*fh->q);
                        fh->ro_hybr[n] = fh->S*fh->q * ro_fh * invden;
        */

        int niter = 0;
        int nk = 0;

      iterstart:
        niter++;
        // FILE *hybrid_terms;
        // if (n == 222){
        //    hybrid_terms = fopen("hybrid_terms.dat", "a");
        //   fprintf(hybrid_terms, "\n %d ", fh->step_MD);
        // }

        fh->summu2[n] = 0;
        fh->m_fhparticle[n] = 0;
        for (d = 0; d < DIM; d++) fh->p_fhparticle[n][d] = 0;

        dvec rsc;
        rsc[0] = r[0] * fh->fh_scale_lin;  ///  =  x_hybr * fh->fh_scale_lin /
                                           ///  fh->hybr_scale_lin;
        rsc[1] = r[1] * fh->fh_scale_lin;
        rsc[2] = r[2] * fh->fh_scale_lin;

        /*        for (int q = 0; q < fh->tssl->GetNpoints(); q++)
                {
                  if (fh->tssl->in_this_cell(q, rsc))
                  {
                    double mu_k = fh->tssl->mu();
                    fh->summu2[n] += mu_k * mu_k;
                    nk++;

                    int k = fh->cells_map[q];
        */
        for (int q = 0; q < 4; q++) {
          double mu_k = fh->cells_data[n].cells_mu[q];
          fh->summu2[n] += mu_k * mu_k;
          nk++;

          int k = fh->cells_data[n].cells[q];

          if (arr[k].natoms_star == 0) {
            printf("empty_cell %d ", k);
            continue;
          }
          fh->m_fhparticle[n] += mu_k * arr[k].avg_mass_fh_particle;
          for (d = 0; d < DIM; d++) {
            fh->p_fhparticle[n][d] +=
                mu_k * arr[k].avg_mass_fh_particle *
                u_interp(fh, k, d);  // arr[k].u_fh[d];
                                     // Maryna special case for v_hybr pure fh
            if (fh->step_MD < 1) {
              fh->v_hybr[n][d] += S * mu_k * u_interp(fh, k, d);
            }
          }
          //        }
        }

        /*		if (!isfinite(norm2(f[n]))) {
                                fprintf (stderr,"\nmd inf   : %d",n);
                                switch (fpclassify(norm2(f[n]))) {
                                        case FP_INFINITE:  fprintf
           (stderr,"infinite");  break; case FP_NAN:       fprintf
           (stderr,"NaN");       break; case FP_ZERO:      fprintf
           (stderr,"zero");      break; case FP_SUBNORMAL: fprintf
           (stderr,"subnormal"); break; case FP_NORMAL:    fprintf
           (stderr,"normal");    break;
                                }
                        }
                        if (!isfinite(dnorm2(v_term))) {
                                fprintf (stderr,"\nhd inf   : %d %g %g
           %g\n",n,v_term[0],v_term[1],v_term[2]); switch
           (fpclassify(dnorm2(v_term))) { case FP_INFINITE:  fprintf
           (stderr,"infinite");  break; case FP_NAN:       fprintf
           (stderr,"NaN");       break; case FP_ZERO:      fprintf
           (stderr,"zero");      break; case FP_SUBNORMAL: fprintf
           (stderr,"subnormal"); break; case FP_NORMAL:    fprintf
           (stderr,"normal");    break;
                                }
                        }
        */

        //		if (nk==0)
        //		    fprintf(stdout, "no cell found! %g %g %g\n", r[0],
        // r[1], r[2]);

        if (nk != 4) {
          //		    fprintf(stdout, "%d neighs - cell %d ", nk,ind);
          fprintf(stdout, "%d  %d %.3f %.3f %.3f   %d neighs \n", fh->step_MD,
                  n, r[0], r[1], r[2], nk);
          nk = 0;

          // if (niter < 10) {
          //   r[0] += 0.0000001;
          //   r[1] += 0.0000002;
          //   r[2] -= 0.0000003;
          // } else {
          //   r[0] += 0.0001;
          //   r[1] -= 0.0002;
          //   r[2] += 0.0003;
          // }
          if (niter < 10) {
            r[0] += (double)rand() / ((double)(RAND_MAX) + 1.0) * 0.0000001;
            r[1] += (double)rand() / ((double)(RAND_MAX) + 1.0) * 0.0000001;
            r[2] += (double)rand() / ((double)(RAND_MAX) + 1.0) * 0.0000001;
          } else {
            r[0] += (double)rand() / ((double)(RAND_MAX) + 1.0) * 0.0001;
            r[1] += (double)rand() / ((double)(RAND_MAX) + 1.0) * 0.0001;
            r[2] += (double)rand() / ((double)(RAND_MAX) + 1.0) * 0.0001;
          }

          for (d = 0; d < DIM; d++)
            if (r[d] >= fh->box[d])
              r[d] -= fh->box[d];
            else if (r[d] < 0)
              r[d] += fh->box[d];

          goto iterstart;
        }

        fh->m_cross_term_prev[n] = fh->m_cross_term[n];
        fh->m_cross_term[n] = 0;
        // Maryna
        // fh->m_cross_term[n] =  10* (1.0 - fh->summu2[n]) * DRNOR() *
        // fh->std_rho / fh->avg_Natoms_percell; fprintf(stderr,"%f
        // ",fh->m_cross_term[n]); Maryna
        if (fh->step_MD > 0) {
          if (fh->m_fhparticle[n] <= 0) {
            fprintf(stdout, "step = %d\n", fh->step_MD);
            fprintf(stdout,
                    "Error: Invalid mass for particle before cross term %d: "
                    "m_fhparticle[%d] = %f\n",
                    n, n, fh->m_fhparticle[n]);
            fprintf(stdout, "m_fhparticle  = %f/n, n = %d/n",
                    fh->m_fhparticle[n], n);
          }
        }
        fh->m_fhparticle[n] += fh->m_cross_term[n];
        fh->m_hybr[n] = ((1.0 - S) * (1.0 - lambda_s)) / invmass[n] +
                        (S + lambda_s - S * lambda_s) * fh->m_fhparticle[n];
        //  fh->m_hybr[n] = (1.0 - S) / invmass[n] + S * fh->m_fhparticle[n];
        //         if (n == 222){
        //   fprintf(hybrid_terms, "%f ", invmass[n]);
        //   fprintf(hybrid_terms, "%f ", fh->m_fhparticle[n]);

        // }

        for (d = 0; d < DIM; d++) {
          fh->p_cross_term_prev[n][d] = fh->p_cross_term[n][d];
          // Maryna
          fh->p_cross_term[n][d] = 0;
          // fh->p_cross_term[n][d] = 10 * (1.0 - fh->summu2[n]) * DRNOR() *
          // fh->std_u * fh->FH_dens / fh->avg_Natoms_percell;
          fh->p_fhparticle[n][d] += fh->p_cross_term[n][d];
          fh->p_hybr[n][d] =
              (1.0 - S) / invmass[n] * v[n][d] + S * fh->p_fhparticle[n][d];
        }

        // fh->summu2[n] = sqrt(fh->summu2[n]); //прибрати, вивести сумарне
        // значення

        /*
                if (outp)
                    printf("  ");
                        if (niter>3)
                            printf("niter = %d\n",niter);
        */

        /*                for (d = 0; d < DIM; d++)
                            vn1[d]  = lg*v[n][d] + dt/mhybr*( (1-S)*f[n][d] -
           S*v_term[d] );

                        dvec diff;
                        dvec_sub(vn1,vn,diff);
                        vn[0] = vn1[0];
                        vn[1] = vn1[1];
                        vn[2] = vn1[2];
                        if (dnorm2(diff)>1e-8)
                            goto iterstart;

                        if (n%100==0)
                            fprintf(stdout,"%d ",niter);
        */

        /*		if (norm2(f[n])>15000*15000)
                            printf("step %d md %g\n", fh->step_MD, norm(f[n]));

                        if (dnorm2(v_term)>15000*15000)
                            printf("step %d          fh %g\n", fh->step_MD,
           dnorm(v_term));
        */

        //	if (outp)
        //  	  printf("\n");

        // FILE *fw;
        // // FILE *hybrid_terms;
        // if (n == 222) fw = fopen("atom223.dat", "a");
        // // if (n == 222) hybrid_terms = fopen("hybrid_terms.dat", "a");
        // FILE *ac;
        // if (n == 222) ac = fopen("atom223acc.dat", "a");

        dvec a_hybr;
        // MAryna
        //  a_hybr[0] = 0;
        //  a_hybr[1] = 0;
        //  a_hybr[2] = 0;

        for (d = 0; d < DIM; d++) {
          //		     if (S<1)    //  we anyway do not use MD forces at
          // S=1
          /*		    if (f[n][d]>1e6)
                                  f[n][d] = 1e6;
                              else
                                  if (f[n][d]<-1e6)
                                      f[n][d] = -1e6;
                                  if (isnan(f[n][d]) || isinf(f[n][d]))
                                          f[n][d] = 0;
          */
          // updating pure MD velocity, by which MD temperature is calculated
          // let pure MD velocities be constant
          //		    v[n][d] = lg*v[n][d] + f[n][d]*dt*invmass[n];

          // acceleration of the hybrid particle (Eq. 13 in formulas_April.pdf)

          // s =1 MARYNA fix
          if (S < 0.99)
          {
            a_hybr[d] = ((1.0 - S) * (1.0 - lambda_s)) * f[n][d] / fh->m_hybr[n];
          }
          else
          {
            a_hybr[d] = 0.0;
            // fh->v_hybr[n][d] = 0;
          }
          // if (S < 0.99)  //  we anyway do not use MD forces at S=1
          //   a_hybr[d] =  (1.0 - S) * f[n][d] / fh->m_hybr[n];
          // else
          //   a_hybr[d] = 0.0;
          //		    v[n][d] = /*lg*/v[n][d] + a_hybr[d]*dt;
          //			a_hybr[d] = 0;

          // if (S >= 2.99) fh->v_hybr[n][d] = 0;
        }

        // if (n == 222) {
        //   fprintf(fw, "step = %d\n", fh->step_MD);
        //   fprintf(fw, "v_hybr = %f %f %f  m_hybr = %f\n", fh->v_hybr[n][0],
        //           fh->v_hybr[n][1], fh->v_hybr[n][2], fh->m_hybr[n]);
        //   fprintf(fw, "f_md = %f %f %f\n", f[n][0], f[n][1], f[n][2]);
        //   fprintf(fw, "x_md = %f %f %f\n", x[n][0], x[n][1], x[n][2]);
        //   // fprintf(fw, "lambda =  %d\n", tcstat[gt].lambda);
        // }
        // if (n == 222) {
        // // Writing only the data values (no text or labels)
        // fprintf(ac, "%d ", fh->step_MD);  // Only the step value
        // fprintf(ac, "%f %f %f ", a_hybr[0], a_hybr[1], a_hybr[2]);
        // fprintf(ac, "%f %f %f ", fh->v_hybr[n][0], fh->v_hybr[n][1],
        // fh->v_hybr[n][2]);  // Only the values for v_hybr and m_hybr

        // }

        // if (n%100==0)
        // printf("%f ",fh->m_hybr[n]);

        // fprintf(stderr,"> %f ",fh->m_hybr[n]);
        // fprintf(stderr,">> %f ",fh->m_cross_term[n]);

        /*		FILE *fw2;
                        if (n==222)
                        fw2 = fopen("atom559_tetr.dat", "a");
        */

        //        for (int q = 0; q < fh->tssl->GetNpoints(); q++)


        /// insert here v_hybr for s >0.99

        if (S >= 0.99){ 
          for (d = 0; d < DIM; d++) {
                fh->v_hybr[n][d] = 0.0;
              }

          for (int q = 0; q < 4; q++) {
            double mu_k = fh->cells_data[n].cells_mu[q];
            int k = fh->cells_data[n].cells[q];
            for (d = 0; d < DIM; d++) {
            fh->v_hybr[n][d] +=
                    mu_k * arr[k].avg_mass_fh_particle * u_interp(fh, k, d);
            }
          }
        }

        for (int q = 0; q < 4; q++) {
          //          if (fh->tssl->in_this_cell(q, rsc))
          {
            //            double mu_k = fh->tssl->mu();
            dvec c_k;
            //            fh->tssl->c(c_k);
            double mu_k = fh->cells_data[n].cells_mu[q];
            c_k[0] = fh->cells_data[n].cells_c[q][0];
            c_k[1] = fh->cells_data[n].cells_c[q][1];
            c_k[2] = fh->cells_data[n].cells_c[q][2];

            //            int k = fh->cells_map[q];
            int k = fh->cells_data[n].cells[q];

            // if (n == 222) {
            //   fprintf(fw, "%f %f = ", fh->arr[k].u_fh[0],
            //           fh->arr_next[k].u_fh[0]);
            //   // fprintf(hybrid_terms, "\n step = %d ", fh->step_MD);
            //   // fprintf(hybrid_terms, "u_fh = %f , u_fh_next = %f ",
            //   fh->arr[k].u_fh[0],
            //   //   fh->arr_next[k].u_fh[0]);
            //   //					printf("%f %f
            //   //%f\n",rsc[0],rsc[1],rsc[2]);
            //   // fh->tssl->PrintFoundSimplex();
            //   // printf("\n");
            // }
            if (arr[k].natoms_star == 0) continue;

            // arr[k].hybr_mass_star += 0.25*fh->m_hybr[n];
            // maryna
            arr[k].hybr_mass_star += mu_k * fh->m_hybr[n];

            if (mu_k > 0) {
              arr[k].hybr_mass_star_without_mu += fh->m_hybr[n];
            }

            if (mu_k >= 0.65 && mu_k < 0.7) {
              arr[k].hybr_mass_star_1 += mu_k * fh->m_hybr[n];
            } else if (mu_k >= 0.6 && mu_k < 0.65) {
              arr[k].hybr_mass_star_2 += mu_k * fh->m_hybr[n];
            } else if (mu_k >= 0.55 && mu_k < 0.6) {
              arr[k].hybr_mass_star_3 += mu_k * fh->m_hybr[n];
            } else if (mu_k >= 0.5 && mu_k < 0.55) {
              arr[k].hybr_mass_star_4 += mu_k * fh->m_hybr[n];
            }

            for (d = 0; d < DIM; d++) {
              if (mu_k > 0) {
                arr[k].hybr_vel_star_without_mu[d] += fh->v_hybr[n][d];
              }
              arr[k].hybr_vel_star[d] += mu_k * fh->v_hybr[n][d];

              if (S < 0.99){
              arr[k].hybr_momentum_star[d] +=
                  mu_k * fh->m_hybr[n] * fh->v_hybr[n][d];  
                }
                else {
                  arr[k].hybr_momentum_star[d] +=
                  mu_k * fh->v_hybr[n][d];  
                }// new stds
            }

            arr[k].tot_mu += mu_k;

            arr[k].mu2m2_sum +=
                mu_k * mu_k * fh->m_hybr[n] * fh->m_hybr[n];  // new stds
            arr[k].mu2m_sum += mu_k * mu_k * fh->m_hybr[n];   // new stds

            double k_c_k_dot_v = c_k[0] * fh->v_hybr[n][0] +
                                 c_k[1] * fh->v_hybr[n][1] +
                                 c_k[2] * fh->v_hybr[n][2];
            k_c_k_dot_v *= fh->fh_scale_lin / fh->hybr_scale_lin;

            // if (n == 222) {
            //   fprintf(fw, "c_k = %f %f %f  k*c_k*v = %f", c_k[0], c_k[1],
            //           c_k[2], k_c_k_dot_v);
            //   fprintf(fw, "  u_fh_interp = %f %f %f\n", u_interp(fh, k, 0),
            //           u_interp(fh, k, 1), u_interp(fh, k, 2));
            //   fprintf(fw, "hybrid terms 1x,2x,1y,2y,1z,2z = ");

            // }

            for (d = 0; d < DIM; d++) {
              /*              if (S >= 2.99)
                            {
                              fh->v_hybr[n][d] +=
                                  mu_k * arr[k].avg_mass_fh_particle *
                 u_interp(fh, k, d); continue;
                            }
              */
              double t = 0;
              // t = 1/ fh->m_hybr[n] * k_c_k_dot_v *
              //     arr[k].avg_mass_fh_particle *
              //     (S*fh->v_hybr[n][d] - (S + lambda_s - S*lambda_s) *
              //     u_interp(fh, k, d) /*arr[k].u_fh[d]*/); - second option

              t = (S + lambda_s - S * lambda_s) / fh->m_hybr[n] * k_c_k_dot_v *
                  arr[k].avg_mass_fh_particle *
                  (fh->v_hybr[n][d] - u_interp(fh, k, d) /*arr[k].u_fh[d]*/);

              if (fabs(t) > 1e4) printf("!%d %f %f ", n, t, fh->v_hybr[n][d]);
              if (d == 1) {
                if (n == 222) {
                  // fprintf(hybrid_terms, "\n step = %d ", fh->step_MD);
                  // fprintf(hybrid_terms, "t_1 = %f ",t);
                  // fprintf(hybrid_terms, "%f ",t);
                  //  fprintf(hybrid_terms, "%f ",fh->m_hybr[n]);
                  //  fprintf(hybrid_terms, "%f ",k_c_k_dot_v);
                  //  fprintf(hybrid_terms, "%f ",arr[k].avg_mass_fh_particle);
                  //  fprintf(hybrid_terms, "%f ",fh->v_hybr[n][d]);
                  //  fprintf(hybrid_terms, "%f ",u_interp(fh, k, 1));
                }
              }
              a_hybr[d] -= t;

              // if (n == 222) {
              //   // fprintf(fw, "%f ", t);
              //   fprintf(th2, "step = %d\n", fh->step_MD);
              //   fprintf(th2, " v_hybr_ near interp[%d] = %f\n", d,
              //   fh->v_hybr[n][d]); fprintf(th2, "CV = %d\n", k); fprintf(th2,
              //   " u_interp[%d] = %f\n", d, u_interp(fh, k, d));
              // }

              // t = 1/ fh->m_hybr[n] * mu_k *
              //     ((arr[k].avg_mass_fh_particle -
              //       arr[k].avg_mass_fh_particle_prev) /
              //          fh->dt_FH *
              //          (S* fh->v_hybr[n][d] -
              //           (S + lambda_s - S*lambda_s) * u_interp(fh, k, d)
              //           /*arr[k].u_fh[d]*/) -
              //      (S + lambda_s- S*lambda_s) *(fh->arr[k].u_fh[d] -
              //      fh->arr_prev[k].u_fh[d]) / fh->dt_FH *
              //          arr[k].avg_mass_fh_particle); - second option

              t = (S + lambda_s - S * lambda_s) / fh->m_hybr[n] * mu_k *
                  ((arr[k].avg_mass_fh_particle -
                    arr[k].avg_mass_fh_particle_prev) /
                       fh->dt_FH *
                       (fh->v_hybr[n][d] -
                        u_interp(fh, k, d) /*arr[k].u_fh[d]*/) -
                   (fh->arr[k].u_fh[d] - fh->arr_prev[k].u_fh[d]) / fh->dt_FH *
                       arr[k].avg_mass_fh_particle);
              if (fabs(t) > 1e4) printf("!!%d %f ", n, t);
              if (d == 1) {
                if (n == 222) {
                  // fprintf(hybrid_terms, "\n step = %d ", fh->step_MD);
                  // fprintf(hybrid_terms, "t_2 = %f ",t);
                  // fprintf(hybrid_terms, "%f ",t);
                  // fprintf(hybrid_terms, "mu = %f ",mu_k);
                  // fprintf(hybrid_terms, "avg_mass_fh_particle = %f
                  // ",arr[k].avg_mass_fh_particle); fprintf(hybrid_terms,
                  // "avg_mass_fh_particle_prev = %f
                  // ",arr[k].avg_mass_fh_particle); fprintf(hybrid_terms,
                  // "u_interp(fh, k, d) = %f ",u_interp(fh, k, d));
                  // fprintf(hybrid_terms, "u_fh = %f ",fh->arr[k].u_fh[d]);
                  // fprintf(hybrid_terms, "u_fh_prev = %f
                  // ",fh->arr_prev[k].u_fh[d]);

                  // fprintf(hybrid_terms, "%f ",mu_k);
                  // fprintf(hybrid_terms, "%f ",arr[k].avg_mass_fh_particle);
                  // fprintf(hybrid_terms, "%f
                  // ",arr[k].avg_mass_fh_particle_prev); fprintf(hybrid_terms,
                  // "%f ",u_interp(fh, k, d)); fprintf(hybrid_terms, "%f
                  // ",fh->arr[k].u_fh[d]); fprintf(hybrid_terms, "%f
                  // ",(fh->arr[k].u_fh[d] - fh->arr_prev[k].u_fh[d]) /
                  // fh->dt_FH); fprintf(hybrid_terms, "%f
                  // ",fh->arr[k].natoms_star);
                }
              }
              a_hybr[d] -= t;

              // if (n == 222) fprintf(fw, "%f ", t);
            }

            // if (n == 222) fprintf(fw, "\n");
          }
        }

        /*		if (n==222)
                                fclose(fw2);
        */
        // if (n == 222) fprintf(fw, "hybrid terms 3x,4x,3y,4y,3z,4z = ");

        for (d = 0; d < DIM; d++) {
          /*          if (S >= 2.99)
                    {
                      fh->v_hybr[n][d] += fh->p_cross_term[n][d];
                      // if (n == 222) {
                      //   fprintf(fw, "v_hybr for s>=0.99 =  %f, d = %d",
             fh->v_hybr[n][d],
                      //           d);
                      // }
                      continue;
                    }
          */
          double t = 0;
          t = (S + lambda_s - S * lambda_s) / fh->m_hybr[n] *
              ((fh->m_cross_term[n] - fh->m_cross_term_prev[n]) / fh->dt_FH *
               fh->v_hybr[n][d]);
          if (fabs(t) > 1e4) printf("!!!%d %f ", n, fh->m_cross_term[n]);
          // if (d == 1){
          //   if (n == 222) {
          //     // fprintf(hybrid_terms, "\n step = %d ", fh->step_MD);
          //     // fprintf(hybrid_terms, "t_3 = %f ",t);
          //     // fprintf(hybrid_terms, "%f ",t);
          //   }
          // }
          a_hybr[d] -= t;

          // if (n == 222) fprintf(fw, "%f ", t);

          t = -(S + lambda_s - S * lambda_s) / fh->m_hybr[n] *
              ((fh->p_cross_term[n][d] - fh->p_cross_term_prev[n][d]) /
               fh->dt_FH);
          if (fabs(t) > 1e4) printf("!!!!%d %f ", n, fh->p_cross_term[n][d]);
          if (d == 1) {
            if (n == 222) {
              // fprintf(hybrid_terms, "\n step = %d ", fh->step_MD);
              // fprintf(hybrid_terms, "t_4 = %f ",t);
              // fprintf(hybrid_terms, "%f ",t);
              // fprintf(hybrid_terms, "%f ",fh->m_hybr[n]);
              // fprintf(hybrid_terms, "%f ",fh->p_cross_term[n][d]);
              // fprintf(hybrid_terms, "%f ",fh->p_cross_term_prev[n][d]);
              // fprintf(hybrid_terms, "%f ",fh->dt_FH);
            }
          }
          a_hybr[d] -= t;
          if (fh->step_MD == 0) {
            fh->a_hybr_prev[n][d] =
                a_hybr[d];  // Initialize to the current acceleration for the
                            // first step - fro verlet
          }

          // if (n == 222) fprintf(fw, "%f ", t);
        }

        // if (n == 222) {
        //   //fprintf(fw, "\n");
        //   fprintf(fw, "\n");
        //   fprintf(fw, "a_hybr = %f %f %f\n", a_hybr[0], a_hybr[1],
        //   a_hybr[2]);
        // }

        // movement of the hybrid particle
        //         if (n == 222) {
        // // Writing only the data values (no text or labels)
        //   // Only the step value
        // fprintf(ac, "%f %f %f", a_hybr[0], a_hybr[1], a_hybr[2]);  // Only
        // the values for v_hybr and m_hybr

        // }
        //                         if (n == 222) {
        //   // Writing only the data values (no text or labels)
        //   fprintf(ac, " before ");  // Only the step value
        //   fprintf(ac, "%f %f %f ", a_hybr[0], a_hybr[1], a_hybr[2]);
        //   fprintf(ac, "%f %f %f ", fh->v_hybr[n][0], fh->v_hybr[n][1],
        //   fh->v_hybr[n][2]);  // Only the values for v_hybr and m_hybr
        //  fprintf(ac, "%f %f %f ", fh->a_hybr_prev[n][0],
        //  fh->a_hybr_prev[n][1], fh->a_hybr_prev[n][2]);
        //   //  fprintf(ac, " /n");
        //   // fprintf(ac, " %d/n", fh->step_MD);
        //   fprintf(ac, " end \n");
        //   }
        for (d = 0; d < DIM; d++)
        {
          if (S >= 0.99)
          {
            double T0 = 0.730094;          // target reduced temperature
            double tau_T = fh->tau_m * dt; // coupling time - was 5 before
            double gammaT = 1.0 / tau_T;   // friction

            double theta = exp(-gammaT * dt);
            double m_md_i = 1.0 / invmass[n]; // MD mass from invmass
            double sigma = sqrt(m_md_i * T0 * (1.0 - theta * theta));
            double beta = fh->fluct_m;

            fh->v_hybr[n][d] += beta * sigma * DRNOR();

            fh->v_hybr[n][d] /= fh->m_hybr[n];
            fh->x_hybr[n][d] = fh->x_hybr[n][d] + fh->v_hybr[n][d] * dt;
          }
          else
          {
            // classical integrator
            //  fh->v_hybr[n][d] += a_hybr[d]*dt;
            // fh->x_hybr[n][d] = fh->x_hybr[n][d] + fh->v_hybr[n][d] * dt;

            // // MARYNA VERLET

            //           if (n == 222) {
            // // Writing only the data values (no text or labels)
            // fprintf(ac, " before after ");  // Only the step value
            // // fprintf(ac, "%f %f %f ", a_hybr[0], a_hybr[1], a_hybr[2]);
            // fprintf(ac, "%f %f %f ", fh->v_hybr[n][0], fh->v_hybr[n][1],
            // fh->v_hybr[n][2]);  // Only the values for v_hybr and m_hybr
            // // fprintf(ac, "%f", dt);
            // // fprintf(ac, "%f %f %f ", fh->a_hybr_prev[n][0],
            // fh->a_hybr_prev[n][1], fh->a_hybr_prev[n][2]);
            // //  fprintf(ac, " /n");

            // }
            // fh->x_hybr[n][d] = fh->x_hybr[n][d] + fh->v_hybr[n][d] * dt + dt
            // * dt * 0.5 * fh->a_hybr_prev[n][d];

            // //fh->v_hybr[n][d] += 0.5 * (a_hybr[d] + fh->a_hybr_prev[n][d]) *
            // dt + S * 0.1 * (2.0 * ((double)rand() / (RAND_MAX + 1.0)) - 1.0)
            // * (0.008) ; fh->v_hybr[n][d] += 0.5 * (a_hybr[d] +
            // fh->a_hybr_prev[n][d]) * dt  ;

            // Old Thermostt Berendsena
            //   if (S < 0.99){
            // fh->p_hybr_m[n][d] = lambda * (1.0 - S) * fh->p_mol_m[n][d] + S *
            // fh->p_fhparticle[n][d];
            // // fh->v_hybr[n][d] = fh->p_hybr_m[n][d] / fh->m_hybr[n]  + S *
            // fh->fluct_m* (2.0 * ((double)rand() / (RAND_MAX + 1.0)) - 1.0) *
            // (0.008)*(1.0 - fh->summu2[n]); fh->v_hybr[n][d] =
            //       fh->p_hybr_m[n][d] / fh->m_hybr[n]
            //       + sqrt(S) * fh->fluct_m * ((1.0 - exp(-2.0 * dt /
            //       fh->tau_m)) * T_0 / fh->m_hybr[n]) * DRNOR();

            //   }

            fh->a_hybr_prev[n][d] = a_hybr[d];

            // === B: kick (use acceleration at r(t))
            fh->v_hybr[n][d] += dt * a_hybr[d];

            // === A: first half drift
            fh->x_hybr[n][d] += 0.5 * dt * fh->v_hybr[n][d];

            if (fabs(fh->v_hybr[n][d]) > 1e3)
              printf("!!!!!%d %f ", n, fh->v_hybr[n][d]);
            //  Maryna therm

              fh->p_hybr_m[n][d] = fh->m_hybr[n] * fh->v_hybr[n][d];
              fh->p_mol_m[n][d] =
                  (fh->p_hybr_m[n][d] - (S + lambda_s - S * lambda_s) * fh->p_fhparticle[n][d]) / ((1.0 - S) * (1.0 - lambda_s));
              double md_part = (fh->p_hybr_m[n][d] - (S + lambda_s - S * lambda_s) * fh->p_fhparticle[n][d]);

              // Thermostat parameters (reduced units, kB = 1)
              double T0 = 0.730094;          // target reduced temperature
              double tau_T = fh->tau_m * dt; // coupling time - was 5 before
              double gammaT = 1.0 / tau_T;   // friction

              double theta = exp(-gammaT * dt);
              double m_md_i = 1.0 / invmass[n]; // MD mass from invmass
              double sigma = sqrt(m_md_i * T0 * (1.0 - theta * theta));
              double beta = fh->fluct_m; // assuming 0 <= S < 1


              // // OU update of MD momentum
              double pmd = fh->p_mol_m[n][d];
              pmd = theta * pmd + sigma * beta * DRNOR(); // DRNOR() = N(0,1)
              fh->p_mol_m[n][d] = pmd;
              // md_part *= theta;

              // Rebuild hybrid momentum and velocity
              fh->p_hybr_m[n][d] = (1.0 - S) * (1.0 - lambda_s) *
              fh->p_mol_m[n][d] + (S + lambda_s - S * lambda_s) *
              fh->p_fhparticle[n][d];
              // fh->p_hybr_m[n][d] =
              //     md_part                      // θ · (1−S) p_MD
              //     + (S + lambda_s - S * lambda_s) * fh->p_fhparticle[n][d] // FH contribution
              //     + beta * sigma * DRNOR();     // hybrid noise
              fh->v_hybr[n][d] = fh->p_hybr_m[n][d] / fh->m_hybr[n];
              fh->x_hybr[n][d] += 0.5 * dt * fh->v_hybr[n][d];
            
          }

          // === A: second half drift
          // NEW INTEGRATOR
          // fh->x_hybr[n][d] += 0.5 * dt * fh->v_hybr[n][d];
          if (fh->x_hybr[n][d] > fh->box[d] * fh->hybr_scale_lin)
            fh->x_hybr[n][d] -= fh->box[d] * fh->hybr_scale_lin;
          else if (fh->x_hybr[n][d] < 0)
            fh->x_hybr[n][d] += fh->box[d] * fh->hybr_scale_lin;

          xprime[n][d] = fh->x_hybr[n][d] / fh->hybr_scale_lin;

          //            xprime[n][d] = x[n][d] + fh->v_hybr[n][d]*dt;
          // run12:            xprime[n][d] = x[n][d] + v[n][d]*dt;
        }

        // fclose(th);
        // fclose(th2);
        // if (n == 222) fclose(fw);
        // if (n == 222) fclose(ac);
        // if (n == 222) fclose(hybrid_terms);
      } else {
        for (d = 0; d < DIM; d++) {
          v[n][d] = 0.0;
          xprime[n][d] = x[n][d];
        }
      }
      // fclose(th);
      // fclose(th2);
    }
    // fclose(th2);
    // fclose(th);
  }

  if (start == 0)  //  the first thread
    fclose(th);

  //	printf ("st%d,a%.4f,k%.4f
  //",fh->step_MD,fh->arr[0].hybr_mass_star,fh->arr[10].hybr_mass_star);

  if (0)
    if (fh->step_MD % 100 == 0) {
      FILE* fw1;
      char fname[32];

      sprintf(fname, "ro_mixed_0.5-0.55_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.5) && (fh->summu2[n] < 0.55))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.55-0.6_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.55) && (fh->summu2[n] < 0.6))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.6-0.65_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.6) && (fh->summu2[n] < 0.65))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.65-0.7_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.65) && (fh->summu2[n] < 0.7))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.7-0.75_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.7) && (fh->summu2[n] < 0.75))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.75-0.8_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.75) && (fh->summu2[n] < 0.8))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.8-0.85_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.8) && (fh->summu2[n] < 0.85))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.85-0.9_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.85) && (fh->summu2[n] < 0.9))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.9-0.95_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.9) && (fh->summu2[n] < 0.95))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);

      sprintf(fname, "ro_mixed_0.95-1_atoms.xvg");
      if ((fw1 = fopen(fname, "a")) == NULL) {
        printf("\n ERROR creating %s for output!\n", fname);
      }

      for (n = start; n < nrend; n++)
        if ((fh->summu2[n] >= 0.95) && (fh->summu2[n] < 1.0))
          fprintf(fw1, "%.4f   %.4f\n", (fh->summu2[n]), fh->ro_hybr[n]);

      fclose(fw1);
    }

  if (fh->write_debug_frames) {
    FILE* f2 = fopen("lastsuccess_fromnewupdate.xyz", "w");
    char an[3] = {'O', 'H', 'H'};
    fprintf(f2, "%d\n", nrend - start + 1);
    fprintf(f2,
            "the frame which has been processed before crash, written in "
            "new_update.cpp\n");
    for (int n = start; n < nrend; n++)
      fprintf(f2, "%c   %g  %g  %g   %g  %g  %g   %g  %g  %g\n", an[n % 3],
              xprime[n][0] * 10.0, xprime[n][1] * 10.0, xprime[n][2] * 10.0,
              v[n][0] * 10.0, v[n][1] * 10.0, v[n][2] * 10.0, f[n][0] * 10.0,
              f[n][1] * 10.0, f[n][2] * 10.0);
    fclose(f2);
  }

  if PAR (cr) gmx_barrier(cr);
}
