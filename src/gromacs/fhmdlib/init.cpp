#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "data_structures.h"
#include "parser.h"
#include "fh_functions.h"
#include "estimate.h"
#include "sfunction.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/hhmd_G/tesselation.h"
#include "macro.h"

#include "macro.h"

void init_ind(rvec x[], int N_atoms, FHMD *fh)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    int        ind;


    /* Collect statistics */
    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        for(int d = 0; d < DIM; d++)
        {
            fh->indv[n][d] = (int)(xn[d]/fh->box[d]*(double)(fh->N_md[d])) + fh->N_shift[d];
        }

        ind = I(fh->indv[n], fh->N);

        fh->ind[n] = ind;
    }
}

int fhmd_init(matrix box, int N_atoms, float mass[], rvec x[], rvec v[], double dt_md, gmx_mtop_t *mtop, t_commrec *cr, FHMD *fh)
// int fhmd_init(matrix box, int N_atoms, float mass[], rvec x[], rvec v[], double dt_md, gmx_mtop_t *mtop, t_commrec *cr, FHMD *fh)

{
    int N_atoms_th = N_atoms;

    FILE *fw;

    if(MASTER(cr))
    {
        char const *fname_in  = "coupling.prm";
        char const *fname_out = "coupling_out.prm";

        /* Initial output */

        printf(MAKE_GREEN "\n  Aston University, Department of Mathematics, Dmitry Nerukh Research Group\n");
        printf("  Queen Mary University of London, School of Engineering and Material Science\n\n");
        printf(MAKE_PURPLE "     Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", FHMD_VERSION);

        /* Default values of FHMD parameters */

        fh->scheme      = 2;
        fh->S           = 0;
        fh->R1          = 0.5;
        fh->R2          = 1;
        fh->Smin        = 0;
        fh->Smax        = 0.5;
        fh->alpha       = 50;
        fh->beta        = 20;
        fh->gamma_x     = 0;
        fh->gamma_u     = 1;
        fh->eps_rho     = 0.05;
        fh->eps_mom     = 0.05;
        fh->S_berendsen = 1;
        fh->N[0]        = 9;
        fh->N[1]        = 9;
        fh->N[2]        = 9;
        fh->N_md[0]     = 5; //fh->N[0];
        fh->N_md[1]     = 5; //fh->N[1];
        fh->N_md[2]     = 5; //fh->N[2];
        fh->FH_EOS      = 1;
        fh->FH_step     = 10;
        fh->FH_equil    = 10000;
        fh->FH_dens     = 600;
        fh->FH_temp     = 298.15;
        fh->FH_blend    = 0.005;
        fh->Noutput     = 100;
        fh->write_debug_frames  =  0;   // VF:
        fh->nsteps_temper = 0;     // VF:
        fh->q           = 0;       // hhmd_G addition
        fh->tau_m       = 1;         //Maryna
        fh->fluct_m     = 1;       //Maryna   

        /* Read FHMD parameters */

        printf(MAKE_GREEN "FHMD: Reading parameters from %s..." RESET_COLOR " ", fname_in);

        int ok = parse_prm(fname_in, fh);

        if(ok == 1) {
            printf(MAKE_GREEN "...OK\n" RESET_COLOR "\n");
            fw = fopen(fname_out, "w");
        } else if(ok == -1) {
            printf(MAKE_RED "\nFHMD: File %s not found. Generating default parameter file...\n" RESET_COLOR "\n", fname_in);
            fw = fopen(fname_in, "w");
        } else {
            printf(MAKE_RED "\nFHMD: ERROR in %s file\n" RESET_COLOR "\n", fname_in);
            exit(2);
        }

        /* Print FHMD parameters to the screen and output file */

        fprintf(fw, "; Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", FHMD_VERSION);

        fprintf(fw, "scheme = %d              ; 0 - Pure MD, 1 - One-way coupling, 2 - Two-way coupling\n\n", fh->scheme);

        fprintf(fw, "S = %g                   ; Parameter S (-1 - fixed sphere, -2 - moving sphere)\n\n", fh->S);
        fprintf(fw, "R1   = %g              ; MD sphere radius for variable S, [0..1]\n", fh->R1);
        fprintf(fw, "R2   = %g              ; FH sphere radius for variable S, [0..1]\n", fh->R2);
        fprintf(fw, "Smin = %g                ; Minimum S for variable S\n", fh->Smin);
        fprintf(fw, "Smax = %g             ; Maximum S for variable S\n\n", fh->Smax);
        fprintf(fw, "tau = %g                   ; tau for berendsen thermostat implementation\n", fh->tau_m);
        fprintf(fw, "magnitude = %g                   ; magnitude of fluctuations for hybrid equations of motion\n", fh->fluct_m);



        switch(fh->scheme)
        {
        case Pure_MD:
            printf(MAKE_RED "FHMD: Starting Pure MD simulation (without MD/FH coupling)\n" RESET_COLOR "\n");
            break;
        case One_Way:
            printf(MAKE_YELLOW "FHMD: One-way MD/FH coupling\n" RESET_COLOR "\n");
            break;
        case Two_Way:
            break;
        }

        if(fh->S >= 0.0)
        {
            printf(MAKE_PURPLE "FHMD: S = %g\n", fh->S);
            fh->S_function = constant_S;
        }
        else
        {
            printf(MAKE_PURPLE "FHMD: S = S(x,y,z) = [%g, %g] with R1 = %g, R2 = %g\n", fh->Smin, fh->Smax, fh->R1, fh->R2);
            fh->S_function = fixed_sphere;

            fh->R1 *= box[0][0]*0.5;
            fh->R2 *= box[0][0]*0.5;
            fh->R12 = fh->R1*fh->R1;
            fh->R22 = fh->R2*fh->R2;
            fh->RS  = (fh->Smax - fh->Smin)/(fh->R2 - fh->R1);
            printf(MAKE_GREEN "FHMD: Absolute values of R [nm]: R1 = %f, R2 = %f\n", fh->R1, fh->R2);

            if(fh->S < -1.5)
            {
                printf(MAKE_PURPLE "FHMD: The MD/FH sphere will follow the protein\n");
                fh->S_function = moving_sphere;
            }
        }

        printf(MAKE_GREEN "FHMD: alpha = %g [nm^2/ps], beta = %g [ps^-1]\n", fh->alpha, fh->beta);
        fprintf(fw, "alpha   = %g           ; Alpha parameter for dx/dt and du/dt equations, nm^2/ps\n", fh->alpha);
        fprintf(fw, "beta    = %g           ; Beta parameter for du/dt equation, ps^-1\n\n", fh->beta);

        if(fh->scheme == Two_Way)
        {
            printf("FHMD: MD dissipator: gamma_x = %g [ps^-1], gamma_u = %g [ps^-1]\n", fh->gamma_x, fh->gamma_u);
            printf("FHMD: FH dissipator: eps_rho = %g [--], eps_mom = %g [--]\n", fh->eps_rho, fh->eps_mom);
            fprintf(fw, "gamma_x = %g             ; Gamma_x parameter (MD density fluctuations dissipator), ps^-1\n", fh->gamma_x);
            fprintf(fw, "gamma_u = %g             ; Gamma_u parameter (MD velocity fluctuations dissipator), ps^-1\n", fh->gamma_u);
            fprintf(fw, "eps_rho = %g             ; Eps_rho parameter (FH density fluctuations dissipator)\n", fh->eps_rho);
            fprintf(fw, "eps_mom = %g             ; Eps_mom parameter (FH momentum fluctuations dissipator)\n\n", fh->eps_mom);
        }

        if(fh->S_berendsen >= 0)
            printf("FHMD: Berendsen thermostat works for S <= %g\n", fh->S_berendsen);
        else
            printf("FHMD: Berendsen thermostat with (1 - S^%g) multiplier\n", -fh->S_berendsen);
        fprintf(fw, "S_berendsen = %g         ; If S_berendsen >= 0, Berendsen thermostat will work for S <= S_berendsen,\n", fh->S_berendsen);
        fprintf(fw, "                        ; otherwise factor (1-S^(-S_berendsen)) will be applied (local thermostat)\n\n");

        for(int d = 0; d < DIM; d++)
        {
            fh->box[d]         = box[d][d];
            fh->box05[d]       = 0.5*fh->box[d];
            fh->protein_com[d] = fh->box05[d];
            fh->N_shift[d]     = (fh->N[d] - fh->N_md[d])/2;
        }

        fh->box_volume = fh->box[0]*fh->box[1]*fh->box[2];

        printf("FHMD: MD box size:  %g x %g x %g [nm]\n", fh->box[0], fh->box[1], fh->box[2]);
        printf("FHMD: FH grid size: %d x %d x %d\n", fh->N[0], fh->N[1], fh->N[2]);
        fprintf(fw, "Nx = %d                  ; Number of FH cells along X axis\n", fh->N[0]);
        fprintf(fw, "Ny = %d                  ; Number of FH cells along Y axis\n", fh->N[1]);
        fprintf(fw, "Nz = %d                  ; Number of FH cells along Z axis\n\n", fh->N[2]);

        if(fh->scheme == Two_Way)
        {
            printf("FHMD: Small-scale MD/FH grid size: %d x %d x %d\n", fh->N_md[0], fh->N_md[1], fh->N_md[2]);
            fprintf(fw, "NxMD = %d                ; Number of small-scale MD-FH cells along X axis\n", fh->N_md[0]);
            fprintf(fw, "NyMD = %d                ; Number of small-scale MD-FH cells along Y axis\n", fh->N_md[1]);
            fprintf(fw, "NzMD = %d                ; Number of small-scale MD-FH cells along Z axis\n\n", fh->N_md[2]);
        }
        else
        {
            for(int d = 0; d < DIM; d++)
            {
                fh->N_shift[d] = 0;
                fh->N_md[d]    = fh->N[d];
            }
        }

        fh->Ntot    = fh->N[0]*fh->N[1]*fh->N[2];
        fh->Ntot_md = fh->N_md[0]*fh->N_md[1]*fh->N_md[2];

        switch(fh->FH_EOS)
        {
        case 0:
            fh->eos = eos_argon;
            printf("FHMD: Equation of state: Liquid Argon (300K)\n");
            break;
        case 1:
            fh->eos = eos_spce;
            printf("FHMD: Equation of state: Rigid SPC/E water\n");
            break;
        default:
            printf(MAKE_RED "\nFHMD: Unknown equation of state (%d) in %s\n" RESET_COLOR "\n", fh->FH_EOS, fname_in);
            exit(18);
        }

        fprintf(fw, "FH_EOS   = %d            ; EOS: 0 - Liquid Argon, 1 - SPC/E water\n", fh->FH_EOS);

        fh->dt_FH = (double)(fh->FH_step)*dt_md;

        printf("FHMD: FH time step dt_FH = %d * dt_MD = %g [ps]\n", fh->FH_step, fh->dt_FH);
        fprintf(fw, "FH_step  = %d           ; FH time step dt_FH = FH_step * dt_MD\n", fh->FH_step);

        if(fh->scheme == Two_Way) fh->FH_equil = 0;
        printf("FHMD: FH equilibration steps: %d\n", fh->FH_equil);
        fprintf(fw, "FH_equil = %d        ; Number of time steps for the FH model equilibration (for 1-way coupling)\n", fh->FH_equil);

        printf("FHMD: FH Density = %g [amu/nm^3], FH Temperature = %g [K]\n", fh->FH_dens, fh->FH_temp);
        fprintf(fw, "FH_dens  = %g          ; FH mean density\n", fh->FH_dens);
        fprintf(fw, "FH_temp  = %g       ; FH mean temperature\n", fh->FH_temp);
        fprintf(fw, "FH_blend = %g        ; FH Blending: -1 - dynamic, or define static blending parameter (0..1)\n\n", fh->FH_blend);

        printf("FHMD: MD/FH arrays will be written every %d MD time steps\n", fh->Noutput);
        fprintf(fw, "Noutput  = %d           ; Write arrays to files every Noutput MD time steps (0 - do not write)\n", fh->Noutput);

        printf(RESET_COLOR "\n");

        fflush(stdout);
    } // if(MASTER(cr))

    fhmd_reset_statistics(fh);

    fh->total_density = 0;
    for(int i = 0; i < N_atoms_th; i++)
        fh->total_density += mass[i];

    /* Broadcast parameters to all threads */

    if(PAR(cr))
    {
        gmx_sumd(1, &fh->total_density, cr);
        gmx_sumi(1, &N_atoms, cr);
        gmx_bcast(sizeof(FHMD), fh, cr);
    }

    fh->total_density /= fh->box_volume;

    if(fh->eos == eos_spce)
    {
        fhmd_find_protein(mtop, N_atoms_th, mass, cr, fh);
        if(fh->protein_N)
            fhmd_find_protein_com(mtop, N_atoms_th, x, mass, cr, fh);
    }
    else
    {
        fh->protein_N    = 0;
        fh->protein_mass = 0;
    }

    if(MASTER(cr))
    {
        printf(MAKE_GREEN "FHMD: Total number of atoms in the box: %d\n", N_atoms);
        printf("FHMD: Total density of the box: %g [amu/nm^3]\n", fh->total_density);

        if(fh->protein_N > 0)
            printf(MAKE_PURPLE "FHMD: Found protein: %d atoms, mass = %g [amu], COM = (%g, %g, %g) [nm]\n",
                    fh->protein_N, fh->protein_mass, fh->protein_com[0], fh->protein_com[1], fh->protein_com[2]);

        printf(RESET_COLOR "\n");
        fflush(stdout);

        fprintf(fw, "\n; You may consider to use FH_dens = %g since this is the total MD density of the box.\n", fh->total_density);
        fprintf(fw, "; NB: Please use spaces before and after '=' in this file, e.g. 'S = 0' (not 'S=0').\n");
        fclose(fw);
    }

    if(fh->scheme == Pure_MD) return 0;     // Start Pure MD simulation

    /* Allocate memory */

    fh->arr  = (FH_arrays*)calloc(fh->Ntot, sizeof(FH_arrays));
    fh->ind  = (int*)calloc(N_atoms, sizeof(int));
    fh->indv = (ivec*)calloc(N_atoms, sizeof(ivec));

    fh->mpi_linear = (double*)malloc(5*fh->Ntot*sizeof(double));   // 5 components: natoms, ro_md, uro_md[3]

    if(fh->arr == NULL || fh->ind == NULL || fh->indv == NULL || fh->mpi_linear == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (array allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    fh->grid.c    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.n    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.h    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.vol  = (double*)malloc(fh->Ntot*sizeof(double));
    fh->grid.ivol = (double*)malloc(fh->Ntot*sizeof(double));
    fh->grid.md   = (FHMD_CELL*)malloc(fh->Ntot*sizeof(FHMD_CELL));

    if(fh->grid.c == NULL || fh->grid.n == NULL || fh->grid.h == NULL || fh->grid.vol == NULL || fh->grid.ivol == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (FH grid allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    fh->stat.avg_rho_md_cell = (double*)calloc(fh->Ntot, sizeof(double));
    fh->stat.avg_rho_fh_cell = (double*)calloc(fh->Ntot, sizeof(double));

    if(fh->stat.avg_rho_md_cell == NULL || fh->stat.avg_rho_fh_cell == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (Statistics allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    /* Create FH grid and initialize FH solver */

    define_FH_grid(cr, fh);

/*
 * hhmd_G initialising
 */ 
/* my modifications for q=1*/
    // fh->fh_scale_lin = fh->alpha;
    // fh->q = fh->fh_scale_lin * fh->fh_scale_lin * fh->fh_scale_lin; //pow(fh->hhmd_G_coeff3, 0.333333);

	// // fh->hybr_scale_vol = 1.0 + fh->S*fh->q;
	// // fh->hybr_scale_lin = pow(fh->hybr_scale_vol, 0.3333333333);


    fh->fh_scale_lin = 1.0;
    fh->q = 1; //pow(fh->hhmd_G_coeff3, 0.333333);

	fh->hybr_scale_vol = 1.0;
	fh->hybr_scale_lin = 1.0;


	fh->nat = N_atoms;

    fh->arr_prev  = (FH_arrays*)calloc(fh->Ntot, sizeof(FH_arrays));
    fh->arr_next  = (FH_arrays*)calloc(fh->Ntot, sizeof(FH_arrays));
    fh->arr_temp  = (FH_arrays*)calloc(fh->Ntot, sizeof(FH_arrays));

    fh->cells_data = (tssl_data*)calloc(N_atoms, sizeof(tssl_data));

	fh->avg_Natoms_percell = ((double) N_atoms)/fh->Ntot;
	fh->x_hybr = (dvec*)calloc(N_atoms, sizeof(dvec));
	fh->ro_hybr = (double*)calloc(N_atoms, sizeof(double));
    fh->v_md = (dvec*)calloc(N_atoms, sizeof(dvec));

	fh->v_hybr = (dvec*)calloc(N_atoms, sizeof(dvec));
    fh->a_hybr_prev = (dvec*)calloc(N_atoms, sizeof(dvec));
    fh->p_hybr_m = (dvec*)calloc(N_atoms, sizeof(dvec));
    fh->p_mol_m = (dvec*)calloc(N_atoms, sizeof(dvec));
	fh->m_fhparticle = (double*)calloc(N_atoms, sizeof(double));
	fh->p_fhparticle = (dvec*)calloc(N_atoms, sizeof(dvec));
	fh->m_hybr = (double*)calloc(N_atoms, sizeof(double));
	fh->p_hybr = (dvec*)calloc(N_atoms, sizeof(dvec));
	fh->m_cross_term = (double*)calloc(N_atoms, sizeof(double));
	fh->p_cross_term = (dvec*)calloc(N_atoms, sizeof(dvec));
	fh->m_cross_term_prev = (double*)calloc(N_atoms, sizeof(double));
	fh->p_cross_term_prev = (dvec*)calloc(N_atoms, sizeof(dvec));
	fh->summu2 = (double*)calloc(N_atoms, sizeof(double));
    for(int n = 0; n < N_atoms; n++) {
		fh->v_md[n][0] = v[n][0];
		fh->v_md[n][1] = v[n][1];
		fh->v_md[n][2] = v[n][2];

		fh->x_hybr[n][0] = x[n][0] * fh->hybr_scale_lin;
		fh->x_hybr[n][1] = x[n][1] * fh->hybr_scale_lin;
		fh->x_hybr[n][2] = x[n][2] * fh->hybr_scale_lin;

		fh->v_hybr[n][0] = v[n][0] * (1.0 - fh->S);// / (fh->hybr_scale_lin);
		fh->v_hybr[n][1] = v[n][1] * (1.0 - fh->S);// / (fh->hybr_scale_lin);
		fh->v_hybr[n][2] = v[n][2] * (1.0 - fh->S);// / (fh->hybr_scale_lin);

        // fh->a_hybr_prev[n][0] = f[n][0] * (1.0 - fh->S) * invmass[n];// / (fh->hybr_scale_lin);
		// fh->a_hybr_prev[n][1] = f[n][1] * (1.0 - fh->S) * invmass[n];// / (fh->hybr_scale_lin);
		// fh->a_hybr_prev[n][2] = f[n][2] * (1.0 - fh->S) * invmass[n];// / (fh->hybr_scale_lin);

//		fh->v_hybr[n][0] = fh->v_hybr[n][1] = fh->v_hybr[n][2] = 0;

    	fh->m_cross_term_prev[n] = 0;
    	fh->p_cross_term_prev[n][0] = fh->p_cross_term_prev[n][1] = fh->p_cross_term_prev[n][2] = 0;
	}


/*	int nX=NX, nY=NY, nZ=NZ;

//	VectorArray<double> X ( nX*nY*nZ, Vector<double>(DIM) );
	double Min[DIM], Max[DIM];
	int np=0;
//	Vector<double> x(3);

	std::vector<Eigen::Vector3d> X ( nX*nY*nZ );
	Eigen::Vector3d x(3);

	for (int i=0; i<nX; i++) {
		for (int j=0; j<nY; j++) {
			for (int k=0; k<nZ; k++) {

				ivec ind;
				ASSIGN_IND(ind, i, j, k);

				x[0] = fh->grid.n[C][0];
				x[1] = fh->grid.n[C][1];
				x[2] = fh->grid.n[C][2];

				if ( k==nZ-1 ) x[2] += fh->box[2];
				if ( j==nY-1 ) x[1] += fh->box[1];
				if ( i==nX-1 ) x[0] += fh->box[0];

				X[np++] = x;

				if ( k==0 ) Min[2] = x[2];
				if ( k==nZ-1 ) Max[2] = x[2];
			}
			if ( j==0 ) Min[1] = x[1];
			if ( j==nY-1 ) Max[1] = x[1];
		}
		if ( i==0 ) Min[0] = x[0];
		if ( i==nX-1 ) Max[0] = x[0];
	}
//printf("create tessel: %f %f\n",X[2][2], Max[1]);
	    fh->tssl = new Tesselation ( X, Min, Max, fh->N );
//	    fh->tssl = tessel_create (fh->grid, fh->N, fh->box);
//printf("tessel created\n");
*/

	int nX=fh->N[0]+1, nY=fh->N[1]+1, nZ=fh->N[2]+1;

	std::vector<Eigen::Vector3d> X ( nX*nY*nZ );
	std::vector<Eigen::Vector3d> Shift ( nX*nY*nZ );

	double Min[DIM], Max[DIM];
	int np=0;

    {
        Eigen::Vector3d x(3), step(3), shift(3);

        step[0] = fh->fh_scale_lin * fh->grid.h[I3(nX-2,0,0,fh->N)][0];
        step[1] = fh->fh_scale_lin * fh->grid.h[I3(0,nY-2,0,fh->N)][1];
        step[2] = fh->fh_scale_lin * fh->grid.h[I3(0,0,nZ-2,fh->N)][2];

        for (int i=0; i<nX; i++) {

            shift[0] = 0.00001;
            if ( i==nX-1 ) shift[0] = -0.00001;

            for (int j=0; j<nY; j++) {

                shift[1] = 0.00001;
                if ( j==nY-1 ) shift[1] = -0.00001;

                for (int k=0; k<nZ; k++) {

                    shift[2] = 0.00001;
                    if ( k==nZ-1 ) shift[2] = -0.00001;
                    Shift[np] = shift;

                    x[0] = fh->fh_scale_lin * fh->grid.n[ I3 ( i==nX-1?i-1:i, j==nY-1?j-1:j, k==nZ-1?k-1:k, fh->N ) ][0];
                    x[1] = fh->fh_scale_lin * fh->grid.n[ I3 ( i==nX-1?i-1:i, j==nY-1?j-1:j, k==nZ-1?k-1:k, fh->N ) ][1];//fhmd.grid.n[I3(i,j,k,fhmd.N)][1];
                    x[2] = fh->fh_scale_lin * fh->grid.n[ I3 ( i==nX-1?i-1:i, j==nY-1?j-1:j, k==nZ-1?k-1:k, fh->N ) ][2];//fhmd.grid.n[I3(i,j,k,fhmd.N)][2];

                    if ( k==nZ-1 ) x[2] += step[2];
                    if ( j==nY-1 ) x[1] += step[1];
                    if ( i==nX-1 ) x[0] += step[0];

                    X[np++] = x;

                    if ( k==0 ) Min[2] = x[2];
                    if ( k==nZ-1 ) Max[2] = x[2];
                }
                if ( j==0 ) Min[1] = x[1];
                if ( j==nY-1 ) Max[1] = x[1];
            }

            if ( i==0 ) Min[0] = x[0];
            if ( i==nX-1 ) Max[0] = x[0];
        }
    }

	int N[3];
	N[0]  = fh->N[0] + 1;
	N[1]  = fh->N[1] + 1;
	N[2]  = fh->N[2] + 1;
//	Tesselation tssl ( X, Min, Max, fh->N );
printf("tssl: %f %f %f %f %f %f\n", Min[0], Min[1], Min[2], Max[0], Max[1], Max[2]);

    int nth = gmx_omp_nthreads_get(emntUpdate);
    fh->tssl = (Tesselation**)calloc(nth, sizeof(Tesselation*));
    for (int i=0; i<nth; i++)
	    fh->tssl[i] = new Tesselation ( X, Min, Max, fh->N );

//printf("tssl: %d, fh: %d\n", fh->tssl->GetNpoints(), fh->Ntot);

    for(int i = 0; i < fh->Ntot; i++) {
    	fh->arr[i].ro_fh = 0;
    	fh->arr[i].u_fh[0] = fh->arr[i].u_fh[1] = fh->arr[i].u_fh[2] = 0;
	}

    fh->tssl_npoints = fh->tssl[0]->GetNpoints();
    fh->cells_map = (int*)calloc(fh->tssl_npoints, sizeof(int));
	for (int i=0; i<N[0]; i++)
		for (int j=0; j<N[1]; j++)
			for (int k=0; k<N[2]; k++)
			{
				int ind = i + j*N[0] + k*N[0]*N[1];
				int ii = i%(fh->N[0]);
				int jj = j%(fh->N[1]);
				int kk = k%(fh->N[2]);  //  beware of ijk order!!
				fh->cells_map[ind] = kk + jj*fh->N[0] + ii*fh->N[0]*fh->N[1];
			}

            FILE *mass_out;
            mass_out = fopen("mass.dat", "a");
            
    for (int n = 0; n < fh->Ntot; n++) {
		fh->arr[n].avg_mass_fh_particle = fh->box_volume * fh->FH_dens / N_atoms;
        // fh->arr[n].avg_mass_fh_particle_prev = fh->arr[n].avg_mass_fh_particle;
        fh->arr[n].vol_star = 0;

        // fprintf(mass_out, "%f ", fh->arr[n].avg_mass_fh_particle);

	}

	// fclose(mass_out);
	for (int q=0; q<fh->tssl_npoints; q++) {
		
		int k = fh->cells_map[q];
//        fh->arr[k].num_star += fh->tssl->GetNSimplices(q);
        fh->arr[k].vol_star += fh->tssl[0]->GetStarVolume(q);
	}

    
    for (int q=0; q<fh->tssl_npoints; q++) {
        int k = fh->cells_map[q];
        fh->arr_prev[k].vol_star = fh->arr[k].vol_star;
        fh->arr_next[k].vol_star = fh->arr[k].vol_star;
    }

//    for (int n = 0; n < fh->Ntot; n++)
//        fprintf(stderr,"%f ",fh->arr[n].vol_star);


//	for (int i=0; i<27; i++)
//		fh->cells[i] = (int) (rand()/(RAND_MAX+1.0) * N_atoms);

/*    for(int k = 0; k < 3; k++)
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
            {
                int ind = 3*(3*i+j)+k;
                printf("%.3f %.3f %.3f ", fh->grid.n[ind][0], fh->grid.n[ind][0], fh->grid.n[ind][0]);
			}
*/			
/*
	FILE *fw1;
    char fname[32];

    sprintf(fname, "mu_map.dat");
    if((fw1 = fopen(fname, "w")) == NULL)  {  printf("\n ERROR creating %s for output!\n", fname);  }

    for(int i = 0; i < 15; i++)
	    for(int j = 0; j < 15; j++)
	        for(int k = 0; k < 15; k++)
	        {
            int ind = 15*(15*i+j)+k;

	        dvec r;
	        r[0] = fh->box[0]*i/15.0 + 0.000001;
	        r[1] = fh->box[1]*j/15.0 + 0.000001;
	        r[2] = fh->box[2]*k/15.0 + 0.000001;
			
			fh->coords[ind][0] = r[0];
			fh->coords[ind][1] = r[1];
			fh->coords[ind][2] = r[2];

			int n2 = 0;
			for (int q=0; q<fh->tssl_npoints; q++)
			    if ( fh->tssl->in_this_cell(q,r) )
			    {
					fh->cells[ind][n2] = q;
					fh->cells_mu[ind][n2] = fh->tssl->mu();
					n2++;
				}	        
			if (n2!=4)
			  printf("not 4 neighs!! ind=%d, n=%d\n",ind,n2);
			fprintf(fw1,"%.4f %.4f %.4f   %.4f %.4f %.4f %.4f\n",r[0],r[1],r[2],
				fh->cells_mu[ind][0],fh->cells_mu[ind][1],fh->cells_mu[ind][2],fh->cells_mu[ind][3]);
	        }
			
	fclose(fw1);		
*/
/*
 * end of hhmd_G initialising
 */

    if(fh->scheme == One_Way)
        FH_init(fh, cr);


   if (PAR(cr))
       gmx_barrier(cr);
       
       

    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->arr[i].Sc = 1;
        for(int d = 0; d < DIM; d++)
        {
            fh->arr[i].Sf[d] = 1;
        }
    }


    if(MASTER(cr))
    {
        if(fh->scheme == One_Way)
        {
            FH_equilibrate(fh);
            printf(MAKE_GREEN "FHMD: Initialization finished. Starting MD/FH solver...\n" RESET_COLOR "\n");
        }
        fflush(stdout);
    }
    for (int n = 0; n < fh->Ntot; n++) {
        // fh->arr[n].avg_mass_fh_particle = fh->box_volume * fh->FH_dens / N_atoms;
        // fh->arr[n].avg_mass_fh_particle_prev = fh->arr[n].avg_mass_fh_particle;
        // fh->arr[n].vol_star = 0;
        fprintf(mass_out, "%f ", fh->arr[n].avg_mass_fh_particle);

    }
    fclose(mass_out);

   return 1;   // Success
}
