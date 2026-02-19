#include "data_structures.h"
#include "macro.h"


#define write_dump(var, name) \
    sprintf(fname, "dump_%s.txt", name); \
    if(fh->step_MD == 0) { \
        fw = fopen(fname, "w"); \
        write_header(fw, fh); \
    } else { \
        fw = fopen(fname, "a"); \
    } \
    fprintf(fw, "\n%d\t", fh->step_MD); \
    for(int i = 0; i < fh->Ntot; i++) fprintf(fw, "%g\t", fh->arr[i].var); \
    fclose(fw);


void write_header(FILE *fw, FHMD *fh)
{
    fprintf(fw, "step\t");

    for(int k = 0; k < NZ; k++)
        for(int j = 0; j < NY; j++)
            for(int i = 0; i < NX; i++)
                fprintf(fw, "cell %d-%d-%d\t", i, j, k);
}


// void fhmd_dump_all(FHMD *fh)
// {
//     FILE *fw2;
//     char fname[32];
//     sprintf(fname, "rho_fh.dat");                 // Tecplot data filename
//     if((fw2 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     fprintf(fw2,"%d  ",fh->step_MD);

//     for (int i=0; i<fh->Ntot; i++)
//       fprintf(fw2,"%.4f ",fh->arr[i].ro_fh);

      

//     fprintf(fw2,"\n");
//     fclose(fw2);

// return;

//     FILE *fw1,/* *fw2 , */ *fw3,*fw4,*fw5,*fw6,*fw7,*fw8,*fw9;
    
// //    char fname[32];

//     sprintf(fname, "rho_md.dat");                 // Tecplot data filename
//     if((fw1 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     sprintf(fname, "rho_fh.dat");                 // Tecplot data filename
//     if((fw2 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     sprintf(fname, "ux_md.dat");                 // Tecplot data filename
//     if((fw3 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     sprintf(fname, "uy_md.dat");                 // Tecplot data filename
//     if((fw4 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     sprintf(fname, "uz_md.dat");                 // Tecplot data filename
//     if((fw5 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     sprintf(fname, "ux_fh.dat");                 // Tecplot data filename
//     if((fw6 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     sprintf(fname, "uy_fh.dat");                 // Tecplot data filename
//     if((fw7 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
//     sprintf(fname, "uz_fh.dat");                 // Tecplot data filename
//     if((fw8 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    
// //    sprintf(fname, "ro_mixed.dat");                 // Tecplot data filename
// //    if((fw9 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    

//     fprintf(fw1,"%d  ",fh->step_MD);
//     fprintf(fw2,"%d  ",fh->step_MD);
//     fprintf(fw3,"%d  ",fh->step_MD);
//     fprintf(fw4,"%d  ",fh->step_MD);
//     fprintf(fw5,"%d  ",fh->step_MD);
//     fprintf(fw6,"%d  ",fh->step_MD);
//     fprintf(fw7,"%d  ",fh->step_MD);
//     fprintf(fw8,"%d  ",fh->step_MD);
// //    fprintf(fw9,"%d  ",fh->step_MD);
//     for (int i=0; i<fh->Ntot; i++) {
//       fprintf(fw1,"%.4f  ",fh->arr[i].ro_md);
//       fprintf(fw2,"%.4f  ",fh->arr[i].ro_fh);
//       fprintf(fw3,"%.5f  ",fh->arr[i].u_md[0]);
//       fprintf(fw4,"%.5f  ",fh->arr[i].u_md[1]);
//       fprintf(fw5,"%.5f  ",fh->arr[i].u_md[2]);
//       fprintf(fw6,"%.5f  ",fh->arr[i].u_fh[0]);
//       fprintf(fw7,"%.5f  ",fh->arr[i].u_fh[1]);
//       fprintf(fw8,"%.5f  ",fh->arr[i].u_fh[2]);
// //      fprintf(fw9,"%.4f  ",fh->arr[i].ro_mixed);
//     }
//     fprintf(fw1,"\n");
//     fprintf(fw2,"\n");
//     fprintf(fw3,"\n");
//     fprintf(fw4,"\n");
//     fprintf(fw5,"\n");
//     fprintf(fw6,"\n");
//     fprintf(fw7,"\n");
//     fprintf(fw8,"\n");
// //    fprintf(fw9,"\n");
    
//     fclose(fw1);
//     fclose(fw2);
//     fclose(fw3);
//     fclose(fw4);
//     fclose(fw5);
//     fclose(fw6);
//     fclose(fw7);
//     fclose(fw8);
// //    fclose(fw9);
// }

void fhmd_dump_all2(FHMD *fh)
{
    FILE *fw1;
    
    char fname[256];

    // sprintf(fname, "grompp_output/hybr_val.dat");                 // Tecplot data filename
    snprintf(fname, sizeof(fname),
         "output/hybr_val_S%.3f_tau%.3f_fluct%.3f.dat",
         fh->S, fh->tau_m, fh->fluct_m);
    if((fw1 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    

    // fprintf(fw1,"step %d\n",fh->step_MD);
    
    for (int i=0; i<fh->Ntot; i++)
    //   fprintf(fw1,"%.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n",fh->arr[i].natoms_star,fh->arr[i].hybr_mass_star,
    //   	fh->arr[i].hybr_vel_star[0],fh->arr[i].hybr_vel_star[1],fh->arr[i].hybr_vel_star[2],fh->arr[i].vol_star);
    fprintf(fw1, "%d %d %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f %.4f  %.4f  %.4f %.4f %.4f  %.4f  %.4f %.4f  %.4f  %.4f  %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
        fh->step_MD,
        i,
        fh->arr[i].natoms_star, 
        fh->arr[i].hybr_mass_star,
        fh->arr[i].hybr_mass_star_without_mu,
        fh->arr[i].hybr_vel_star[0], 
        fh->arr[i].hybr_vel_star[1], 
        fh->arr[i].hybr_vel_star[2], 

        fh->arr[i].hybr_momentum_star[0], 
        fh->arr[i].hybr_momentum_star[1], 
        fh->arr[i].hybr_momentum_star[2],


        fh->arr[i].hybr_vel_star_without_mu[0], 
        fh->arr[i].hybr_vel_star_without_mu[1], 
        fh->arr[i].hybr_vel_star_without_mu[2], 
        fh->arr[i].vol_star,
        fh->arr[i].u_fh[0],  // u_fh[0]
        fh->arr[i].u_fh[1],  // u_fh[1]
        fh->arr[i].u_fh[2],
        fh->arr[i].ro_fh,
        fh->arr[i].tot_mu,

        fh->arr[i].mu2m2_sum, 
        fh->arr[i].mu2m_sum, 


        fh->arr[i].hybr_mass_star_1,
        fh->arr[i].hybr_mass_star_2,
        fh->arr[i].hybr_mass_star_3,
        fh->arr[i].hybr_mass_star_4); 
        // fprintf(fw1,"%.4f ", 4.0*fh->arr[i].hybr_mass_star/fh->arr[i].vol_star);
        // fprintf(fw1, "\n");
    
	fclose(fw1);

/*    for (int i=0; i<fh->Ntot; i++) {
        fprintf(fw1,"%.4f  ",fh->arr[i].ro_mixed);
    }
*/
/*	for (int k=0; k<10; k++)
      fprintf(fw1,"%.4f  ", fh->ro_hybr_avg[k]);
	for (int k=0; k<10; k++)
      fprintf(fw1,"%.4f  ", fh->ro_hybr_sd[k]);
*/
//	for (int i=0; i<270; i++)
//        fprintf(fw1,"%.4f  ", fh->ro_hybr[fh->cells[i]]);
//	for (int k=0; k<10; k++) {
int cells[27] = 
// cell centers:
{1780,4275,3079,276,878,4244,1198,1887,1457,2821,1374,2190,1738,2347,2768,4055,48,146,3516,2610,3617,3008,2166,658,221,2840,3591};
// cell corners:
//{187,210,315,354,599,739,741,988,1156,2240,2259,2275,2414,2477,2780,2782,4439,2530,3231,3313,94,1018,4134,2725,477,1825,1571};
//  cell corners 2:
//{94,187,210,315,354,477,599,739,741,988,1018,1156,1571,1825,2240,2259,2275,2414,2477,2530,2725,2780,2782,3231,3313,4134,4439};

/*
int k=0;
	{
		for (int i=0+k*27; i<27+k*27; i++)
			cells[i]--;
		double ro_hybr_avg = 0;
		for (int i=0+k*27; i<27+k*27; i++)
			ro_hybr_avg += fh->ro_hybr[fh->cells[i]];
		ro_hybr_avg /= 27.0;
//        fprintf(fw1,"%.4f  ", ro_hybr_avg);

		double ro_hybr_sd = 0;
		for (int i=0+k*27; i<27+k*27; i++)
			ro_hybr_sd += (ro_hybr_avg - fh->ro_hybr[cells[i]])*(ro_hybr_avg - fh->ro_hybr[cells[i]]);
//        fprintf(fw1,"%.4f  ", ro_hybr_sd);
		ro_hybr_sd = sqrt( ro_hybr_sd / (27.0 - 1.0) );
        fprintf(fw1,"%.4f  ", ro_hybr_sd);
		}
    fprintf(fw1,"\n");
    
    fclose(fw1);
    
    sprintf(fname, "ro_mixed_all.dat");                 // Tecplot data filename
    if((fw1 = fopen(fname, "a")) == NULL)    {        printf("\n ERROR creating %s for output!\n", fname);        return;    }    

    fprintf(fw1,"%d  ",fh->step_MD);
	for (int i=0; i<27; i++)
        fprintf(fw1,"%.4f  ", fh->ro_hybr[fh->cells[i]]);
    fprintf(fw1,"\n");
    
    fclose(fw1);
*/}
void fhmd_dump_all(FHMD *fh)
{
    FILE *fw;
    char  fname[64];

    write_dump(ro_md,   "ro_md");
    write_dump(uro_md[0],   "uro_md_x");
    write_dump(uro_md[1],   "uro_md_y");
    write_dump(uro_md[2],   "uro_md_z");
    write_dump(ro_fh,   "ro_fh");
    write_dump(u_md[0], "u_md_X");
    write_dump(u_md[1], "u_md_Y");
    write_dump(u_md[2], "u_md_Z");
    write_dump(u_fh[0], "u_fh_X");
    write_dump(u_fh[1], "u_fh_Y");
    write_dump(u_fh[2], "u_fh_Z");


}
void fhmd_dump_all_prev(FHMD *fh)
{
    FILE *fw;
    char  fname[64];

    // write_dump(ro_md,   "ro_md");
    write_dump(ro_fh,   "ro_fh");
    // write_dump(u_md[0], "u_md_X");
    // write_dump(u_md[1], "u_md_Y");
    // write_dump(u_md[2], "u_md_Z");
    write_dump(u_fh[0], "u_fh_X");
    write_dump(u_fh[1], "u_fh_Y");
    write_dump(u_fh[2], "u_fh_Z");

// #ifdef FHMD_DEBUG
//     write_dump(f_fh[0],       "f_fh_X");
//     write_dump(f_fh[1],       "f_fh_Y");
//     write_dump(f_fh[2],       "f_fh_Z");
//     write_dump(alpha_term[0], "alpha_term_X");
//     write_dump(alpha_term[1], "alpha_term_Y");
//     write_dump(alpha_term[2], "alpha_term_Z");
//     write_dump(beta_term[0],  "beta_term_X");
//     write_dump(beta_term[1],  "beta_term_Y");
//     write_dump(beta_term[2],  "beta_term_Z");
//     write_dump(grad_ro[0],    "grad_ro_X");
//     write_dump(grad_ro[1],    "grad_ro_Y");
//     write_dump(grad_ro[2],    "grad_ro_Z");

//     write_dump(ro_prime,      "ro_prime");
//     write_dump(ro_star,       "ro_star");
//     write_dump(m_prime[0],    "m_prime_X");
//     write_dump(m_prime[1],    "m_prime_Y");
//     write_dump(m_prime[2],    "m_prime_Z");
//     write_dump(m_star[0],     "m_star_X");
//     write_dump(m_star[1],     "m_star_Y");
//     write_dump(m_star[2],     "m_star_Z");

//     write_dump(S,             "S");
//     write_dump(Sf[0],         "Sf_X");
//     write_dump(Sf[1],         "Sf_Y");
//     write_dump(Sf[2],         "Sf_Z");
// #endif
}


void writePoint(FILE *fout, FHMD *fh, const char *line, char ch)
{
    int c = 0;
    ivec ind;

    fprintf(fout, "%s", line);

    for(int k = 0; k <= NZ; k++) {
        for(int j = 0; j <= NY; j++) {
            for(int i = 0; i <= NX; i++) {
                ASSIGN_IND(ind, i - (int)(i/NX), j - (int)(j/NY), k - (int)(k/NZ));
                     if(ch == 'X' && i < NX)  fprintf(fout, "%e ", fh->grid.n[C][0]);
                else if(ch == 'X' && i == NX) fprintf(fout, "%e ", fh->grid.n[C][0] + fh->grid.h[C][0]);
                else if(ch == 'Y' && j < NY)  fprintf(fout, "%e ", fh->grid.n[C][1]);
                else if(ch == 'Y' && j == NY) fprintf(fout, "%e ", fh->grid.n[C][1] + fh->grid.h[C][1]);
                else if(ch == 'Z' && k < NZ)  fprintf(fout, "%e ", fh->grid.n[C][2]);
                else if(ch == 'Z' && k == NZ) fprintf(fout, "%e ", fh->grid.n[C][2] + fh->grid.h[C][2]);
                if(c++ == 10) {fprintf(fout, "\n"); c = 0;}
            }
        }
    }
}


void writeCons(FILE *fout, FHMD *fh, const char *line, char ch)
{
    int c = 0;
    ivec ind;

    fprintf(fout, "%s", line);

    for(int k = 0; k < NZ; k++) {
        for(int j = 0; j < NY; j++) {
            for(int i = 0; i < NX; i++) {
                ASSIGN_IND(ind, i, j, k);
                     if(ch == 'U') fprintf(fout, "%e ", fh->arr[C].u_fh[0]);
                else if(ch == 'V') fprintf(fout, "%e ", fh->arr[C].u_fh[1]);
                else if(ch == 'W') fprintf(fout, "%e ", fh->arr[C].u_fh[2]);
                else if(ch == 'R') fprintf(fout, "%e ", fh->arr[C].ro_fh);
                if(c++ == 10) {fprintf(fout, "\n"); c = 0;}
            }
        }
    }
}


// TODO: This function is a copy/paste from the previous code -- should be rewritten completely
void fhmd_write_tecplot_data(FHMD *fh, int step, double time)
{
    FILE *fout;                                                     // Tecplot data file

    static int zc = 0;

    char fname[32];

    sprintf(fname, "tecplot/data_%7.7d.dat", step);                 // Tecplot data filename

    if((fout = fopen(fname, "w")) == NULL)                          // Open file
        printf("\n ERROR creating %s for output!\n", fname);

    zc++;

    fprintf(fout, "TITLE=\"OUT\"\nVARIABLES=\"X\",\"Y\",\"Z\",\"U_tilde\",\"V_tilde\",\"W_tilde\",\"RHO_tilde\"\n");
    fprintf(fout, "ZONE T=\"%f\", I=%d, J=%d, K=%d, DATAPACKING=BLOCK\nVARLOCATION=([4-7]=CELLCENTERED)\n", time, NX+1, NY+1, NZ+1);
    fprintf(fout, "STRANDID=%d\nSOLUTIONTIME=%g", zc, time);

    writePoint(fout, fh, "\n\n# X:\n", 'X');                        // Write X
    writePoint(fout, fh, "\n\n# Y:\n", 'Y');                        // Write Y
    writePoint(fout, fh, "\n\n# Z:\n", 'Z');                        // Write Z

    writeCons(fout, fh, "\n\n# U:\n",   'U');                       // Write U
    writeCons(fout, fh, "\n\n# V:\n",   'V');                       // Write V
    writeCons(fout, fh, "\n\n# W:\n",   'W');                       // Write W
    writeCons(fout, fh, "\n\n# RHO:\n", 'R');                       // Write RHO

    fclose(fout);
}


#ifdef FHMD_PARAVIEW
void write_paraview_data(HHMD *hhmd) {

    const int ln = 5;       // Numbers per one line

    FILE *fout;
    char fname[32];
    static int zc = 0;
    int lc = 1;

    int nx = hhmd->cellNum.x;
    int ny = hhmd->cellNum.y;
    int nz = hhmd->cellNum.z;

    sprintf(fname, "paraview_%d.vtk", zc);

    if((fout = fopen(fname, "w")) == NULL)
        printf("\n ERROR creating %s for output!\n", fname);

    zc++;

    fprintf(fout, "# vtk DataFile Version 3.0\nData Output\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS %d %d %d\n", nx + 1, ny + 1, nz + 1);
    fprintf(fout, "POINTS %d float\n", (nx + 1)*(ny + 1)*(nz + 1));

    for(int k = 0; k < nz + 1; k++) {
        for(int j = 0; j < ny + 1; j++) {
            for(int i = 0; i < nx + 1; i++) {
                fprintf(fout, "%f %f %f\n", hhmd->grid.n[i].x, hhmd->grid.n[j].y, hhmd->grid.n[k].z);
            }
        }
    }

    fprintf(fout, "CELL_DATA %d\n", nx*ny*nz);
    fprintf(fout, "FIELD FieldData %d\n", 1);
    fprintf(fout, "Density %d %d float\n", 1, nx*ny*nz);

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                fprintf(fout, "%f", hhmd->arr->ro_fh[i][j][k]);
                if(lc++ >= ln) {
                    fprintf(fout, "\n");
                    lc = 1;
                } else {
                    fprintf(fout, " ");
                }
            }
        }
    }

    fprintf(fout, "VECTORS Velocity float\n");

    lc = 1;

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                fprintf(fout, "%f %f %f", hhmd->arr->u_fh[i][j][k].x, hhmd->arr->u_fh[i][j][k].y, hhmd->arr->u_fh[i][j][k].z);
                if(lc++ >= ln) {
                    fprintf(fout, "\n");
                    lc = 1;
                } else {
                    fprintf(fout, "   ");
                }
            }
        }
    }

    fprintf(fout, "POINT_DATA %d\n", (nx + 1)*(ny + 1)*(nz + 1));

    fclose(fout);

}
#endif


