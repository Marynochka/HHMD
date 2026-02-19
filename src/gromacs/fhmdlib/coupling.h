#ifndef FHMD_COUPLING_H_
#define FHMD_COUPLING_H_

void fhmd_update_MD_in_FH(rvec x[], rvec v[], float mass[], rvec f[], int N_atoms, FHMD *fh);
void fhmd_sum_arrays(t_commrec *cr, FHMD *fh);
void fhmd_update_MD_in_FH2(rvec x[], rvec v[], float mass[], rvec f[], int N_atoms, FHMD *fh);
void fhmd_sum_arrays2(t_commrec *cr, FHMD *fh);
//void fhmd_calculate_MDFH_terms(FHMD *fh);
void fhmd_calculate_MDFH_terms(rvec x[], rvec v[], float mass[], rvec f[], int N_atoms, FHMD *fh, int th);

#endif /* FHMD_COUPLING_H_ */
