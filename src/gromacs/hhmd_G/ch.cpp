/* ch.c : numerical functions for hull computation */

/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */




#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <string.h>


#include "hull.h"

static short check_overshoot_f=0;

double logb(double x) {  /* on SGI machines: returns floor of log base 2 */
	if (x<=0) return -1e38;
	return log(x)/log(2.0);
}

#define SWAP(X,a,b) {X t; t = a; a = b; b = t;}

#define check_overshoot(x)							\
	{if (check_overshoot_f && ((x)>9e15))		\
		printf("overshot exact arithmetic");}			\


#define DELIFT 0

h_Coord ClarksonHull::Vec_dot(h_point x, h_point y) {
	int i;
	h_Coord sum = 0;
	for (i=0;i<rdim;i++) sum += x[i] * y[i];
	return sum;
}

h_Coord ClarksonHull::Vec_dot_pdim(h_point x, h_point y) {
	int i;
	h_Coord sum = 0;
	for (i=0;i<pdim;i++) sum += x[i] * y[i];
	return sum;
}

h_Coord ClarksonHull::Norm2(h_point x) {
	int i;
	h_Coord sum = 0;
	for (i=0;i<rdim;i++) sum += x[i] * x[i];
	return sum;
}

void ClarksonHull::Ax_plus_y(h_Coord a, h_point x, h_point y) {
	int i;
	for (i=0;i<rdim;i++) {
		*y++ += a * *x++;
	}
}

void ClarksonHull::Ax_plus_y_test(h_Coord a, h_point x, h_point y) {
	int i;
	for (i=0;i<rdim;i++) {
		check_overshoot(*y + a * *x);
		*y++ += a * *x++;
	}
}

void ClarksonHull::Vec_scale_test(int n, h_Coord a, h_Coord *x)
{
	h_Coord *xx = x,
		*xend = xx + n	;

	if (check_overshoot_f && ((a)>9e15))
		printf("overshot exact arithmetic");

	while (xx!=xend) {
		*xx *= a;
		if (check_overshoot_f && ((*xx)>9e15))
			printf("overshot exact arithmetic");
		xx++;
	}
}


#define VA(x) ((x)->vecs+rdim)
#define VB(x) ((x)->vecs)




/* tables for runtime stats */
//static int A[100]={0}, B[100] ={0}, C[100] = {0}, D[100] ={0};

double ClarksonHull::sc(h_basis_s *v,h_simplex *s, int k, int j) {
/* amount by which to scale up vector, for reduce_inner */

	double		labound;

	if (j<10) {
		labound = logb(v->sqa)/2;  //!!!!!!!!!!!1
		sc_max_scale = exact_bits - labound - 0.66*(k-2)-1  -0;
		if (sc_max_scale<1) {
			printf("overshot exact arithmetic");
			sc_max_scale = 1;

		}

		if (j==0) {
			int	i;
			h_neighbor *sni;
			h_basis_s *snib;

			sc_ldetbound = 0;

			sc_Sb = 0;
			for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
				snib = sni->basis;
				sc_Sb += snib->sqb;
				sc_ldetbound += logb(snib->sqb)/2 + 1;
				sc_ldetbound -= snib->lscale;
			}
		}
	}
	if (sc_ldetbound - v->lscale + logb(v->sqb)/2 + 1 < 0) {
		return 0;					
	} else {
		sc_lscale = (int)logb(2*sc_Sb/(v->sqb + v->sqa*b_err_min))/2;	
		if (sc_lscale > sc_max_scale) {
			sc_lscale = (int)sc_max_scale;
		} else if (sc_lscale<0) sc_lscale = 0;
		v->lscale += sc_lscale;
		return (sc_lscale<20) ? 1<<sc_lscale : ldexp(1.0,sc_lscale);
	}
}


int ClarksonHull::reduce_inner(h_basis_s *v, h_simplex *s, int k) {

	h_point	va = v->vecs+rdim,
		vb = v->vecs;
	int	i,j;
	double	dd;//,
//		ldetbound = 0,
//		Sb = 0;
	double	scale;
	h_basis_s	*snibv;
	h_neighbor *sni;

	v->sqa = v->sqb = Norm2(vb);
	if (k<=1) {
		memcpy(vb,va,h_basis_vec_size);
		return 1;
	}

	for (j=0;j<250;j++) {

		memcpy(vb,va,h_basis_vec_size);
		for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
			snibv = sni->basis;
			dd = -Vec_dot(snibv->vecs,vb)/ snibv->sqb;
			Ax_plus_y( dd, snibv->vecs+rdim, vb);		
		}
		v->sqb = Norm2(vb);		
		v->sqa = Norm2(va);
		
		if (2*v->sqb >= v->sqa) {/*B[j]++;*/ return 1;}

		Vec_scale_test(rdim,scale = sc(v,s,k,j),va);
		
		for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
			snibv = sni->basis;
			dd = -Vec_dot(VB(snibv),va)/snibv->sqb;	
			dd = floor(dd+0.5);
			Ax_plus_y_test( dd, VA(snibv), va);
		}		
	}
	return 0;
}

#define trans(z,p,q) {int i; for (i=0;i<pdim;i++) z[i+rdim] = z[i] = p[i] - q[i];}
#define lift(z,s) {z[2*rdim-1] =z[rdim-1]= ldexp(Vec_dot_pdim(z,z), -DELIFT);}
				/*not scaling lift to 2^-DELIFT */



int ClarksonHull::reduce(h_basis_s **v, h_point p, h_simplex *s, int k) {

	h_point	z;
	h_point	tt = s->neigh[0].vert;

	if (!*v) {
		*v = h_basis_s_list ? h_basis_s_list : new_block_h_basis_s(1);
		assert(*v);
		h_basis_s_list = (*v)->next;
		(*v)->ref_count = 1;
	}
	else (*v)->lscale = 0;
	z = (*v)->vecs;
	if (p==infinity) memcpy(*v,infinity_basis,h_basis_s_size);
	else {
		trans(z,p,tt); 
		lift(z,s);
	}
	return reduce_inner(*v,s,k);
}


void ClarksonHull::get_basis_sede(h_simplex *s) {

	int	k=1;
	h_neighbor *sn = s->neigh+1,
		 *sn0 = s->neigh;

	if (sn0->vert == infinity && cdim >1) {
		{h_neighbor t; t = *sn0; *sn0 = *sn; *sn = t;}
		if ((sn0->basis) && --(sn0->basis)->ref_count == 0) {
			memset(sn0->basis,0,h_basis_s_size);
			(sn0->basis)->next = h_basis_s_list;
			h_basis_s_list = sn0->basis;
		}
		sn0->basis = NULL;
		sn0->basis = tt_basisp;
		tt_basisp->ref_count++;
	} else {
		if (!sn0->basis) {
			sn0->basis = tt_basisp;
			tt_basisp->ref_count++;
		} else while (k < cdim && sn->basis) {k++;sn++;}
	}
	while (k<cdim) {
		if ((sn->basis) && --(sn->basis)->ref_count == 0) {
			memset(sn->basis,0,h_basis_s_size);
			(sn->basis)->next = h_basis_s_list;
			h_basis_s_list = sn->basis;
		}
		sn->basis = NULL;
		reduce(&sn->basis,sn->vert,s,k);
		k++; sn++;
	}
}


int ClarksonHull::out_of_flat(h_simplex *root, h_point p) {

	if (!out_p_neigh.basis) out_p_neigh.basis = (h_basis_s*)memory(h_basis_s_size);

	out_p_neigh.vert = p;
	cdim++;
	root->neigh[cdim-1].vert = root->peak.vert;
	if ((root->neigh[cdim-1].basis) && --(root->neigh[cdim-1].basis)->ref_count == 0) {
		memset(root->neigh[cdim-1].basis,0,h_basis_s_size);
		(root->neigh[cdim-1].basis)->next = h_basis_s_list;
		h_basis_s_list = root->neigh[cdim-1].basis;
	}
	root->neigh[cdim-1].basis = NULL;

	get_basis_sede(root);
	if (root->neigh[0].vert == infinity) return 1;
	reduce(&out_p_neigh.basis,p,root,cdim);
	if (out_p_neigh.basis->sqa != 0) return 1;
	cdim--;
	return 0;
}


void ClarksonHull::get_normal_sede(h_simplex *s) {

	h_neighbor *rn;
	int i,j;

	get_basis_sede(s);
	if (rdim==3 && cdim==3) {
		h_point	c,
			a = VB(s->neigh[1].basis),
			b = VB(s->neigh[2].basis);
		NEWLRC(h_basis_s,s->normal);
		c = VB(s->normal);
		c[0] = a[1]*b[2] - a[2]*b[1];
		c[1] = a[2]*b[0] - a[0]*b[2];
		c[2] = a[0]*b[1] - a[1]*b[0];
		s->normal->sqb = Norm2(c);
		for (i=cdim+1,rn = ch_root->neigh+cdim-1; i; i--, rn--) {
			for (j = 0; j<cdim && rn->vert != s->neigh[j].vert;j++);
			if (j<cdim) continue;
			if (rn->vert==infinity) {
				if (c[2] > -b_err_min) continue;
			} else  if (!sees(rn->vert,s)) continue;
			c[0] = -c[0]; c[1] = -c[1]; c[2] = -c[2];
			break;
		}
		return;
	}	
		
	for (i=cdim+1,rn = ch_root->neigh+cdim-1; i; i--, rn--) {
		for (j = 0; j<cdim && rn->vert != s->neigh[j].vert;j++);
		if (j<cdim) continue;
		reduce(&s->normal,rn->vert,s,cdim);
		if (s->normal->sqb != 0) break;
	}

}


int ClarksonHull::sees(h_site p, h_simplex *s) {

	h_point	tt,zz;
	double	dd,dds;
	int i;


	if (!see_b) {see_b = (h_basis_s*)memory(h_basis_s_size);assert(see_b);}
	else see_b->lscale = 0;
	zz = VB(see_b);
	if (cdim==0) return 0;
	if (!s->normal) {
		get_normal_sede(s);
		for (i=0;i<cdim;i++) NULLIFY(h_basis_s,s->neigh[i].basis);
	}
	tt = s->neigh[0].vert;

	if (p==infinity) memcpy(see_b,infinity_basis,h_basis_s_size);
	else {trans(zz,p,tt); lift(zz,s);}

	for (i=0;i<3;i++) {
		dd = Vec_dot(zz,s->normal->vecs);	
		if (dd == 0.0) {
			return 0;
		} 
		dds = dd*dd/s->normal->sqb/Norm2(zz);
		if (dds > b_err_min_sq) return (dd<0);
		get_basis_sede(s);
		reduce_inner(see_b,s,cdim);
	}
	return 0;
}



#define swap_points(a,b) {point t; t=a; a=b; b=t;}

h_simplex *ClarksonHull::build_convex_hull(short dim) {

/*
	get_s		returns next site each call;
			hull construction stops when NULL returned;
	site_numm	returns number of site when given site;
	dim		dimension of point set;
	vdd		if (vdd) then return Delaunay triangulation


*/

	h_simplex *s, *root;

	if (!Huge) Huge = DBL_MAX;

	cdim = 0;
	pdim = dim;

	exact_bits = int(DBL_MANT_DIG*log((double)FLT_RADIX)/log(2.0));
	b_err_min = float(DBL_EPSILON*H_MAXDIM*(1<<H_MAXDIM)*H_MAXDIM*3.01);
	b_err_min_sq = b_err_min * b_err_min;

	rdim = pdim+1;

	h_point_size = h_site_size = sizeof(h_Coord)*pdim;
	h_basis_vec_size = sizeof(h_Coord)*rdim;
	h_basis_s_size = sizeof(h_basis_s)+ (2*rdim-1)*sizeof(h_Coord);
	h_simplex_size = sizeof(h_simplex) + (rdim-1)*sizeof(h_neighbor);
	h_Tree_size = sizeof(h_Tree);
	h_fg_size = sizeof(h_fg);

	root = NULL;

	p = infinity;
	infinity_basis = h_basis_s_list ? h_basis_s_list : new_block_h_basis_s(1);
	assert(infinity_basis);
	h_basis_s_list = infinity_basis->next;
	infinity_basis->ref_count = 1;

	infinity_basis->vecs[2*rdim-1] = infinity_basis->vecs[rdim-1] = 1;
	infinity_basis->sqa = infinity_basis->sqb = 1;

 	root = h_simplex_list ? h_simplex_list : new_block_h_simplex(1);
	assert(root);
 	h_simplex_list = root->next;

	ch_root = root;

 	s = h_simplex_list ? h_simplex_list : new_block_h_simplex(1);
	assert(s);
 	h_simplex_list = s->next;
	memcpy(s,root,h_simplex_size);
	{
		int i;
		h_neighbor *mrsn;

		for (i=-1,mrsn=root->neigh-1;i<cdim;i++,mrsn++)
			inc_ref(basis_s, mrsn->basis);
	}

	root->peak.vert = p;
	root->peak.simp = s;
	s->peak.simp = root;

	buildhull(root);
	return root;
}


void ClarksonHull::free_hull_storage(void) {

	free_h_basis_s_storage();
	free_h_simplex_storage();
	free_h_Tree_storage();
	free_h_fg_storage();
}

